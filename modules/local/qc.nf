
process FASTQC {
    tag { meta.id }
    label 'qc'
    publishDir "${params.outdir}/qc/fastqc_${fqtype}/${meta.id}", mode: "${params.publish_dir_mode}"

    input:
        tuple val(meta), path(fastq), val(fqtype)
    output:
        path("${meta.id}*.html"), emit: html
        path ("${meta.id}*.zip"), emit: zip

    script:
    """
    fastqc \
      $fastq \
      -t $task.cpus \
      -o .
    """

    stub:
    """
    touch ${meta.id}_fastqc.html ${meta.id}_fastqc.zip
    """
}

process FASTQ_SCREEN {
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(fastq), path(conf)
    output:
        path("${meta.id}*_screen.*"), emit: screen

    script:
    """
    fastq_screen -c $conf --threads $task.cpus $fastq
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${meta.id}.trimmed_screen.\${EXT}
    done
    """
}

process PRESEQ {
    """
    Calls preseq c_curve and lc_extrap, and calls bin/parse_preseq_log.py to get statistics from the log.
    """
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(bam)

    output:
        tuple path("*.c_curve"), path("*.preseq"), path("*.preseqlog"), path("*.txt"), emit: preseq_files

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    preseq c_curve -B -o ${prefix}.c_curve $bam
    preseq lc_extrap -B -D -o ${prefix}.preseq $bam -seed 12345 -v -l 100000000000 2> ${prefix}.preseqlog
    parse_preseq_log.py ${prefix}.preseqlog > ${prefix}.preseqlog.nrf.txt
    """

    stub:
    """
    touch ${meta.id}.c_curve ${meta.id}.preseq ${meta.id}.preseqlog ${meta.id}.preseqlog.nrf.txt
    """
}

process PHANTOM_PEAKS { // https://github.com/kundajelab/phantompeakqualtools
    // TODO: set tmpdir as lscratch if available https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/InitialChIPseqQC.snakefile#L504
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        path("${meta.id}.ppqt.pdf"), emit: pdf
        path("${meta.id}.ppqt"), emit: ppqt

    script: // TODO: for PE, just use first read of each pair
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RUN_SPP=\$(which run_spp.R)
    Rscript \$RUN_SPP -c=${bam} -savp=${prefix}.ppqt.pdf -out=${prefix}.ppqt
    """

    stub:
    """
    touch ${meta.id}.ppqt.pdf ${meta.id}.ppqt
    """
}

process DEDUPLICATE {
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(bam), path(chrom_sizes)

    output:
        tuple val(meta), path("${meta.id}.TagAlign.bed"), emit: tag_align
        tuple val(meta), path("${bam.baseName}.dedup.bam"), emit: bam
        path("${bam.baseName}.dedup.bam.flagstat"), emit: flagstat
        path("${bam.baseName}.dedup.bam.idxstat"), emit: idxstat

    script:
    """
    macs2 filterdup -i ${bam} -g ${params.align.effective_genome_size} --keep-dup="auto" -o TmpTagAlign1.bed
    awk -F"\\t" -v OFS="\\t" '{{if (\$2>0 && \$3>0) {{print}}}}' TmpTagAlign1.bed > TmpTagAlign2.bed
    awk -F"\\t" -v OFS="\\t" '{{print \$1,1,\$2}}' ${chrom_sizes} | sort -k1,1 -k2,2n > GenomeFile.bed
    bedtools intersect -wa -f 1.0 -a TmpTagAlign2.bed -b GenomeFile.bed > ${meta.id}.TagAlign.bed
    bedtools bedtobam -i ${meta.id}.TagAlign.bed -g ${chrom_sizes} | samtools sort -@ ${task.cpus} -o ${bam.baseName}.dedup.bam
    samtools index ${bam.baseName}.dedup.bam
    samtools flagstat ${bam.baseName}.dedup.bam > ${bam.baseName}.dedup.bam.flagstat
    samtools idxstats ${bam.baseName}.dedup.bam > ${bam.baseName}.dedup.bam.idxstat
    """

    stub:
    """
    touch ${meta.id}.TagAlign.bed ${bam.baseName}.dedup.bam ${bam.baseName}.dedup.bam.flagstat ${bam.baseName}.dedup.bam.idxstat
    """
}

process NGSQC_GEN {
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(bed), path(chrom_sizes)

    output:
        path("NGSQC_report.txt"), emit: report

    script:
    """
    ngsqc -v ${bed} ${chrom_sizes}
    """

    stub:
    """
    touch NGSQC_report.txt
    """
}

process DEEPTOOLS_BAMCOV {
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(bam), path(bai)
        path(ppqt)

    output:
        val(meta.id), emit: meta_id
        path("${meta.id}.bw"), emit: bigwig

    script: // https://deeptools.readthedocs.io/en/2.1.0/content/tools/bamCoverage.html
    """
    frag_len=\$(cut -f 2 ${ppqt})
    bamCoverage \
      --bam ${bam} \
      -o ${meta.id}.bw \
      --binSize ${params.deeptools.bin_size} \
      --smoothLength ${params.deeptools.smooth_length} \
      --ignoreForNormalization ${params.deeptools.excluded_chroms} \
      --numberOfProcessors ${task.cpus} \
      --normalizeUsing ${params.deeptools.normalize_using} \
      --effectiveGenomeSize ${params.align.effective_genome_size} \
      --extendReads \$frag_len
    """

    stub:
    """
    touch ${meta.id}.bw
    """

}
process DEEPTOOLS_BIGWIG_SUM {
    label 'qc'

    input:
        val(meta_ids)
        path(bigwigs)

    output:
        path("bigWigSum.npz")

    script:
    """
    multiBigwigSummary bins \
      -b ${bigwigs} \
      --labels ${meta_ids} \
      -o bigWigSum.npz
    """

    stub:
    """
    echo "${bigwigs}" > bigWigSum.npz
    """
}

/*
process PLOT_NGSQC {
    // TODO refactor bin/ngsqc_plot.py for simplicity

}


process MULTIQC {
    tag { meta.id }
}

 // TODO -- come back to this later. hard to deal with on biowulf and long-running. have to copy entire db to lscratch on biowulf.
process KRAKEN_SE {
    tag { meta.id }
}
*/
