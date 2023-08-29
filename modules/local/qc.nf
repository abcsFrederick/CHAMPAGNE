
process FASTQC {
    tag { sample_id }
    publishDir "$params.outdir/qc/fastqc_${fqtype}/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(fastq), val(fqtype)
    output:
        path("${sample_id}*.html"), emit: html

    script:
    """
    fastqc \
        $fastq \
        -t $task.cpus \
        -o .
    """

    stub:
    """
    touch ${sample_id}_fastqc.html
    """
}

process FASTQ_SCREEN {
    tag { sample_id }
    publishDir "${params.outdir}/qc/fastq_screen/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(fastq), path(conf)
    output:
        path("${sample_id}*_screen.*"), emit: screen

    script:
    """
    fastq_screen -c $conf --threads $task.cpus $fastq
    """

    stub:
    """
    for EXT in html txt png; do
        touch ${sample_id}.trimmed_screen.\${EXT}
    done
    """
}

process ALIGN_BLACKLIST {
    tag { sample_id }
    publishDir "${params.outdir}/qc/align/${sample_id}", mode: "$params.filePublishMode"

    input:
        tuple val(sample_id), path(fastq)
        path(blacklist_files)

    output:
        tuple val(sample_id), path("${sample_id}.no_blacklist.fastq.gz"), emit: reads

    script: // TODO use samtools -f4 for single-end and -f12 for paired to get unmapped reads https://broadinstitute.github.io/picard/explain-flags.html
    """
    bwa mem -t $task.cpus $params.align.blacklist $fastq |\
    samtools view -@ $task.cpus -f4 -b |\
    samtools bam2fq |\
    pigz -p $task.cpus > ${sample_id}.no_blacklist.fastq.gz
    """

    stub:
    """
    touch ${sample_id}.no_blacklist.fastq.gz
    """
}

process ALIGN_GENOME {
    tag { sample_id }
    publishDir "${params.outdir}/qc/align/${sample_id}", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
        path(reference_files)

    output:
        tuple val(sample_id), path("${sample_id}.aligned.filtered.bam"), emit: bam

    script:
    """
    bwa mem -t ${task.cpus} ${params.align.genome} $fastq |\
    samtools sort -@ ${task.cpus} |\
    samtools view -b -q ${params.align.min_quality} -o ${sample_id}.aligned.filtered.bam
    """

    stub:
    """
    touch ${sample_id}.aligned.filtered.bam
    """

}

process PRESEQ {
    """
    Calls preseq c_curve and lc_extrap, and calls bin/parse_preseq_log.py to get statistics from the log.
    """
    tag { sample_id }
    publishDir "${params.outdir}/qc/preseq/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple path("*.c_curve"), path("*.preseq"), path("*.preseqlog"), path("*.txt"), emit: preseq_files

    script:
    """
    preseq c_curve -B -o ${sample_id}.c_curve $bam
    preseq lc_extrap -B -D -o ${sample_id}.preseq $bam -seed 12345 -v -l 100000000000 2> ${sample_id}.preseqlog
    parse_preseq_log.py ${sample_id}.preseqlog > ${sample_id}.preseqlog.nrf.txt
    """

    stub:
    """
    touch ${sample_id}.c_curve ${sample_id}.preseq ${sample_id}.preseqlog ${sample_id}.preseqlog.nrf.txt
    """
}

process PHANTOM_PEAKS { // https://github.com/kundajelab/phantompeakqualtools
    // TODO: set tmpdir as lscratch if available https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/InitialChIPseqQC.snakefile#L504
    tag { sample_id }
    publishDir "${params.outdir}/qc/ppqt/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(bam), path(bai)

    output:
        path("${sample_id}.ppqt.pdf"), emit: pdf
        path("${sample_id}.ppqt"), emit: ppqt

    script: // TODO: for PE, just use first read of each pair
    """
    RUN_SPP=\$(which run_spp.R)
    Rscript \$RUN_SPP -c=${bam} -savp=${sample_id}.ppqt.pdf -out=${sample_id}.ppqt
    """

    stub:
    """
    touch ${sample_id}.ppqt.pdf ${sample_id}.ppqt
    """
}

process DEDUPLICATE {
    tag { sample_id }
    publishDir "${params.outdir}/qc/dedup/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(bam), path(chrom_sizes)

    output:
        tuple val(sample_id), path("${sample_id}.TagAlign.bed"), emit: tag_align
        tuple val(sample_id), path("${bam.baseName}.dedup.bam"), emit: bam
        path("${bam.baseName}.dedup.bam.flagstat"), emit: flagstat
        path("${bam.baseName}.dedup.bam.idxstat"), emit: idxstat

    script:
    """
    macs2 filterdup -i ${bam} -g ${params.align.effective_genome_size} --keep-dup="auto" -o TmpTagAlign1.bed
    awk -F"\\t" -v OFS="\\t" '{{if (\$2>0 && \$3>0) {{print}}}}' TmpTagAlign1.bed > TmpTagAlign2.bed
    awk -F"\\t" -v OFS="\\t" '{{print \$1,1,\$2}}' ${chrom_sizes} | sort -k1,1 -k2,2n > GenomeFile.bed
    bedtools intersect -wa -f 1.0 -a TmpTagAlign2.bed -b GenomeFile.bed > ${sample_id}.TagAlign.bed
    bedtools bedtobam -i ${sample_id}.TagAlign.bed -g ${chrom_sizes} | samtools sort -@ ${task.cpus} -o ${bam.baseName}.dedup.bam
    samtools index ${bam.baseName}.dedup.bam
    samtools flagstat ${bam.baseName}.dedup.bam > ${bam.baseName}.dedup.bam.flagstat
    samtools idxstats ${bam.baseName}.dedup.bam > ${bam.baseName}.dedup.bam.idxstat
    """

    stub:
    """
    touch ${sample_id}.TagAlign.bed ${bam.baseName}.dedup.bam ${bam.baseName}.dedup.bam.flagstat ${bam.baseName}.dedup.bam.idxstat
    """
}

process INDEX_BAM {
    tag { sample_id }
    publishDir "${params.outdir}/qc/align/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(bam)

    output:
        tuple val(sample_id), path("*.bam"), path("*.bai"), emit: bam
        path("*.idxstats"), emit: idxstats

    script:
    """
    samtools sort -@ ${task.cpus} -o ${sample_id}.sorted.bam $bam
    samtools index ${sample_id}.sorted.bam   # creates ${sample_id}.sorted.bam.bai
    samtools idxstats ${sample_id}.sorted.bam > ${sample_id}.sorted.bam.idxstats
    """

    stub:
    """
    touch ${sample_id}.sorted.bam ${sample_id}.sorted.bam.bai ${sample_id}.sorted.bam.idxstats
    """
}

process NGSQC_GEN {
    tag { sample_id }
    publishDir "${params.outdir}/qc/ngsqc/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(bed), path(chrom_sizes)

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
    tag { sample_id }
    publishDir "${params.outdir}/qc/deeptools_bigwigs/${sample_id}", mode: "${params.filePublishMode}"

    input:
        tuple val(sample_id), path(bam), path(bai)
        path(ppqt)

    output:
        val(sample_id), emit: sample_id
        path("${sample_id}.bw"), emit: bigwig

    script: // https://deeptools.readthedocs.io/en/2.1.0/content/tools/bamCoverage.html
    """
    frag_len=\$(cut -f 2 ${ppqt})
    bamCoverage \
      --bam ${bam} \
      -o ${sample_id}.bw \
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
    touch ${sample_id}.bw
    """

}
process DEEPTOOLS_BIGWIG_SUM {
    tag { sample_id }
    publishDir "${params.outdir}/qc/deeptools_qc", mode: "${params.filePublishMode}"

    input:
        path(bigwigs)

    output:
        path("bigWigSum.npz")

    script:
    """
    multiBigwigSummary bins \
      -b ${bigwigs}
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
    tag { sample_id }
}

 // TODO -- come back to this later. hard to deal with on biowulf and long-running. have to copy entire db to lscratch on biowulf.
process KRAKEN_SE {
    tag { sample_id }
}
*/
