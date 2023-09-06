
process FASTQC {
    tag { meta.id }
    label 'qc'
    publishDir "${params.outdir}/qc/fastqc_${fqtype}/${meta.id}", mode: "${params.publish_dir_mode}"

    input:
        tuple val(meta), path(fastq), val(fqtype)
    output:
        path("${fastq.getBaseName(2)}*.html"), emit: html
        path ("${fastq.getBaseName(2)}*.zip"), emit: zip

    script:
    """
    fastqc \
      $fastq \
      -t $task.cpus \
      -o .
    """

    stub:
    """
    touch ${fastq.getBaseName(2)}_fastqc.html ${fastq.getBaseName(2)}_fastqc.zip
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
    Calls preseq c_curve and lc_extrap
    """
    tag { meta.id }
    label 'qc'

    // preseq is known to fail inexplicably, especially on small datasets.
    // https://github.com/nf-core/methylseq/issues/161
    errorStrategy 'ignore'

    input:
        tuple val(meta), path(bam)

    output:
        path("*.c_curve"), emit: c_curve
        path("*.lc_extrap.txt"), emit: preseq
        tuple val(meta), path("*.preseq.log"), emit: log

    script:
    // TODO handle paired: https://github.com/nf-core/rnaseq/blob/3bec2331cac2b5ff88a1dc71a21fab6529b57a0f/modules/nf-core/preseq/lcextrap/main.nf#L25
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    preseq c_curve -B -o ${prefix}.c_curve ${bam}
    preseq lc_extrap -B -D -o ${prefix}.lc_extrap.txt ${bam} -seed 12345 -v -l 100000000000 2> ${prefix}.preseq.log
    """

    stub:
    """
    touch ${meta.id}.c_curve ${meta.id}.preseq ${meta.id}.preseqlog ${meta.id}.preseqlog.nrf.txt
    """
}

process PARSE_PRESEQ_LOG {
    """
    Calls bin/parse_preseq_log.py to get NRF statistics from the preseq log.
    """
    input:
        tuple val(meta), path(log)
    output:
        path("*nrf.txt"), emit: nrf
    script:
    """
    parse_preseq_log.py ${prefix}.preseq.log > ${prefix}.preseq.nrf.txt
    """

}

process PHANTOM_PEAKS {
    // TODO: set tmpdir as lscratch if available https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/InitialChIPseqQC.snakefile#L504
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(bam), path(bai)

    output:
        path("${meta.id}.ppqt.pdf"), emit: pdf
        path("${meta.id}.spp.out"), emit: spp
        tuple val(meta), path("${meta.id}.fraglen.txt"), emit: fraglen

    script: // TODO: for PE, just use first read of each pair
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    RUN_SPP=\$(which run_spp.R)
    Rscript \$RUN_SPP -c=${bam} -savp=${prefix}.ppqt.pdf -out=${prefix}.spp.out
    frag_len=`cut -f 3 ${prefix}.spp.out | sed 's/,.*//g'`
    echo \$frag_len > "${meta.id}.fraglen.txt"
    """

    stub:
    """
    touch ${meta.id}.ppqt.pdf ${meta.id}.spp.out
    """
}

process PPQT_PROCESS { // refactor of https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/InitialChIPseqQC.snakefile#L513-L541

    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(fraglen)
    output:
        path("${fraglen.baseName}.process.txt"), emit: fraglen

    script:
    """
    #!/usr/bin/env python

    import warnings
    with open("${fraglen}", 'r') as infile:
        fragment_length = int(infile.read().strip())
    min_frag_len = ${params.min_fragment_length}
    if fragment_length < min_frag_len:
        warnings.warn(f"The estimated fragment length was {fragment_length}. Using default of {min_frag_len} instead.")
        fragment_length = min_frag_len
    with open("${fraglen.baseName}.process.txt", 'w') as outfile:
        outfile.write(str(fragment_length))
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
        tuple path("${bam.baseName}.dedup.bam.flagstat"), path("${bam.baseName}.dedup.bam.idxstat"), emit: flagstat

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

process NGSQC_GEN { // TODO segfault
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

process QC_STATS {
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(raw_fastq)
        path(align_flagstat)
        tuple path(dedup_flagstat), path(idxstat)
        path(preseq_nrf)
        path(ppqt_spp)
        path(ppqt_fraglen)


    output:
        path("${meta.id}.qc_stats.txt")

    script:
    // TODO: handle paired reads
    def outfile = "${meta.id}.qc_stats.txt"
    """
    touch ${outfile}
    # Number of reads
    zcat ${raw_fastq} | wc -l | filterMetrics.py ${meta.id} tnreads >> ${outfile}
    # Number of mapped reads
    grep 'mapped (' ${align_flagstat} | awk '{{print \$1,\$3}}' | filterMetrics.py ${meta.id} mnreads >> ${outfile}
    # Number of uniquely mapped reads
    grep 'mapped (' ${dedup_flagstat} | awk '{{print \$1,\$3}}' | filterMetrics.py ${meta.id} unreads >> ${outfile}
    # NRF, PCB1, PCB2
    cat ${preseq_nrf} | filterMetrics.py ${meta.id} nrf >> ${outfile}
    # NSC, RSC, Qtag
    awk '{{print \$(NF-2),\$(NF-1),\$NF}}' ${ppqt_spp} | filterMetrics.py ${meta.id} ppqt >> ${outfile}
    # Fragment Length
    fragLen=\$(cat ${ppqt_fraglen})
    echo "${meta.id}\tFragmentLength\t\$fragLen" >> ${outfile}
    """

}

process QC_TABLE {
    label 'qc'
    input:
        path(qc_stats)

    output:
        path("qc_table.txt"), emit: txt

    script:
    """
    cat ${qc_stats.join(' ')} | createtable.py > qc_table.txt
    """

}

process MULTIQC {
    label 'qc'

    input:
        path(multiqc_conf)

        path(fastqc_raw)
        path(fastqc_trimmed)
        path(fastq_screen)
        path(dedup)
        path(phantom_peaks)
        path(qc_table)
        path(plot_fingerprint_matrix)
        path(plot_fingerprint_metrics)
        path(plot_corr)
        path(plot_pca)
        path(plot_profile)

    output:
        path('multiqc_report.html'), emit: html

    script:
    """
    multiqc \\
        -c ${multiqc_conf} \\
        .
    """

    stub:
    """
    touch multiqc_report.html
    """
}

/*
process PLOT_NGSQC {
    // TODO refactor bin/ngsqc_plot.py for simplicity

}

 // TODO -- come back to this later. hard to deal with on biowulf and long-running. have to copy entire db to lscratch on biowulf.
process KRAKEN_SE {
    tag { meta.id }
}
*/
