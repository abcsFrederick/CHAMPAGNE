
process FASTQC {
    tag { sample_id }
    publishDir "$params.outdir/qc/$sample_id/fastqc_$fqtype", mode: "$params.filePublishMode"

    input:
        tuple val(sample_id), path(fastq), val(fqtype)
    output:
        path("${sample_id}*.html")

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
    publishDir "$params.outdir/qc/$sample_id/fastq_screen", mode: "$params.filePublishMode"

    input:
        tuple val(sample_id), path(fastq), path(conf)
    output:
        path("${sample_id}*_screen.*")

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
    publishDir "$params.outdir/qc/$sample_id/align/", mode: "$params.filePublishMode"

    input:
        tuple val(sample_id), path(fastq)
        path(blacklist_files)

    output:
        tuple val(sample_id), path("${sample_id}.no_blacklist.fastq.gz")

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
    publishDir "${params.outdir}/qc/${sample_id}/align/", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
        path(reference_files)

    output:
        path("${sample_id}.aligned.filtered.bam")

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

/* // TODO -- hard to deal with on biowulf. come back to this later. have to copy entire db to lscratch on biowulf.
process KRAKEN_SE {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: "$params.filePublishMode"

}
*/
/*

process PRESEQ {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: "$params.filePublishMode"

    input:
        path(bam)

    script:
    """
    preseq c_curve -B -o ${sample_id}.ccurve $bam
    """
    stub:
    """
    touch ${sample_id}.ccurve
    """
}

process MACS2 {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: "$params.filePublishMode"

}

process NGSQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: "$params.filePublishMode"

    stub:
    """
    touch NGSQC_report.txt
    """

}

process MULTIQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: "$params.filePublishMode"

}
*/
