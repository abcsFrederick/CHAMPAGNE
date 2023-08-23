
process FASTQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
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
    touch ${sample_id}_fastqc.html ${sample_id}.cutadapt_fastqc.html
    """
}

process FASTQ_SCREEN {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
        path(conf)
    output:
        path("${sample_id}_screen.*")

    script:
    """
    fastq_screen -c $conf --threads $task.cpus $fastq
    """
}
/* // TODO -- hard to deal with on biowulf. come back to this later. have to copy entire db to lscratch on biowulf.
process KRAKEN_SE {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

}
*/
/*
process BWA_SE {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

}

process PRESEQ {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

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
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

}

process NGSQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

    stub:
    """
    touch NGSQC_report.txt
    """

}

process MULTIQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

}
*/
