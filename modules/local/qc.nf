
process FASTQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc/", mode: 'copy'

    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("${sample_id}*.html")

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
