process FASTQC {
    tag { sample_id }
    publishDir "$params.outdir/$sample_id/qc", mode: 'symlink'
    container 'nciccbr/ccrgb_qctools:latest'

    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("fastqc")

    script:
    """
    fastqc \
        $fastq \
        -t $task.cpus
        -o fastqc
    """
}
