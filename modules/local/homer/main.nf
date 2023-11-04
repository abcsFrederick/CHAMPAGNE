process HOMER_MOTIFS {
    tag "${sample}.${tool}"
    label 'peaks'
    label 'process_medium'

    container 'nciccbr/ccbr_homer_4.11:v1'

    input:
        tuple val(sample), val(tool), path(bed), path(genome_fasta)
        val(de_novo) // true or false
        path(motif_db)

    output:
        tuple val(sample), val(tool), path("${sample}.${tool}/*"), emit: motifs

    script:
    def args = de_novo ? "" : " -nomotif "
    """
    findMotifsGenome.pl \\
        ${bed} \\
        ${genome_fasta} \\
        ${sample}.${tool}/ \\
        -mknown ${motif_db} \\
        -size given \\
        -p ${task.cpus} \\
        -len 8,10 \\
        -dumpFasta \\
        -preparsedDir preparsed \\
        ${args}
    """

    stub:
    """
    mkdir ${sample}.${tool}/
    touch ${sample}.${tool}/${sample}.${tool}
    """
}
