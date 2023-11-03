process HOMER_MOTIFS {
    tag "${sample}_${tool}"
    label 'peaks'
    label 'process_medium'

    container 'nciccbr/ccbr_homer_4.11:v1'

    input:
        tuple val(sample), val(tool), path(bed), path(genome_fasta)
        val(de_novo) // true or false

    output:
        tuple val(sample), val(tool), path("${sample}.${tool}/*")

    script:
    def args = de_novo ? "" : " -nomotif "
    """
    findMotifsGenome.pl \\
        ${bed} \\
        ${genome_fasta} \\
        ${sample}.${tool}/ \\
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
