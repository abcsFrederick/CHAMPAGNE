process HOMER_MOTIFS {
    tag { sample_tool }
    label 'peaks'
    label 'process_medium'

    container 'nciccbr/ccbr_homer_4.11:v1'

    input:
        tuple val(sample_tool), val(sample_basename), path(bed)
        val(de_novo) // true or false

    output:
        tuple val(sample_tool), path("${sample_tool}/*")

    script:
    def args = de_novo ? "" : " -nomotif "
    """
    findMotifsGenome.pl \\
        ${bed} \\
        ${params.genome} \\
        ${sample_tool} \\
        -size given \\
        -p ${task.cpus} \\
        -len 8,10 \\
        -dumpFasta \\
        -preparsedDir preparsed \\
        ${args}
    """

    stub:
    """
    mkdir ${sample_tool}/
    touch ${sample_tool}/${sample_tool}.txt
    """
}
