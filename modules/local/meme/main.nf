process MEME_AME {
    tag "${meta.id}.${meta.group}"
    label 'peaks'
    label 'process_medium'

    container 'nciccbr:ccbr_meme_5.5.4:v1'

    input:
        tuple val(meta), path(homer_outdir)
        path(motifs)

    output:
        tuple val(meta), path("${meta.id}.${meta.group}_ame/*")

    script:
    """
    ame \\
        --o ${meta.id}.${meta.group}_ame/
        --noseq \\
        --control ${homer_outdir}/background.fa \\
        --seed 20231103 \\
        --verbose 1 \\
        ${homer_outdir}/target.fa
        ${motifs}
    """

    stub:
    """
    mkdir ${meta.id}.${meta.group}_ame/
    touch ${meta.id}.${meta.group}_ame/${meta.id}.${meta.group}.blank.txt
    """
}
