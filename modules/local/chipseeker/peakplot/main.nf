process CHIPSEEKER_PEAKPLOT {
    tag "${meta.id}.${group}"
    label 'peaks'
    label 'process_medium'

    container 'nciccbr/ccbr_chipseeker:1.1.2'

    input:
        tuple val(meta), path(bed), val(group)

    output:
        tuple val(meta), path("*.png"), emit: plots

    script:
    """
    chipseeker_peakplot.R \\
        --peak ${bed} \\
        --outfile-prefix ${meta.id}.${group} \\
        --genome ${params.genome}
    """

    stub:
    """
    for ftype in annotated.txt summary.txt genelist.txt annotation.Rds .png
    do
        touch ${meta.id}.${group}.\${ftype}
    done
    """
}
