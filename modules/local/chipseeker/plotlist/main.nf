process CHIPSEEKER_PLOTLIST {
    tag "${meta.id}.${meta.group}"
    label 'peaks'
    label 'process_medium'

    container 'nciccbr/ccbr_chipseeker:v1.0.1'

    input:
        path(rds)

    output:
        path("*.png"), emit: plots

    script:
    """
    chipseeker_plotlist.R \\
        --annotations ${rds.join(' ')} \\
        --outfile plot_anno_bar.png
    """

    stub:
    """
    touch plot_anno_bar.png
    """
}
