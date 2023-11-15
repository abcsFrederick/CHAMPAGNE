process CHIPSEEKER_ANNOTATE {
    tag "${meta.id}.${meta.group}"
    label 'peaks'
    label 'process_medium'

    container 'nciccbr/ccbr_atacseq:v1.10'

    input:
        tuple val(meta), path(bed)

    output:
        tuple val(meta), path("*.annotated.txt"), path("*.summary.txt"), path("*.genelist.txt"), emit: annot

    script:
    """
    annotate_peaks.R \\
        --narrowpeak ${bed} \\
        --annotated ${meta.id}.${meta.group}.annotated.txt \\
        --atypefreq ${meta.id}.${meta.group}.summary.txt \\
        --genelist ${meta.id}.${meta.group}.genelist.txt \\
        --genome ${params.genome} \\
        --uptss 2000 \\
        --downtss 2000 \\
        --toppromoterpeaks 1000
    """

    stub:
    """
    for ftype in annotated summary genelist
    do
        touch ${meta.id}.${meta.group}.\${ftype}.txt
    done
    """
}
