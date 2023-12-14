process MANORM_PAIRWISE {
    """
    adapted from: https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/ChIPseq.snakefile#L709-L785
    """
    tag "${meta1.id}.${meta2.id}"
    label 'process_single'

    container ''

    input:
        tuple val(meta1), val(meta2), path(tagalign1), path(tagalign2), path(peak1), path(peak2)

    output:
        path("manorm_output/*"), emit: dir

    script:
    """
    manorm \\
        --p1 ${peak1} \\
        --p2 ${peak2} \\
        --r1 ${tagalign2} \\
        --r2 ${tagalign2} \\
        --s1 ${meta1.fraglen} \\
        --s2 ${meta2.fraglen} \\
        -o manorm_output/ \\
        --name1 ${meta1.group} \\
        --name2 ${meta2.group}
    """

    stub:
    """
    mkdir manorm_output/
    touch manorm_output/blank.txt
    """
}
