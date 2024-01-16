process SORT_BED {
    tag { meta.id }
    label 'process_single'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v5'

    input:
    tuple val(meta), path(bed)

    output:
    tuple val(meta), path("*.sorted.bed"), emit: bed
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ' -k 1,1 -k2,2n '
    def outfile = "${bed.baseName}.sorted.bed"
    """
    sort \\
        ${args} \\
        ${bed} > ${outfile}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(echo \$(sort --version | head -n 1 | sed 's/.* //'))
    END_VERSIONS
    """

    stub:
    def outfile = "${bed.baseName}.sorted.bed"
    """
    touch ${outfile}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sort: \$(echo \$(sort --version | head -n 1 | sed 's/.* //'))
    END_VERSIONS
    """
}
