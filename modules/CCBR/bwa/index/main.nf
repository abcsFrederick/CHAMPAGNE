process BWA_INDEX {
    tag { fasta }
    label 'process_single'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v5'

    input:
        tuple val(meta), path(fasta)

    output:
        tuple val(meta), path(bwa) , emit: index
        path "versions.yml"        , emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    def args   = task.ext.args ?: ''
    """
    mkdir bwa
    bwa \\
        index \\
        ${args} \\
        -p bwa/${prefix} \\
        ${fasta}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${fasta.baseName}"
    """
    mkdir bwa
    for ext in amb ann bwt pac sa; do
        touch bwa/${prefix}.\$ext
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
    END_VERSIONS
    """
}
