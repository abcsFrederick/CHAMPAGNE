process CAT_CAT {
    tag { meta.id }
    label 'process_low'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v5'

    input:
    tuple val(meta), path(files_in)

    output:
    tuple val(meta), path("${prefix}"), emit: file_out
    path "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def file_list = files_in.sort({ a, b -> a.baseName <=> b.baseName }).collect{ it.toString() }

    // | input     | output     | command1 | command2 |
    // |-----------|------------|----------|----------|
    // | gzipped   | gzipped    | cat      |          |
    // | ungzipped | ungzipped  | cat      |          |
    // | gzipped   | ungzipped  | zcat     |          |
    // | ungzipped | gzipped    | cat      | pigz     |

    // Use input file ending as default
    prefix   = task.ext.prefix ?: "${meta.id}${file_list[0].substring(file_list[0].lastIndexOf('.'))}"
    out_zip  = prefix.endsWith('.gz')
    in_zip   = file_list[0].endsWith('.gz')
    command1 = (in_zip && !out_zip) ? 'zcat' : 'cat'
    command2 = (!in_zip && out_zip) ? "| pigz -c -p $task.cpus $args2" : ''
    """
    $command1 \\
        $args \\
        ${file_list.join(' ')} \\
        $command2 \\
        > ${prefix}

    echo "Concatenated file order: ${file_list.join(' ')}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """

    stub:
    def file_list = files_in.collect { it.toString() }
    prefix   = task.ext.prefix ?: "${meta.id}${file_list[0].substring(file_list[0].lastIndexOf('.'))}"
    """
    touch $prefix

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pigz: \$( pigz --version 2>&1 | sed 's/pigz //g' )
    END_VERSIONS
    """
}
