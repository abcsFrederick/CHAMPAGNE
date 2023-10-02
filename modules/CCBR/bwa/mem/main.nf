process BWA_MEM {
    tag { meta.id }
    label 'process_higher'

    container = 'nciccbr/ccbr_ubuntu_base_20.04:v5'

    input:
        tuple val(meta), path(fastq)
        path(index_files)
        val(index_name)

    output:
        tuple val(meta), path("*.bam"), emit: bam
        path  "versions.yml"          , emit: versions

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # current working directory is a tmpdir when 'scratch' is set
    TMP=tmp/
    mkdir \$TMP
    trap 'rm -rf "\$TMP"' EXIT

    bwa mem -t ${task.cpus} ${index_name} ${fastq} -o \$TMP/align.bam

    samtools sort \\
      -@ ${task.cpus} \\
      -m 2G \\
      -T \$TMP \\
      \$TMP/align.bam > ${prefix}.bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
