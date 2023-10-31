
process CUSTOM_COUNTFASTQ {
    tag { meta.id }
    label 'process_single'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v6.1'

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("*.txt"), emit: count

    script:
    def txt_filename = "${meta.baseName}.txt"
    template 'count-fastq.py'

    stub:
    """
    count=-1
    echo \$count
    """
}
