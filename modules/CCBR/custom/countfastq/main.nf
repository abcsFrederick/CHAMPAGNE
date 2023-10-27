
process CUSTOM_COUNTFASTQ {
    tag { meta.id }
    label 'process_single'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v5'

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), env(count), emit: count

    script:
    """
    count=`zcat ${fastq} | grep "^@" | wc -l`
    echo \$count
    """

    stub:
    """
    count=-1
    echo \$count
    """
}
