
process CUSTOM_COUNTFASTQ {
    tag { meta.id }
    label 'process_single'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v6.1'

    input:
        tuple val(meta), path(fastq)

    output:
        tuple val(meta), path("*.txt"), emit: count
        path('versions.yml'),           emit: versions

    when:
        task.ext.when == null || task.ext.when

    script:
    template 'count-fastq.py'

    stub:
    """
    count=-1
    echo \$count > ${meta.id}.count.txt
    touch versions.yml
    """
}
