process BAM_TO_BED {
    tag { meta.id }

    container "${params.containers.base}"

    input:
        tuple val(meta), path(bam)

    output:
        tuple val(meta), path(bed), emit: bed

    script:
    """
    bedtools bamtobed -i ${bam} > ${bam.baseName}.bed
    """
    stub:
    """
    touch ${bam.baseName}.bed
    """
}
