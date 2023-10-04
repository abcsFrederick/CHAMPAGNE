process GTF2BED {
    tag { gtf }
    label 'process_single'

    container 'quay.io/biocontainers/bedops:2.4.41--h4ac6f70_1'

    input:
        path(gtf)

    output:
        path("${gtf.baseName}.bed")

    script:
    """
    cat ${gtf} | gtf2bed > ${gtf.baseName}.bed
    """
    stub:
    """
    touch ${gtf.baseName}.bed
    """
}
process SPLIT_REF_CHROMS {
    tag { fasta }
    label 'process_single'
    container "${params.containers.base}"

    input:
        path(fasta)

    output:
        path("${fasta.baseName}.chrom.sizes"), emit: chrom_sizes
        path("chroms/*")                     , emit: chrom_files

    script:
    """
    splitRef.py ${fasta} ${fasta.baseName}.chrom.sizes chroms
    """
    stub:
    """
    touch ${fasta.baseName}.chrom.sizes
    """
}
