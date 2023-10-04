process GTF2BED {
    tag { gtf }
    label 'process_single'

    container 'nciccbr/ccbr_ubuntu_base_20.04:v5'

    input:
        path(gtf)

    output:
        path("${gtf.baseName}.bed"), emit: bed

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

process CREATE_GENOME_CONFIG {
    label 'process_single'
    container "${params.containers.base}"

    input:
        tuple val(meta_ref), path(reference_index)
        tuple val(meta_bl), path(blacklist_index)
        path(chrom_sizes)
        path(chrom_files)
        path(gene_info)
        val(effective_genome_size)

    output:
        path('genome/')

    script:
    def genome_name = 'custom_genome' // TODO allow user to set this
    """
    #!/usr/bin/env python
    import os
    import pprint

    os.mkdir('genome/')
    genome = dict(reference_index = "${reference_index}",
                  blacklist_index = "${blacklist_index}",
                  effective_genome_size = "${effective_genome_size}",
                  chrom_sizes = "${chrom_sizes}",
                  gene_info = "${gene_info}",
                  chromosomes_dir = "${chrom_files}"
    )
    print(genome)
    with open('genome/${genome_name}.config', 'w') as conf_file:
        conf_file.write("'${genome_name}' {")
        for k, v in genome:
            conf_file.write(f"\t{k} = \${{params.index_dir}}/{v}")
        conf_file.write("}")
    """
}
