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
        path("chroms/")                      , emit: chrom_dir

    script:
    """
    splitRef.py ${fasta} ${fasta.baseName}.chrom.sizes chroms
    """
    stub:
    """
    touch ${fasta.baseName}.chrom.sizes
    mkdir -p chroms/
    touch chroms/chr1.fa
    """
}

process WRITE_GENOME_CONFIG {
    label 'process_single'
    container "${params.containers.base}"

    input:
        tuple val(meta_ref), path(reference_index)
        tuple val(meta_bl), path(blacklist_index)
        path(chrom_sizes)
        path(chrom_dir)
        path(gene_info)
        val(effective_genome_size)

    output:
        path("*.config"), emit: conf
        path("custom_genome/*"), emit: files

    script:
    def genome_name = 'custom_genome'
    """
    #!/usr/bin/env python
    import os
    import pprint
    import shutil
    os.makedirs("${genome_name}/")
    for subdir, filelist in (('reference/', "${reference_index}"), ('blacklist', "${blacklist_index}")):
        dirpath = f"${{genome_name}}/{subdir}"
        os.mkdir(dirpath)
        for file in filelist.split():
            shutil.copy(file, dirpath)
    shutil.copytree("${chrom_dir}", '${genome_name}/chroms/')
    for file in ("${chrom_sizes}", "${gene_info}"):
        shutil.copy(file, "${genome_name}/")

    genome = dict(reference_index = "\${params.index_dir}/${genome_name}/reference/*",
                    blacklist_index = "\${params.index_dir}/${genome_name}/blacklist/*",
                    effective_genome_size = "${effective_genome_size}",
                    chrom_sizes = "\${params.index_dir}/${genome_name}/${chrom_sizes}",
                    gene_info = "\${params.index_dir}/${genome_name}/${gene_info}",
                    chromosomes_dir = "\${params.index_dir}/${genome_name}/chroms/*"
    )
    pprint.pprint(genome)
    with open('${genome_name}.config', 'w') as conf_file:
        conf_file.write("'${genome_name}' {")
        for k, v in genome.items():
            conf_file.write(f"\\t{k} = {v}\\n")
        conf_file.write("}")
    """
}
