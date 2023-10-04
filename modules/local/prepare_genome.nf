process EMIT_META {
    input:
        path(file)
    output:
        tuple val(file.baseName), path(file)
    script:
    """
    echo $file
    """
}

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

    publishDir = [
        path: { "${params.outdir}/genome" },
        mode: "copy"
    ]
    input:
        tuple val(meta_ref), path(reference_index)
        tuple val(meta_bl), path(blacklist_index)
        path(chrom_sizes)
        path(chrom_dir)
        path(gene_info)
        val(effective_genome_size)

    output:
        path("*.config"), emit: conf
        path("custom_genome/"), emit: files // TODO can't use genome_name variable here, nextflow thinks it's null??

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

    genome = dict(reference_index = '"\${params.index_dir}/${genome_name}/reference/*"',
                  blacklist_index = '"\${params.index_dir}/${genome_name}/blacklist/*"',
                  chromosomes_dir = '"\${params.index_dir}/${genome_name}/chroms/"',
                  chrom_sizes = '"\${params.index_dir}/${genome_name}/${chrom_sizes}"',
                  gene_info = '"\${params.index_dir}/${genome_name}/${gene_info}"',
                  effective_genome_size = "${effective_genome_size}"
    )
    pprint.pprint(genome)
    with open('${genome_name}.config', 'w') as conf_file:
        head = ["params {\\n",
                '\\tindex_dir = "${params.outdir}/genome/"\\n',
                "\\tgenomes {\\n"
                "\\t\\t'${genome_name}' {\\n"]
        conf_file.writelines(head)
        for k, v in genome.items():
            conf_file.write(f'\\t\\t\\t{k} = {v}\\n')
        tail = ["\\t\\t}\\n",
                "\\t}\\n",
                "}\\n"]
        conf_file.writelines(tail)
    """
}
