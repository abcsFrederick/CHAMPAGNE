workflow PREPARE_GENOME {
    main:
        ch_blacklist_name = params.genomes[ params.genome ].blacklist ? Channel.value(params.genomes[ params.genome ].blacklist) : "${params.genome}.blacklist"

        if (params.fasta && params.gtf) {
            ch_fasta = file(params.fasta)
            ch_gtf = file(params.gtf)
        } else {
            Channel.fromPath(params.genomes[ params.genome ].blacklist_files, checkIfExists: true)
                .collect()
                .set{ ch_blacklist_files }
            Channel.fromPath(params.genomes[ params.genome ].reference_files, checkIfExists: true)
                .collect()
                .set{ ch_reference_files }
            Channel.fromPath(params.genomes[ params.genome ].chrom_sizes, checkIfExists: true)
                .set{ ch_chrom_sizes }
            Channel.fromPath("${params.genomes[ params.genome ].chromosomes_dir}", type: 'dir', checkIfExists: true)
                .set{ ch_chrom_files }
            Channel.fromPath(params.genomes[ params.genome ].gene_info,
                         checkIfExists: true)
                .set{ ch_gene_info }
            Channel.value(params.genomes[ params.genome ].effective_genome_size)
                .set{ ch_egs}
        }

    emit:
        blacklist_files = ch_blacklist_files
        blacklist_name = ch_blacklist_name
        reference_files = ch_reference_files
        chrom_sizes = ch_chrom_sizes
        chrom_files = ch_chrom_files
        gene_info = ch_gene_info
        effective_genome_size = ch_egs
}
