include { BWA_INDEX as BWA_INDEX_BL
          BWA_INDEX as BWA_INDEX_REF } from "../../modules/CCBR/bwa/index"
include { KHMER_UNIQUEKMERS          } from '../../modules/nf-core/khmer/uniquekmers'
include { SPLIT_REF_CHROMS           } from '../../modules/local/prepare_genome.nf'
include { GTF2BED                    } from '../../modules/local/prepare_genome.nf'

workflow PREPARE_GENOME {
    main:
        if (params.genome_fasta && params.genes_gtf) {
            ch_fasta = file(params.genome_fasta)
            ch_blacklist_input = file(params.blacklist)
            if (ch_blacklist_input.endsWith('.fa') || ch_blacklist_input.endsWith('.fna') || ch_blacklist_input.endsWith('.fasta')) {
                ch_blacklist_fasta = ch_blacklist_input
            } else { // TODO optionally convert bed to fasta
                error "The blacklist file must be in fasta nucleotide format. \n\tBlacklist: ${params.blacklist}"
            }
            ch_blacklist_index = BWA_INDEX_BL(ch_blacklist_fasta).index
            ch_reference_index = BWA_INDEX_REF(ch_fasta).index
            KHMER_UNIQUEKMERS(ch_fasta, params.read_length)
            ch_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
            ch_gene_bed = GTF2BED ( file(params.gtf) ).bed
            //ch_chrom_sizes = CUSTOM_GETCHROMSIZES ( ch_fasta ).sizes
            SPLIT_REF_CHROMS(ch_fasta)
            ch_chrom_sizes = SPLIT_REF_CHROMS.out.chrom_sizes
            ch_chrom_files = SPLIT_REF_CHROMS.out.chrom_files
            ch_bwa_index = BWA_INDEX ( ch_fasta ).index
            ch_bowtie2_index = BOWTIE2_BUILD ( ch_fasta ).index

        } else {
            Channel.fromPath(params.genomes[ params.genome ].blacklist_files, checkIfExists: true)
                .collect()
                .set{ ch_blacklist_index }
            Channel.fromPath(params.genomes[ params.genome ].reference_files, checkIfExists: true)
                .collect()
                .set{ ch_reference_index }
            Channel.fromPath(params.genomes[ params.genome ].chrom_sizes, checkIfExists: true)
                .set{ ch_chrom_sizes }
            Channel.fromPath("${params.genomes[ params.genome ].chromosomes_dir}", type: 'dir', checkIfExists: true)
                .set{ ch_chrom_files }
            Channel.fromPath(params.genomes[ params.genome ].gene_info,
                         checkIfExists: true)
                .set{ ch_gene_info }
            Channel.value(params.genomes[ params.genome ].effective_genome_size)
                .set{ ch_gsize}
        }

    emit:
        blacklist_index = ch_blacklist_index
        reference_index = ch_reference_index
        chrom_sizes = ch_chrom_sizes
        chrom_files = ch_chrom_files
        gene_info = ch_gene_info
        effective_genome_size = ch_gsize
}
