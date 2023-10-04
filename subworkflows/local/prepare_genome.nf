include { BWA_INDEX as BWA_INDEX_BL
          BWA_INDEX as BWA_INDEX_REF } from "../../modules/CCBR/bwa/index"
include { KHMER_UNIQUEKMERS          } from '../../modules/nf-core/khmer/uniquekmers'
include { BEDTOOLS_GETFASTA          } from '../../modules/nf-core/bedtools/getfasta/main'
include { SPLIT_REF_CHROMS           } from '../../modules/local/prepare_genome.nf'
include { GTF2BED                    } from '../../modules/local/prepare_genome.nf'
include { WRITE_GENOME_CONFIG       } from "../../modules/local/prepare_genome.nf"

workflow PREPARE_GENOME {
    main:
        if (params.genome_fasta && params.genes_gtf) {
            ch_fasta = file(params.genome_fasta)
            ch_blacklist_input = file(params.blacklist)
            if (ch_blacklist_input.endsWith('.bed')) {
                BEDTOOLS_GETFASTA(ch_blacklist_input, ch_fasta)
                ch_blacklist_fasta = BEDTOOLS_GETFASTA.out.fasta
            } else {
                ch_blacklist_fasta = ch_blacklist_input
            }
            ch_blacklist_index = BWA_INDEX_BL([['id': 'blacklist'], ch_blacklist_fasta]).index
            ch_reference_index = BWA_INDEX_REF([['id': 'genome'], ch_fasta]).index
            KHMER_UNIQUEKMERS(ch_fasta, params.read_length)
            ch_gsize = KHMER_UNIQUEKMERS.out.kmers.map { it.text.trim() }
            ch_gene_info = GTF2BED ( file(params.genes_gtf) ).bed
            SPLIT_REF_CHROMS(ch_fasta)
            ch_chrom_sizes = SPLIT_REF_CHROMS.out.chrom_sizes
            ch_chrom_dir = SPLIT_REF_CHROMS.out.chrom_dir

            // save reference files
            if (params.save_reference) {
                //file("${params.outdir}/genome/").mkdirs()
                WRITE_GENOME_CONFIG(
                    ch_reference_index,
                    ch_blacklist_index,
                    ch_chrom_sizes,
                    ch_chrom_dir,
                    ch_gene_info,
                    ch_gsize
                )

            }

        } else {
            Channel.fromPath(params.genomes[ params.genome ].blacklist_index, checkIfExists: true)
                .collect()
                .set{ ch_blacklist_index }
            Channel.fromPath(params.genomes[ params.genome ].reference_index, checkIfExists: true)
                .collect()
                .set{ ch_reference_index }
            Channel.fromPath(params.genomes[ params.genome ].chrom_sizes, checkIfExists: true)
                .set{ ch_chrom_sizes }
            Channel.fromPath("${params.genomes[ params.genome ].chromosomes_dir}", type: 'dir', checkIfExists: true)
                .set{ ch_chrom_dir }
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
        chrom_dir = ch_chrom_dir
        gene_info = ch_gene_info
        effective_genome_size = ch_gsize
}
