
include { CALC_GENOME_FRAC
          MACS_BROAD
          MACS_NARROW
          SICER
          CONVERT_SICER
          GEM
          FRACTION_IN_PEAKS
          CONCAT_FRIPS
          PLOT_FRIP
          JACCARD_INDEX
          CONCAT_JACCARD
          PLOT_JACCARD      } from "../../modules/local/peaks.nf"


workflow CALL_PEAKS {
    take:
        chrom_sizes
        deduped_tagalign
        deduped_bam
        frag_lengths

    main:
        // peak calling
        genome_frac = CALC_GENOME_FRAC(chrom_sizes)
        // create channel with [ meta, chip_tag, input_tag, fraglen, genome_frac]
        deduped_tagalign
            .combine(deduped_tagalign)
            .map {
                meta1, tag1, meta2, tag2 ->
                    meta1.control == meta2.id ? [ meta1, tag1, tag2 ]: null
            }
            .join(frag_lengths)
            .combine(genome_frac)
            .set { ch_tagalign }

        deduped_tagalign
            .combine(deduped_tagalign)
            .map {
                meta1, tag1, meta2, tag2 ->
                    meta1.control == meta2.id ? [ meta1, tag1, tag2 ]: null
            }
            .combine(Channel.fromPath(params.gem_read_dists, checkIfExists: true))
            .combine(chrom_sizes)
            .combine(Channel.fromPath("${params.genomes[ params.genome ].chromosomes_dir}", type: 'dir', checkIfExists: true))
            .set { ch_gem }

        ch_tagalign | MACS_BROAD
        ch_tagalign | MACS_NARROW
        ch_tagalign | SICER | CONVERT_SICER
        GEM(ch_gem)

        CONVERT_SICER.out.peak
            .mix(GEM.out.peak,
                 MACS_BROAD.out.peak,
                 MACS_NARROW.out.peak
                ).set{ ch_peaks }


        // Create Channel with meta, deduped bam, peak file, peak-calling tool, and chrom sizes fasta
        deduped_bam.cross(ch_peaks)
            .map{ it ->
               it.flatten()
            }
            .map{  meta1, bam, bai, meta2, peak, tool ->
                meta1 == meta2 ? [ meta1, bam, bai, peak, tool ] : null
            }
            .combine(chrom_sizes)
            .set{ ch_bam_peaks }
        ch_bam_peaks | FRACTION_IN_PEAKS
        FRACTION_IN_PEAKS.out.collect() | CONCAT_FRIPS | PLOT_FRIP

        ch_peaks
            .combine(ch_peaks) // jaccard index on all-vs-all samples & peak-calling tools
            .combine(chrom_sizes)
            .set{ pairwise_peaks }
        pairwise_peaks | JACCARD_INDEX
        JACCARD_INDEX.out.collect() | CONCAT_JACCARD | PLOT_JACCARD

    emit:
        ch_bam_peaks
}
