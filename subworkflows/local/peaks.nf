
include { CALC_GENOME_FRAC
          MACS_BROAD
          MACS_NARROW
          SICER
          CONVERT_SICER
          GEM
          FILTER_GEM
          FRACTION_IN_PEAKS
          CONCAT_FRIPS
          PLOT_FRIP
          JACCARD_INDEX
          CONCAT_JACCARD
          PLOT_JACCARD
          GET_PEAK_META
          CONCAT_PEAK_META
          PLOT_PEAK_WIDTHS } from "../../modules/local/peaks.nf"
include { BAM_TO_BED       } from "../../modules/local/bedtools.nf"
include { CONSENSUS_PEAKS  } from "../../modules/local/consensus_peaks"
include { HOMER_MOTIFS     } from "../../modules/local/homer"
include { MEME_AME         } from "../../modules/local/meme"

workflow CALL_PEAKS {
    take:
        chrom_sizes
        chrom_dir
        deduped_tagalign
        deduped_bam
        frag_lengths
        effective_genome_size
        genome_fasta

    main:
        genome_frac = CALC_GENOME_FRAC(chrom_sizes, effective_genome_size)

        // create channel with [ meta, chip_tag, input_tag, format ]
        deduped_tagalign
            .combine(deduped_tagalign)
            .map {
                meta1, tag1, format1, meta2, tag2, format2 ->
                    meta1.control == meta2.id && format1 == format2 ? [ meta1, tag1, tag2, format1 ]: null
            }
            .set{ ch_tagalign }

        // create macs channel with [ meta, chip_tag, input_tag, format, fraglen, genome_frac]
        ch_tagalign
            .join(frag_lengths)
            .combine(genome_frac)
            .combine(effective_genome_size)
            .set { ch_macs }

        // create sicer channel containing only bed files with [ meta, chip_tag, input_tag, fraglen, genome_frac]
        deduped_tagalign
            .branch{ meta, tag, format ->
                bam: format == 'BAMPE'
                    return(tuple(meta, tag))
                bed: format == 'BED'
                    return(tuple(meta, tag))
            }.set{ tag_split }
        tag_split.bam | BAM_TO_BED
        BAM_TO_BED.out.bed.set{ tag_split_converted }
        tag_split.bed.mix(tag_split_converted).set{ tag_all_bed }
        tag_all_bed
            .combine(tag_all_bed)
            .map {
                meta1, tag1, meta2, tag2 ->
                    meta1.control == meta2.id ? [ meta1, tag1, tag2 ]: null
            }
            .join(frag_lengths)
            .combine(genome_frac)
            .set { ch_sicer }


        // create gem channel with [ meta, chip_tag, input_tag, format, read_dists, chrom_sizes, chrom_dir, effective_genome_sizes ]
        ch_tagalign
            .combine(Channel.fromPath(params.gem.read_dists, checkIfExists: true))
            .combine(chrom_sizes)
            .combine(chrom_dir)
            .combine(effective_genome_size)
            .set { ch_gem }

        ch_peaks = Channel.empty()
        if (params.run.macs_broad) {
            ch_macs | MACS_BROAD
            ch_peaks = ch_peaks.mix(MACS_BROAD.out.peak)
        }
        if (params.run.macs_narrow) {
            ch_macs | MACS_NARROW
            ch_peaks = ch_peaks.mix(MACS_NARROW.out.peak)
        }
        if (params.run.sicer) {
            ch_sicer | SICER | CONVERT_SICER
            ch_peaks = ch_peaks.mix(CONVERT_SICER.out.peak)
        }
        if (params.run.gem) {
            ch_gem | GEM
            GEM.out.peak
                .combine(chrom_sizes) | FILTER_GEM
            ch_peaks = ch_peaks.mix(FILTER_GEM.out.peak)
        }

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
            .map{ meta1, peak1, tool1, meta2, peak2, tool2 ->
                (meta1 != meta2 || tool1 != tool2) ? [ meta1, peak1, tool1, meta2, peak2, tool2 ] : null
            }
            .combine(chrom_sizes) | JACCARD_INDEX
        JACCARD_INDEX.out.collect() | CONCAT_JACCARD | PLOT_JACCARD

        ch_bam_peaks | GET_PEAK_META
        GET_PEAK_META.out.collect() | CONCAT_PEAK_META | PLOT_PEAK_WIDTHS

        ch_plots = PLOT_FRIP.out
            .mix(PLOT_JACCARD.out)
            .mix(PLOT_PEAK_WIDTHS.out)

        // consensus peak calling on replicates
        ch_peaks
            .map{ meta, bed, tool ->
                [ "${meta.sample_basename}_${tool}", meta, bed, tool,]
            }
            .groupTuple(by: 0) // group by the sample_basename and peak-calling tool
            .set{ peak_reps }
        // assert that sample_basenames & tools match
        peak_reps.subscribe { basename_tool, metas, beds, tools ->
            assert metas.collect{ it.sample_basename }.toSet().size() == 1
            assert tools.toSet().size() == 1
        }
        peak_reps
            .map { basename_tool, metas, beds, tools ->
                [ [id: metas[0].sample_basename, group: tools[0]], beds ]
            }
            .set{
                peaks_grouped
            }
        peaks_grouped | CONSENSUS_PEAKS

        HOMER_MOTIFS( CONSENSUS_PEAKS.out.peaks.combine(genome_fasta),
                      params.homer.de_novo,
                      file(params.homer.jaspar_db, checkIfExists: true)
                    )
        if (params.genomes[ params.genome ].meme_motifs) {
            MEME_AME( peaks_groupped,
                      file(params.genomes[ params.genome ].meme_motifs, checkIfExists: true)
                    )
        }

    emit:
        peaks = ch_bam_peaks
        plots = ch_plots
}
