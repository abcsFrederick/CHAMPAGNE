
include { CALC_GENOME_FRAC  } from "../../modules/local/peaks.nf"
include { MACS_BROAD        } from "../../modules/local/peaks.nf"
include { MACS_NARROW       } from "../../modules/local/peaks.nf"
include { SICER             } from "../../modules/local/peaks.nf"
include { CONVERT_SICER     } from "../../modules/local/peaks.nf"
include { GEM               } from "../../modules/local/peaks.nf"
include { FRACTION_IN_PEAKS } from "../../modules/local/peaks.nf"
include { CONCAT_FRIPS      } from "../../modules/local/peaks.nf"


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

        MACS_BROAD.out.peak.set{ macs_broad_peaks }
        MACS_NARROW.out.peak.set{ macs_narrow_peaks }
        CONVERT_SICER.out.peak.set{ sicer_peaks }
        GEM.out.peak.set{ gem_peaks }
        sicer_peaks.mix(gem_peaks,
            macs_broad_peaks, macs_narrow_peaks).set{ ch_peaks }

        // Create Channel with meta, deduped bam, peak file, peak-calling tool, and chrom sizes fasta
        deduped_bam.cross(ch_peaks)
            .map{ it ->
               it.flatten()
            }
            .map{  meta1, bam, bai, meta2, peak, tool ->
                [ meta1, bam, bai, peak, tool ]
            }
            .combine(chrom_sizes)
            .set{ ch_bam_peaks }
        ch_bam_peaks | FRACTION_IN_PEAKS
        FRACTION_IN_PEAKS.out.collect() | CONCAT_FRIPS

    emit:
        ch_bam_peaks
}
