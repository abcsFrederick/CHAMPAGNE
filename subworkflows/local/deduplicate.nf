include { MACS2_DEDUP; PICARD_DEDUP      } from "../../modules/local/deduplicate.nf"
include { SAMTOOLS_INDEX as INDEX_SINGLE
          SAMTOOLS_INDEX as INDEX_PAIRED } from "../../modules/local/samtools_index.nf"
workflow DEDUPLICATE {
    take:
        aligned_bam
        chrom_sizes

    main:

        // branch reads: single end to macs2 filterdup, paired end to picard markduplicates
        aligned_bam.branch { meta, bam ->
            single: meta.single_end
            paired: !meta.single_end

        }.set{ bam }
        bam.single.combine(chrom_sizes) | MACS2_DEDUP
        MACS2_DEDUP.out.bam | INDEX_SINGLE
        bam.paired | PICARD_DEDUP
        PICARD_DEDUP.out.bam | INDEX_PAIRED
        INDEX_PAIRED.out.bam
            .mix(INDEX_SINGLE.out.bam)
            .set{ ch_bam_bai }

        // mix single bed and paired bam channel for peak callers, with variable to designate format
        INDEX_PAIRED.out.bam.combine(Channel.value("BAMPE"))
            .mix(MACS2_DEDUP.out.bed.combine(Channel.value("NO_BAI")).combine(Channel.value("TagAlign")))
            .set{ mixed_tagAlign }
        mixed_tagAlign | view


    emit:
        bam = ch_bam_bai
        tag_align = mixed_tagAlign
}
