include { BWA_MEM           } from "../../modules/CCBR/bwa/mem"
include { FILTER_QUALITY    } from "../../modules/local/align.nf"
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/samtools/flagstat/main'
include { SAMTOOLS_SORT     } from '../../modules/CCBR/samtools/sort/main' // TODO use ccbr samtools/sort

workflow ALIGN_GENOME {

    take:
        reads
        reference

    main:
        BWA_MEM(reads, reference)
        FILTER_QUALITY( BWA_MEM.out.bam )
        SAMTOOLS_SORT( FILTER_QUALITY.out.bam )
        SAMTOOLS_FLAGSTAT( SAMTOOLS_SORT.out.bam )

        ch_versions = Channel.empty().mix(
            BWA_MEM.out.versions,
            SAMTOOLS_FLAGSTAT.out.versions
        )

    emit:
        bam = FILTER_QUALITY.out.bam
        flagstat = SAMTOOLS_FLAGSTAT.out.flagstat
        versions = ch_versions
}
