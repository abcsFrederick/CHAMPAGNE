
include { BAM_COVERAGE             } from "../../modules/local/deeptools.nf"
include { BIGWIG_SUM               } from "../../modules/local/deeptools.nf"
include { BED_PROTEIN_CODING       } from "../../modules/local/deeptools.nf"
include { COMPUTE_MATRIX           } from "../../modules/local/deeptools.nf"
include { PLOT_FINGERPRINT         } from "../../modules/local/deeptools.nf"
include { PLOT_CORRELATION         } from "../../modules/local/deeptools.nf"
include { PLOT_PCA                 } from "../../modules/local/deeptools.nf"
include { PLOT_HEATMAP             } from "../../modules/local/deeptools.nf"
include { PLOT_PROFILE             } from "../../modules/local/deeptools.nf"

workflow DEEPTOOLS {
    take:
        deduped_bam
        frag_lengths
        effective_genome_size
        gene_info

    main:

        deduped_bam.join(frag_lengths).combine(effective_genome_size) | BAM_COVERAGE
        BAM_COVERAGE.out.bigwig.collect().set{ bigwig_list }
        BIGWIG_SUM(bigwig_list)
        BIGWIG_SUM.out.array.combine(Channel.from('heatmap', 'scatterplot')) | PLOT_CORRELATION
        BIGWIG_SUM.out.array | PLOT_PCA

        // Create channel: [ meta, [ ip_bam, control_bam ] [ ip_bai, control_bai ] ]
        deduped_bam
            .combine(deduped_bam)
            .map {
                meta1, bam1, bai1, meta2, bam2, bai2 ->
                    meta1.control == meta2.id ? [ meta1, [ bam1, bam2 ], [ bai1, bai2 ] ] : null
            }
            .set { ch_ip_ctrl_bam_bai }
        ch_ip_ctrl_bam_bai | PLOT_FINGERPRINT
        gene_info | BED_PROTEIN_CODING
        COMPUTE_MATRIX(bigwig_list,
                       BED_PROTEIN_CODING.out.bed.combine(Channel.from('metagene','TSS'))
        )
        PLOT_HEATMAP(COMPUTE_MATRIX.out.mat)
        PLOT_PROFILE(COMPUTE_MATRIX.out.mat)

        // Create channel: [ meta, ip_bw, control_bw ]
        BAM_COVERAGE.out.meta
            .merge(BAM_COVERAGE.out.bigwig)
            .set { bigwigs }
        bigwigs
            .combine(bigwigs)
            .map {
                meta1, bw1, meta2, bw2 ->
                    meta1.control == meta2.id ? [ meta1, bw1, bw2 ] : null
            }
            .set { ch_ip_ctrl_bigwig }

    emit:
        bigwig              = ch_ip_ctrl_bigwig
        fingerprint_matrix  = PLOT_FINGERPRINT.out.matrix
        fingerprint_metrics = PLOT_FINGERPRINT.out.metrics
        corr                = PLOT_CORRELATION.out.tab
        pca                 = PLOT_PCA.out.tab
        profile             = PLOT_PROFILE.out.tab
}
