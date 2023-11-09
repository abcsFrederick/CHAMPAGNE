include { BAM_COVERAGE
          NORMALIZE_INPUT
          BIGWIG_SUM
          BED_PROTEIN_CODING
          COMPUTE_MATRIX
          PLOT_FINGERPRINT
          PLOT_CORRELATION
          PLOT_PCA
          PLOT_HEATMAP
          PLOT_PROFILE        } from "../../../modules/local/deeptools.nf"

workflow DEEPTOOLS {
    take:
        deduped_bam
        frag_lengths
        effective_genome_size
        gene_info

    main:

        deduped_bam.join(frag_lengths).combine(effective_genome_size) | BAM_COVERAGE
        BAM_COVERAGE.out.bigwig
            .set { bigwigs }
        bigwigs.map{meta, bigwig -> bigwig}.collect().set{ bigwig_list }
        // Create channel: [ meta, ip_bw, control_bw ]
        bigwigs
            .combine(bigwigs)
            .map {
                meta1, bw1, meta2, bw2 ->
                    meta1.control == meta2.id ? [ meta1, bw1, bw2 ] : null
            }
            .set { ch_ip_ctrl_bigwig }
        ch_ip_ctrl_bigwig | NORMALIZE_INPUT
        NORMALIZE_INPUT.out.bigwig.map{ meta, bigwig -> bigwig }.collect().set{ bigwigs_norm }

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
        COMPUTE_MATRIX(bigwigs_norm,
                       BED_PROTEIN_CODING.out.bed.combine(Channel.from('metagene','TSS'))
        )
        PLOT_HEATMAP(COMPUTE_MATRIX.out.mat)
        PLOT_PROFILE(COMPUTE_MATRIX.out.mat)

    emit:
        bigwig              = ch_ip_ctrl_bigwig
        fingerprint_matrix  = PLOT_FINGERPRINT.out.matrix
        fingerprint_metrics = PLOT_FINGERPRINT.out.metrics
        corr                = PLOT_CORRELATION.out.tab
        pca                 = PLOT_PCA.out.tab
        profile             = PLOT_PROFILE.out.tab
}
