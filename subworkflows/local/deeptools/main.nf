include { BAM_COVERAGE
          NORMALIZE_INPUT
          BIGWIG_SUM as BIGWIG_SUM_NORM
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

        // Create channel: [ meta, ip_bw, control_bw ]
        bigwigs
            .combine(bigwigs)
            .map {
                meta1, bw1, meta2, bw2 ->
                    meta1.control == meta2.id ? [ meta1, bw1, bw2 ] : null
            }
            .set { ch_ip_ctrl_bigwig }
        ch_ip_ctrl_bigwig | NORMALIZE_INPUT
        NORMALIZE_INPUT.out.bigwig.map{ meta, bigwig -> bigwig }.set{ bigwigs_norm }
        BIGWIG_SUM_NORM(bigwigs_norm.collect())
        BIGWIG_SUM_NORM.out.array.combine(Channel.of('heatmap', 'scatterplot')) | PLOT_CORRELATION
        BIGWIG_SUM_NORM.out.array | PLOT_PCA

        // group raw bigwigs by sample basename to group replicates & inputs together
        bigwigs_raw = ch_ip_ctrl_bigwig
            .map{ meta, sample_bw, control_bw ->
                [ [ id: meta.sample_basename ], sample_bw ]
            }
            .concat(ch_ip_ctrl_bigwig
                .map{ meta, sample_bw, control_bw ->
                    [ [ id: meta.sample_basename], control_bw ]
                }
            )
        gene_info | BED_PROTEIN_CODING
        beds = BED_PROTEIN_CODING.out.bed_prot.mix(BED_PROTEIN_CODING.out.bed_all)

        // create plots with:
        //    - raw and normalized bigwigs,
        //    - protein coding and all genes
        //    - metagene or TSS
        ch_all_bigwigs = Channel.value([ id: 'norm' ])
            .combine(bigwigs_norm)
            .mix(bigwigs_raw)
            .groupTuple()
            .combine(beds)
            .combine(Channel.of('metagene','TSS'))
        COMPUTE_MATRIX(ch_all_bigwigs)
        PLOT_HEATMAP(COMPUTE_MATRIX.out.mat)
        PLOT_PROFILE(COMPUTE_MATRIX.out.mat)

        ch_controls = deduped_bam
            .combine(deduped_bam)
            .map {
                meta1, bam1, bai1, meta2, bam2, bai2 ->
                    meta1.control == meta2.id ? [ [id: meta1.sample_basename], bam2, bai2 ] : null
            }
        ch_samples = deduped_bam
            .combine(deduped_bam)
            .map {
                meta1, bam1, bai1, meta2, bam2, bai2 ->
                    meta1.control == meta2.id ? [ [id: meta1.sample_basename], bam1, bai1 ] : null
            }
        ch_controls.mix(ch_samples)
            .groupTuple()
            .set { ch_ip_ctrl_bam_bai }
        ch_ip_ctrl_bam_bai | PLOT_FINGERPRINT

    emit:
        fingerprint_matrix  = PLOT_FINGERPRINT.out.matrix
        fingerprint_metrics = PLOT_FINGERPRINT.out.metrics
        corr                = PLOT_CORRELATION.out.tab
        pca                 = PLOT_PCA.out.tab
        profile             = PLOT_PROFILE.out.tab
}
