
process DEEPTOOLS_BAMCOV {
    tag { meta.id }
    label 'qc'

    input:
        tuple val(meta), path(bam), path(bai)
        path(ppqt)

    output:
        val(meta.id), emit: meta_id
        path("${meta.id}.bw"), emit: bigwig

    script: // https://deeptools.readthedocs.io/en/2.1.0/content/tools/bamCoverage.html
    """
    frag_len=\$(cut -f 2 ${ppqt})
    bamCoverage \
      --bam ${bam} \
      -o ${meta.id}.bw \
      --binSize ${params.deeptools.bin_size} \
      --smoothLength ${params.deeptools.smooth_length} \
      --ignoreForNormalization ${params.deeptools.excluded_chroms} \
      --numberOfProcessors ${task.cpus} \
      --normalizeUsing ${params.deeptools.normalize_using} \
      --effectiveGenomeSize ${params.align.effective_genome_size} \
      --extendReads \$frag_len
    """

    stub:
    """
    touch ${meta.id}.bw
    """

}
process DEEPTOOLS_BIGWIG_SUM {
    label 'qc'

    input:
        val(meta_ids)
        path(bigwigs)

    output:
        path("bigWigSum.npz")

    script:
    """
    multiBigwigSummary bins \
      -b ${bigwigs} \
      --labels ${meta_ids} \
      -o bigWigSum.npz
    """

    stub:
    """
    echo "${bigwigs}" > bigWigSum.npz
    """
}

process DEEPTOOLS_PLOTS {
    label 'qc'

    input:
        path(array)

    output:
        tuple path("*.pdf"), path("*.png"), emit: plots

    script:
    """
    plotCorrelation \
      -in ${array} \
      -o spearman_heatmap.pdf \
      -c 'spearman' \
      -p 'heatmap' \
      --skipZeros \
      --removeOutliers \
      --plotNumbers
    plotCorrelation \
      -in ${array} \
      -o spearman_scatterplot.pdf \
      -c 'spearman' \
      -p 'scatterplot' \
      --skipZeros \
      --removeOutliers
    plotPCA \
      -in ${array} \
      -o pca.pdf
    plotCorrelation \
      -in ${array} \
      -o spearman_heatmap_mqc.png \
      -c 'spearman' \
      -p 'heatmap' \
      --skipZeros \
      --removeOutliers \
      --plotNumbers
    """

    stub:
    """
    touch spearman_heatmap.pdf spearman_scatterplot.pdf pca.pdf spearman_heatmap_mqc.png
    """
}

process DEEPTOOLS_FINGERPRINT { // TODO aggregate metadata and bamfiles.
    label: 'qc'

    input:
        tuple val(meta), path(bams)

    output:
        path("quality_metrics.tsv"), emit: metrics
        path("fingerprint.pdf"), emit: plot

    script:
    // TODO bams should be space-delimited list of bam files. labels should be antibodies + inputs
    // TODO -e only for single end. should it be hardcoded here? above it's taken from ppqt.
    """
    plotFingerprint \
      -b ${bams} \
      --labels ${labels} \
      -p ${task.cpus} \
      --skipZeros \
      --outQualityMetrics quality.tsv \
       --plotFile fingerprint.pdf \
       -e 200
    """
}
