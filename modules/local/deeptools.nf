
process BAM_COVERAGE {
    tag { meta.id }
    label 'qc'
    label 'deeptools'

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
process BIGWIG_SUM {
    label 'qc'
    label 'deeptools'

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

process ARRAY_PLOTS {
    label 'qc'
    label 'deeptools'

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

process FINGERPRINT {
    label 'qc'
    label 'deeptools'

    input:
      tuple val(meta), path(bams), path(bais)

    output:
      tuple val(meta), path("*.pdf")          , emit: pdf
      tuple val(meta), path("*.matrix.txt")      , emit: matrix
      tuple val(meta), path("*.qcmetrics.txt"), emit: metrics

    script:
    // TODO handle extendReads for single vs paired https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/modules/nf-core/modules/deeptools/plotfingerprint/main.nf
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    plotFingerprint \\
      --bamfiles ${bams.join(' ')} \\
      --plotFile ${prefix}.plotFingerprint.plot.pdf \\
      --outRawCounts ${prefix}.plotFingerprint.matrix.txt \\
      --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
      --numberOfProcessors ${task.cpus} \\
      --skipZeros \\
      --extendReads 200
    """

}

/*
TODO: process deeptools_genes only on protein coding for now https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/InitialChIPseqQC.snakefile#L704-L743
*/
