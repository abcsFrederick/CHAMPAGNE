
process BAM_COVERAGE {
    tag { meta.id }
    label 'qc'
    label 'deeptools'

    input:
        tuple val(meta), path(bam), path(bai)
        path(ppqt)

    output:
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
        path(bigwigs)

    output:
        path("bigWigSum.npz"), emit: array

    script:
    """
    multiBigwigSummary bins \\
      -b ${bigwigs.join(' ')} \\
      --smartLabels \\
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
    tuple val(meta), path("*.mat.txt")      , emit: matrix
    tuple val(meta), path("*.qcmetrics.txt"), emit: metrics

  script:
  // TODO handle extendReads for single vs paired https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/modules/nf-core/modules/deeptools/plotfingerprint/main.nf
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  plotFingerprint \\
    --bamfiles ${bams.join(' ')} \\
    --plotFile ${prefix}.plotFingerprint.plot.pdf \\
    --outRawCounts ${prefix}.plotFingerprint.mat.txt \\
    --outQualityMetrics ${prefix}.plotFingerprint.qcmetrics.txt \\
    --numberOfProcessors ${task.cpus} \\
    --skipZeros \\
    --smartLabels \\
    --extendReads 200
  """

  stub:
  """
  for ext in plot.pdf mat.txt qcmetrics.txt; do
    touch ${meta.id}.plotFingerprint.\$ext
  done
  """

}

process BED_PROTEIN_CODING {
    label 'qc'
    label 'deeptools'

    input:
      path(bed)

    output:
      path("*.bed"), emit: bed

    script:
    """
    grep --line-buffered 'protein_coding' ${bed} | awk -v OFS='\t' -F'\t' '{{print \$1, \$2, \$3, \$5, \".\", \$4}}' > ${params.align.genome}.protein_coding.bed
    """

    stub:
    """
    touch ${params.align.genome}.protein_coding.bed
    """
}

// TODO: best way to have one COMPUTE_MATRIX process with different arguments for metagene vs TSS?
process COMPUTE_MATRIX_METAGENE {
  label 'qc'
  label 'deeptools'

  input:
    path(bigwigs)
    path(bed)

  output:
    path("*.mat.gz"), emit: matrix

  script:
  """
  computeMatrix scale-regions \\
    -S ${bigwigs.join(' ')} \\
    -R ${bed} \\
    -p ${task.cpus} \\
    -o metagene.mat.gz \\
    --upstream 1000 \\
    --regionBodyLength 2000 \\
    --downstream 1000 \\
    --skipZeros \\
    --smartLabels
  """

  stub:
  """
  touch metagene.mat.gz
  """
}

process COMPUTE_MATRIX_TSS {
  label 'qc'
  label 'deeptools'

  input:
    path(bigwigs)
    path(bed)

  output:
    path("*.mat.gz"), emit: matrix

  script:
  """
  computeMatrix reference-point \\
    -S ${bigwigs.join(' ')} \\
    -R ${bed} \\
    -p ${task.cpus} \\
    -o TSS.mat.gz \\
    --referencePoint TSS \\
    --upstream 3000 \\
    --downstream 3000 \\
    --skipZeros \\
    --smartLabels
  """

  stub:
  """
  touch TSS.mat.gz
  """
}
process PLOT_HEATMAP {
  label 'qc'
  label 'deeptools'

  input:
    path(mat)

  output:
    path("*.pdf")

  script:
  // set colorMap to "BuGn" if "metagene" in matrix filename, otherwise use "BuPu"
  def color_map = mat.baseName.contains('metagene') ? 'BuGn' : 'BuPu'
  """
  plotHeatmap \\
    -m ${mat} \\
    -out ${mat.baseName}.heatmap.pdf \\
    --colorMap ${color_map} \\
    --yAxisLabel 'average RPGC' \\
    --regionsLabel 'genes' \\
    --legendLocation 'none'"
  """

  stub:
  """
  touch ${mat.baseName}.heatmap.pdf
  """
}

process PLOT_PROFILE {
  label 'qc'
  label 'deeptools'

  input:
    path(mat)

  output:
    path("*.pdf")

  script: // TODO set plotType for SE vs paired
  def legend_loc = mat.baseName.contains('metagene') ? 'upper-right' : 'upper-left'
  """
  plotProfile \\
    -m ${mat} \\
    -out ${mat.baseName}.lineplot.pdf \\
    --plotHeight 15 \\
    --plotWidth 15 \\
    --perGroup \\
    --yAxisLabel 'average RPGC' \\
    --plotType 'se' \\
    --legendLocation ${legend_loc}"
  """

  stub:
  """
  touch ${mat.baseName}.lineplot.pdf
  """
}
