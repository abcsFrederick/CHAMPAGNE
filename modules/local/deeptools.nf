
process BAM_COVERAGE {
    tag { meta.id }
    label 'qc'
    label 'deeptools'
    label 'process_higher'

    input:
        tuple val(meta), path(bam), path(bai)
        tuple val(meta), val(fraglen)

    output:
        val(meta), emit: meta
        path("${meta.id}.bw"), emit: bigwig

    script: // https://deeptools.readthedocs.io/en/2.1.0/content/tools/bamCoverage.html
    """
    bamCoverage \
      --bam ${bam} \
      -o ${meta.id}.bw \
      --binSize ${params.deeptools.bin_size} \
      --smoothLength ${params.deeptools.smooth_length} \
      --ignoreForNormalization ${params.deeptools.excluded_chroms} \
      --numberOfProcessors ${task.cpus} \
      --normalizeUsing ${params.deeptools.normalize_using} \
      --effectiveGenomeSize ${params.align.effective_genome_size} \
      --extendReads ${fraglen}
    """

    stub:
    """
    touch ${meta.id}.bw
    """

}
process BIGWIG_SUM {
    label 'qc'
    label 'deeptools'
    label 'process_high'

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

process PLOT_CORRELATION {
    label 'qc'
    label 'deeptools'

    input:
        tuple path(array), val(plottype)

    output:
        path("*.pdf"), emit: pdf
        path("*.tab"), emit: tab

    script:
    def args = (plottype == 'heatmap')? '--plotNumbers': ''
    // TODO throw error if plottype != either'heatmap' or 'scatterplot'
    """
    plotCorrelation \\
      -in ${array} \\
      -o ${array.baseName}.spearman_${plottype}.pdf \\
      --outFileCorMatrix ${array.baseName}.spearman_${plottype}.tab \\
      -c 'spearman' \\
      -p '${plottype}' \\
      --skipZeros \\
      --removeOutliers ${args}
    """

    stub:
    """
    touch ${array.baseName}.spearman_${plottype}.pdf ${array.baseName}.spearman_${plottype}.tab
    """
}

process PLOT_PCA {
    label 'qc'
    label 'deeptools'

    input:
        path(array)

    output:
        path("*.pdf"), emit: pdf
        path("*.tab"), emit: tab

    script:
    """
    plotPCA \\
      -in ${array} \\
      -o ${array.baseName}.pca.pdf \\
      --outFileNameData ${array.baseName}.plotPCA.tab
    """

    stub:
    """
    touch ${array.baseName}.pca.pdf ${array.baseName}.plotPCA.tab
    """
}

process PLOT_FINGERPRINT {
  label 'qc'
  label 'deeptools'
  label 'process_higher'

  input:
    tuple val(meta), path(bams), path(bais)

  output:
    path("*.pdf")          , emit: pdf
    path("*.mat.txt")      , emit: matrix
    path("*.qcmetrics.txt"), emit: metrics

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

    input:
      path(bed)

    output:
      path("*.bed"), emit: bed

    script:
    """
    grep --line-buffered 'protein_coding' ${bed} | awk -v OFS='\t' -F'\t' '{{print \$1, \$2, \$3, \$5, \".\", \$4}}' > ${params.genome}.protein_coding.bed
    """

    stub:
    """
    touch ${params.genome}.protein_coding.bed
    """
}

process COMPUTE_MATRIX {
  label 'qc'
  label 'deeptools'
  label 'process_higher'

  input:
    path(bigwigs)
    tuple path(bed), val(mattype)

  output:
    path("*.mat.gz"), emit: mat

  script:
  def cmd = null
  def args = null
  if (mattype == 'TSS') {
    cmd = 'reference-point'
    args = { ['--referencePoint TSS',
              '--upstream 3000',
              '--downstream 3000'
              ].join(' ').trim()
            }
  } else if (mattype == 'metagene') {
    cmd = 'scale-regions'
    args = { ['--upstream 1000',
              '--regionBodyLength 2000',
              '--downstream 1000'
              ].join(' ').trim()
            }
  } else {
    error "Invalid matrix type: ${mattype}"
  }
  """
  echo "$mattype" > file.txt
  computeMatrix ${cmd} \\
    -S ${bigwigs.join(' ')} \\
    -R ${bed} \\
    -p ${task.cpus} \\
    -o ${mattype}.mat.gz \\
    --skipZeros \\
    --smartLabels \\
    ${args}
  """

  stub:
  """
  touch ${mattype}.mat.gz
  """
}
process PLOT_HEATMAP {
  label 'qc'
  label 'deeptools'

  input:
    path(mat)

  output:
    path("*.pdf"), emit: pdf

  script:
  // sets colorMap to "BuGn" if "metagene" in matrix filename, otherwise use "BuPu"
  def color_map = mat.baseName.contains('metagene') ? 'BuGn' : 'BuPu'
  """
  plotHeatmap \\
    -m ${mat} \\
    -out ${mat.baseName}.heatmap.pdf \\
    --colorMap ${color_map} \\
    --yAxisLabel 'average RPGC' \\
    --regionsLabel 'genes' \\
    --legendLocation 'none'
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
    path("*.pdf"), emit: pdf
    path("*.tab"), emit: tab

  script: // TODO set plotType for SE vs paired
  def legend_loc = mat.baseName.contains('metagene') ? 'upper-right' : 'upper-left'
  """
  plotProfile \\
    -m ${mat} \\
    --outFileName ${mat.baseName}.lineplot.pdf \\
    --outFileNameData ${mat.baseName}.plotProfile.tab \\
    --plotHeight 15 \\
    --plotWidth 15 \\
    --perGroup \\
    --yAxisLabel 'average RPGC' \\
    --plotType 'se' \\
    --legendLocation ${legend_loc}
  """

  stub:
  """
  touch ${mat.baseName}.lineplot.pdf ${mat.baseName}.plotProfile.tab
  """
}

process NORMALIZE_INPUT {
  label 'qc'
  label 'deeptools'
  label 'process_higher'

  input:
    tuple val(meta), path(chip), path(input)

  output:
    tuple val(meta), path("*.norm.bw"), emit: bigwig

  script:
  """
  bigwigCompare \\
  --binSize ${params.deeptools.bin_size} \\
  --outFileName ${meta.id}.norm.bw \\
  --outFileFormat 'bigwig' \\
  --bigwig1 ${chip} \\
  --bigwig2 ${input} \\
  --operation 'subtract' \\
  --skipNonCoveredRegions \\
  --numberOfProcessors ${task.cpus}
  """

  stub:
  """
  touch ${meta.id}.norm.bw
  """
}
