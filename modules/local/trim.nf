// https://github.com/CCBR/ChiP-seek/blob/9ba449e4855f9710e86f2db7c1d9560de634b3f1/workflow/rules/align.smk#L21
// https://github.com/nf-core/ampliseq/blob/dev/subworkflows/local/cutadapt_workflow.nf
// https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/modules/nf-core/modules/trimgalore/main.nf
process TRIM_SE {
  tag { meta.id }
  label 'qc'
  label 'process_high'

  input:
    tuple val(meta), path(fastq)

  output:
    tuple val(meta), path("*.fastq.gz")

  script:
  def prefix = task.ext.prefix ?: "${meta.id}"
  if (meta.single_end) {
    """
    nseqs_raw=\$(zgrep "^@" ${fastq} | wc -l)
    echo "\$nseqs_raw in ${fastq}"
    cutadapt \
      --nextseq-trim=2 \
      --trim-n \
      -n 5 -O 5 \
      -q ${params.cutadapt.leadingquality},${params.cutadapt.trailingquality} \
      -m ${params.cutadapt.minlen} \
      -b file:${params.cutadapt.adapters} \
      -j $task.cpus \
      $fastq |\
    pigz -p ${task.cpus} > ${prefix}.trimmed.fastq.gz
    nseqs_trimmed=\$(zgrep "^@" ${prefix}.trimmed.fastq.gz | wc -l)
    echo "\$nseqs_trimmed in ${prefix}.trimmed.fastq.gz"
    """
  } else {
    "paired end reads are not supported yet"
  }

  stub:
  """
  touch ${meta.id}.trimmed.fastq.gz
  """
}
