// https://github.com/CCBR/ChiP-seek/blob/9ba449e4855f9710e86f2db7c1d9560de634b3f1/workflow/rules/align.smk#L21
// https://github.com/nf-core/ampliseq/blob/dev/subworkflows/local/cutadapt_workflow.nf
process TRIM_SE {
  tag { sample_id }
  publishDir "$params.outdir/qc/$sample_id/trimmed", mode: "$params.filePublishMode"

  input:
    tuple val(sample_id), path(fastq)

  output:
    tuple val(sample_id), path("${sample_id}.fastq.gz")

  script:
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
  pigz -p ${task.cpus} > ${sample_id}.trimmed.fastq.gz
  nseqs_trimmed=\$(zgrep "^@" ${sample_id}.trimmed.fastq.gz | wc -l)
  echo "\$nseqs_trimmed in${sample_id}.trimmed.fastq.gz"
  """

  stub:
  """
  touch ${sample_id}.fastq.gz
  """
}
