log.info """\
         CHAMPAGNE 🍾
         =============
         NF version   : $nextflow.version
         runName      : $workflow.runName
         username     : $workflow.userName
         configs      : $workflow.configFiles
         profile      : $workflow.profile
         cmd line     : $workflow.commandLine
         start time   : $workflow.start
         projectDir   : $workflow.projectDir
         launchDir    : $workflow.launchDir
         workDir      : $workflow.workDir
         homeDir      : $workflow.homeDir
         reads        : ${params.reads}
         """
         .stripIndent()

include { TRIM_SE } from "./modules/local/trim.nf"
include { FASTQC as FASTQC_RAW } from "./modules/local/qc.nf"
include { FASTQC as FASTQC_TRIMMED } from "./modules/local/qc.nf"
include { FASTQ_SCREEN } from "./modules/local/qc.nf"
include { ALIGN_BLACKLIST } from "./modules/local/qc.nf"
include { ALIGN_GENOME } from "./modules/local/qc.nf"
include { INDEX_BAM } from "./modules/local/qc.nf"
include { PRESEQ } from "./modules/local/qc.nf"
include { PHANTOM_PEAKS } from "./modules/local/qc.nf"

workflow {
  raw_fastqs = Channel
                    .fromPath(params.reads)
                    .map { file -> tuple(file.simpleName, file) }
  raw_fastqs.combine(Channel.value("raw")) | FASTQC_RAW
  raw_fastqs | TRIM_SE
  trimmed_fastqs = TRIM_SE.out
  trimmed_fastqs.combine(Channel.value("trimmed")) | FASTQC_TRIMMED
  //trimmed_fastqs = TRIM_SE.out | map { it -> it[1] } | collect // obtain list of just fastqs without sample_ids
  trimmed_fastqs.combine(Channel.fromPath(params.fastq_screen.conf)) | FASTQ_SCREEN
  blacklist_files = Channel
                    .fromPath("${params.align.index_dir}${params.align.blacklist}*")
                    .collect()
  ALIGN_BLACKLIST(trimmed_fastqs, blacklist_files)
  reference_files = Channel
                    .fromPath("${params.align.index_dir}${params.align.genome}*")
                    .collect()
  ALIGN_GENOME(ALIGN_BLACKLIST.out, reference_files)
  //PRESEQ(ALIGN_GENOME.out)
  INDEX_BAM(ALIGN_GENOME.out)
  PHANTOM_PEAKS(ALIGN_GENOME.out)
}
