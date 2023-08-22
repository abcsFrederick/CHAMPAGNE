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

workflow convert2fasta {
  Channel.fromPath(params.input) | any2fasta | view
}

workflow {
  raw_fastqs = Channel
                    .fromPath(params.reads)
                    .map { file -> tuple(file.simpleName, file) }
  raw_fastqs | FASTQC_RAW
  raw_fastqs | TRIM_SE | FASTQC_TRIMMED
}
