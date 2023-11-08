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
         input        : ${params.input}
         genome       : ${params.genome}
         """
         .stripIndent()

// SUBWORKFLOWS

include { INPUT_CHECK              } from './subworkflows/local/input_check.nf'
include { PREPARE_GENOME           } from './subworkflows/local/prepare_genome.nf'
include { FILTER_BLACKLIST         } from './subworkflows/CCBR/filter_blacklist'
include { ALIGN_GENOME             } from "./subworkflows/local/align.nf"
include { DEDUPLICATE              } from "./subworkflows/local/deduplicate.nf"
include { QC                       } from './subworkflows/local/qc.nf'
include { CALL_PEAKS               } from './subworkflows/local/peaks.nf'


// MODULES
include { CUTADAPT                 } from "./modules/CCBR/cutadapt"
include { PHANTOM_PEAKS            } from "./modules/local/qc.nf"
include { PPQT_PROCESS
          MULTIQC                  } from "./modules/local/qc.nf"
include { NORMALIZE_INPUT          } from "./modules/local/deeptools.nf"

workflow.onComplete {
    if (!workflow.stubRun && !workflow.commandLine.contains('-preview')) {
        println "Running spooker"
        def message = Utils.spooker(workflow)
        if (message) {
            println message
        }
    }
}

workflow MAKE_REFERENCE {
    PREPARE_GENOME()
}

// MAIN WORKFLOW
workflow {
    CHIPSEQ()
}

workflow CHIPSEQ {
    INPUT_CHECK(file(params.input), params.seq_center)
    INPUT_CHECK.out.reads.set { raw_fastqs }
    raw_fastqs | CUTADAPT
    CUTADAPT.out.reads.set{ trimmed_fastqs }

    PREPARE_GENOME()
    chrom_sizes = PREPARE_GENOME.out.chrom_sizes

    effective_genome_size = PREPARE_GENOME.out.effective_genome_size
    FILTER_BLACKLIST(trimmed_fastqs, PREPARE_GENOME.out.blacklist_index)
    ALIGN_GENOME(FILTER_BLACKLIST.out.reads, PREPARE_GENOME.out.reference_index)
    ALIGN_GENOME.out.bam.set{ aligned_bam }

    DEDUPLICATE(aligned_bam, chrom_sizes, effective_genome_size)
    DEDUPLICATE.out.bam.set{ deduped_bam }
    DEDUPLICATE.out.tag_align.set{ deduped_tagalign }

    deduped_bam | PHANTOM_PEAKS
    PHANTOM_PEAKS.out.fraglen | PPQT_PROCESS
    PPQT_PROCESS.out.fraglen.set {frag_lengths }

    ch_multiqc = Channel.of()
    if (params.run.qc) {
        QC(raw_fastqs, trimmed_fastqs, FILTER_BLACKLIST.out.n_surviving_reads,
           aligned_bam, ALIGN_GENOME.out.aligned_flagstat, ALIGN_GENOME.out.filtered_flagstat,
           deduped_bam, DEDUPLICATE.out.flagstat,
           PHANTOM_PEAKS.out.spp, frag_lengths,
           PREPARE_GENOME.out.gene_info,
           effective_genome_size
           )
        QC.out.bigwigs.set{ ch_ip_ctrl_bigwig }

        if (params.run.normalize_input) {
            ch_ip_ctrl_bigwig | NORMALIZE_INPUT
        }

        ch_multiqc = ch_multiqc.mix(QC.out.multiqc_input)
    }
    if (params.run.call_peaks && [params.run.macs_broad, params.run.macs_narrow, params.run.gem, params.run.sicer].any()) {
        CALL_PEAKS(chrom_sizes,
                   PREPARE_GENOME.out.chrom_dir,
                   deduped_tagalign,
                   deduped_bam,
                   frag_lengths,
                   effective_genome_size,
                   PREPARE_GENOME.out.fasta
                   )
        ch_multiqc = ch_multiqc.mix(CALL_PEAKS.out.plots)
    }

    MULTIQC(
        file(params.multiqc.config, checkIfExists: true),
        file(params.multiqc.logo, checkIfExists: true),
        ch_multiqc.collect()
    )

}
