include { CHECK_CONTRASTS   } from "../../../modules/local/check_contrasts/"
include { PREP_DIFFBIND     } from "../../../modules/local/diffbind/prep/"
include { RMARKDOWNNOTEBOOK } from '../../../modules/nf-core/rmarkdownnotebook/'

workflow DIFF {
    take:
        bam_peaks
        samplesheet_file
        contrasts_file

    main:

        CHECK_CONTRASTS(samplesheet_file, contrasts_file)
            .csv
            .flatten()
            .splitCsv( header: true, sep: ',' )
            .map{ it ->
                meta = get_contrast_meta(it)
                [ meta.sample_basename, [group: meta.group, contrast: meta.contrast] ]
            }
            .unique()
            .set{ contrasts }
        bam_peaks
            .combine( contrasts )
            .map{ sample_basename1, peak_meta, bam, bai, peak, ctrl_bam, ctrl_bai, sample_basename2, con_meta ->
                sample_basename1 == sample_basename2 ? [ peak_meta + con_meta, bam, bai, peak, ctrl_bam, ctrl_bai ] : null
            }
            .unique()
            .set{ ch_peaks_contrasts }
        ch_peaks_contrasts | PREP_DIFFBIND

        ch_peaks_contrasts
            .map{ meta, bam, bai, peak, ctrl_bam, ctrl_bai ->
                [ [contrast: meta.contrast, tool: meta.tool], [bam, bai, peak, ctrl_bam, ctrl_bai] ]
            }
            .groupTuple()
            .map{ meta, files ->
                [ meta, files.flatten().unique() ]
            }
            .set{ data_files }

        PREP_DIFFBIND.out.csv
            .collectFile(storeDir: "${params.outdir}/diffbind/contrasts") { meta, row ->
                [ "${meta.contrast}.${meta.tool}.csv", row ]
            }
            .map{ contrast_file ->
                meta = [:]
                meta_list = contrast_file.baseName.tokenize('.')
                meta.contrast = meta_list[0]
                meta.tool = meta_list[1]
                [ meta, contrast_file ]
            }
            .cross(data_files)
            .map{ contrast, files ->
                meta = contrast[0]
                meta.csvfile = contrast[1]
                file_list = [contrast[1], files[1]].flatten().unique()
                [ meta, file_list ] // meta, csv, filelist
            }
            .set{contrast_files}
        contrast_files
            .combine(Channel.fromPath(file(params.diffbind.report, checkIfExists: true)))
            .map{ meta, files, rmarkdown ->
                [ meta, rmarkdown ]
            }
            .set{ch_rmarkdown}
        contrast_files.map{ meta, files -> files}.flatten().unique().collect().set{file_list}
        RMARKDOWNNOTEBOOK( ch_rmarkdown, [], file_list )

    emit:
        diff_peaks = bam_peaks

}

def get_contrast_meta(LinkedHashMap row) {
    def meta = [:]
    meta.id              = row.sample
    meta.sample_basename = row.sample_basename
    meta.rep             = row.rep
    meta.single_end      = row.single_end.toBoolean()
    meta.antibody        = row.antibody
    meta.control         = row.control
    meta.group           = row.group
    meta.contrast        = row.contrast

    return meta
}
