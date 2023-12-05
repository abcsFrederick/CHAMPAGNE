include { CHECK_CONTRASTS          } from "../../../modules/local/check_contrasts/"
include { PREP_DIFFBIND            } from "../../../modules/local/diffbind/prep/"

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
            .set{ contrasts }
        bam_peaks.map{ meta, bam, bai, peak, tool ->
            [ meta.sample_basename, meta + [tool: tool], bam, bai, peak ]
            }
            .combine( contrasts )
            .map{ sample_basename1, peak_meta, bam, bai, peak, sample_basename2, con_meta ->
                sample_basename1 == sample_basename2 ? [ peak_meta + con_meta, bam, bai, peak ] : null
            }
            .unique()
            .set{ ch_peaks_contrasts }
        ch_peaks_contrasts
            .map{ meta, bam, bai, peak ->
                [ "${meta.contrast}.${meta.tool}", meta, bam, bai, peak ]
            }
            .groupTuple()
            .map{ it -> // drop meta.contrast
              it[1..-1]
            }
            .view()
            | PREP_DIFFBIND
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
