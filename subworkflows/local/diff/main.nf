include { CHECK_CONTRASTS          } from "../../../modules/local/check_contrasts/"

workflow DIFF {
    take:
        peaks
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
        peaks.map{ meta, bed ->
            [ meta.sample_basename, meta, bed]
            }
            .cross( contrasts )
            .map{ it.flatten() }
            .map{ sample_basename1, peak_meta, bed, sample_basename2, con_meta ->
                assert sample_basename1 == sample_basename2
                meta = peak_meta + con_meta
                [ meta, bed ]
            }
            .set{ ch_peaks_contrasts }

    emit:
        diff_peaks = peaks

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
