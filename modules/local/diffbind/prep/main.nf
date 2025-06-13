
process CONCAT_CSV {
    tag { "${csv_filename}" }
    label 'process_low'

    container "${params.containers_base}"

    input:
        tuple val(csv_filename), val(csv_rows)

    output:
        path(csv_filename), emit: csv

    script:
    def csv_header = ['SampleID', "Replicate", 'Condition', 'bamReads', "ControlID", "bamControl", 'Peaks', 'PeakCaller'].join(',')
    def csv_rows_joined = csv_rows.join(System.lineSeparator)
    """
    echo -e "${csv_header}" > ${csv_filename}
    echo -e "${csv_rows_joined}" >> ${csv_filename}
    """

    stub:
    """
    touch ${csv_filename}
    """
}
