#!/usr/bin/env python3

"""
source: https://github.com/nf-core/chipseq/blob/51eba00b32885c4d0bec60db3cb0a45eb61e34c5/bin/check_samplesheet.py
"""

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/chipseq samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:
    sample,rep,fastq_1,fastq_2,antibody,control
    SPT5_T0_REP1,1,SRR1822153_1.fastq.gz,SRR1822153_2.fastq.gz,SPT5,SPT5_INPUT_REP1
    SPT5_T0_REP2,2,SRR1822154_1.fastq.gz,SRR1822154_2.fastq.gz,SPT5,SPT5_INPUT_REP2
    SPT5_INPUT_REP1,1,SRR5204809_Spt5-ChIP_Input1_SacCer_ChIP-Seq_ss100k_R1.fastq.gz,SRR5204809_Spt5-ChIP_Input1_SacCer_ChIP-Seq_ss100k_R2.fastq.gz,,
    SPT5_INPUT_REP2,2,SRR5204810_Spt5-ChIP_Input2_SacCer_ChIP-Seq_ss100k_R1.fastq.gz,SRR5204810_Spt5-ChIP_Input2_SacCer_ChIP-Seq_ss100k_R2.fastq.gz,,
    For an example see:
    https://raw.githubusercontent.com/nf-core/test-datasets/chipseq/samplesheet/v2.0/samplesheet_test.csv
    """

    sample_mapping_dict = {}
    with open(file_in, "r", encoding="utf-8-sig") as fin:
        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "rep", "fastq_1", "fastq_2", "antibody", "control"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print(
                f"ERROR: Please check samplesheet header -> {','.join(header)} != {','.join(HEADER)}"
            )
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(lspl) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in lspl if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(
                        MIN_COLS
                    ),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, rep, fastq_1, fastq_2, antibody, control = lspl[:]
            print(lspl)
            if sample.find(" ") != -1:
                print(
                    f"WARNING: Spaces have been replaced by underscores for sample: {sample}"
                )
                sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastQ file extension
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line",
                            line,
                        )

            ## Check antibody and control columns have valid values
            if antibody:
                if antibody.find(" ") != -1:
                    print(
                        f"WARNING: Spaces have been replaced by underscores for antibody: {antibody}"
                    )
                    antibody = antibody.replace(" ", "_")
                if not control:
                    print_error(
                        "Both antibody and control columns must be specified!",
                        "Line",
                        line,
                    )
            if control:
                if control.find(" ") != -1:
                    print(
                        f"WARNING: Spaces have been replaced by underscores for control: {control}"
                    )
                    control = control.replace(" ", "_")
                if not antibody:
                    print_error(
                        "Both antibody and control columns must be specified!",
                        "Line",
                        line,
                    )

            ## Auto-detect paired-end/single-end
            if not sample or not fastq_1:
                print_error("Invalid combination of columns provided!", "Line", line)
            is_single = str(int(bool(fastq_1 and not fastq_2)))
            # prepare sample info
            sample_basename = sample.rstrip(rep)
            sample_info = [
                sample_basename,
                rep,
                is_single,
                fastq_1,
                fastq_2,
                antibody,
                control,
            ]
            print(sample_info)

            ## Create sample mapping dictionary = {sample: [[ single_end, fastq_1, fastq_2, antibody, control ]]}
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            fout.write(
                ",".join(
                    [
                        "sample",
                        "sample_basename",
                        "rep",
                        "single_end",
                        "fastq_1",
                        "fastq_2",
                        "antibody",
                        "control",
                    ]
                )
                + "\n"
            )
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
                if not all(
                    x[0] == sample_mapping_dict[sample][0][0]
                    for x in sample_mapping_dict[sample]
                ):
                    print_error(
                        f"Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end!",
                        "Sample",
                        sample,
                    )

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    control = val[-1]
                    if control and control not in sample_mapping_dict.keys():
                        print_error(
                            f"Control identifier has to match does a provided sample identifier!",
                            "Control",
                            control,
                        )
                    plus_T = (
                        f"_T{idx+1}" if len(sample_mapping_dict[sample]) > 1 else ""
                    )  # do not append _T{idx} if not needed
                    fout.write(",".join([f"{sample}{plus_T}"] + val) + "\n")
    else:
        print_error(f"No entries to process!", "Samplesheet: {file_in}")


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT)


if __name__ == "__main__":
    sys.exit(main())
