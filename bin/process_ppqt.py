#!/usr/bin/env python

"""
Process fragment length output from PPQT.
Refactor of https://github.com/CCBR/Pipeliner/blob/86c6ccaa3d58381a0ffd696bbf9c047e4f991f9e/Rules/InitialChIPseqQC.snakefile#L513-L541
"""

import sys
import warnings


def main(args):
    fraglen_filename = sys.argv[1]
    min_frag_len = int(sys.argv[2])
    print(get_fragment_length(fraglen_filename, min_frag_len))


def get_fragment_length(fraglen_filename, min_frag_len):
    with open(fraglen_filename, "r") as infile:
        fragment_length = int(infile.read().strip())
    if fragment_length < min_frag_len:
        warnings.warn(
            f"The estimated fragment length was {fragment_length}. Using default of {min_frag_len} instead."
        )
        fragment_length = min_frag_len
    return fragment_length


if __name__ == "__main__":
    main(sys.argv)
