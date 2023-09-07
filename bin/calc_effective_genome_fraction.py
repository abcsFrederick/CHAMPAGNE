#!/usr/bin/env python3

"""
Usage:
    calc_effective_genome_fraction.py <chrom_sizes_file> <effective_genome_size> <excluded_chroms>

Example:
    calc_effective_genome_fraction.py hg38.fa.sizes 2700000000 chrM chrX chrY
"""

import sys


def calc_effective_genome_fraction(
    chrom_sizes_filename, effective_genome_size, excluded_chroms=""
):
    with open(chrom_sizes_filename, "r") as infile:
        # creates dictionary with { chromosome: length }
        chrom_lengths = {line.split()[0]: int(line.split()[1]) for line in infile}
    chrom_len_sum = sum(
        chrom_lengths[chrom]
        for chrom in chrom_lengths
        if "_" not in chrom and chrom not in excluded_chroms
    )
    return round(chrom_len_sum / effective_genome_size, 3)


def main(args):
    chrom_sizes_filename = args[1]
    effective_genome_size = int(args[2])
    excluded_chroms = set(args[3:]) if len(args) > 3 else set()
    print(
        calc_effective_genome_fraction(
            chrom_sizes_filename, effective_genome_size, excluded_chroms
        )
    )


if __name__ == "__main__":
    main(sys.argv)
