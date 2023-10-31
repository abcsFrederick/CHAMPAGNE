#!/usr/bin/env python
import Bio.SeqIO
import gzip
import os


def main():
    count = 0
    for fastq_filename in "${fastq}".split():
        with gzip.open(fastq_filename, "rt") as file_handle:
            n_seqs = sum(1 for rec in Bio.SeqIO.parse(file_handle, "fastq"))
        count += n_seqs
    with open("${meta.id}.count.txt", "w") as out_file:
        out_file.write(str(count))
    return count


if __name__ == "__main__":
    print(main())
