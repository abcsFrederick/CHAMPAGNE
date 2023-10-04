#!/usr/bin/env python

from __future__ import print_function
import sys


def formatSequencelength(seq, stringlen):
    fseq = ""
    for i in range(len(seq)):
        index = i + 1
        if index % 80 == 0:
            fseq += "{}{}".format(seq[i], "\n")
        else:
            fseq += seq[i]
    return fseq


def parsed(filename):
    fh = open(filename, "r")
    sequence = ""
    chrom = ""
    seqindex = 0
    seqlen = 0
    for line in fh:
        line = line.strip()
        if line.startswith(">") and sequence != "":
            yield chrom, formatSequencelength(sequence, seqlen), len(sequence)
            chrom = line.split(" ")[0]
            sequence = ""
        elif line.startswith(">"):
            chrom = line.split(" ")[0]
        else:
            seqindex += 1
            sequence += line
            if seqindex == 1:
                seqlen = len(line)
    else:
        # formatSequencelength(sequence, seqlen)
        yield chrom, formatSequencelength(sequence, seqlen), len(sequence)
    fh.close()


def main():
    file = sys.argv[1]
    chromsizesfh = open("chrom.sizes", "w")

    for chrom, seq, chromsize in parsed(file):
        chromsizesfh.write("{}\t{}\n".format(chrom.replace(">", ""), chromsize))
        outfilename = chrom.replace(">", "") + ".fa"
        outfh = open(outfilename, "w")
        print("{}\n".format(chrom))
        outfh.write("{}\n{}\n".format(chrom, seq.rstrip()))
        outfh.close()

    chromsizesfh.close()


if __name__ == "__main__":
    main()
