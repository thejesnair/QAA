#!/usr/bin/env python

import bioinfo
import gzip
import argparse
import matplotlib.pyplot as plt

def get_args():
    '''Takes inputs for paired reads (adapter and quality trimmed) and outputs a plot of qual distr.'''
    parser = argparse.ArgumentParser(description="Plot trimmed (adapter and quality) paired reads on same plot.")
    parser.add_argument("-R1", "--read1", help="read1 filename", type=str)
    parser.add_argument("-R2", "--read2", help="read2 filename", type=str)
    parser.add_argument("-o", "--outputpng", help="name for output png", type=str)
    parser.add_argument("-b", "--binsize", type=int, default=5, help="Bin size for histogram (default: 5)")
    parser.add_argument("-t", "--title", help="title on plot", type=str)


    return parser.parse_args()

args = get_args()
R1 = args.read1
R2 = args.read2
outputpng = args.outputpng
title = args.title
binsize = args.binsize

def read_lengths(file: str) -> list:
    """Return list of read lengthss from fastq file"""
    length = []
    with gzip.open(file, "rt", encoding="utf-8") as f:
    #with open(file, "r") as f:
        for i, line in enumerate(f):
            if i % 4 == 1: #seq line, bc seq second line
                line = line.strip()
                length.append(len(line.strip()))
    return length


R1_distr = read_lengths(R1)
R2_distr = read_lengths(R2)

def get_bincounts(read_lengths, binsize):
    counts = {}
    for read_len in read_lengths:
        bins = (read_len // binsize) * binsize
        counts[bins] = counts.get(bins, 0) +1
    return counts

R1_counts = get_bincounts(R1_distr, binsize)
R2_counts = get_bincounts(R2_distr, binsize)



def plot_hist(R1_counts: dict, R2_counts: dict, outputpng: str, title:str):
    """plot histogram of R1 and R2 read distr lens"""
    all_bins = sorted(R1_counts.keys() | R2_counts.keys()) #can also use set()
    R1_vals = [R1_counts.get(bins, 0) for bins in all_bins]
    R2_vals = [R2_counts.get(bins, 0) for bins in all_bins]

    x = range(len(all_bins))
    width = 0.4

    plt.bar([i - width/2 for i in x], R1_vals, width=width, label="R1")
    plt.bar([i + width/2 for i in x], R2_vals, width=width, label="R2")

    plt.xticks(x, all_bins, rotation=45)
    plt.xlabel("Read length (bp)")
    plt.ylabel("Read count")
    plt.title(title)
    plt.legend()
    plt.tight_layout()
    plt.savefig(outputpng)
    plt.close()

plot_hist(R1_counts, R2_counts, outputpng, title)