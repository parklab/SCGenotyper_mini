#!/usr/bin/env python

"""
Top-level executible for SCGentyper_mini

Usage: ./run_mini.py -h

Author: Maxwell A. Sherman
contact: maxwell_sherman@hms.harvard.edu
Last revision: 18 June 2018
Version: 0.0.1
"""
import argparse
import pysam
import os
import pandas as pd
import numpy as np
import sys

from scGTmini import scGTmini

def gt_formatter(gt, nref, nalt, hap, filt, bal, p, ll):
    if gt == 'RR':
        GT = "0|0"
    elif gt == 'AA':
        GT = '1|1'
    elif gt == 'RB':
        GT = '0/.'
    elif gt == 'Ukn':
        GT = "./."
    else:
        if hap == 'PAT':
            GT = "1|0"
        else:
            GT = "0|1"

    if np.isnan(bal):
        bal = "."
        logp = "."
        ll_int = ["." for i in ll]
    else:
        bal = round(bal, 2)
        logp = int(round(-10*np.log10(p)))
        ll_int = np.array([int(round(-10*l)) for l in ll])
        ll_int = ll_int - ll_int.min()

    fmt = "{}:{},{}:{}:{}:{}:{},{},{},{}".format(GT, nref, nalt, filt, bal, logp, *ll_int) 

    return fmt

def vcf_writer(d, keys, bam_names, fout=None):
    vcf = ''
    for line in open("vcf_header.txt", 'r'):
        vcf += line

    vcf += "##python_command=" + " ".join(sys.argv) + "\n"
    vcf += "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + bam_names) + "\n"
    fmt = "GT:AD:FT:BAL:PRB:PL"

    for key in keys:
        l_fmt = [gt_formatter(*l_gt) for l_gt in d[key]]
        l_flt = ["{}:{}".format(l_gt[0], l_gt[4]) for l_gt in d[key]]

        # Set FILTER:
        if np.all(["NOSNP" in fmt for fmt in l_fmt]):
            FILTER = "NOSNP"
        elif ('RA:PASS' in l_flt) or ('AA:PASS' in l_flt):
            FILTER = "PASS"
        # elif np.all(["PASS" in fmt for fmt in l_fmt]):
        #     FILTER = "PASS"
        elif np.all(["LOWQUAL" in fmt for fmt in l_fmt]):
            FILTER = "LOWQUAL"
        elif np.all(["NODP" in fmt for fmt in l_fmt]):
            FILTER = "NODP"
        elif np.all(["ARTIF" in fmt for fmt in l_fmt]):
            FILTER = "ARTIF"
        else:
            FILTER = '.'
        vcf += "\t".join([key, ".", FILTER, ".", fmt] + l_fmt) + "\n"

    if fout:
        with open(fout, 'w') as f:
            f.write(vcf)
    else:
        print(vcf)

def mk_bam_names(bam_list):
    bam_names = []
    for i, bam in enumerate(bam_list):
        s = bam.split('/')
        idx = len(s) - 1
        s_tmp = s[idx]
        while s_tmp in bam_names:
            idx -= 1
            if idx < 0:
                print("WARNING: {} is a duplicate bam file".format(bam))
                s_tmp = s_tmp + "_{}".format(i)

            else:
                s_tmp = s[idx] + '/' + s_tmp

        bam_names.append(s_tmp)

    return bam_names

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-b', '--bam', nargs='+', default=[], type=str, help='bam(s) to genotype')
    parser.add_argument('-B', '--fbam', type=str, help='file containing paths to bam files to analyze. Overrides -b.')
    parser.add_argument('-r', '--pos', required=True, type=str, help='tab-separated file of positions to genotype with columns CHR, POS, REF, ALT')
    parser.add_argument('-v', '--vcf', required=True, type=str, help='vcf file of HAPLOTYPE-PHASED germline SNPs')
    parser.add_argument('-o', '--fout', type=str, help='output file (vcf format)')

    parser.add_argument('-l', '--window', default=10000, type=int, help='window size for kernel')
    parser.add_argument('--p-art', default=0.125,type=float, help='probability of amplification artifact')
    parser.add_argument('--p-err', default=0.003, type=float, help='probability of sequencer artifact')
    parser.add_argument('--DP-min', default=0, type=int, help='Minimum depth to consider a training site')
    parser.add_argument('--mapq', default=5, type=int, help='Minimum mapq to consider a read')

    args = parser.parse_args()

    if not args.bam and not args.fbam:
        message = "the following arguments are required: -b/--bam OR -B/--fbam"
        parser.error(message)

    return args

def run():
    args = parse_args()

    l_bams = args.bam
    if args.fbam:
        l_bams = [line.strip() for line in open(args.fbam, 'r')]

    bam_names = mk_bam_names(l_bams)
    # bam_names = [bam.split('/')[-1].split('.bam')[0] for bam in l_bams]

    dtype = {'CHR': str, 'POS': int, 'REF': str, 'ALT': str}
    pos = pd.read_table(args.pos, names=['CHR', 'POS', 'REF', 'ALT'], dtype=dtype)
    keys = ["{}\t{}\t.\t{}\t{}".format(row.CHR, row.POS, row.REF, row.ALT) for i, row in pos.iterrows()]
    d = {key: [] for key in keys}

    l_gt = [scGTmini(b, args.vcf) for b in l_bams]
    for i, row in pos.iterrows():
        print("chr{}:{}".format(row.CHR, row.POS))
        key = "{}\t{}\t.\t{}\t{}".format(row.CHR, row.POS, row.REF, row.ALT) 
        for scgt in l_gt:
            gt = scgt.caller(row.CHR, row.POS, row.REF, row.ALT, l=args.window, p_art=args.p_art, p_err=args.p_err, DP_min=args.DP_min, drop_pos=True, mapq=args.mapq)
            d[key].append(gt)

    vcf_writer(d, keys, bam_names, fout=args.fout)

if __name__ == "__main__":
    run()
