#!/usr/bin/env python

import sys
import fnmatch
import os
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f",
                  "--SV_calls",
                  dest="SV_calls",
                  help="SV calls file"
                  )

options, args = parser.parse_args()


if not options.SV_calls:
    parser.error('No input file provided')

variants = options.SV_calls

base_name = (os.path.splitext(variants)[0])
out_base = base_name.split('_')[0]

pattern = out_base + '*.cnv'

big_window='/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data'
samll_window='/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data/w500'

def select_cnvfile(sv_size, bp1, bp2):
    if float(sv_size) <= 1000:
        start = bp1 - 100000
        if start < 0:
            start = 0
        stop  = bp2 + 100000
        tick = 100000

        files = os.listdir(samll_window)
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                name = "data/w500/" + name
                return[name, start, stop, tick]
    else:
        start = bp1 - 1000000
        if start < 0:
            start = 0
        stop  = bp2 + 1000000
        tick = 1000000

        files = os.listdir(big_window)
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                name = "data/" + name
                return[name, start, stop, tick]

def print_R_command(SV_calls_file):
    """Generates command that can be used to plot CNVs with https://github.com/nriddiford/cnvPlotteR"""
    i=1
    with open(SV_calls_file, 'U') as calls:
        for l in calls:
            parts = l.rstrip().split('\t')

            if i == 1 and parts[0] == 'event':
                continue
            elif parts[10] != 'somatic_tumour' or parts[10] != 'somatic_normal':
                continue

            event, source, type, chrom1, bp1, chrom2, bp2 = parts[0:7]
            bp1 = int(bp1)
            bp2 = int(bp2)
            # size = parts[10] change for new genotype col
            size = parts[11]
            # if chrom1 == chrom2 and type != 'INV':
            if chrom1 == chrom2:

                fields = select_cnvfile(size, bp1, bp2)
                cnv_file = fields[0]
                start = fields[1]
                end = fields[2]
                tick = fields[3]
                print("SV event: %s, type: %s, size: %s") % (event, type, size)
                print("regionPlot(cnv_file=\"%s\", from=%s, to=%s, bp1=%s,bp2=%s,chrom=\"%s\", tick=%s, title=\"%sKb %s on %s\")") % (cnv_file, start, end, bp1, bp2, chrom1, tick, size, type, chrom1)

print_R_command(variants)
