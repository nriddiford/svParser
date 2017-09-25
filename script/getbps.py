#!/usr/bin/env python

# this script reads each reannotated_SVs file and appends coordinates for bps to a new file 'all_bps_cleaned'
# This has to be done in a second pass (i.e. not in sv2gene --reannotate) in case of modification of SV types etc - sv2gene.pl will ignore previously annotated calls

import sys

from optparse import OptionParser
import os

parser = OptionParser()

parser.add_option("-f",
                  "--SV_calls",
                  dest="SV_calls",
                  help="SV calls file"
                  )

options, args = parser.parse_args()

if not options.SV_calls:
    parser.error('No input file provided')

base_name = (os.path.splitext(options.SV_calls)[0])
sample = base_name.split('_')[0]

event = ''
firstline = True
with open(options.SV_calls, 'U') as sv_calls_file, open('all_bps_cleaned.txt', 'a+') as all_bps_out:
    for l in sv_calls_file:
        parts = l.rstrip().split('\t')

        if firstline:
            firstline = False
            continue

        if event == parts[0]:
            continue

        event = parts[0]

        if parts[17] == '-' or parts[18] == '-':
            continue

        if parts[17] == 'intergenic':
            bp1_feature = 'intergenic'
            bp1_gene = 'intergenic'
        else:
            bp1_gene, bp1_feature = parts[17].split(',')

        if parts[18] == 'intergenic':
            bp2_feature = 'intergenic'
            bp2_gene = 'intergenic'

        else:
            bp2_gene, bp2_feature = parts[18].split(',')


        bp2_gene = bp2_gene.strip('" ')
        bp2_feature = bp2_feature.strip('" ')
        bp1_gene = bp1_gene.strip('" ')
        bp1_feature = bp1_feature.strip('" ')


        bp1_line = [event, 'bp1', sample, parts[3], parts[4], bp1_gene, bp1_feature, parts[2], parts[10]]
        bp2_line = [event, 'bp2', sample, parts[5], parts[6], bp2_gene, bp2_feature, parts[2], parts[10]]

        all_bps_out.write('\t'.join(bp1_line) + '\n')
        all_bps_out.write('\t'.join(bp2_line) + '\n')
