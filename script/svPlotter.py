import sys
import fnmatch
import os
from optparse import OptionParser
import pandas as pd
import ntpath


def select_cnvfile(sv_size, bp1, bp2, options):
    big_window = '/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data'
    small_window = '/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data/w500'

    sample = ntpath.basename(options.SV_calls).split("_")[0]

    # base_name = (os.path.splitext(options.SV_calls)[0])
    # out_base = base_name.split('_')[0]
    pattern = sample + '*.cnv'

    if float(sv_size) <= 1000:
        start = bp1 - 100000
        if start < 0:
            start = 0
        stop = bp2 + 100000
        tick = 100000

        files = os.listdir(small_window)
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                name = "data/w500/" + name
                return [name, start, stop, tick]
    else:
        start = bp1 - 1000000
        if start < 0:
            start = 0
        stop = bp2 + 1000000
        tick = 1000000

        files = os.listdir(big_window)
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                name = "data/" + name
                return [name, start, stop, tick]


def print_R_command(options):
    """Generates R command that can be used to plot CNVs with https://github.com/nriddiford/cnvPlotteR"""
    df = pd.read_csv(options.SV_calls, delimiter="\t")
    for idx, row in df.iterrows():
        event, source, type, chrom1, bp1, chrom2, bp2, genotype, size, true = row[['event', 'source', 'type', 'chromosome1', 'bp1', 'chromosome2', 'bp2', 'genotype', 'length(Kb)', 'T/F']]
        if genotype != 'somatic_tumour': continue
        if options.true_only and true == 'F': continue

        if chrom1 == chrom2:
            cnv_file, start, end, tick = select_cnvfile(size, bp1, bp2, options)
            print("SV event: %s, type: %s, size: %s") % (event, type, size)
            print("regionPlot(cnv_file=\"%s\", from=%s, to=%s, bp1=%s,bp2=%s,chrom=\"%s\", tick=%s, title=\"%sKb %s on %s\")") % (cnv_file, start, end, bp1, bp2, chrom1, tick, size, type, chrom1)

def main():
    parser = OptionParser()

    parser.add_option("-f", "--SV_calls", dest="SV_calls", help="SV calls file")
    parser.add_option("-t", "--true_only", dest="true_only", action="store_true", help="Don't print calls marked as FP")

    options, args = parser.parse_args()

    if options.SV_calls is None:
        parser.print_help()
        print()
    else:
        try:
            print_R_command(options)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())
