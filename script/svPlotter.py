import sys
import fnmatch
import os
from optparse import OptionParser
import pandas as pd
import ntpath
from collections import defaultdict


def select_cnvfile(sv_size, bp1, bp2, options):
    big_window = options.big_window
    small_window = options.small_window


    sample = ntpath.basename(options.variants).split("_")[0]
    pattern = sample + '*.cnv'
    dir = "data/"
    file_name = None

    if bp1 is None:
        bp1 = 0
    if float(sv_size) <= 300:
        start = bp1 - 50000
        if start < 0:
            start = 0
        stop = bp2 + 50000
        tick = 50000

        files = os.listdir(small_window)
        dir = os.path.join(dir, "w500")

    else:
        start = bp1 - 1000000
        if start < 0:
            start = 0
        stop = bp2 + 1000000
        tick = 1000000

        files = os.listdir(big_window)

    for f in files:
        # if f.endswith(".cnv") and fnmatch.fnmatch(f, pattern):
        if f.endswith(".cnv") and sample == f.split(".")[0]:
            file_name = os.path.join(dir, f)

    return file_name, start, stop, tick


def print_R_command(options):
    """Generates R command that can be used to plot CNVs with https://github.com/nriddiford/cnvPlotteR"""
    df = pd.read_csv(options.variants, delimiter="\t", index_col=False)
    for idx, row in df.iterrows():
        event, source, svtype, chrom1, bp1, chrom2, bp2, genotype, size, true = row[['event', 'source', 'type', 'chromosome1', 'bp1', 'chromosome2', 'bp2', 'genotype', 'length(Kb)', 'status']]

        if genotype != 'somatic_tumour': continue
        if options.true_only and true == 'F': continue

        # if chrom1 == chrom2:
        #     print(select_cnvfile(size, bp1, bp2, options), event)

        if chrom1 == chrom2:
            cnv_file, start, end, tick = select_cnvfile(size, bp1, bp2, options)
            if cnv_file:
                print("SV event: %s, type: %s, size: %s") % (event, svtype, size)
                print("regionPlot(cnv_file=\"%s\", from=%s, to=%s, bp1=%s,bp2=%s,chrom=\"%s\", tick=%s, title=\"%sKb %s on %s\")") % (cnv_file, start, end, bp1, bp2, chrom1, tick, size, svtype, chrom1)
            else:
                print("Can't find CNV file")


def write_notebook(options):
    """Write R commands to an R markdown notebook"""
    sample = ntpath.basename(options.variants).split('_')[0]
    nb_out = sample + '_variants.Rmd'

    print("Writing variants to file %s" % (nb_out))

    with open(nb_out, 'w') as nb:
        df = pd.read_csv(options.variants, delimiter="\t", index_col=False)
        df = df[df['status']!='F']
        var_count = len(df['event'].unique())

        nb.write(write_header(sample))
        nb.write(write_summary(sample, var_count))

        seen_event = defaultdict(lambda: defaultdict(int))

        for idx, row in df.iterrows():
            event, source, type, chrom1, bp1, chrom2, bp2, genotype, size, true = row[['event', 'source', 'type', 'chromosome1', 'bp1', 'chromosome2', 'bp2', 'genotype', 'length(Kb)', 'status']]

            if genotype != 'somatic_tumour': continue

            if seen_event[event][type] >= 1:
                continue
            seen_event[event][type] += 1
            nb.write(write_title(row))

            if chrom1 == chrom2:
                cnv_file, start, end, tick = select_cnvfile(size, bp1, bp2, options)
                cnv_file = os.path.join("/Users/Nick_curie/Desktop/script_test/cnvPlotteR", cnv_file)
                nb.write("```{r}\n")
                r_cmd = ("lightPlot(cnv_file=\"%s\", position = \'%s:%s-%s\', title=\"%sKb %s on %s\")") % (cnv_file, chrom1, bp1, bp2, size, type, chrom1)
                nb.write(r_cmd + "\n")
                nb.write("```\n")



def write_summary(sample, n):
    summary = """
# Summary

* {0} events in sample {1}

***
""".format(n, sample)
    return summary


def write_title(r):
    t1 = """
# Event {0} - {1} Kb {2} on {3}
""".format(r['event'], r['length(Kb)'], r['type'], r['chromosome1'])
    gene_count = len(r['affected_genes'].split(','))

    t2 = """
* Allele frequencty: {0}
* Bp1 locus: {1}
* Bp2 locus: {2}
* {3} genes affected
""".format(r['allele_frequency'], r['bp1_locus'], r['bp2_locus'], gene_count)

    if not (pd.isna(r['notes'])):
        t3 = """
* {0}
""".format(r['notes'])
        line = t1 + t2 + t3
    else:
        line = t1 + t2
    return line


def write_header(sample):
    title = """
---
title: "Variant calls for {0}"
""".format(sample)
    header = """

date: "`r format(Sys.time(), '%B %e, %Y')`"
output:
  html_notebook:
    theme: flatly
    toc: yes
    toc_float: yes
    number_sections: true
    code_folding: hide
    df_print: paged
  html_document:
    theme: flatly
    toc: yes
    toc_float: yes
    number_sections: true
    code_folding: hide
    df_print: paged
---

<style>
 body .main-container {
    max-width: 1600px !important;
    margin-left: auto;
    margin-right: auto;
  }
</style>


```{r setup, include=FALSE, cache=TRUE}
source("~/Desktop/Notebooks/R/utils.R")
library(cnvPlotteR)
```
"""
    return title + header

def write_table(line):
    l1 = """
```{r}
"""
    l2 = """
inTable({0})
```
""".format(line)
    return l1 + l2


def main():
    parser = OptionParser()

    parser.add_option("-f", "--variants", dest="variants", help="SV calls file")
    parser.add_option("-t", "--true_only", dest="true_only", action="store_true", help="Don't print calls marked as FP")
    parser.add_option("--cnv_large", dest="big_window", action="store", help="CN file with large window")
    parser.add_option("--cnv_small", dest="small_window", action="store", help="CN file with small window")

    parser.add_option("-n", "--notebook", dest="notebook", action="store_true", help="Write out a notebook")
    parser.set_defaults(big_window='/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data',
                        small_window='/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data/w500')

    options, args = parser.parse_args()

    if options.variants is None:
        parser.print_help()
    else:
        try:
            if options.notebook:
                write_notebook(options)
            else:
                print_R_command(options)
        except IOError as err:
            sys.stderr.write("IOError " + str(err) + "\n")
            return


if __name__ == "__main__":
    sys.exit(main())
