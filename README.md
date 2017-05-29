# svParser

Parse VCF output from structural variant callers LUMPY and DELLY (also supports input from novoBreak).

Run without options or with `--help` or `-h` to print usage statement

```
usage: sv_parse.1.0.pl [-h] [-v FILE] [-p] [-t STR] [-i STR] [-d] [-f key=val] [-c STR]

svParser
author: Nick Riddiford (nick.riddiford@curie.fr)
version: v1.0
description: Browse vcf output from several SV callers LUMPY, DELLY and novobreak

arguments:
  -h, --help            show this help message and exit
  -v FILE, --vcf
                        VCF input [required]
  -p, --print
                        print filtered vcf and summary file to './filtered'
  -t STRING, --type
                        specify input source [default: guess from input]
                        -l = LUMPY
                        -d = DELLY
                        -n = novobreak
  -i STRING, --id
                        breakpoint id to inspect
  -d, --dump            cycle through breakpoints
  -c STRING, --chromosome
                        limit search to chromosome and/or region (e.g. X:10000-20000)
                        can be used in conjunction with -d
  -f KEY=VAL, --filter
                        filters to apply:
                        -f su=INT [number of tumour reads supporting var]
                        -f dp=INT [minimum depth for both tumour normal at variant site]
                        -f rdr=FLOAT [supporting reads/tumour depth - a value of 1 would mean all reads support variant]
                        -f sq=INT [phred-scaled variant likelihood]
                        -f chr=1 [only show chromosomes 2L 2R 3L 3R 4 X Y. Only use for Drosophila]
                        -f, -f a = apply default filters [ -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 ]
```

# Parsing structural variants from VCF files called by LUMPY DELLY and novoBreak

## Summarise variants

Print a summary of variants called in VCF file to see the number of DELS/DUPS/INV/TRA called  

#### Examples

* Read vcf file from lumpy (or delly / novobreak) and see summary of variants called:

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -t l`

* See summary of variants that passed default filters (only appropriate for Drosophila):

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -t l -f a`

* See summary of variants that have >= 4 reads supporting call in tumour (-f su=4)

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -t l -f su=4`


## Get variant (-i)

* Investigate a specific variant (by ID):

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -i 4706`


## Browse variants (-d)

Browse variants using the `-d` flag. This cycles through each line of the VCF file, printing an easy-to-read breakdown for each variant.  
Press any key to move to the next variant, or `q` to quit

#### Examples

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -d`

* Browse all variants with su>=5 on X chromosome (-c X -f su=5):

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -f su=5 -d -c X`

* Browse all variants within a specific window on X chromosome (-c X:3000000-3500000):

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -d -c X:3000000-3500000`

* Browse all variants that passed read depth filter (dp=20) filter within a specific window on X chromosome:

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -f dp=20 -d -c X:3000000-3500000`


## Print variants (-p)

Write a new VCF file (`filtered/input.filtered.vcf`) and a summary txt file containing useful info (`filtered/summary/input.filtered.summary.txt`).  
To be used in conjunction with filter options.

* Write all variants that passed filter (-p) to 'HUM-7.tagged.SC.lumpy.gt_all.filtered.vcf'

`perl script/sv_parse.1.0.pl -v data/lumpy/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -p`


# Combining calls from multiple sources

* Generate `filtered/input.filtered.vcf` and `filtered/summary/input.filtered.summary.txt` files for calls made by delly, lumpy and novoBreak

```
for file in data/lumpy/*.vcf
do
  perl script/sv_parse.1.0.pl -v $file -f a -t l -p
done
```

```
for file in data/delly/*.vcf
do
  perl script/sv_parse.1.0.pl -v $file -f a -t d -p
done
```

```
for file in data/novobreak/*.vcf
do
  perl script/sv_parse.1.0.pl -v $file -f a -t n -p
done
```

* Run [merge_vcf](https://github.com/ljdursi/mergevcf/tree/master/mergevcf) on all filtered files. (Run from `filtered` directory):

```
cd filtered
perl ../script/merge_vcf.pl
```

* Run svMerger.pl for all files in sample group. E.g. `input1.delly.filtered.summary.txt` `input1.lumpy.filtered.summary.txt` `input1.novobreak.filtered.summary.txt`

```
cd filtered/summary
perl ../../script/svMerger.pl -f input1.*
```


# To do
- [x] User-controlled filter params
- [x] Write output for specified params (rather than for post-filter)
- [ ] Custom sorting of vcf summary (in `summarise_variants`), e.g. by SQ value
