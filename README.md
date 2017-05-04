# SV_parser

Parse VCF output from structural variant callers Lumpy and Delly (also supports input from novobreak).

Run without options or with `--help` or `-h` to print usage statement

```
svParser
author: Nick Riddiford (nick.riddiford@curie.fr)
version: v1.0
description: Browse vcf output from several SV callers LUMPY, DELLY and novobreak

arguments:
  -h, --help            show this help message and exit
  -v FILE, --vcf FILE
                        VCF input [required]
  -o PATH, --output PATH
                        path to write filtered file to
  -t STRING, --type STRING
                        specify input source [default: guess from input]
                        -l = LUMPY
                        -d = DELLY
                        -n = novobreak
  -i STRING, --id STRING
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

# Examples


## Summarise variants

* Read vcf file from lumpy or delly (or novobreak) and see summary of variants called:

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -t l`

* See summary of variants that passed default filters (only appropriate for Drosophila):

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -t l -f a`

* See summary of variants that have >= 4 reads in tumour (-f 4)

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -t l -f su=4`


## Get variant (-i)

* Investigate a specific variant (by ID):

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -i 4706`


## Browse variants (-d)

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -d`

* Browse all variants with su>=5 on X chromosome (-c -f su=5):

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f su=5 -d -c X`

* Browse all variants within a specific window on X chromosome:

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -d -c X:3000000-3500000`

* Browse all variants that passed read depth filter (dp=20) filter within a specific window on X chromosome:

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f dp=20 -d -c X:3000000-3500000`


## Print variants (-p)

* Write all variants that passed filter (-o) to 'HUM-7.tagged.SC.lumpy.gt_all.filtered.vcf'

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -o .`


# To do
- [x] User-controlled filter params
- [x] Write output for specified params (rather than for post-filter)
- [ ] Custom sorting of vcf summary (in `summarise_variants`), e.g. by SQ value
