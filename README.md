# SV_parser

Parse VCF output from structural variant callers Lumpy and Delly.

Run without options or with `--help` or `-h` to print usage statement

```
********** sv_parse.1.0.pl ***********
Usage: sv_parse.1.0.pl [options]
  --vcf = VCF file for parsing
  --id = extract information for a given variant
  --dump = cycle through all variants (can be combined with both -f and -c)
  --filter = apply filters and mark filtered variants
  --print = write out variants that pass filters
  --chromosome = used in conjunction with --dump will cycle though variants on chromosome speciified in -c
  --help

Examples:
  Print to screen all vars on chromosome '2L': perl script/sv_parse.1.0.pl -v [file.vcf] -d -c 2L
  Print to screen all vars on chromosome '2L' within range 100000-200000: perl script/sv_parse.1.0.pl -v [file.vcf] -d -c 2L:100000-200000
  Print to screen all vars that pass the fitleres on chromosome '2L' within range 100000-200000: perl script/sv_parse.1.0.pl -v [file.vcf] -f -d -c 2L:100000-200000
  Filter vars and write to file: perl script/sv_parse.1.0.pl -v [file.vcf] -f -p

Nick Riddiford 2017
```


# Examples 


## Summarise variants

* Read vcf file from lumpy or delly and see summary of variants called: 

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf`


## Get variant (-i)

#### Investigate a specific variant (by ID):

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -i 4706`

#### Investigate a specific variant (by ID) and see if it's filtered:

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -i 4706 -f`


## Browse variants (-d)

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -d`

#### Browse all variants that passed filter on X chromosome (-c -f): 
 
`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -d -c X`

#### Browse all variants within a speicifc window on X chromosome: 

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -d -c X:20000-100000`

#### Browse all variants that passed filter within a speicifc window on X chromosome: 

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -d -f -c X:20000-100000`


## Print variants (-p)

#### Print all variants that passed filter (-p)

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -p`