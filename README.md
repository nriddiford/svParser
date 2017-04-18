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
  --output = write out variants that pass filters to specified dir
  --chromosome = used in conjunction with --dump will cycle though variants on chromosome speciified in -c
  --help

Examples:
o Browse all variants that passed filter within a speicifc window on X chromosome:
->  perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -d -c X:3000000-3500000
o Filter vars and write to file in cwd:
->  perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -o .

Nick Riddiford 2017
```

# Examples 


## Summarise variants

* Read vcf file from lumpy or delly and see summary of variants called: 

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf`

* See summary of variants that passed filter: 

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f`


## Get variant (-i)

* Investigate a specific variant (by ID):

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -i 4706`


## Browse variants (-d)

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -d`

* Browse all variants that passed filter on X chromosome (-c -f): 
 
`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -d -c X`

* Browse all variants within a speicifc window on X chromosome: 

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -d -c X:3000000-3500000`

* Browse all variants that passed filter within a speicifc window on X chromosome: 

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -d -c X:3000000-3500000`


## Print variants (-p)

* Print all variants that passed filter (-o)

`perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -o .`


# To do
- [ ] User-controlled filter params 
- [ ] Custom sorting of vcf summary (in `summarise_variants`), e.g. by SQ value
- [ ] Write output for specified params (rather than for post-filter)
- [ ] Specify file listing chromosomes to filter for

