# SV_parser

Parse VCF output from structural variant callers Lumpy and Delly.

Run without options or with `--help` or `-h` to print usage statement

```
********** SV_parser ***********
Usage: script/sv_parse.1.0.pl [options]
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