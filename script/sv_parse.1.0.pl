#!/usr/bin/perl

use strict;
use warnings;

use 5.18.2;

use FindBin qw($Bin);
use FindBin '$Script';

use File::Spec;
use lib File::Spec->catdir($FindBin::Bin, '..', 'bin/');

use SV_parser;

use feature qw/ say /;
use Data::Dumper;
use Getopt::Long qw/ GetOptions /;

use File::Basename;
use File::Path qw(make_path);

my $vcf_file;
my $help;
my $id;
my $dump;
my $chromosome;
my $type = "guess";
my $print;

my %filters;

GetOptions( 'vcf=s'           =>    \$vcf_file,
            'type=s'          =>    \$type,
            'id=s'            =>    \$id,
            'dump'            =>    \$dump,
            'filter:s'        =>    \%filters,
            'print'           =>    \$print,
            'chromosome=s'    =>    \$chromosome,
            'help'            =>    \$help
    ) or die usage();

if ($help) { exit usage() }

if (not $vcf_file) {
   exit usage();
}


my ($filtered_out, $summary_out);
if ($print){
  $filtered_out = "$Bin/../filtered/";
  eval { make_path($filtered_out) };

  if ($@) {
    print "Couldn't create '$filtered_out': $@";
  }

  $summary_out = "$Bin/../filtered/summary/";
  eval { make_path($summary_out) };

  if ($@) {
    print "Couldn't create '$summary_out': $@";
  }
}

my $filter = 0;

if ( scalar keys %filters > 0 ){
  print "\n";
  if ( exists $filters{'a'} ){
    say "Running in filter mode, using all default filters:";
    say " o Read support >= 4";
    say " o Read depth (in both tumor and normal) > 10";
    say " o Read support / depth > 0.1";
    say " o SQ quality > 10";
    say " o Chromosomes 2L 2R 3L 3R 4 X Y";

    %filters = ("su"  =>  4,
                "dp"  =>  10,
                "rdr" =>  0.1,
                "sq"  =>  10,
                "chr" =>  1
               );
    $filter = 1;

  }
  elsif ( $filters{'su'} or $filters{'dp'} or $filters{'rdr'} or $filters{'sq'} or $filters{'chr'} ) {
    say "Running in filter mode, using custom filters:";
    say " o Read support >= $filters{'su'}" if $filters{'su'};
    say " o Read depth (in both tumor and normal) > $filters{'dp'}" if $filters{'dp'};
    say " o Read support / depth > $filters{'rdr'}" if $filters{'rdr'};
    say " o SQ quality > $filters{'sq'}" if $filters{'sq'};
    say " o Chromosomes 2L 2R 3L 3R 4 X Y" if $filters{'chr'};
    $filter = 1;
  }
  else {
    my $illegals = join(",", keys %filters);
    say "Illegal filter option used: '$illegals'. Please specify filters to run with (or use '-f or -f a' to run all defaults)";
    say "Filter options available:";
    say " o Read support: su=INT";
    say " o Read depth: dp=INT";
    say " o Read support / depth: rdr=FLOAT";
    say " o SQ quality: sq=INT";
    say " o Chromosomes 2L 2R 3L 3R 4 X Y";
    die "Please check filter specification\n";
     }
}

my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

print "\n";

# Retun SV and info hashes
my ( $SVs, $info, $filtered_vars ) = SV_parser::typer($vcf_file, $type, %filters);

# Print all info for specified id
SV_parser::summarise_variants( $SVs, $filter, $chromosome ) unless $id or $dump;

# Print all info for specified id
SV_parser::get_variant( $id, $SVs, $info, $filter ) if $id;

# Dump all variants to screen
SV_parser::dump_variants( $SVs, $info, $filter, $chromosome ) if $dump;

# Write out variants passing filters
SV_parser::print_variants ( $SVs, $filtered_vars, $name, $filtered_out ) if $print;

# Write out some useful info to txt file
SV_parser::write_summary ( $SVs, $name, $summary_out, $type ) if $print;

sub usage {
  print
"
usage: $Script [-h] [-v FILE] [-p] [-t STR] [-i STR] [-d] [-f key=val] [-c STR]

svParser
author: Nick Riddiford (nick.riddiford\@curie.fr)
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
"
}
