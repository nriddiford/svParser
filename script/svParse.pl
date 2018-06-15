#!/usr/bin/perl

use strict;
use warnings;

use FindBin qw/ $Bin /;
use FindBin qw/ $Script /;
use File::Spec;
use lib File::Spec->catdir($FindBin::Bin, '..', 'bin/');

use svParser;

use feature qw/ say /;
use Data::Dumper;
use File::Basename;
use File::Path qw/ make_path / ;
use File::Slurp;

use Getopt::Long qw/ GetOptions /;

my $vcf_file;
my $help;
my $id;
my $dump;
my $chromosome;
my $method = "guess";
my $print;
my $exclude;
my %filters;
my $PON_print = 5;
my %opt;
my $true_positives;

GetOptions( 'vcf=s'             =>    \$vcf_file,
            'method=s'          =>    \$method,
            'id=s'              =>    \$id,
            'dump'              =>    \$dump,
            'filter:s'          =>    \%filters,
            'print'             =>    \$print,
            'normalPrint=s'     =>    \$PON_print,
            'chromosome=s'      =>    \$chromosome,
            'exclude=s'         =>    \$exclude,
            'true_calls=s'      =>    \$true_positives,
            'help'              =>    \$help
) or die usage();

if ($help) { exit usage() }

if (not $vcf_file) {
   exit usage();
}

my @keys;
if (-e 'chroms.txt'){
  say "\nReading chromosomes file from 'chroms.txt'";
  @keys = read_file('chroms.txt', chomp=>1);
}
else {
  say "\nFiltering for Drosophila chroms: 2L 2R 3L 3R 4 X Y ";
  @keys = qw / 2L 2R 3L 3R 4 X Y /;
}

my @exclude_regions;
if ($exclude){
  my @exclude = read_file($exclude, chomp=>1);
  my %chroms;
  @chroms{@keys} = ();
  # Should speed things up a bit by removing any chromosomes in exclude (.bed) file that aren't in our list of chroms we want to look at (@keys)
  foreach(@exclude){
    my $chrom = (split)[0];
    next unless exists $chroms{$chrom};
    push @exclude_regions, $_;
  }
  undef @exclude;
  $filters{'e'} = 1;
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

my $filter_switch= 0;

if ( scalar keys %filters > 0 ){
  print "\n";
  if ( exists $filters{'a'} ){
    $filter_switch = 1;
    my $filter_ref = allFilters($exclude, \@keys, \%filters);
  }
  elsif ( $filters{'su'} or $filters{'dp'} or $filters{'rdr'} or $filters{'sq'} or $filters{'chr'} or $filters{'g'} or $filters{'e'} or $filters{'n'} or $filters{'s'} ) {
    say "Running in filter mode, using custom filters:";
    say " o Read support >= $filters{'su'}" if $filters{'su'};
    say " o Read depth (in both tumour and normal) > $filters{'dp'}" if $filters{'dp'};
    say " o Read support / depth > $filters{'rdr'}" if $filters{'rdr'};
    say " o SQ quality > $filters{'sq'}" if $filters{'sq'};
    say " o Chromosomes: " . join(' ', @keys) if $filters{'chr'};
    say " o Running in germline mode" if $filters{'g'};
    say " o Running in somatic TUMOUR mode" if $filters{'s'};
    say " o Running in somatic NORMAL mode" if $filters{'n'};
    say " o Excluding calls overlapping: $exclude" if $filters{'e'};
    $filter_switch = 1;
  }
  else {
    my $illegals = join(",", keys %filters);
    say "Illegal filter option used: '$illegals'. Please specify filters to run with (or use '-f or -f a' to run all defaults)";
    say "Filter options available:";
    say " o Read support: su=INT";
    say " o Read depth: dp=INT";
    say " o Read support / depth: rdr=FLOAT";
    say " o SQ quality: sq=INT";
    say " o Chromosomes: " . join(' ', @keys);
    say " o Germline only: g=1";
    say " o Exclude calls in bed file: e [bed file]";
    die "Please check filter specification\n";
     }
}

my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

print "\n";

# Retun SV and info hashes
my ( $SVs, $info, $filtered_vars, $call_lookup) = svParser::typer( $vcf_file, $method, \@exclude_regions, \@keys, \%filters );

testCalls($true_positives, $SVs, $info, $filter_switch, $PON_print, $call_lookup) if $true_positives;

if ($method ne 'snp') {
  # Print summary to screen
  svParser::summarise_variants( $SVs, $filter_switch, $chromosome ) unless $id or $dump or $true_positives;

  # Print all info for specified id
  svParser::get_variant( $id, $SVs, $info, $filter_switch, $PON_print) if $id and not $true_positives;

  # Dump all variants to screen
  svParser::dump_variants( $SVs, $info, $filter_switch, $chromosome, $method, $PON_print ) if $dump and not $true_positives;

  # Write out variants passing filters
  # Write out some useful info to txt file
  if ( $filters{'g'} ){
    print "Printing germline events\n";
    svParser::print_variants( $SVs, $filtered_vars, $name, $filtered_out, 1 ) if $print;
    svParser::write_summary( $SVs, $name, $summary_out, $method, 1) if $print;
  }
  else {
    svParser::print_variants( $SVs, $filtered_vars, $name, $filtered_out, 0 ) if $print;
    svParser::write_summary( $SVs, $name, $summary_out, $method, 0 ) if $print;
  }
}

if ($method eq 'snp') {
  # Dump all variants to screen
  svParser::dump_variants( $SVs, $info, $filter_switch, $chromosome, $method) if $dump;

  if ($print or $id ){
    die "Print and get variants not supported for SNP data\n";
  }
}


sub testCalls {
  my ($truePositves, $SVs, $info, $filter_switch, $PON_print, $call_lookup) = @_;
  open my $tps, '<', $true_positives or die $!;
  my %truth;
  my %calls = %{ $call_lookup };
  my $found = 0;
  my $not_found = 0;
  my $total = 0;
  while(<$tps>){
    chomp;
    my ($sample, $caller, $type, $c1, $bp1, $c2, $bp2) = split;
    my $key = $c1 . "_" . $bp1 . "_" . $c2 . "_" . $bp2;
    $truth{$key} = 1;

    if (exists $calls{$key}){
      say "o True positve found in filtered data id: $calls{$key}";
      $found++;
    }
    else{
      say "!! True positve NOT found in filtered data!";
      say "  - $_";
      $not_found++;
    }
    $total++;
  }

  print "$found/$total true positves retained in filtered data\n";
}


sub allFilters {
  my ($exclude_file, $chroms, $f) = @_;
  my %filters = %{$f};
  say "Running in filter mode, using all default filters:";
  say " o Read support >= 4";
  say " o Read depth (in both tumour and normal) > 10";
  say " o Read support / depth > 0.1";
  say " o SQ quality > 10";
  say " o Chromosomes: " . join(' ', @{$chroms} );
  say " o Excldung calls in regions: $exclude_file";

  $filters{'su'} = 4;
  return(\%filters)
}


sub usage {
  print
"
usage: perl $Script -v variant_file.vcf [options]

svParser
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: $VERSION
description: Browse vcf output from several SV callers LUMPY, DELLY and novobreak

arguments:
  -h, --help                show this help message and exit
  -v file, --vcf            VCF input    [required]
  -p, --print               print filtered vcf and summary file to './filtered'
  -m str, --method            specify input source [default: guess from input]
                                -l = LUMPY
                                -d = DELLY
                                -n = novobreak
                                -snp = snp/snv

  -i str, --id              breakpoint id to inspect
  -d, --dump                print each variant called to screen, press any key to advance, q to exit

  -n, --normalPrint         maximum number of PON members to print to screen when using -d or -i [default: 5]
  -c str, --chromosome      limit search to chromosome and/or region (e.g. X:10000-20000)
                            can be used in conjunction with -d

  -e file, --exclude        path to .bed file containing regions to exlcude

  -f key=val, --filter
                            filters to apply:
                            -f su=INT    [number of tumour reads supporting var]
                            -f dp=INT    [minimum depth for both tumour normal at variant site]
                            -f rdr=FLOAT [fraction of supporting reads/tumour depth - a value of 1 would mean all reads support variant]
                            -f sq=INT    [phred-scaled variant likelihood]
                            -f chr=1     [only show chromosomes in 'chroms.txt'. [Default use Drosophila chroms: 2L 2R 3L 3R 4 X Y]
                            -f s=1       [only keep somatic tumour events]
                            -f, -f a     [apply default filters: -f su=4 -f dp=10 -f rdr=0.1 -f sq=10 -f chr=1 ]
"
}
