#!/usr/bin/perl
use strict;
use warnings;

use Data::Dumper;
use feature qw/ say /;

use Getopt::Long qw/ GetOptions /;

my $help;

GetOptions( 'help'          =>    \$help
          ) or die usage();

if ($help) { exit usage() }

my %vcfs;
my %types;

opendir my $dir, "." or die "Cannot open directory: $!";
my @files = readdir $dir;
closedir $dir;

foreach (@files){
  chomp;

  unless (/.vcf$/){
    next;
  }

  my @parts = split(/\./, $_);
  my ($id) = $parts[0];

  if (/(lumpy|delly|novoBreak)/){
    push @{$types{$id}}, $1;
  }

  push @{$vcfs{$id}}, $_;
}

foreach (keys %vcfs){
  my $files = join(' ', @{$vcfs{$_}} );

  my $types = join(",",  @{$types{$_}} );
  my $merged_vcf = "$_\_merged_SVs.vcf";

  print "Running mergevcf $files -l $types -o $merged_vcf\n";
  system("mergevcf $files -l $types -o $merged_vcf");
}

sub usage {
  print
"
usage: perl merge_vcf.pl [-h]

svMerger
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v0.1
description: Read all VCF files in working directory and merge files sharing a common id (filename: `id\.*.vcf`)

arguments:
  -h, --help            show this help message and exit
"
}
