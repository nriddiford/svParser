#!/usr/bin/perl
use strict;
use warnings;
use autodie;

use feature qw/ say /;

use Data::Dumper;

my (%samples, %chroms, %genes, %gene_data, %features);

my $in_file = $ARGV[0];
open my $in, '<', $in_file;

while(<$in>){
  chomp;
  my ($sample, $chrom, $bp, $gene, $feature) = (split);
  $samples{$sample}++;
  $chroms{$chrom}++;
  $genes{$gene}++ unless $gene eq 'intergenic';
  push @{$gene_data{$gene}} , $sample;
  $features{$feature}++ unless $feature eq 'intergenic';
}

my $top_samples = 0;
print "Top 10 samples for structural variant count:\n";
for ( sort { $samples{$b} <=> $samples{$a} } keys %samples ){
  say "$_: " . ($samples{$_}/2);
  $top_samples++;
  last if $top_samples == 10;
}

print "\n";

print "Structural variant calls per chromosome for all samples:\n";
say "$_: " . $chroms{$_} for sort { $chroms{$b} <=> $chroms{$a} } keys %chroms;
print "\n";

my $top_genes = 0;
print "10 most hit genes accross all samples\n";
for ( sort { $genes{$b} <=> $genes{$a} } keys %genes ){
  say "$_: " . "$genes{$_} " . join(", ", @{$gene_data{$_}});
  $top_genes++;
  last if $top_genes == 10;
}

print "\n";

my $top_features = 0;
print "10 most hit features accross all samples\n";
for ( sort { $features{$b} <=> $features{$a} } keys %features ){
  print join("\t", $_, $features{$_}) . "\n";
  $top_features++;
  last if $top_features == 10;
}

print "\n";
