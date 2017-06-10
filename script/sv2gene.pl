#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Basename;

use Data::Printer;

use feature qw/ say /;

# my $bed_file = "/Users/Nick_curie/Documents/Curie/Data/Genomes/Dmel_v6.12/Features/dmel-all-r6.12.gtf";

my $bed_file = $ARGV[0];

open my $bed_in, '<', $bed_file;

my %genes;
while(<$bed_in>){
  chomp;
  next unless (split)[2] eq 'gene';
  my ($chrom, $start, $stop, $name) = (split)[0,3,4,11];
  ($name) = $name =~ /"(.*)\";/;
  push @{$genes{$chrom}{$name}} , $start, $stop;
}

my $in = $ARGV[1];
open my $SV_in, '<', $in;

my ( $name, $extention ) = split(/\.([^.]+)$/, basename($in), 2);

my ($sample) = split(/_/, $name, 3);

open my $annotated_svs, '>', $sample . ".anno_SVs.txt";

open my $genes_out, '>>', 'all_genes.txt';

open my $bp_out, '>>', 'all_bps.txt';

my $call = 1;
my %hits;

while(<$SV_in>){
  chomp;
  if (/source/){
    print $annotated_svs join("\t", $_, "bp1 gene", "bp2 gene", "affected genes", "notes") . "\n";
    next;
  }
  my ($type, $chrom1, $sv_start, $chrom2, $sv_end, $length) = (split)[1..5, 9];

  my @hit_genes;
  my $hit_bp1 = "-";
  my $hit_bp2 = "-";

  for my $gene ( sort { $genes{$chrom1}{$a}[0] <=> $genes{$chrom1}{$b}[0] } keys %{$genes{$chrom1}} ){

    my ($gene_start, $gene_stop) = @{$genes{$chrom1}{$gene}};

    # for DELS and DUPs, save all genes with both start and stop within BP1 and BP2
    if ( ( $gene_start >= $sv_start and $gene_stop <= $sv_end ) and ( $type eq 'DEL' or $type eq 'DUP' ) ) {
      push @hit_genes, $gene;
      $hits{$gene}++;
    }

    # save gene containing BP1
    if ( $sv_start >= $gene_start and $sv_start <= $gene_stop ) {
      push @hit_genes, $gene unless $hits{$gene};
      $hits{$gene}++;
      $hit_bp1 = $gene;
      print $bp_out join ("\t", $gene, $sample, $sv_start) . "\n";
    }

    # save gene containing BP2
    if ( $sv_end >= $gene_start and $sv_end <= $gene_stop ) {
      push @hit_genes, $gene unless $hits{$gene};
      $hits{$gene}++;
      $hit_bp2 = $gene;
      print $bp_out join ("\t", $gene, $sample, $sv_end) . "\n";
    }

  }

  print $genes_out $_ . "\n" foreach @hit_genes;

  my $affected_genes = scalar @hit_genes;
  my $joined_genes = join(", ", @hit_genes);

  if ($type eq 'DEL' or $type eq 'DUP'){
    say "SV $call: $length kb $type affecting $affected_genes genes: $joined_genes ";
  }
  else {
    say "SV $call: $length kb $type with break points in $sv_start and $sv_end: $joined_genes";
  }

  if (scalar(@hit_genes) == 0){
    push @hit_genes, "-";
  }

  my $joined_genes2print = join(", ", @hit_genes);

  print $annotated_svs join("\t", $_, $hit_bp1, $hit_bp2, $joined_genes2print, " ") . "\n";
  $call++;
}
