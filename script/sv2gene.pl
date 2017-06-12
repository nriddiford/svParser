#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Basename;

use Data::Printer;

use Data::Dumper;

use feature qw/ say /;

my $bed_file = $ARGV[0];

open my $bed_in, '<', $bed_file;

my (%genes, %features);
while(<$bed_in>){
  chomp;
  my ($chrom, $feature, $start, $stop, $gene) = (split)[0,2,3,4,11];
  ($gene) = $gene =~ /\"(.*)\";/;

  # Take longest feature for each gene (i.e. longest exon for gene 'X')
  if ( (not exists $features{$chrom}{$gene}{$feature}) or ( ($features{$chrom}{$gene}{$feature}[1] - $features{$chrom}{$gene}{$feature}[0]) <= ($stop - $start)) ){
    $features{$chrom}{$gene}{$feature} = [$start, $stop];
  }
  $genes{$chrom}{$gene} = [$start, $stop] if $feature eq 'gene';
}

my $in = $ARGV[1];
open my $SV_in, '<', $in;

my ( $name, $extention ) = split(/\.([^.]+)$/, basename($in), 2);

my ($sample) = split(/_/, $name, 3);

open my $annotated_svs, '>', $sample . ".anno_SVs.txt";

open my $genes_out, '>>', 'all_genes.txt';

open my $bp_out, '>>', 'all_bps.txt';

my $call = 1;

while(<$SV_in>){
  chomp;
  if (/source/){
    print $annotated_svs join("\t", $_, "bp1 feature", "bp2 feature", "affected genes", "notes") . "\n";
    next;
  }
  my ($type, $chrom1, $bp1, $chrom2, $bp2, $length) = (split)[1..5, 9];

  my (%hits, $hits);
  my @hit_genes;
  my $hit_genes;
  my $hit_bp1 = "intergenic";
  my $hit_bp2 = "intergenic";

  if ($type eq "DEL" or $type eq "DUP"){
    ($hit_bp1, $hit_genes, $hits) = getbps($chrom1, $bp1, $hit_bp1, \@hit_genes, \%hits);
    ($hit_genes, $hits)           = getgenes($chrom1, $bp1, $bp2, $hit_genes, $hits);
    ($hit_bp2, $hit_genes, $hits) = getbps($chrom2, $bp2, $hit_bp2, $hit_genes, $hits);
    @hit_genes = @{ $hit_genes };
    %hits = %{$hits};
  }
  else {
    ($hit_bp1, $hit_genes, $hits) = getbps($chrom1, $bp1, $hit_bp1, \@hit_genes, \%hits);
    ($hit_bp2, $hit_genes, $hits) = getbps($chrom2, $bp2, $hit_bp2, $hit_genes, $hits);
    @hit_genes = @{ $hit_genes };
    %hits = %{$hits};
  }

  print $genes_out $_ . "\n" foreach @hit_genes;

  my $affected_genes = scalar @hit_genes;
  my $joined_genes = join(", ", @hit_genes);

  if ($type eq 'DEL' or $type eq 'DUP'){
    say "SV $call: $length kb $type affecting $affected_genes genes: $joined_genes. Bp1: $hit_bp1 Bp2: $hit_bp2 ";
  }
  else {
    say "SV $call: $length kb $type with break points in $chrom1\:$bp1 ($hit_bp1) and $chrom2\:$bp2 ($hit_bp2)";
  }

  if (scalar(@hit_genes) == 0){
    push @hit_genes, "-";
  }

  my $joined_genes2print = join(", ", @hit_genes);

  print $annotated_svs join("\t", $_, $hit_bp1, $hit_bp2, $joined_genes2print, " ") . "\n";
  $call++;
}

sub getgenes {

  my ($chrom1, $bp1, $bp2, $hit_genes, $hits) = @_;

  my @hit_genes = @{ $hit_genes };
  my %hits = %{$hits};

  for my $gene ( sort { $genes{$chrom1}{$a}[0] <=> $genes{$chrom1}{$b}[0] } keys %{$genes{$chrom1}} ){
    my ($gene_start, $gene_stop) = @{$genes{$chrom1}{$gene}};
    # for DELS and DUPs, save all genes with both start and stop within BP1 and BP2
    if ( $gene_start >= $bp1 and $gene_stop <= $bp2 ) {
      push @hit_genes, $gene;
      $hits{$gene}++;
    }
  }
  return (\@hit_genes, \%hits);
}

sub getbps {

  my ($chrom, $bp, $hit_bp, $hit_genes, $hits) = @_;

  my @hit_genes = @{ $hit_genes };
  my %hits = %{$hits};
  my %smallest_hit_feature;

  for my $gene ( sort { $genes{$chrom}{$a}[0] <=> $genes{$chrom}{$b}[0] } keys %{$genes{$chrom}} ){
    for my $feature (reverse sort keys %{$features{$chrom}{$gene}}){

      my ($feature_start, $feature_stop) = @{$features{$chrom}{$gene}{$feature}};
      my $length = ($feature_stop - $feature_start);

      # if ($gene eq 'CR32773' ){
      #   print join(",", $feature, $length) . "\n";
      # }

      # if breakpoint contained in feature
      if ( $bp >= $feature_start and $bp <= $feature_stop ) {
        # save gene containing BP
        push @hit_genes, $gene unless $hits{$gene};
        $hits{$gene}++;
        # find smallest feature and save
        if ( (not exists $smallest_hit_feature{$gene}) or ($smallest_hit_feature{$gene} > $length ) ){
            $smallest_hit_feature{$gene} = $length;
            $hit_bp = "$gene, $feature";
            print $bp_out join ("\t", $gene, $feature, $sample, $chrom, $bp) . "\n";
        }
      }
    }

  }
  return ($hit_bp, \@hit_genes, \%hits);
}
