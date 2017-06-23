#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Basename;
use Data::Printer;
use Data::Dumper;
use feature qw/ say /;

my (%transcript_length, %genes, %features);

my ($sample, $annotated_svs, $genes_out, $bp_out);

make_gene_hash($ARGV[0]);

annotate_SVs($ARGV[1]);

sub make_gene_hash {
  my $bed_file = shift;

  open my $bed_in, '<', $bed_file;
  my %exons;

  while(<$bed_in>){
    chomp;
    my ($chrom, $feature, $start, $stop, $gene) = (split)[0,2,3,4,11];
    ($gene) = $gene =~ /\"(.*)\";/;

    if ($feature eq 'gene'){
      $genes{$chrom}{$gene} = [$start, $stop];
    }

    else {
      my $transcript = (split)[13];
      ($transcript) = $transcript =~ /\"(.*)\";/;
      my $feature_length = ($stop - $start);
      if ($feature eq 'exon'){
        $exons{$transcript}{$feature}++;
        $feature = $feature . "_" . $exons{$transcript}{$feature};
      }
      $transcript_length{$chrom}{$gene}{$transcript} += $feature_length;
      $features{$chrom}{$gene}{$transcript}{$feature} = [$start, $stop, $feature_length];
    }
  }

}

sub annotate_SVs {
  my $in = shift;

  open my $SV_in, '<', $in;

  my ( $name, $extention ) = split(/\.([^.]+)$/, basename($in), 2);
  ($sample) = split(/_/, $name, 3);
  open $annotated_svs, '>', $sample . ".annotated_SVs.txt";
  open $genes_out, '>>', 'all_genes.txt';
  open $bp_out, '>>', 'all_bps.txt';

  my $call = 1;

  while(<$SV_in>){
    chomp;
    if (/source/){
      print $annotated_svs join("\t", $_, "bp1 locus", "bp2 locus", "affected genes", "notes") . "\n";
      next;
    }
    my ($event, $source, $type, $chrom1, $bp1, $chrom2, $bp2, $length) = (split)[0..6, 10];

    my (%hits, $hits);
    my @hit_genes;
    my $hit_genes;
    my $hit_bp1 = "intergenic";
    my $hit_bp2 = "intergenic";

    if ($type eq "DEL" or $type eq "DUP"){
      ($hit_bp1, $hit_genes, $hits) = getbps('bp1', $event, $type, $chrom1, $bp1, $hit_bp1, \@hit_genes, \%hits);
      ($hit_genes, $hits)           = getgenes($chrom1, $bp1, $bp2, $hit_genes, $hits);
      ($hit_bp2, $hit_genes, $hits) = getbps('bp2', $event, $type, $chrom2, $bp2, $hit_bp2, $hit_genes, $hits);
      @hit_genes = @{ $hit_genes };
      %hits = %{ $hits };
    }
    else {
      ($hit_bp1, $hit_genes, $hits) = getbps('bp1', $event, $type, $chrom1, $bp1, $hit_bp1, \@hit_genes, \%hits);
      ($hit_bp2, $hit_genes, $hits) = getbps('bp2', $event, $type, $chrom2, $bp2, $hit_bp2, $hit_genes, $hits);
      @hit_genes = @{ $hit_genes };
      %hits = %{ $hits };
    }

    print $genes_out $_ . "\n" foreach @hit_genes;

    my $affected_genes = scalar @hit_genes;
    my $joined_genes = join(", ", @hit_genes);

    if ($type eq 'DEL' or $type eq 'DUP'){
      say "SV $call: $length kb $type affecting $affected_genes genes: $joined_genes Bp1: $hit_bp1 Bp2: $hit_bp2 ";
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
  my ($bp_id, $event, $type, $chrom, $bp, $hit_bp, $hit_genes, $hits) = @_;

  my @hit_genes = @{ $hit_genes };
  my %hits = %{$hits};
  my %smallest_hit_feature;
  my $bp_feature = "intergenic";
  my $bp_gene = "intergenic";


  for my $gene ( sort { $genes{$chrom}{$a}[0] <=> $genes{$chrom}{$b}[0] } keys %{$genes{$chrom}} ){
    # Smallest features are last (and will then replace larger overlapping features i.e. exon over CDS)
    my %smallest_hit_feature;

    for my $transcript ( sort keys %{$transcript_length{$chrom}{$gene}}){

    for my $feature ( sort { $features{$chrom}{$gene}{$transcript}{$b}[2] <=> $features{$chrom}{$gene}{$transcript}{$a}[2] } keys %{$features{$chrom}{$gene}{$transcript}}){

      # print join(",", $gene, $transcript, $transcript_length{$chrom}{$gene}{$transcript}, $feature, $features{$chrom}{$gene}{$transcript}{$feature}[2]  ) . "\n";

      my ($feature_start, $feature_stop, $length) = @{$features{$chrom}{$gene}{$transcript}{$feature}};

      $feature = "intron" if $feature eq 'gene';
      $feature = "intron" if $feature eq 'mRNA';

      # if ($gene eq 'bun' ){
      #   print join(",", $gene, $feature, $length) . "\n";
      # }

      # if breakpoint contained in feature
      if ( $bp >= $feature_start and $bp <= $feature_stop ) {
        # save gene containing BP
        push @hit_genes, $gene unless $hits{$gene};
        $hits{$gene}++;

        # take smallest feature that is hit accross all transcript
        if ( (not exists $smallest_hit_feature{$gene}) or ($smallest_hit_feature{$gene} > ($feature_stop - $feature_start)) ){
          $smallest_hit_feature{$gene} = ($feature_stop - $feature_start);
          $bp_feature = $feature;
          $bp_gene = $gene;
          $hit_bp = "$gene, $feature";
        }

      }
    }
  }

}
  print $bp_out join ("\t", $event, $bp_id, $sample, $chrom, $bp, $bp_gene, $bp_feature, $type) . "\n";

  return ($hit_bp, \@hit_genes, \%hits);
}
