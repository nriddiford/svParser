#!/usr/bin/perl
use strict;
use warnings;
use autodie;
use File::Basename;
use Data::Dumper;
use Data::Printer;

use feature qw/ say /;
use FindBin '$Script';
use Getopt::Long qw/ GetOptions /;

my $sv_calls;
my $help;
my $features;
my $reannotate;
my $blacklist;
my $whitelist;
my $somatic = 0;
my $mark = 0;

GetOptions( 'infile=s'         =>    \$sv_calls,
            'features=s'       =>    \$features,
            'somatic'          =>    \$somatic,
            'mark'             =>    \$mark,
            're-annotate'      =>    \$reannotate,
            'blacklist=s'      =>    \$blacklist,
            'whitelist=s'      =>    \$whitelist,
            'help'             =>    \$help

    ) or die usage();

if ($help) { exit usage() }

if (not $sv_calls and not $features ){
  exit usage();
}

if ($somatic){
  say "Running in somatic mode - will annotate all germline events as 'NA'";
}
if ($mark){
  say "Running in mark mode. Somatic dels and dups with abs(log2(FC) < 0.26) will be marked as F in T/F column";
}

my (%false_positives, %true_positives);

if ($blacklist){
  say "Marking calls in $blacklist as false positives";
  open my $blacklist_file, '<', $blacklist;
  while(<$blacklist_file>){
    chomp;
    $false_positives{$_}++;
  }
}

if ($whitelist){
  say "Using whitelist $whitelist";
  open my $whitelist_file, '<', $whitelist;
  while(<$whitelist_file>){
    chomp;
    my @line = split(//);
    my $lookup = (split(""))[0];
    # shouldn't include caller (?) 12.1.18
    #$lookup = $lookup . "_" . $line[2];
    my @cols = @line[1..$#line];
    $true_positives{$lookup} = [@cols];
  }
}

my (%transcript_length, %genes, %features);
my ($sample, $annotated_svs, $genes_out, $bp_out);

make_gene_hash($features);
annotate_SVs($sv_calls);

sub make_gene_hash {
  my $bed_file = shift;
  open my $bed_in, '<', $bed_file;
  my %exons;
  while(<$bed_in>){
    chomp;
    my ($chrom, $feature, $start, $stop, $id, $gene) = (split)[0,2,3,4,9,11];
    ($gene) = $gene =~ /\"(.*)\";/;
    $gene =~ s/\"//g;
    ($id) = $id =~ /\"(.*)\";/;
    if ($feature eq 'gene'){
      $genes{$chrom}{$gene} = [$start, $stop, $id];
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
  ($sample) = (split(/_/, $name))[0];

  if ($reannotate){
    open $annotated_svs, '>', $sample . "_reannotated_SVs.txt";
    print "Reannotating SV calls from $sample\n";
    open $genes_out, '>>', 'all_genes_filtered.txt';
    open $bp_out, '>>', 'all_bps_filtered.txt';
  }
  else{
    open $annotated_svs, '>', $sample . "_annotated_SVs.txt";
    print "Annotating SV calls from $sample\n";
    open $genes_out, '>>', 'all_genes.txt';
    open $bp_out, '>>', 'all_bps.txt';
  }

  my $call = 1;
  # This allows files to be parsed irrespective of line ending. Not recommened for large files
  local $/ = undef;
  my $content = <$SV_in>;
  my @lines = split /\r\n|\n|\r/, $content;

  my (%events, %record, @header);

  my @col_order = qw/ event source type chromosome1 bp1 chromosome2 bp2 split_reads disc_reads genotype id length(Kb) position consensus microhomology configuration allele_frequency mechanism log2(cnv) bp1_locus	bp2_locus	affected_genes status	notes /;

  for my $line (@lines){
    chomp;
    my @cells = split(/\t/, $line);
    if ($line =~ /event/){
      push @header, @cells;
      print $annotated_svs join("\t", @col_order) . "\n";
      next;
    }

    @record{@header} = @cells;

    # foreach my $item (@header){
    #   print $record{$item} . "\t";
    # }
    # print "\n";

    # print join("\t", @record{@col_order} ) . "\n";
    my $event = $record{event};
    my $source = $record{source};
    my $type = $record{type};
    my $chrom1 = $record{chromosome1};
    my $bp1 = $record{bp1};
    my $sr = $record{split_reads};
    my $pe = $record{disc_reads};
    my $chrom2 = $record{chromosome2};
    my $bp2 = $record{bp2};
    my $genotype = $record{genotype};
    my $af = $record{allele_frequency};
    my $length = $record{'length(Kb)'};
    my $cnv = $record{'log2(cnv)'};

    my $precise;
    if ($sr eq '-' and $pe eq '-'){
      $precise = 'imprecise';
    }
    else{
      $precise = 'precise';
    }

    # my ($event, $source, $type, $chrom1, $bp1, $chrom2, $bp2, $sr, $pe, $genotype, $id, $length, undef, undef, undef, undef, $af, $mech, $cnv, $bp1_locus,	$bp2_locus,	$affected_genes, $status,	$notes) = @record{@header};
    # Check to see if the SV has already been annotated - print and skip if next
    if ($genotype !~ 'somatic_tumour' and $somatic){
      print $annotated_svs join("\t", $line, "NA", "NA", "-", "", "") . "\n" if not $reannotate;
      print $annotated_svs "$line\n" if $reannotate;
      next;
    }
    if ($record{status} and $record{status} eq "F"){
      foreach my $col (@col_order){
        if (exists($record{$col})){
          print $annotated_svs "$record{$col}\t";
        }
        else{
          print $annotated_svs "\t";
        }
      }
      print $annotated_svs "\n";
      next;
    }
    # if ( $cells[18] and $cells[18] ne ' ' and $cells[18] ne '-' and $reannotate ){
    #   print $annotated_svs "$_\n";
    #   next if $events{$sample}{$event}++;
    #   my ($bp1_locus, $bp2_locus, $affected_genes, undef, undef) = @cells[18..22];
    #   my @genes;
    #
    #   if ($affected_genes =~ /\,/){
    #     push @genes, split(', ', $affected_genes);
    #   }
    #   else{
    #     push @genes, $affected_genes;
    #   }
    #   print $genes_out join("\t", $event, $sample, $genotype, $type, $chrom1, $_) . "\n" foreach @genes;
    #
    #   $bp1_locus =~ s/"//g;
    #   $bp2_locus =~ s/"//g;
    #
    #   my $bp1_feature = 'intergenic';
    #   my $bp2_feature = 'intergenic';
    #   my $bp1_gene    = 'intergenic';
    #   my $bp2_gene    = 'intergenic';
    #
    #   for my $bptype ('bp1', 'bp2'){
    #     my ($bp1_gene, $bp1_feature) = split(", ", $bp1_locus);
    #     $bp1_feature = 'intergenic' if not length $bp1_feature;
    #     my ($bp2_gene, $bp2_feature) = split(", ", $bp2_locus);
    #     $bp2_feature = 'intergenic' if not length $bp2_feature;
    #
    #     print $bp_out join("\t", $event, $bptype, $sample, $genotype, $chrom1, $bp1, $bp1_gene, $bp1_feature, $chrom2, $bp2, $bp2_gene, $bp2_feature, $type, $length, $af) . "\n" if $bptype eq 'bp1';
    #     print $bp_out join("\t", $event, $bptype, $sample, $genotype, $chrom2, $bp2, $bp2_gene, $bp2_feature, $chrom1, $bp1, $bp1_gene, $bp1_feature, $type, $length, $af) . "\n" if $bptype eq 'bp2';
    #   }
    #   next;
    # }

    # if it hasn't been annotated, trim off blank cells and proceed

    # elsif ( $reannotate or $whitelist ){
    #   no warnings;
    #   $line = join("\t", @cells[0..17]);
    # }

    my (%hits, $hits);
    my @hit_genes;
    my $hit_genes;
    my ($bp1_id, $bp2_id);
    my $hit_bp1 = "intergenic";
    my $hit_bp2 = "intergenic";

    if ( ($type =~ 'DEL|DUP' or $type =~ 'BND') and ($cnv eq '-' or abs($cnv) > 0.3) ){
      ($hit_bp1, $hit_genes, $hits) = getbps('bp1', $event, $type, $genotype, $chrom1, $bp1, $hit_bp1, $length, \@hit_genes, \%hits);
      ($hit_genes, $hits)           = getgenes($chrom1, $bp1, $bp2, $hit_genes, $hits);
      ($hit_bp2, $hit_genes, $hits) = getbps('bp2', $event, $type, $genotype, $chrom2, $bp2, $hit_bp2, $length, $hit_genes, $hits);
      @hit_genes = @{ $hit_genes };
      %hits = %{ $hits };
    }
    else {
      ($hit_bp1, $hit_genes, $hits) = getbps('bp1', $event, $type, $genotype, $chrom1, $bp1, $hit_bp1, $length, \@hit_genes, \%hits);
      ($hit_bp2, $hit_genes, $hits) = getbps('bp2', $event, $type, $genotype, $chrom2, $bp2, $hit_bp2, $length, $hit_genes, $hits);
      @hit_genes = @{ $hit_genes };
      %hits = %{ $hits };
    }

    # Does this really append to 'all_genes_filtered'? 12.12.17
    print $genes_out join("\t", $event, $sample, $genotype, $type, $record{allele_frequency}, $chrom1, $_) . "\n" foreach @hit_genes;
    for my $bptype ('bp1', 'bp2'){
      my ($bp1_gene, $bp1_feature) = split(", ", $hit_bp1);
      $bp1_feature = 'intergenic' if not length $bp1_feature;
      my ($bp2_gene, $bp2_feature) = split(", ", $hit_bp2);
      $bp2_feature = 'intergenic' if not length $bp2_feature;

      # per sample/event/type?
      if(not $events{$event}{$bptype}{$chrom1}{$bp1}{$chrom2}{$bp2}++){
        print $bp_out join("\t", $event, $bptype, $sample, $genotype, $chrom1, $bp1, $bp1_gene, $bp1_feature, $chrom2, $bp2, $bp2_gene, $bp2_feature, $type, $length, $record{allele_frequency}, $precise) . "\n" if $bptype eq 'bp1';
        print $bp_out join("\t", $event, $bptype, $sample, $genotype, $chrom2, $bp2, $bp2_gene, $bp2_feature, $chrom1, $bp1, $bp1_gene, $bp1_feature, $type, $length, $record{allele_frequency}, $precise) . "\n" if $bptype eq 'bp2';
      }
    }

    my $affected_genes = scalar @hit_genes;
    my $joined_genes = join(", ", @hit_genes);

    if ($genotype eq 'somatic_tumour') {
      if ($type =~ 'DEL|DUP'){
        say "SV $call: $length kb $type affecting $affected_genes genes: Bp1: $hit_bp1 Bp2: $hit_bp2 ";
      }
      else {
        say "SV $call: $length kb $type with break points in $chrom1\:$bp1 ($hit_bp1) and $chrom2\:$bp2 ($hit_bp2)";
      }
    }

    if (scalar(@hit_genes) == 0){
      push @hit_genes, "-";
    }

    my $joined_genes2print = join(", ", @hit_genes);

    # If blacklist specified, check to see if location is blacklisted (this is done by adding 'F' to T/F col, and running the file through 'clean.py')
    # If location is blacklisted, then carry over the 'F' tag
    my $whitelookup = join("_", $sample, $chrom1, $bp1, $chrom2, $bp2 );
    my $blacklookup = join("_", $sample, $chrom1, $bp1, $chrom2, $bp2);

    if ($whitelist and exists $true_positives{$whitelookup}){
      say "* Annotating call from whitelist: $whitelookup";
      my @cols = @{$true_positives{$whitelookup}};
      print $annotated_svs join("\t", $event, @cols[1..$#cols]) . "\n";
      $call++;
      next;
    }

    # if ($mark and $cnv ne '-') {
    #   # 4.10.18 # Need to find a sweet spot here !!
    #   # setting a threshold of log2(fc)*2 < 0.58 (fc == 1.5) on autosomes allows some spurious calls through...
    #
    #   my $mult = 1;
    #   if ($chrom1 eq 'X' or $chrom1 eq 'Y'){
    #     $mult = 1;
    #   }
    #   # If called by a CN approach then mark as FP unless > log2(1.5)
    #   if ( (abs($cnv)*$mult < 0.58) and $sr eq '-' and $pe eq '-' ){
    #     # say "Mark as FP: abs($cnv)*$mult < 0.58 with no sr";
    #     print $annotated_svs join("\t", $_, $hit_bp1, $hit_bp2, $joined_genes2print, "F", "Low FC in $type called by CN") . "\n";
    #     $call++;
    #     next;
    #   }
    #   # Unless there's read support, in which case only mark if < log2(1.2)
    #   elsif ( (abs($cnv)*$mult < 0.26) and ($type eq 'DEL' or $type eq 'DUP') and $length > .8 ) {
    #     # say "Mark as FP: abs($cnv)*$mult < 0.26 with sr";
    #     print $annotated_svs join("\t", $_, $hit_bp1, $hit_bp2, $joined_genes2print, "F", "Low FC in $type") . "\n";
    #     $call++;
    #     next;
    #   }
    # }

    if ($blacklist and exists $false_positives{$blacklookup}){
      say "* Marking blacklisted call as FP: $blacklookup";
      print $annotated_svs join("\t", $line, $hit_bp1, $hit_bp2, $joined_genes2print, "F") . "\n";
      $call++;
      # next;
    }

    else {
      # event	source	type	chromosome1	bp1	chromosome2	bp2	split_reads	disc_reads	genotype	id	length(Kb)	position	consensus	microhomology	configuration	allele_frequency	log2(cnv)	bp1_locus	bp2_locus	affected_genes	T/F	notes
      #
      # my $event = $record{event};
      # my $source = $record{source};
      # my $type = $record{type};
      # my $chrom1 = $record{chromosome1};
      # my $bp1 = $record{bp1};
      # my $sr = $record{split_reads};
      # my $pe = $record{disc_reads};
      # my $chrom2 = $record{chromosome2};
      # my $bp2 = $record{bp2};
      # my $genotype = $record{genotype};
      # my $af = $record{allele_frequency};
      $record{bp1_locus} = $hit_bp1;
      $record{bp2_locus} = $hit_bp2;
      $record{affected_genes} = $joined_genes2print;

      foreach my $col (@col_order){
        if (defined($record{$col})){
          # $record{$col} = $record{$col} =~ s/^\-$//;
          print $annotated_svs "$record{$col}\t";
        }
        else{
          print $annotated_svs "\t";
        }
      }
      print $annotated_svs "\n";

      # print $annotated_svs join("\t", @record{@col_order}, " ") . "\n";


      # print $annotated_svs join("\t", $line, $hit_bp1, $hit_bp2, $joined_genes2print, " ") . "\n";
      $call++;
    }

  }
}

sub getgenes {
  my ($chrom1, $bp1, $bp2, $hit_genes, $hits) = @_;

  my @hit_genes = @{ $hit_genes };
  my %hits = %{ $hits };

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
  my ($bp_id, $event, $type, $genotype, $chrom, $bp, $hit_bp, $length, $hit_genes, $hits) = @_;

  my @hit_genes = @{ $hit_genes };
  my %hits = %{$hits};
  my %smallest_hit_feature;
  my $bp_feature = "intergenic";
  my $bp_gene = "intergenic";
  my $bp_gene_id = "-";

  for my $gene ( sort { $genes{$chrom}{$a}[0] <=> $genes{$chrom}{$b}[0] } keys %{$genes{$chrom}} ){
    # Smallest features are last (and will then replace larger overlapping features i.e. exon over CDS)
    my %smallest_hit_feature;

    my $gene_id = $genes{$chrom}{$gene}[2];

    for my $transcript ( sort keys %{$transcript_length{$chrom}{$gene}} ){

    for my $feature ( sort { $features{$chrom}{$gene}{$transcript}{$b}[2] <=> $features{$chrom}{$gene}{$transcript}{$a}[2] } keys %{$features{$chrom}{$gene}{$transcript}} ){

      # print join(",", $gene, $transcript, $transcript_length{$chrom}{$gene}{$transcript}, $feature, $features{$chrom}{$gene}{$transcript}{$feature}[2]  ) . "\n";

      my ($feature_start, $feature_stop, $length) = @{$features{$chrom}{$gene}{$transcript}{$feature}};

      $feature = "intron" if $feature eq 'gene';
      $feature = "intron" if $feature eq 'mRNA';

      # if breakpoint contained in feature
      if ( $bp >= $feature_start and $bp <= $feature_stop ) {
        # save gene containing BP
        push @hit_genes, $gene unless $hits{$gene};
        $hits{$gene}++;

        # take smallest feature that is hit accross all transcript
        if ( (not exists $smallest_hit_feature{$gene}) or ($smallest_hit_feature{$gene} > ($feature_stop - $feature_start)) ){
          $smallest_hit_feature{$gene} = ($feature_stop - $feature_start);
          $bp_feature = $feature;
          $feature = 'exon' if $feature eq 'CDS'; # Adapted from snv2gene
          $bp_gene = $gene;
          $bp_gene_id = $gene_id;
          $hit_bp = $gene . ", " . $feature;
        }

      }
    }
  }

}
  # if ( $reannotate ){
  # Why did all_bps.txt not have the breakpoint id? 12.12.17
  # print $bp_out join("\t", $event, $bp_id, $sample, $genotype, $chrom, $bp, $bp_gene, $bp_feature, $type, $length) . "\n";
  # }
  # else {
  # print $bp_out join ("\t", $event, $sample, $chrom, $bp, $bp_gene, $bp_feature, $type, $length) . "\n";
  # }

  return ($hit_bp, \@hit_genes, \%hits);
}

sub usage {
  print
"
usage: $Script -i variants -f features.gtf -b blacklist -w whitelist

sv2gene
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v1.0
description: Annotate breakpoints of structural variants in svParser summary file

arguments:
  -h, --help            show this help message and exit
  -i, --infile          variant calls file (as produced by svParser)[required]
  -f, --features        features file to annotate from (must be in .gtf format)
  -s, --somatic         switch on somatic mode - only annotate somatic events
  -r, --reannotate      run in re-annotate mode and
  -b, --blacklist       specify blacklist to exclude variants
  -w, --whitelist       specify whitelist to proect variants
"
}
