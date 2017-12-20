#!/usr/bin/perl
use strict;
use warnings;
use feature qw/ say /;
use Data::Printer;
use Data::Dumper;
use autodie;

use File::Basename;
use FindBin qw($Bin);
use FindBin '$Script';

exit usage() unless scalar @ARGV == 4;

my @keys = qw / 2L 2R 3L 3R 4 X Y /;
my %chrom_filt;

$chrom_filt{$_} = 1 for (@keys);

my ($vars_in, $allele_bal, $fusions_in, $read_support_filter) = @ARGV;

my ($var_ref, $clustered_events_ref) = extractVars($vars_in);

# my $allele_bal = shift;
my $al_ref = extractAlleles($allele_bal);

# my $fusions_in = shift;
my ($line_ref) = extractFusions($fusions_in, $var_ref, $al_ref, $clustered_events_ref);

my @name_fields = split( /\./, basename($vars_in) );
my $dir = "$Bin/../filtered/summary/";

# my $read_support_filter = shift;

say $read_support_filter;

open my $out, '>', "$dir" . $name_fields[0] . ".meerkat.filtered.summary.txt";
print $out join("\t", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2", "split_reads", "disc_reads", 'genotype', "id", "length(Kb)", "position", "consensus|type", "microhomology", "configuration", "allele_frequency", "mechanism|log2(cnv)") . "\n";

my @lines = @{ $line_ref };
print $out "$_\n" foreach @lines;


sub extractVars {
  my $in = shift;
  open my $vars, '<', $in;
  my %vars;
  my %clustered_events;
  while(<$vars>){
    chomp;
    my @parts = split(/\t/);
    my ($type, $mechanism, $id, $PE, $SR, $chr) = @parts[0..5];
    my @clusters = (split(/\//, $id));

    for my $i ( 0 .. $#clusters ){
      $clustered_events{$clusters[$i]} = $id if scalar @clusters > 1;
      foreach my $item (@parts) {
        if ( scalar (split /\//, $item) > 1 ){
          # say "$item can be split";
          push @{$vars{$clusters[$i]}}, (split /\//, $item, scalar @clusters)[$i];
        }
        else {
            # say "$item can't be split";
          push @{$vars{$clusters[$i]}}, $item;
        }
      }
    }

  }
  # p(%vars);
  return(\%vars, \%clustered_events );
}


sub extractAlleles {
  my $in = shift;
  open my $ab, '<', $in;
  my %abs;

  while(<$ab>){
    chomp;
    my @parts = split(/\t/, $_);
    my $id = $parts[2];
    my @clusters = (split(/\//, $id));

    for my $i ( 0 .. $#clusters ){
      foreach my $item (@parts) {
        if ( scalar (split /\//, $item) > 1 ){
          push @{$abs{$clusters[$i]}}, (split /\//, $item, scalar @clusters)[$i];
        }
        else {
          push @{$abs{$clusters[$i]}}, $item;
        }
      }
    }
  }
  # p(%abs);
  return(\%abs);
}


sub extractFusions {
  my ($in, $var_ref, $al_ref, $clustered_events_ref ) = @_;
  open my $fusions, '<', $in;
  my @cols = qw(type1	type2	type3	chrA	posA	oriA	geneA	exon_intronA	chrB	posB	oriB	geneB	exon_intronB	event_type	mechanism	event_id	disc_pair	split_read	homology	partners);
  my @lines;

  my %variants = %{ $var_ref };
  my %alleles  = %{ $al_ref };
  my %clusters = %{ $clustered_events_ref };

  while(<$fusions>){
    chomp;
    my @parts = split(/\t/);

    # say "$_\t$cols[$_]=$parts[$_]" for 0..$#cols;
    # say "---";

    my ($config) = join('_', grep { /\S/ } @parts[0..2]);
    my $event = $parts[13];
    ($event) = 'DEL' if $event =~ /del/;
    ($event) = 'DEL_INV' if $event =~ /del_inv/;
    ($event) = 'INV' if $event =~ /^inv/;
    ($event) = 'TANDUP' if $event =~ /dup/;
    ($event) = 'INS' if $event =~ /^ins/;
    ($event) = 'TRA' if $event =~ /^trans/;

    my ( @var_parts ) = @{ $variants{$parts[15]}};

    # say "$_ $var_parts[$_]" for 0..$#var_parts;
    # say "---";

    my ($chr1, $bp1, $chr2, $bp2) = ($parts[3], $parts[4], $parts[8], $parts[9]);

    my $lookup = "$chr1:$bp1-$bp2";
    $lookup = "$chr1:" . abs($bp1) . " " . "$chr2:" . abs($bp2) if $event eq 'TRA';

    my (@ab_parts) = @{ $alleles{$parts[15]} };
    my $ab_group =  $ab_parts[-1];
    my ( $disc_reads, $split_reads, $concordant_reads_bp1, $concordant_reads_bp2 ) = split("_", $ab_group);
    $disc_reads = $disc_reads =~ /(\d+)/;

    $disc_reads = 0 unless $disc_reads;
    $split_reads = 0 unless $split_reads;
    $concordant_reads_bp1 = 0 unless $concordant_reads_bp1;
    $concordant_reads_bp2 = 0 unless $concordant_reads_bp2;

    # A: number of full length discordant read pairs
    # B: number of discordant read pairs from partially mapped reads (clipped reads)
    # C: number of concordant read pairs at the first breakpoint
    # D: number of concordant read pairs at the second breakpoint

    my $total_disc_reads = $disc_reads + $split_reads;

    # Filter out if doesn't pass $read_support_filter
    my $tum_support = $parts[16] + $parts[17];
    if ( $tum_support < $read_support_filter) {
      say "Filtering on read support: $tum_support < $read_support_filter";
      next;
    }

    my $bp1_reads = $total_disc_reads + $concordant_reads_bp1 + 0.001;
    my $bp2_reads = $total_disc_reads + $concordant_reads_bp2 + 0.001;

    my $bp1_freq = sprintf("%.2f", $total_disc_reads/$bp1_reads);
    my $bp2_freq = sprintf("%.2f", $total_disc_reads/$bp2_reads);
    my $allele_frequency = sprintf("%.2f", ($bp1_freq + $bp2_freq)/2);

    # say "---";
    # say "disc:$disc_reads";
    # say "split:$split_reads";
    # say "conBP1:$concordant_reads_bp1";
    # say "conBP2:$concordant_reads_bp2";
    # say "freqBP1: $bp1_freq";
    # say "freqBP2: $bp2_freq";

    # say "$_ $ab_parts[$_]" for 0..$#ab_parts;
    # say "---";

    my $length = $var_parts[8];

    if ($event eq 'TRA'){
      $length = 'NA';
      $bp1 = abs($bp1);
      $bp2 = abs($bp2)
    }

    $length = sprintf("%.1f", abs($length)/1000) unless $event eq 'TRA';

    # Skip var unless both chroms is fully assembled
    next unless (exists $chrom_filt{$parts[3]} and exists $chrom_filt{$parts[8]});

    next if $allele_frequency < 0.1;

    my $mehanism;
    if ( $clusters{$parts[15]} ){
      my @linked_events = split /\//, $clusters{$parts[15]};
      my $linked = join(', ', @linked_events);
      $mehanism = $var_parts[1] . " ($linked)";
    }
    else{
      $mehanism = $var_parts[1];
    }

    push @lines, join("\t", "Meerkat",                          # source
                            $event,                             # type
                            $chr1,                              # chrom1
                            $bp1,                               # bp1
                            $chr2,                              # chrom2
                            $bp2,                               # bp2
                            $parts[17],                         # split reads
                            $parts[16],                         # paired reads
                            'somatic_tumour',                   # tag as somatic tumour
                            $parts[15],                         # id
                            $length,                            # length
                            $lookup,                            # IGV
                            $parts[13],                         # misc1 (type)
                            $parts[18],                         # microhomology
                            $config,                            # configuration
                            $allele_frequency,                  # allele frequency
                            $mehanism);                         # misc (mechanism)

  }
  return(\@lines);
}


sub usage {
  print
"
usage: $Script [variants file] [allele_frequency file] [fusions file]

sv2gene
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v.0.1
description: Extract breakpoint info for Meerkat data

example:
  perl $Script sample.variants sample.variants_af sample.fusions
"
}
