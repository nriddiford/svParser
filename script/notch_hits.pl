use strict;
use warnings;
use autodie;
use Data::Dumper;
use feature qw/ say /;
use Data::Printer;

open my $out, '>', 'Notch_hits.txt';
open my $new_genes, '>', 'all_genes_cleaned.txt';

# This allows files to be parsed irrespective of line ending. Not recommened for large files
local $/ = undef;
my $content = <>;
my @lines = split /\r\n|\n|\r/, $content;

my $N_count = 0;
my $N_DEL_count = 0;
my $sample_count = 0;
my %events;
my %seen;
my %samples;
my %hit_samples;
my %not_hit;

my %hit_genes;
my %notchSVs;
foreach(@lines){
  chomp;
  print $out "$_\n" if /sample/;
  my @cells = split(/\t/);

  next if $events{$cells[0]}{$cells[1]}++;
  $samples{$cells[0]}++ unless /sample/;
  my @affected_genes = split(/,/, $cells[20]);
  my $line = $_;

  my ($sample, $event, $caler, $type, $chrom1, $bp1, $chrom2, $bp2) = @cells[0..7];
  next unless $chrom1 and $chrom2 eq 'X';
  my ($nstart, $nstop, $buff) = (3134870, 3172221,200000);
  if($type eq 'INV' or $type eq 'BND' or $type eq 'DUP'){
    if ( ( ($bp1 >= $nstart and $bp1 <= $nstop) or ($bp2 >= $nstart and $bp2 <= $nstop) ) or (($bp1 <= $nstart and $bp1 >= $nstart-$buff) and ($bp2 >= $nstop and $bp2 <= $nstop+$buff)) ){
      # print join("\t", $sample, $type, $bp1, $bp2) . "\n";
      $hit_samples{$sample}++;
      print $out "$line\n";
    }
  }

  elsif($type eq 'TRA'){
    if ( ($bp1 >= $nstart and $bp1 <= $nstop) or ($bp2 >= $nstart and $bp2 <= $nstop) ){
      # print join("\t", $sample, $type, $bp1, $bp2) . "\n";
      $hit_samples{$sample}++;
      print $out "$line\n";
    }
  }


  for(@affected_genes){
    s/\"//g;
    s/ //g;
    next if /-/;
    print $new_genes "$event\t$sample\t$type\t$chrom1\t$_" . "\n";
    $hit_genes{$_}++;
    if (/^N$/){
      if ($type eq 'DEL' or $type eq 'DUP'){
        print $out "$line\n" if $type eq 'DEL';
        $N_DEL_count++ if $type eq 'DEL';
      }
      $hit_samples{$sample}++;
    }
  }
  $not_hit{$sample}++ unless $hit_samples{$sample};
}

delete $not_hit{$_} for keys %hit_samples;

say "No hits in $_" for keys %not_hit;

$sample_count++ for keys %samples;
$N_count++ for keys %hit_samples;

say "$sample_count samples";
say "$N_count samples with SVs affecting Notch";
say "$N_DEL_count samples with DELS affecting Notch";


for my $gene (sort keys %hit_genes){
  say "$gene\t$hit_genes{$gene}" if $hit_genes{$gene} > 5;
}

# event   source  type    chromosome1     bp1     chromosome2     bp2     split reads     pe reads        id      length(Kb)      position        consensus|type  microhomology   configuration   allele_frequency        mechanism|log2(cnv)     bp1_locus       bp2_locus       affected_genes  T/F     notes
# 2       lumpy   DEL     X       3117762 X       3177110 2       2       4349    59.3    X:3117762-3177110       -       -       -       0.19    -       kirre, intron   dnc, intron
#     kirre, N, CG18508, Fcp3C, CG3939, dnc
