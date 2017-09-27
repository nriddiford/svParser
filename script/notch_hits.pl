use strict;
use warnings;
use autodie;
use Data::Dumper;
use feature qw/ say /;

open my $out, '>', 'Notch_hits.txt';

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

foreach(@lines){
  chomp;
  print $out "$_\n" if /sample/;
  my @cells = split(/\t/);

  next if $events{$cells[0]}{$cells[1]}++;
  $samples{$cells[0]}++ unless /sample/;
  my @affected_genes = split(/,/, $cells[20]);
  my $line = $_;

  for(@affected_genes){
    s/\"//g;
    s/ //g;
    if (/^N$/){
      print $out "$line\n";
      my $sample = $cells[0];
      $hit_samples{$sample}++;
      say "$cells[0] has no hits in Notch" unless $hit_samples{$sample};

      unless ($seen{$cells[0]}++){
        $N_count++;
        $N_DEL_count++ if $cells[3] eq 'DEL';
      }
    }
  }

}

$sample_count++ for keys %samples;
say "$sample_count samples";
say "$N_count samples with SVs affecting Notch";
say "$N_DEL_count samples with DELS affecting Notch";
