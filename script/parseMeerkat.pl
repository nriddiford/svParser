#!/usr/bin/perl
use strict;
use warnings;
use feature qw/ say /;
use Data::Printer;
use Data::Dumper;
use autodie;

use File::Basename;
use FindBin qw($Bin);

# open my $vars, '<', '/Users/Nick_curie/Desktop/script_test/svParser/data/A573R25.RG.somatic_g.variants';
# open my $fusions, '<', '/Users/Nick_curie/Desktop/script_test/svParser/data/A573R25.RG.somatic_g.fusions';
# open my $var_hads, '<', '/Users/Nick_curie/Desktop/script_test/svParser/data/Meerkat_vars_heads.txt';

my $vars_in = shift;
my $var_ref = extractVars($vars_in);

my $fusions_in = shift;
my ($line_ref) = extractFusions($fusions_in);

my @name_fields = split( /\./, basename($vars_in) );
my $dir = "$Bin/../filtered/summary/";

open my $out, '>', "$dir" . $name_fields[0] . ".meerkat.filtered.summary.txt";

print $out join("\t", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2", "split reads", "pe reads", "id", "length(Kb)", "position", "consensus", "microhomology", "configuration", "read_depth_ratio", "misc") . "\n";

my @lines = @{$line_ref};
print $out "$_\n" foreach @lines;

my %vars;
my %fusions;

sub extractVars {
  my $in = shift;
  open my $vars, '<', $in;
  while(<$vars>){
    chomp;
    my @parts = split(/\t/);
    my ($type, $mechanism, $id, $PE, $SR, $chr) = @parts[0..5];

    $vars{$id} = [@parts];
  }
  return(\%vars);
}

sub extractFusions {
  my $in = shift;
  open my $fusions, '<', $in;
  my @cols = qw(type1	type2	type3	chrA	posA	oriA	geneA	exon_intronA	chrB	posB	oriB	geneB	exon_intronB	event_type	mechanism	event_id	disc_pair	split_read	homology	partners);
  my @lines;
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

    my (@var_parts ) = @{ $vars{$parts[15]}};

    # say "$_ $var_parts[$_]" for 0..$#var_parts;
    # say "---";

    my $size = $var_parts[8];
    $size = sprintf("%.1f", abs($size)/1000);

    $size = 'NA' if $event eq 'TRA';

    push @lines, join("\t", "Meerkat", $event, $parts[3], $parts[4], $parts[8], $parts[9], $parts[17], $parts[16], $parts[15], $size, "$parts[3]:$parts[4]-$parts[9]", "-", $parts[18], $config, $parts[13], "-", $var_parts[1] );
    # print join("\t", "Meerkat", $event, $parts[3], $parts[4], $parts[8], $parts[9], $parts[17], $parts[16], $parts[15], $size, "$parts[3]:$parts[4]-$parts[9]", $parts[18], $config, $parts[13], $var_parts[1] ) . "\n";

  }
  return(\@lines);
}


# my %var_heads;
#
# getHeads($var_hads);
# # print Dumper \%var_heads;
#
# sub getHeads {
#   my $in = shift;
#   while(<$in>){
#     chomp;
#     my $type = (split)[0];
#     my @parts = grep { /\S/ } split(/\t/,$_); # ignore missing cells
#     $var_heads{$type} = [ @parts[1..$#parts] ];
#   }
#   return(\%var_heads);
# }


# source	type	chromosome1	bp1	chromosome2	bp2	split reads	pe reads	id	length(Kb)	position	consensus_microhomology_length	configuration	read_depth_ratio	read_depth_evidence
# source	type	chromosome1	bp1	chromosome2	bp2	split_reads	pe_reads	id	length(Kb)	position	consensus_microhomology_length	configuration	mechanism	read_depth_evidence
