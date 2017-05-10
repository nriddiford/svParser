#!/usr/bin/perl
use strict;
use warnings;

use File::Basename;

use feature qw/ say /;

use Data::Printer;
use Data::Dumper;

die "Must list at least 1 file\n" unless @ARGV > 0;

my @files = @ARGV;

my @parts = split(/\./, $files[0]);

open my $out, '>', $parts[0] . "_merged_SVs.tsv" or die $!;

my %SVs;
my %seen;
my $header;
foreach (@files){
  open my $in, '<', $_ or die $!;

  while(<$in>){
    chomp;
    if ($. == 1){
      $header = $_;
      next;
    }
    my ($source, $chromosome1, $bp1, $chromosome2, $bp2) = (split)[ 0,2,3,4,5 ];

    push @{$SVs{$chromosome1}{$bp1}{$bp2}{$chromosome2}}, $_;

  }
}

print $out "$header\n";

for my $chr1 (sort keys %SVs){
  for my $bp1 (sort { $a <=> $b } keys $SVs{$chr1}){
    for my $chromosome2 (sort keys $SVs{$chr1}{$bp1}){
      for my $bp2 (sort keys $SVs{$chr1}{$bp1}{$chromosome2} ){
        print $out "$_\n" foreach @{$SVs{$chr1}{$bp1}{$chromosome2}{$bp2}};
      }
    }
  }
}

__DATA__
source	type	chromosome	bp1	bp2	split reads	pe reads	id	length(Kb)	position	consensus	microhomology length	configuration
lumpy	DEL	2L	2L:5890205	2L:5890720	5	11	287	0.5	2L:5890205-5890720	-	-	-
lumpy	BND	2L	2L:22220720	2L:22255744	0	7	1062_1	35.0	2L:22220720 2L:22255744	-	-	N]2L:22255744]
lumpy	INV	3L	3L:15568694	3L:15568866	0	5	3083	0.2	3L:15568694-15568866	-	-	-
lumpy	DEL	3R	3R:14006281	3R:14008253	7	10	3810	2.0	3R:14006281-14008253	-	-	-
lumpy	BND	3R	3R:32060908	3R:32061196	0	16	4665_1	0.3	3R:32060908 3R:32061196	-	-	[3R:32061196[N
lumpy	BND	3R	3R:32066206	3R:32068392	11	26	4666_1	2.2	3R:32066206 3R:32068392	-	-	N]3R:32068392]
lumpy	DEL	X	X:4574313	X:4576607	4	5	4781	2.3	X:4574313-4576607	-	-	-
lumpy	BND	X	X:11414158	X:11531996	0	5	4910_1	117.8	X:11414158 X:11531996	-	-	N]X:11531996]
