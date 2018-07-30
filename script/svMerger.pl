#!/usr/bin/perl
use strict;
use warnings;

use 5.18.2;

use File::Basename;
use File::Path qw(make_path);
use FindBin qw($Bin);
use feature qw/ say /;
use Data::Dumper;
use Getopt::Long qw/ GetOptions /;
use Cwd;
use autodie;

my @files;
my $help;
my $output_dir;

# foo=s{1,} indicates one or more values
GetOptions( 'files=s{1,}'   =>   \@files,
            'output_dir=s'  =>   \$output_dir,
            'help'          =>   \$help
          ) or die usage();

if ($help) { exit usage() }

if (@files == 0){
  say "Exiting. Must list at least 1 file";
  exit usage();
}

eval { make_path($output_dir) };
if ($@) {
  print "Couldn't create $output_dir: $@";
}

my ($name, $extention) = split(/\.([^.]+)$/, basename($files[0]), 2);
$name = (split /\./, $name)[0];
my $outfile =  $name . '_merged_SVs.txt';
my $outpath = File::Spec->catdir( $output_dir, $outfile );
open my $out, '>', $outpath;
say "Writing merged files to: '$outpath'";
say "Merging files: ";

my %SVs;
my $header;
foreach (@files){
  say "* $_";
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
  for my $bp1 (sort { $a <=> $b } keys %{ $SVs{$chr1} }){
    for my $chromosome2 (sort keys %{ $SVs{$chr1}{$bp1} }){
      for my $bp2 (sort keys %{ $SVs{$chr1}{$bp1}{$chromosome2} }){
        print $out "$_\n" foreach @{ $SVs{$chr1}{$bp1}{$chromosome2}{$bp2} };
      }
    }
  }
}

sub usage {
  print
"
usage: perl svMerger.pl [-h] [-f FILE]

svMerger
author: Nick Riddiford (nick.riddiford\@curie.fr)
version: v0.1
description: Merge summary output of SV calls from svParser

arguments:
  -h, --help            show this help message and exit
  -f FILE, --file
                        All summary files to be merged for a sample [required]
"
}
