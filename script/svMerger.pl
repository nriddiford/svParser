#!/usr/bin/perl
use strict;
use warnings;

use 5.18.2;

use File::Basename;
use File::Path qw(make_path);
use FindBin qw($Bin);
use feature qw/ say /;
use Getopt::Long qw/ GetOptions /;
use Cwd;
use Data::Printer;
use autodie;

my @files;
my $help;
my $output_dir;

# foo=s{1,} indicates one or more values
GetOptions( 'files=s{1,}'   =>   \@files,
            'output_dir=s'  =>   \$output_dir,
            'help'          =>   \$help
          ) or die usage();

if (scalar @files == 0) { exit usage() }
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
    my ($source, $chromosome1, $bp1, $chromosome2, $bp2) = (split)[ 0,2..5 ];
    push @{$SVs{$chromosome1}{$bp1}{$chromosome2}{$bp2}}, $_;
  }
}

print $out "$header\n";

my %chr1_seen;
my %chr2_seen;
for my $chr1 (sort keys %SVs){
  for my $bp1 (sort { $a <=> $b } keys %{ $SVs{$chr1} }){
    my $b1_key = "$chr1\_$bp1";
    $chr1_seen{$b1_key}++;
    for my $chr2 (sort keys %{ $SVs{$chr1}{$bp1} }){
      for my $bp2 (sort { $a <=> $b } %{ $SVs{$chr1}{$bp1}{$chr2} }){
        my $b2_key = "$chr2\_$bp2";
        if( $chr1_seen{$b2_key} ){
          say "Duplicated entry: $chr1:$bp1    $chr2:$bp2";
          p(@{ $SVs{$chr1}{$bp1}{$chr2}{$bp2} });
          next;
        }
        print $out "$_\n" foreach @{ $SVs{$chr1}{$bp1}{$chr2}{$bp2} };
      }
    }
  }
}

sub usage {
  print "usage: perl svMerger.pl -f file1 .. fileN\n";
}
