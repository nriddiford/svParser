#!/usr/bin/env perl

use strict;
use warnings;

use autodie;
use Data::Printer;
use File::Basename;
use File::Path qw(make_path);
use feature qw/ say /;
use Cwd;
use FindBin qw($Bin);
use FindBin qw/ $Script /;
use Getopt::Long qw/ GetOptions /;

my $variants;
my $distance = 50;
my $output_dir = getcwd;
my $help;

GetOptions( 'variants=s'     =>   \$variants,
            'distance=s'     =>   \$distance,
            'output_dir=s'   =>   \$output_dir,
            'help'           =>   \$help
) or die usage();

if (not $variants) {
   exit usage();
}

eval { make_path($output_dir) };
if ($@) {
  print "Couldn't create $output_dir: $@";
}

my ($name, $extention) = split(/\.([^.]+)$/, basename($variants), 2);
$name = (split /\./, $name)[0];
my $outfile =  $name . "_clustered_SVs.txt";
my $outpath = File::Spec->catdir( $output_dir, $outfile );

open my $out, '>', $outpath;
say "Writing clustered sv file to: " . $outpath;

my $sample_count;

run($variants, $distance, $name);

sub run {
    my ($file, $distance, $name) = @_;

    say "Clustering variants +/- $distance bps";
    open my $fh, '<', $file;

    my @header = split(/\t/, scalar <$fh>);
    chomp(@header);

    my @events = ([]);

    while (<$fh>) {
      chomp;
      $sample_count++;
      my %event;
      @event{ @header } = split(/\t/);
      # $events[-1] last element of array

      unless (@{ $events[-1] }) {
          push @{ $events[-1] }, \%event;
          next;
      }

      if (same_event($events[-1][-1], \%event, $distance)) {
          push @{ $events[-1] }, \%event;
          next;
      }
      push @events, [ \%event ];
    }

    print $out join("\t", event => @header), "\n";

    say $sample_count . " calls were grouped into " . scalar @events . " events for " . $name;

    for my $i (1 .. @events) {
        for my $ev (@{ $events[$i - 1] }) {
            print $out join("\t", $i, @{$ev}{@header}), "\n";
        }

    }
}

sub same_event {
    my ($x, $y, $threshold) = @_;
    return if $x->{genotype} =~ 'germline';
    return if $x->{chromosome1} ne $y->{chromosome1};
    return if abs($x->{bp1} - $y->{bp1}) > $threshold;
    return if abs($x->{bp2} - $y->{bp2}) > $threshold;
    return 1;
}

sub usage {
  print "usage: perl $Script -v variants -d distance -o output_dir\n";
}
