#!/usr/bin/perl

use strict;
use warnings;
use feature qw/ say /;
use autodie;

use Data::Dumper;
use Data::Printer;

my @files = @ARGV;

p(@files);

open my $all, '>', 'all_samples.txt';

open my $file1, '<', $ARGV[0];
my $header = <$file1>;
close $file1;

print $all "sample" . "\t" . $header;

my %seen;
my @lines;
for (@files) {
  say;
  open my $in, '<', $_ or die "$!";
  my ($sample) = split(/\_/);
  say "Sample" if $sample eq 'A573R25';

    while(<$in>){
      chomp;
      say if $sample eq 'A573R25';

      next if /^event/;
      push @lines, "$sample\t$_" unless $seen{$_}++;
    }
}

print $all "$_\n" foreach @lines;
