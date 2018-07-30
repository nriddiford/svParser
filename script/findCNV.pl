#!/usr/bin/perl
use feature qw/ say /;
use Data::Printer;
use autodie;
use File::Basename;
use FindBin qw/ $Script /;

use Getopt::Long qw/ GetOptions /;

my $cnv_file;
my $output_dir;
my $help;
my @variants_files;

GetOptions( 'cnvs=s'        =>   \$cnv_file,
            'variants=s{1,}'   =>    \@variants_files,
            'output_dir=s'  =>   \$output_dir,
            'help'          =>   \$help
) or die usage();

if (not $cnv_file) {
   exit usage();
}

sub getFreec {
  my ($file, $c1, $bp1, $bp2) = @_;
  open my $in, '<', $file;
  my ($name, $extention) = split(/\.([^.]+)$/, basename($file), 2);
  say $name;
  my %cnvs;
  say "Getting CNV";
  my ($r, $i) = (0,0);

  while(<$in>){
    chomp;
    next if /Chromosome/;
    next if /^\s+/;

    my ($chr, $start, $ratio, $medratio, $cn) = (split)[0..4];
    if($c1 eq $chr and $bp1 <= $start and $bp2 >= $start){
      $r += $ratio;
      say $ratio;
      $i++;
    }
  }
  return($r, $i);
}

sub readCNV {
  my ($file) = shift;
  open my $in, '<', $file;
  my ($name, $extention) = split(/\.([^.]+)$/, basename($file), 2);
  my %cnvs;
  say "Getting CNVs from '$name'";

  while(<$in>){
    chomp;
    next if /chromosome/;
    next if /^\s+/;
    my ($chr, $start, $end, $test, $ref, $p, $log2) = (split)[0..6];
    $chr =~ s/\"//g;
    next if $log2 eq 'NA' or $log2 =~ /inf/i;
    $cnvs{$chr}{$start} = [$end, $test, $ref, $p, $log2];
  }
  return(\%cnvs);
}


sub findCNV {
  my($d, $c1, $bp1, $bp2) = @_;

  my %cnvs = %{ $d };
  my ($r, $i) = (0, 0);
  for my $st (keys %{$cnvs{$c1}}){
    my ($end, $test, $ref, $p, $log2) = @{$cnvs{$c1}{$st}};
    if($bp1 <= $st and $bp2 >= $end){
      $r += $log2;
      $i++;
    }
  }
  return($r, $i);
}

sub readParser {
  my ($d, $variants_file) = @_;
  my ($name, $extention) = split(/\.([^.]+)$/, basename($variants_file), 2);
  my $outfile =  $name . '.cnv.txt';
  my $outpath = File::Spec->catdir( dirname($variants_file), $outfile );
  open my $in, '<', $variants_file;
  open my $out, '>', $outpath or die $!;
  my @header = qw/ source type chromosome1 bp1 chromosome2 bp2 split_reads disc_reads genotype id length(Kb) position consensus microhomology configuration allele_frequency log2(cnv) /;
  while(<$in>){
    chomp;
    next if /^source/;
    my @parts = split;
    my ($chr, $bp1, $bp2, $genotype) = @parts[2,3,5,8];
    my ($r, $i, $av_r) = (0, 0, 0);

    if($genotype =~ 'somatic'){
      ($r, $i) = findCNV($d, $chr, $bp1, $bp2);
      $av_r = sprintf("%.2f", $r/$i);
    }
    splice(@parts, 16, 1, $av_r);
    print $out join("\t", @parts) . "\n";
    # say "Average log2 FC: $av_r over $i windows";
  }
}

my $d = readCNV($cnv_file);
say "Read file '$cnv_file'";
foreach(@variants_files){
  readParser($d, $_)
}

sub usage {
  print "usage: perl $Script -c cnvs -o output_dir\n";
}
