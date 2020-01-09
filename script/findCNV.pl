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
my $locus;

GetOptions( 'cnvs=s'           =>   \$cnv_file,
            'variants=s{1,}'   =>   \@variants_files,
            'locus=s'          =>   \$locus,
            'output_dir=s'     =>   \$output_dir,
            'help'             =>   \$help
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

sub read_CNV {
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
  my ($r, $i) = (0, 1);
  for my $st (keys %{$cnvs{$c1}}){
    my ($end, $test, $ref, $p, $log2) = @{$cnvs{$c1}{$st}};
    if($bp1 <= $st and $bp2 >= $end){
      $r += $log2;
      $i++;
    }
  }
  return($r, $i);
}

sub read_parser {
  my ($d, $variants_file) = @_;
  my ($name, $extention) = split(/\.([^.]+)$/, basename($variants_file), 2);
  my $outfile =  $name . '.cnv.txt';
  my $outpath = File::Spec->catdir( dirname($variants_file), $outfile );
  open my $in, '<', "$variants_file";
  open my $out, '>', $outpath or die $!;
  my @header = qw/ source type chromosome1 bp1 chromosome2 bp2 split_reads disc_reads genotype id length(Kb) position consensus microhomology configuration allele_frequency mechanism log2(cnv) status	notes /;
  print $out join("\t", @header) . "\n";

  while(<$in>){
    chomp;
    next if /^source/;
    my @parts = split(/\t/);
    my ($type, $chr1, $bp1, $chr2, $bp2, $genotype) = @parts[1..5,8];

    my ($r, $i, $av_r) = (0, 0, 0);

    if( $chr1 eq $chr2 and $genotype =~ 'somatic' ){
      ($r, $i) = findCNV($d, $chr1, $bp1, $bp2);
      $av_r = sprintf("%.2f", $r/$i);
    }
    splice(@parts, 17, 1, $av_r);
    print $out join("\t", @parts) . "\n";
  }

}


sub get_locus {
  my $locus = shift;
  die "Error parsing the specified region.\nPlease specify chromosome regions using the folloing format:\tchrom:start-stop\n" if ($locus !~ /-/);

  my ($chromosome, $query_region) = split(/:/, $locus);
  my ($query_start, $query_stop) = split(/-/, $query_region);
  die "Error: End of window must be larger than start:\tchrom:start-stop\n" if ($query_stop - $query_start < 0);

  $query_start =~ s/,//g;
  $query_stop =~ s/,//g;
  say "Searching within region: '$chromosome:$query_start-$query_stop'";
  return($chromosome, $query_start, $query_stop);
}

sub get_depth{
  my ($d, $l) = @_;
  my ($chromosome, $query_start, $query_stop) = get_locus($l);
  ($r, $i) = findCNV($d, $chromosome, $query_start, $query_stop);
  $av_r = sprintf("%.2f", $r/$i);
  say "Average log2(FC) over region '$l': $av_r";
}

my $d = read_CNV($cnv_file);
say "Read file '$cnv_file'";
if($locus){
  get_depth($d, $locus);
  exit(1);
}

foreach(@variants_files){
  read_parser($d, $_)
}

sub usage {
  print "usage: perl $Script -c cnvs -o output_dir\n";
}
