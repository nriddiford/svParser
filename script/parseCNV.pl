#!/usr/bin/perl
use strict;
use warnings;
use feature qw/ say /;
use Data::Printer;
use Data::Dumper;
use autodie;

use File::Basename;
use FindBin qw($Bin);
use FindBin qw/ $Script /;

use Getopt::Long qw/ GetOptions /;

my $cnvs;
my $output_dir;
my $help;

GetOptions( 'cnvs=s'        =>   \$cnvs,
            'output_dir=s'  =>   \$output_dir,
            'help'          =>   \$help
) or die usage();

if (not $cnvs and $output_dir) {
   exit usage();
}

my $cnv_ref = extractVars($cnvs);

my @name_fields = split( /\_/, basename($cnvs) );

my $outfile =  $name_fields[0] . ".cnvseq.filtered.summary.txt";
my $outpath = File::Spec->catdir( $output_dir, $outfile );

open my $out, '>', $outpath;

my @header = qw/ source type chromosome1 bp1 chromosome2 bp2 split_reads disc_reads genotype id length(Kb) position consensus microhomology configuration allele_frequency log2(cnv) /;

print $out join("\t", @header) . "\n";

my @lines = @{$cnv_ref};
print $out "$_\n" foreach @lines;


sub extractVars {
  my $in = shift;
  open my $annotated_cnvs, '<', $in;
  my @cnv;
  while(<$annotated_cnvs>){
    chomp;
    my ($chrom, $start, $stop, $type, $length, $fc) = split;

    push @cnv, join("\t", "CNV-Seq",               # source
                            $type,                 # type
                            $chrom,                # chrom1
                            $start,                # bp1
                            $chrom,                # chrom2
                            $stop,                 # bp2
                            '-',                   # split reads
                            '-',                   # paired reads
                            'somatic_tumour',      # genotype
                            '-',                   # id
                            $length,               # length
                            "$chrom:$start-$stop", # IGV
                            '-',                   # misc1 (type)
                            '-',                   # microhomology
                            '-',                   # configuration
                            '-',                   # allele frequency
                            $fc);                  # misc (cnv)

  }
  return(\@cnv);
}

sub usage {
  print "usage: perl $Script -c cnvs -o output_dir\n";
}
