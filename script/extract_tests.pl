#!/usr/bin/perl
use feature qw/ say /;
use Data::Printer;
use autodie;

sub getData {
  my ($file) = @_;
  open my $in, '<', $file;
  my %tps;
  say "Getting data";
  while(<$in>){
    chomp;
    next if /sample/;
    next if /^\s+/;
    my ($sample, $id, $source, $type, $ch1, $bp1, $ch2, $bp2, $genotype) = (split)[0..7, 10];
    next unless $genotype eq 'somatic_tumour';
    push @{$tps{$sample}{$source}{$id}}, $type, $ch1, $bp1, $ch2, $bp2;
  }
  return(\%tps);
}

sub printFiles{
  my ($tps) = @_;
  my %data = %{ $tps };
  say "Printing data";
  for my $s (keys %data){
    for my $c (keys $data{$s}){
      open my $out, '>', join("_", $s, $c, 'somatic.txt');
      for my $id (keys $data{$s}{$c}){
        my @parts = @{$data{$s}{$c}{$id}};
        print $out join("\t", $s, $c, @parts) . "\n";
      }
    }
  }
}

my $tps = getData($ARGV[0]);
printFiles($tps);



__DATA__
sample	event	source	type	chromosome1	bp1	chromosome2	bp2	split_reads	disc_reads	genotype	id	length(Kb)	position	consensus|type	microhomology	configuration	allele_frequency	mechanism|log2(cnv)	bp1_locus	bp2_locus	affected_genes	T/F	notes

A373R11	1	lumpy	DEL	2L	73837	2L	73896	7	0	germline_recurrent	6	0.1	2L:73837-73896	-	-	-	0.39	-	galectin, exon_2	galectin, exon_2	galectin
A373R11	2	lumpy	DEL	2L	116597	2L	116755	2	0	germline_recurrent	7	0.2	2L:116597-116755	-	-	-	0.28	-	ND-15, intron	ND-15, intron	ND-15
A373R11	3	lumpy	DEL	2L	263351	2L	263445	9	0	germline_recurrent	9	0.1	2L:263351-263445	-	-	-	0.48	-	CG17075, 3UTR	CG17075, 3UTR	CG3645, CG17075
A373R11	4	lumpy	INV	2L	496301	2L	496361	7	3	somatic_tumour	19	0.1	2L:496301-496361	-	-	-	0.28	-	ush, intron	ush, intron	ush
A373R11	5	lumpy	DUP	2L	605109	2L	606021	1	11	germline_recurrent	1	0.9	2L:605109-606021	-	-	-	0.43	-	intergenic	CR43648, ncRNA	CR43648
A373R11	6	lumpy	DEL	2L	795313	2L	795369	12	0	somatic_tumour	29	0.1	2L:795313-795369	-	-	-	0.33	-	intergenic	intergenic	-
A373R11	7	lumpy	DEL	2L	892451	2L	892509	4	0	germline_recurrent	30	0.1	2L:892451-892509	-	-	-	0.44
