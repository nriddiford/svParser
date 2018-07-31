#!/usr/bin/perl
use strict;
use warnings;
use feature qw/ say /;
use autodie;
use File::Basename;
use File::Path qw/ make_path / ;
use Data::Printer;
use Getopt::Long qw/ GetOptions /;
use FindBin qw/ $Script /;

my $variants;
my $output_dir;
my $blacklist;
my $whitelist;
my $help;

GetOptions( 'variants=s'    =>   \$variants,
            'output_dir=s'  =>   \$output_dir,
            'blacklist=s'   =>   \$blacklist,
            'whitelist=s'   =>   \$whitelist,
            'help'          =>   \$help
) or die usage();

if ($help) { exit usage() }

if (not $variants or not $output_dir) {
   exit usage();
}

sub get_vars {
  my ($variants) = @_;

  open my $sv_file, '<', $variants;
  # This allows files to be parsed irrespective of line ending. Not recommened for large files
  my @lines;
  {
  local $/ = undef;
  my $content = <$sv_file>;
  @lines = split /\r\n|\n|\r/, $content;
  }
  return(\@lines);
}


sub mark_events{
  my ($outfile, $output_dir, $tp_fh, $fp_fh, $sample, $l, $tp, $fp) = @_;

  my @lines       = @{ $l };
  my %false_calls = %{ $fp };
  my %true_calls  = %{ $tp };
  my %events;

  my $outpath = File::Spec->catdir( $output_dir, $outfile );
  open my $cleaned_file, '>', $outpath;

  my $filtered_calls = 0;
  my $whitelisted_calls = 0;
  my $new_fp = 0;
  my $new_tp = 0;

  for (@lines){
    chomp;
    my @cells = split(/\t/);

    if ($cells[0] eq 'event'){
      print $cleaned_file $_ . "\n";
      next;
    }

    my $match = join("_", $sample, @cells[3..6]);

    my $notes;
    if (@cells >= 21){
      $notes = $cells[21];
    }
    $notes = '-' if !defined $notes || $notes eq '';

    if ($notes eq 'F') {
      $filtered_calls++;
      if (not $false_calls{$match}){
        say "* Adding new false call to blacklist: $match";
        $new_fp++;
        print $fp_fh $match . "\n";
      }
    }
    elsif ($false_calls{$match}) {
      say "* Unmarked false call: $match";
      $filtered_calls++;
    }

    elsif ($notes eq 'T') {
      if ( not $true_calls{$match} ){
        say "* Added new call to whitelist: $match";
        print $tp_fh $match . "\n";
        $new_tp++;
        $whitelisted_calls++;
      }
      else{
        $whitelisted_calls++;
      }
    }

    if ($notes ne 'F'){
      print $cleaned_file $_ . "\n";
    }
  }

  if ($filtered_calls>0){
    printf "Removed %s calls marked as false positives from %s [%s new]\n", $filtered_calls, $sample, $new_fp;
  }
  if ($whitelisted_calls>0){
    printf "%s calls in whitelist in %s [%s new]\n", $whitelisted_calls, $sample, $new_tp;
  }
}


sub get_TP {
  my $tp_file = shift;
  open my $true_positives, '<', $tp_file if -e $tp_file;
  my %true_calls;
  while(<$true_positives>){
    chomp;
    my $id = (split)[0];
    $true_calls{$id}++;
  }
  return(\%true_calls);
}


sub get_FP {
  my $fp_file = shift;
  open my $false_positives, '<', $fp_file if -e $fp_file;
  my %false_calls;
  while(<$false_positives>){
    chomp;
    $false_calls{$_}++;
  }
  return(\%false_calls);
}


sub makeFiles {
  my ($vars, $out_dir, $true_calls_file, $false_calls_file, $blacklist, $whitelist) = @_;

  my ( $name, $extention ) = split(/\.([^.]+)$/, basename($vars), 2);
  my ($sample) = (split /_/, $name)[0];
  my $outfile = $sample . "_" . "cleaned_SVs.txt";
  my $tp_out = File::Spec->catdir( $out_dir, $true_calls_file );
  my $fp_out = File::Spec->catdir( $out_dir, $false_calls_file );

  my $fp_fh;
  my $tp_fh;
  if($blacklist){
    say "Ammending fp calls to $fp_out";
    open $fp_fh, '>>', $fp_out;
  }
  else{
    say "Writing fp calls to $fp_out";
    open $fp_fh, '>', $fp_out;
  }

  if($whitelist){
    say "Ammending tp calls to $tp_out";
    open $tp_fh, '>>', $tp_out;
  }
  else{
    say "Writing tp calls to $tp_out";
    open $tp_fh, '>', $tp_out;
  }

  return($outfile, $tp_fh, $fp_fh, $sample)
}

my $tp = {'dummy' => 1};
my $fp = {'dummy' => 1};

my ($true_calls_file, $false_calls_file);
if ($whitelist){
  $tp = get_TP($whitelist);
  $true_calls_file = basename($whitelist);
}
else{
  $true_calls_file = 'all_samples_whitelist.txt';
}

if ($blacklist){
  $fp = get_TP($blacklist);
  $false_calls_file = basename($blacklist);
}
else{
  $false_calls_file = 'all_samples_blacklist.txt';
}

my ($out_file, $tp_fh, $fp_fh, $sample) = makeFiles($variants, $output_dir, $true_calls_file, $false_calls_file, $blacklist, $whitelist );
my ($l) = get_vars( $variants );
mark_events($out_file, $output_dir, $tp_fh, $fp_fh, $sample, $l, $tp, $fp);

sub usage {
  print
"
usage: $Script -h -v variants -o out_dir

arguments:
  -h, --help            show this help message and exit
  -v, --variants        variant calls file (as produced by svParser)[required]
  -o, --output_dir      output directory
  -b, --blacklist       specify blacklist to exclude variants
  -w, --whitelist       specify whitelist to proect variants
"
}
