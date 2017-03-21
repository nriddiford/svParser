#!/urs/bin/perl

use strict;
use warnings;

use FindBin;
use FindBin '$Script';

use File::Spec;
use lib File::Spec->catdir($FindBin::Bin, '..', 'bin/SV_parser');

use SV_parser;

use feature qw/ say /;
use Data::Dumper;
use Getopt::Long qw/ GetOptions /;

my $vcf_file; 
my $help;
my $id;
my $dump;
my $filter;
my $chromosome;
my $print;

# Should add score threshold option
GetOptions( 'vcf=s'	        	=>		\$vcf_file,
			'id=s'				=>		\$id,
			'dump'				=>		\$dump,
			'filter'			=>		\$filter,
			'print'				=>		\$print,
			'chromosome=s'		=>		\$chromosome,
			'help'              =>      \$help
	  ) or die usage();

if ($help) { exit usage() } 

if (not $vcf_file) {
	 exit usage();
} 


# Retun SV and info hashes 
my ( $SVs, $info, $filtered_vars ) = SV_parser::typer($vcf_file);

# Print all infor for specified id


SV_parser::summarise_variants( $SVs, $filter ) unless $id or $dump;

# Print all infor for specified id

SV_parser::get_variant( $id, $SVs, $info, $filter ) if $id;

# Dump all variants to screen
SV_parser::dump_variants( $SVs, $info, $filter, $chromosome ) if $dump;

SV_parser::print_variants ( $SVs, $filtered_vars ) if $print;

sub usage {
	say "********** $Script ***********";
    say "Usage: $Script [options]";
	say "  --vcf = VCF file for parsing";
	say "  --id = extract information for a given variant";
	say "  --dump = cycle through all variants (can be combined with both -f and -c)";
	say "  --filter = apply filters and mark filtered variants";
	say "  --print = write out variants that pass filters";
	say "  --chromosome = used in conjunction with --dump will cycle though variants on chromosome speciified in -c";
	say "  --help\n";
	
	say "Examples: ";
	say "  Print to screen all vars on chromosome '2L': perl $0 -v [file.vcf] -d -c 2L";
	say "  Print to screen all vars on chromosome '2L' within range 100000-200000: perl $0 -v [file.vcf] -d -c 2L:100000-200000";
	say "  Print to screen all vars that pass the fitleres on chromosome '2L' within range 100000-200000: perl $0 -v [file.vcf] -f -d -c 2L:100000-200000";
	say "  Filter vars and write to file: perl $0 -v [file.vcf] -f -p\n";	
	say "Nick Riddiford 2017";
}