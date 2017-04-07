#!/urs/bin/perl

use strict;
use warnings;

use SV_parser;

use 5.18.2;

use FindBin;
use FindBin '$Script';

use File::Spec;
use lib File::Spec->catdir($FindBin::Bin, '..', 'bin/');

use feature qw/ say /;
use Data::Dumper;
use Getopt::Long qw/ GetOptions /;

use File::Basename;

my $vcf_file; 
my $help;
my $id;
my $dump;
my $filter;
my $chromosome;

#
my $tumour_read_support = 4;

# Change to write to cwd by default (-o . is messy and confusing)
my $output_dir;

GetOptions( 'vcf=s'	        	=>		\$vcf_file,
			'id=s'				=>		\$id,
			'dump'				=>		\$dump,
			'filter'			=>		\$filter,
            'output_dir=s'     	=>      \$output_dir,
			'chromosome=s'		=>		\$chromosome,
			'help'              =>      \$help
	  ) or die usage();

if ($help) { exit usage() } 

if (not $vcf_file) {
	 exit usage();
} 

# Need to make sure this is stable
$output_dir =~ s!/*$!/! if $output_dir; # Add a trailing slash

my ($name, $extention) = split(/\.([^.]+)$/, basename($vcf_file), 2);

# Retun SV and info hashes 
my ( $SVs, $info, $filtered_vars ) = SV_parser::typer($vcf_file);

# Print all info for specified id

SV_parser::summarise_variants( $SVs, $filter, $chromosome ) unless $id or $dump;

# Print all infor for specified id

SV_parser::get_variant( $id, $SVs, $info, $filter ) if $id;

# Dump all variants to screen
SV_parser::dump_variants( $SVs, $info, $filter, $chromosome ) if $dump;

SV_parser::print_variants ( $SVs, $filtered_vars, $name, $output_dir ) if $output_dir;

sub usage {
	say "********** $Script ***********";
    say "Usage: $Script [options]";
	say "  --vcf = VCF file for parsing";
	say "  --id = extract information for a given variant";
	say "  --dump = cycle through all variants (can be combined with both -f and -c)";
	say "  --filter = apply filters and mark filtered variants";
	say "  --output = write out variants that pass filters to specified dir";
	say "  --chromosome = used in conjunction with --dump will cycle though variants on chromosome speciified in -c";
	say "  --help\n";
	
	say "Examples: ";
	say "o Browse all variants that passed filter within a speicifc window on X chromosome:"; 

	say "->  perl script/sv_parse.1.0.pl -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -d -c X:3000000-3500000";
	say "o Filter vars and write to file in cwd:";
	say "->  perl $0 -v data/HUM-7.tagged.SC.lumpy.gt_all.vcf -f -o .\n";	
	say "Nick Riddiford 2017";
}