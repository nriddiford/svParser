#!/urs/bin/perl

package SV_parser;
use strict;
use warnings;

use feature qw/ say /;
use Data::Dumper;

sub typer {
	my $file = shift;
	my $type;
	if (`grep "source=LUMPY" $file`){
		say "Recognised $file as Lumpy input";
		$type = 'lumpy';
		parse($file, $type);
	}
	
	elsif (`grep "DELLY" $file`){
		say "Recognised $file as Delly input";
		$type = 'delly';
		parse($file, $type);
	}
	# elsif (`grep "VarScan" $file`){
	# 	say "Recognised $file as VarScan2 input";
	# 	$type = 'varscan2';
	# 	parse($file, $type);
	# }
	
	else {
		die "This VCF is not from lumpy or delly. Abort";
	}
}

sub parse {
	my ($file, $type) = @_;
	open my $in, '<', $file or die $!;

	my @headers;
	
	my (%SVs, %info, %filtered_SVs);
	
	my ($tumour_name, $control_name);
	
	my %format_long;

	my %info_long;
	my $filter_count;

	my @samples;
	
	while(<$in>){
		chomp;

		if (/^#{2}/){
	 		push @headers, $_;
			$filtered_SVs{$.} = $_;

			if (/##FORMAT/){
				
				my ($format_long) = $_ =~ /\"(.*?)\"/;
				my ($available_format_info) = $_ =~ /ID=(.*?),/;
				$format_long{$available_format_info} = $format_long;
			}

			if (/##INFO/) {
				my ($info_long) = $_ =~ /\"(.*?)\"/;
				my ($available_info) = $_ =~ /ID=(.*?),/;
				$info_long{$available_info} = $info_long;
			}
			next;
		}

		if (/^#{1}/){
	 		push @headers, $_;
			$filtered_SVs{$.} = $_;
			my @split = split;
			push @samples, $_ foreach @split[9..$#split];
 			
			$tumour_name = $samples[0];
			$control_name = $samples[1];
			next;
		}

		my @fields = split;
				
	    my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, @sample_info) = @fields;

 		my %sample_parts;

 		push @{$sample_parts{$samples[$_]}}, split(/:/, $sample_info[$_]) for 0..$#samples;
		
		my @tumour_parts 	= split(/:/, $sample_info[0]);
	    my @normal_parts 	= split(/:/, $sample_info[1]);
						
		my @format 		 	= split(/:/, $format_block);
		my @info_parts		= split(/;/, $info_block);
		
		my %sample_info;
		
		for my $sample (@samples){

			for my $info (0..$#format){
				$sample_info{$id}{$sample}{$format[$info]} = $sample_parts{$sample}[$info];
			}
			
		}
								
 		my @normals = @samples[1..$#samples];
		
		my @filter_reasons;
		
		###################
		# Genotype filter #
		###################
		
		# Filter if ANY of the controls are NOT 'GT = 0/0' (hom ref) 
		
		for my $normal (@normals){
			if ($sample_info{$id}{$normal}{'GT'} eq '1/1' or $sample_info{$id}{$normal}{'GT'} eq '0/1'){
				push @filter_reasons, "$normal\_not_homo_ref=" . $sample_info{$id}{$normal}{'GT'};
			}
		}
		
		my %information;
		
					
		foreach(@info_parts){
			
			my ($info_key, $info_value);

			if (/=/){
				($info_key, $info_value) = $_ =~ /(.*)=(.*)/;
			}

			else {
				($info_key) = $_ =~ /(.*)/;
				$info_value = "TRUE";
			}
			$information{$id}{$info_key} = $info_value;
		}
				
	    my ($SV_type) = $info_block =~ /SVTYPE=(.*?);/;
				
		my ($SV_length, $chr2, $stop, $t_SR, $t_PE, $filters);
		
		if ($type eq 'lumpy'){
			( $SV_length, $chr2, $stop, $t_SR, $t_PE, $filters ) = lumpy( $id, $info_block, $SV_type, $alt, $start, \%sample_info, $tumour_name, $control_name, \@samples, \@normals, \@filter_reasons );
		}
		
		elsif ($type eq 'delly'){
			( $SV_length, $chr2, $stop, $t_SR, $t_PE, $filters ) = delly( $info_block, $start, $SV_type, \@filter_reasons );
		}
		
		$filters = chrom_filter( $chr, $chr2, $filters );
				
		$SVs{$id} = [ @fields[0..10], $SV_type, $SV_length, $stop, $chr2, $t_SR, $t_PE, $filters, \@samples ];
		
		$info{$id} = [ [@format], [%format_long], [%info_long], [@tumour_parts], [@normal_parts], [%information], [%sample_info] ];
		
		if (scalar @{$filters} == 0){
			$filtered_SVs{$.} = $_;
		}

			
	}	
	return (\%SVs, \%info, \%filtered_SVs);
}

sub print_variants {
	open my $out, '>', 'filtered_vars.vcf' or die $!;
	
	my ( $SVs, $filtered_SVs ) = @_;
	
	my %filtered_SVs = %{ $filtered_SVs };
	my $sv_count = 0;
	for (sort {$a <=> $b} keys %filtered_SVs){
		my $line = $filtered_SVs{$_};
				
		if ($line =~ /^#/){
			print $out $line . "\n"; 
		}
		else {
			$sv_count++;
			my @cols = split("\t", $line);
			print $out join("\t", @cols[0..5], "PASS", @cols[7..$#cols]) . "\n";
		}
	}	
	say "$sv_count variants passed all filters";
}
		
sub lumpy {
	my ( $id, $info_block, $SV_type, $alt, $start, $sample_info, $tumour, $control, $samples, $normals, $filters ) = @_;
	
	my @filter_reasons = @{ $filters };
	my @normals = @{ $normals };
	
	my @samples = @{ $samples };
	
	my %sample_info = %{ $sample_info };
		
		
	my ($SV_length) = $info_block =~ /SVLEN=(.*?);/;
	
			
		my ($t_PE, $t_SR, $all_c_PE, $all_c_SR, $c_PE, $c_SR) = (0,0,0,0,0,0);
		
		########################
		# Read support filters #
		########################
		
		# Get read support for tumour
		$t_PE = $sample_info{$id}{$tumour}{'PE'};
		$t_SR = $sample_info{$id}{$tumour}{'SR'};
		
		my $tumour_read_support = ( $t_PE + $t_SR );
		
		# Create temp pseudo counts to avoid illegal division by 0
		my $sc_tumour_read_support = $tumour_read_support + 0.001;
		
		

		# Anything with less than 4 supporting reads is filtered
		if ( $tumour_read_support < 4 ){
			push @filter_reasons, 'tumour_reads<4=' . $tumour_read_support;

		}
			
		# for precise variants:
		if ($info_block !~ /IMPRECISE;/){
						
			# We want to be strict, so include all controls used for genotyping (and sum read support)
			for my $normal (@normals){
				$sample_info{$id}{$normal}{'PE'} eq '.' ? $sample_info{$id}{$normal}{'PE'} = '0' : $all_c_PE += $sample_info{$id}{$normal}{'PE'};
				$sample_info{$id}{$normal}{'SR'} eq '.' ? $sample_info{$id}{$normal}{'SR'} = '0' : $all_c_SR += $sample_info{$id}{$normal}{'SR'};			
			}
	
			my $all_control_read_support = ( $all_c_PE + $all_c_SR );
						
			# Filter if # tumour reads supporting var is less than 5 * control reads
			# Or if there are more than 2 control reads
			if ( $all_control_read_support > 0 ){
				push @filter_reasons, 'precise_var_with_normal_read_support=' . $all_control_read_support;
			}
		}
		
		# for imprecise variants:
		if ($info_block =~ /IMPRECISE;/){
						
			# Filter if # tumour reads supporting var is less than 5 * control reads
			# Or if there are more than 2 control reads
			
			# Get read support for direct control
			$c_PE =  $sample_info{$id}{$control}{'PE'};
			$c_SR =  $sample_info{$id}{$control}{'SR'};
	
			my $direct_control_read_support = ( $c_PE + $c_SR );
			my $sc_direct_control_read_support = $direct_control_read_support + 0.001;
			
			if ( $direct_control_read_support > 1 ){
				push @filter_reasons, "imprecise_var_with_$control\_read_support>1=" . $direct_control_read_support;
			}
			
		}
		
		##################
		# Quality filter #
		##################
		
		$sample_info{$id}{$tumour}{'SQ'} = 0 if $sample_info{$id}{$tumour}{'SQ'} eq '.';

		if ( $sample_info{$id}{$tumour}{'SQ'} <= 10 ){
			push @filter_reasons, "SQ=" . $sample_info{$id}{$tumour}{'SQ'};
		}

		
		######################
		# Read depth filters #
		######################
		
		my $t_DP =  $sample_info{$id}{$tumour}{'DP'};
		my $c_DP =  $sample_info{$id}{$control}{'DP'};
			
		# Flag if either control or tumour has depth < 10 at site
		if ( $t_DP <= 10 ){
			push @filter_reasons, 'tumour_depth<10=' . $t_DP;
		}

		if ( $c_DP <= 10 ){
			push @filter_reasons, 'control_depth<10=' . $c_DP;
		}
						
		my ($chr2, $stop) = 0,0;

		if ($SV_type eq 'BND'){
			$chr2 = $alt =~ s/[\[\]N]//g;
			($chr2, $stop) = $alt =~ /(.+)\:(\d+)/;
			$SV_length = $stop - $start;
		}
		else {
		    ($stop) = $info_block =~ /;END=(.*?);/;
			
		}
		
	return ($SV_length, $chr2, $stop, $t_SR, $t_PE, \@filter_reasons);
}

sub delly {
	my ($info_block, $start, $SV_type, $filters) = @_;
	
	my @filter_reasons = @{ $filters };
	
    my ($stop) = $info_block =~ /;END=(.*?);/;
		
	my ($SV_length) = ($stop - $start);
	
		my ($t_SR, $t_PE) = (0,0);
		
		if ($info_block =~ /;SR=(\d+);/){
	   		$t_SR = $1;
	   	}
	   
	   if ($info_block =~ /;PE=(\d+);/){
	   		$t_PE = $1;
	   }
 	  		
		if ($start > $stop){
			die "Start bigger than stop - shouldn't be here!!!!";
		}
		
		my ($chr2) = 0;
				
		if ($SV_type eq 'TRA'){
			($chr2) = $info_block =~ /CHR2=(.*?);/;
		}
		
	return ($SV_length, $chr2, $stop, $t_SR, $t_PE, \@filter_reasons );
}

sub summarise_variants {
	my ($SVs, $filter_flag) = @_;
	
	my ($dels, $dups, $trans, $invs, $filtered) = (0,0,0,0);
	
	my %support_by_chrom;
	
	my $read_support;
	
	my %filtered_sv;	
	
	for (keys %{ $SVs } ){
   	   
       my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $filters ) = @{ $SVs->{$_} };
	   
	   my @filter_reasons = @{ $filters };
	   	   	   
		foreach (@filter_reasons){
			my ($reason) = $_ =~ /(.+)=/;
			$filtered_sv{$reason}++ ;
		}
	   	
		if ( scalar @filter_reasons > 0 ){
			$filtered++;
			next if $filter_flag;
		}
			
	   $read_support = ( $SR + $PE );
	   $support_by_chrom{$id} = [ $read_support, $sv_type, $chr ];
	  	   	   	   
	   $dels++ if $sv_type eq 'DEL';
	   $dups++ if $sv_type eq 'DUP';
	   $trans++ if $sv_type eq 'BND' or $sv_type eq 'TRA';
	   $invs++ if $sv_type eq 'INV';
	  
	}
	
	print "\n";
	if ($filter_flag){
		say "Running in filter mode: $filtered calls filtered out:";
		say " - $_: $filtered_sv{$_}" for sort {$filtered_sv{$b} <=> $filtered_sv{$a} } keys %filtered_sv;
		print "\n";		
	}

	say "$dels deletions";
	say "$dups duplications";
	say "$trans translocations";
	say "$invs inversions";
	
	my $top_count = 0;
	print "\nTop 5 SVs by read count:\n";
	for ( sort { $support_by_chrom{$b}[0] <=> $support_by_chrom{$a}[0] } keys %support_by_chrom ){
		$top_count++;
		print join("\n", "ID: $_", "TYPE: $support_by_chrom{$_}[1]", "CHROM: $support_by_chrom{$_}[2]", "READS: $support_by_chrom{$_}[0]\n");
		print "\n";
		last if $top_count >= 5;
	}	
}

sub get_variant {
	
	my ($id_lookup, $SVs, $info, $filter_flag) = @_;
	
	if (not $info->{$id_lookup}){
		say "Couldn't find any variant with ID: '$id_lookup' in file. Abort";
		exit;
	}
	
	my (@format) 		= @{ $info->{$id_lookup}->[0]};
	my (%format_long) 	= @{ $info->{$id_lookup}->[1]}; #
	my (%info_long)		= @{ $info->{$id_lookup}->[2]};
	
	my (@tumour_parts) 	= @{ $info->{$id_lookup}->[3]};
	my (@normal_parts) 	= @{ $info->{$id_lookup}->[4]};
	
	my (%information)	= @{ $info->{$id_lookup}->[5]};
	my (%sample_info)	= @{ $info->{$id_lookup}->[6]};
					
	my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $filters, $samples ) = @{ $SVs->{$id_lookup} };
    
	my @filter_reasons = @{ $filters };
		
	my @samples = @{ $samples };
	
	if (scalar @filter_reasons > 0 and $filter_flag){
	say "\n______________________________________________";	
	say "Variant '$id' was filtered for the following reasons:";
	say "* $_" foreach @filter_reasons;
	say "______________________________________________\n";
	}
	
	say "ID:     $id";
	say "TYPE:   $sv_type";
	$chr2 ? say "CHROM1: $chr" : say "CHROM:  $chr";
	say "CHROM2: $chr2" if $chr2;	
	say "START:  $start";
	say "STOP:   $stop";
	$chr2 eq $chr ? say "IGV:	$chr:$start-$stop" : say "IGV:	$chr:$start";
	say "LENGTH: $SV_length";
	say "PE:	$PE";
	say "SR:	$SR";
	say "QUAL:   $quality_score";
	say "FILT:   $filt";
	say "REF:    $ref";
	say "ALT:    $alt";
	
	
	say "__________________________________________________________________________________________________________________";
	printf "%-20s", "INFO";
	printf "%-20s", $_ for @samples;
	printf "%-s\n", "EXPLAINER";
	say "__________________________________________________________________________________________________________________";
	
	foreach my $format_block (@format){
		printf "%-20s", "$format_block";
		foreach (@samples){
			printf "%-20s", "$sample_info{$id_lookup}{$_}{$format_block}";
		}
		printf "%-s", "$format_long{$format_block}";
		print "\n";
	}

	say "____________________________________________________________________________________";
	printf "%-20s %-20s %-s\n", "INFO", "VALUE", "EXPLAINER";
	say "____________________________________________________________________________________";

	for (sort keys %{$information{$id_lookup}}){
		printf "%-20s %-20s %-s\n", $_, $information{$id_lookup}{$_}, $info_long{$_};
	}
	say "____________________________________________________________________________________";
		
}

sub dump_variants {
	
	my ( $SVs, $info, $filter_flag, $chromosome ) = @_;
			
	say "Enter any key to start cycling through calls or enter 'q' to exit";
	say "Running in filter mode - not displaying filtered calls" if $filter_flag;
	
	for (sort { @{ $SVs->{$a}}[0] cmp @{ $SVs->{$b}}[0] or
				@{ $SVs->{$a}}[1] <=> @{ $SVs->{$b}}[1] 
			}  keys %{ $SVs } ){
		
		my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $filters, $samples ) = @{ $SVs->{$_} };
		
		if ($chromosome){
			next if $chr ne $chromosome;
		}
		
		my (@format) 		= @{ $info->{$_}->[0]};
		my (%format_long) 	= @{ $info->{$_}->[1]}; #
		my (%info_long)		= @{ $info->{$_}->[2]};
	
		my (@tumour_parts) 	= @{ $info->{$_}->[3]};
		my (@normal_parts) 	= @{ $info->{$_}->[4]};
	
		my (%information)	= @{ $info->{$_}->[5]};
		my (%sample_info)	= @{ $info->{$_}->[6]};
	
		
		my @filter_reasons = @{ $filters };
				
		my @samples = @{ $samples };
		
		if (scalar @filter_reasons > 0 ){
			next if $filter_flag;
		}
		
		my $next_line = <>;
		say "Displaying info for variant '$id'. Enter any key to go to the next variant or type 'q' to exit\n";
		
		if ( $next_line ){
			chomp($next_line);
			exit if $next_line eq 'q';
			
			if (scalar @filter_reasons > 0 ){
				say "______________________________________________";	
				say "Variant '$id' was filtered for the following reasons:";
				say "* $_" foreach @filter_reasons;
				say "______________________________________________\n";
			}
	
			printf "%-10s %-s\n",  			"ID:", 		$id;
			printf "%-10s %-s\n", 			"TYPE:", 	$sv_type;
			$chr2 ? printf "%-10s %-s\n",  	"CHROM1:", 	$chr : printf "%-10s %-s\n",  "CHROM:",  $chr;
			printf "%-10s %-s\n", 			"CHROM2:", 	$chr2 if $chr2;	
			printf "%-10s %-s\n",   		"START:", 	$start;
			printf "%-10s %-s\n",   		"STOP:",  	$stop;
			$chr2 ? printf "%-10s %-s\n",   "IGV:",   	"$chr:$start" : printf "%-10s %-s\n", "IGV:", "$chr:$start-$stop";
			printf "%-10s %-s\n", 			"LENGTH:", 	$SV_length;
			printf "%-10s %-s\n",   		"PE:",    	$PE;
			printf "%-10s %-s\n",   		"SR:",		$SR;
			printf "%-10s %-s\n",   		"QUAL:",   	$quality_score;
			printf "%-10s %-s\n",   		"FILT:",   	$filt;
			printf "%-10s %-s\n",   		"REF:",   	$ref;
			printf "%-10s %-s\n",   		"ALT:",   	$alt;
			
			say "__________________________________________________________________________________________________________________";
			printf "%-20s", 				"INFO";
			printf "%-20s", 				$_ for @samples;
			printf "%-s\n", 				"EXPLAINER";
			say "__________________________________________________________________________________________________________________";
	
			foreach my $format_block (@format){
				printf "%-20s", 			$format_block;
				foreach (@samples){
					printf "%-20s", 		$sample_info{$id}{$_}{$format_block};
				}
				printf "%-s", 				$format_long{$format_block};
				print "\n";
			}

			say "____________________________________________________________________________________";
			printf "%-20s %-20s %-s\n", "INFO", "VALUE", "EXPLAINER";
			say "____________________________________________________________________________________";

			for (sort keys %{$information{$id}}){
				printf "%-20s %-20s %-s\n", $_, $information{$id}{$_}, $info_long{$_};
			}
			say "____________________________________________________________________________________";
	
		}
	}
	
}

sub chrom_filter {
	my ( $chr, $chr2, $filters ) = @_;
	
	my @keys = qw / 2L 2R 3L 3R 4 X Y /;
	my %chrom_filt;

	$chrom_filt{$_} = 1 for (@keys);
	
	my @filter_reasons = @{ $filters };
		 
	if ($chr2 eq '0'){
		$chr2 = $chr;
	}
	
	if ( not $chrom_filt{$chr} ){
		push @filter_reasons, 'chrom1=' . $chr;
	}
		
	elsif ( not $chrom_filt{$chr2} ){
		push @filter_reasons, 'chrom2=' . $chr2;
	}
	
	return (\@filter_reasons);

}

1;