#!/usr/bin/perl

package svParser;
use strict;
use warnings;
use autodie;

use feature qw/ say /;
use Data::Dumper;
use Data::Printer;

sub typer {
  my ($file, $type, $exclude_regions, $filters ) = @_;

  if ( $type eq 'l' ){
    say "Specified $file as Lumpy output";
    $type = 'lumpy';
    parse($file, $type, $exclude_regions, $filters);
  }

  elsif ( $type eq 'd' ){
    say "Specified $file as Delly output";
    $type = 'delly';
    parse($file, $type, $exclude_regions, $filters);
  }

  elsif ( $type eq 'n' ){
    say "Specified $file as novoBreak output";
    $type = 'novobreak';
    parse($file, $type, $exclude_regions, $filters);
  }
  elsif ($type eq 'snp'){
    say "Forcing parsing of $file";
    $type = 'snp';
    parse($file, $type, $exclude_regions, $filters);
  }

  elsif ( $type eq 'guess' ){

    if ( `grep "source=LUMPY" $file` ){
      say "Recognised $file as Lumpy input";
      $type = 'lumpy';
      parse($file, $type, $exclude_regions, $filters);
    }
    elsif ( `grep "DELLY" $file` ){
      say "Recognised $file as Delly input";
      $type = 'delly';
      parse($file, $type, $exclude_regions, $filters);
    }
    elsif ( `grep "bamsurgeon spike-in" $file` ){
      say "Recognised $file as novoBreak input";
      $type = 'novobreak';
      parse($file, $type, $exclude_regions, $filters);
    }
    else {
      die "This VCF can not be parsed. Try specfiying type '-t' explicitly. See -h for details. Abort";
    }
  }

  else {
    die "This VCF can not be parsed. Try specfiying type '-t' explicitly. See -h for details. Abort";
  }
}

sub parse {
  my ($file, $type, $exclude_regions, $filter_flags ) = @_;
  open my $in, '<', $file or die $!;

  my @headers;

  my %filter_flags = %{ $filter_flags };

  my (%SVs, %info, %filtered_SVs);
  my ($tumour_name, $control_name);
  my %format_long;
  my %info_long;
  my $filter_count;
  my @samples;
  my $replacement_id = 1;

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

    # NovoBreak doesn't assign ids
    if ($id eq 'N' or $id eq '.'){
      $id = $replacement_id++;
    }

    my %sample_parts;

    push @{$sample_parts{$samples[$_]}}, split(/:/, $sample_info[$_]) for 0..$#samples;

    my @tumour_parts   = split(/:/, $sample_info[0]);
    my @normal_parts   = split(/:/, $sample_info[1]) if @samples > 1; # In case there are no control samples...

    my @format        = split(/:/, $format_block);
    my @info_parts    = split(/;/, $info_block);

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

    if ( $filter_flags{'g'} ){

      #############################
      # Filter OUT somatic events #
      #############################

      # Filter if ANY of the panel of normals (but not direct control) are NOT 'GT = 0/0' (hom ref)
      for my $normal (@normals[1..$#normals]){
        if ( $sample_info{$id}{$normal}{'GT'} eq '1/1' or $sample_info{$id}{$normal}{'GT'} eq '0/1' ){
          push @filter_reasons, "$normal\_not_homo_ref=" . $sample_info{$id}{$normal}{'GT'};
        }
      }
      # If tumour is het/hom alt (0/1, 1/1), but direct control is hom ref(0/0), filter as somatic
      if ( ( $sample_info{$id}{$tumour_name}{'GT'} eq '0/1' or $sample_info{$id}{$tumour_name}{'GT'} eq '1/1' ) and $sample_info{$id}{$control_name}{'GT'} eq '0/0' ){
          push @filter_reasons, "$tumour_name\_somatic_event=" . $sample_info{$id}{$tumour_name}{'GT'};
        }
      elsif ( $sample_info{$id}{$tumour_name}{'GT'} eq '0/0' and $sample_info{$id}{$control_name}{'GT'} ne '0/0' ){
          push @filter_reasons, "$control_name\_exclusive_normal_event=" . $sample_info{$id}{$control_name}{'GT'};
        }

    }

    else {
      # Filter if ANY of the controls are NOT 'GT = 0/0' (hom ref)

      for my $normal (@normals){
        if ( $sample_info{$id}{$normal}{'GT'} eq '1/1' or $sample_info{$id}{$normal}{'GT'} eq '0/1' ){
          push @filter_reasons, "$normal\_not_homo_ref=" . $sample_info{$id}{$normal}{'GT'};
        }
      }
    }


    ###################

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

    my ($SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $filter_list);

    if ($type eq 'lumpy'){
      ( $SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $filter_list ) = lumpy( $id, $chr, $info_block, $SV_type, $alt, $start, \%sample_info, $tumour_name, $control_name, \@samples, \@normals, \@filter_reasons, \%filter_flags );
    }

    elsif ($type eq 'delly'){
      ( $SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, $filter_list ) = delly( $id, $info_block, $start, $SV_type, \@filter_reasons, \%filter_flags, $tumour_name, \%sample_info );
    }

    elsif ($type eq 'novobreak'){
      @samples = qw/tumour normal/;
      my ( $sample_info_novo, $format_novo, $format_long_novo );
      ( $SV_length, $chr2, $stop, $t_PE, $t_SR, $filter_list, $sample_info_novo, $format_novo, $format_long_novo ) = novobreak( $id, $info_block, $start, $SV_type, \@filter_reasons, \@sample_info, \%filter_flags );
      %sample_info = %{ $sample_info_novo };
      $ab = "-";
      @format = @{ $format_novo };
      %format_long = %{ $format_long_novo };
    }

    elsif ($type eq 'snp'){
      $chr2 = $chr;
      $filter_list = \@filter_reasons;
    }

    if ( $filter_flags{'chr'} ){
      $filter_list = chrom_filter( $chr, $chr2, $filter_list );
    }

    if ( $filter_flags{'e'} ){
      ###################
      # Region exclude ##
      ###################

      $filter_list = region_exclude_filter($chr, $start, $chr2, $stop, $exclude_regions, \@filter_reasons);
    }

    $SV_length = abs($SV_length);
    $SVs{$id} = [ @fields[0..10], $SV_type, $SV_length, $stop, $chr2, $t_SR, $t_PE, $ab, $filter_list, \@samples ];

    $info{$id} = [ [@format], [%format_long], [%info_long], [@tumour_parts], [@normal_parts], [%information], [%sample_info] ];

    if (scalar @{$filter_list} == 0){
      # say "$id passes all filters";
      $filtered_SVs{$.} = $_;
    }

  }
  return (\%SVs, \%info, \%filtered_SVs);
}

sub novobreak {
  my ( $id, $info_block, $start, $SV_type, $filters, $info, $filter_flags) = @_;

  my @filter_reasons = @{ $filters };

  my %filter_flags = %{ $filter_flags };

  my @info = @{ $info };

  my %sample_info;
  my %format_long;

  my $bp1_SR = $info[9];  # high qual split reads bp1
  my $bp1_PE = $info[26]; # PE reads bp1
  my $bp2_SR = $info[19]; # high qual split reads bp2
  my $bp2_PE = $info[28]; # PE reads bp2

  my $t_SR = $bp1_SR + $bp2_SR;
  my $t_PE = $bp1_PE + $bp2_PE;

  $t_SR = $t_SR/2;
  $t_PE = $t_PE/2;

  $t_SR = int($t_SR + 0.5);
  $t_PE = int($t_PE + 0.5);

  # my $t_PE = 0; # Don't believe PE read support!

  my $tumour_read_support = ( $t_PE + $t_SR );

  ########################
  # Read support filters #
  ########################


  if (exists $filter_flags{'su'}){
    my $filtered_on_reads = read_support_filter($tumour_read_support, $filter_flags{'su'}, \@filter_reasons);
    @filter_reasons = @{$filtered_on_reads};
  }

  my $n_bp1_SR = $info[14]; # high qual split reads bp1
  my $n_bp1_PE = $info[27]; # PE reads bp1
  my $n_bp2_SR = $info[24]; # high qual split reads bp2
  my $n_bp2_PE = $info[29]; # PE reads bp2

  my $all_control_read_support = ( $n_bp1_SR + $n_bp1_PE + $n_bp2_SR + $n_bp2_PE )/2;

  # Filter if there are more than 1 control reads
  if ( $all_control_read_support > 1 ){
    push @filter_reasons, "precise var with read support in a normal sample=" . $all_control_read_support;
  }

  ######################
  # Read depth filters #
  ######################

  my ($c_DP_1, $c_DP_2) = @info[6,16];
  my ($t_DP_1, $t_DP_2) = @info[11,21];

  my $c_DP = ($c_DP_1 + $c_DP_2)/2;
  my $t_DP = ($t_DP_1 + $t_DP_2)/2;

  if ( exists $filter_flags{'dp'} and $c_DP <= $filter_flags{'dp'} ){
    push @filter_reasons, 'control_depth<' . $filter_flags{'dp'} . '=' . $c_DP;
  }

  if ( exists $filter_flags{'dp'} and $t_DP <= $filter_flags{'dp'} ){
    push @filter_reasons, 'tumour_depth<' . $filter_flags{'dp'} . '=' . $t_DP;
  }

  my @tumour_parts = @info[6..10, 16..20, 26, 28];
  my @normal_parts = @info[11..15, 21..25, 27, 29];

  my @short_format = (
  "Breakpoint 1 depth",
  "Breakpoint 1 split reads",
  "Breakpoint 1 quality score",
  "Breakpoint 1 high quality split reads",
  "Breakpoint 1 high quality quality score",
  "Breakpoint 2 depth",
  "Breakpoint 2 split reads",
  "Breakpoint 2 quality score",
  "Breakpoint 2 high quality split reads",
  "Breakpoint 2 high quality quality score",
  "Breakpoint 1 discordant reads",
  "Breakpoint 2 discordant reads"
  );

  my @format = qw / DP1 SR1 Q1 HCSR1 HCQ1 DP2 SR2 Q2 HCSR2 HCQ2 PE1 PE2 /;

  $format_long{$format[$_]} = $short_format[$_] for 0..$#format;

  $sample_info{$id}{'tumour'}{$format[$_]} = $tumour_parts[$_] for 0..$#format;
  $sample_info{$id}{'normal'}{$format[$_]} = $normal_parts[$_] for 0..$#format;

  my ($stop) = $info_block =~ /;END=(.*?);/;

  my ($SV_length) = ($stop - $start);

  my ($chr2) = $info_block =~ /CHR2=(.*?);/;

  # if ($start > $stop){
  #   my $old_start = $start;
  #   my $old_stop = $stop;
  #   $start = $old_stop;
  #   $stop = $old_start;
  # }

  return ($SV_length, $chr2, $stop, $t_PE, $t_SR, \@filter_reasons, \%sample_info, \@format, \%format_long );
}


sub lumpy {
  my ( $id, $chr, $info_block, $SV_type, $alt, $start, $sample_info, $tumour, $control, $samples, $normals, $filters, $filter_flags ) = @_;

  my @filter_reasons = @{ $filters };
  my @normals = @{ $normals };

  my %filter_flags = %{ $filter_flags };

  my @samples = @{ $samples };

  my %sample_info = %{ $sample_info };

  my ($SV_length) = $info_block =~ /SVLEN=(.*?);/;

  # switch tumour normal
  if ($filter_flags{'g'}){
    my $tum2norm = $control;
    $control = $tumour;
    $tumour = $tum2norm;
  }

  my ($t_PE, $t_SR, $all_c_PE, $all_c_SR, $c_PE, $c_SR) = (0,0,0,0,0,0);

  # Get allele balance
  my $ab = $sample_info{$id}{$tumour}{'AB'};

  ########################
  # Read support filters #
  ########################

  # Get read support for tumour
  $t_PE = $sample_info{$id}{$tumour}{'PE'};
  $t_SR = $sample_info{$id}{$tumour}{'SR'};

  my $tumour_read_support = ( $t_PE + $t_SR );

  # Create temp pseudo counts to avoid illegal division by 0
  my $pc_tumour_read_support = $tumour_read_support + 0.001;

  my $pc_direct_control_read_support = 0;

  if (exists $filter_flags{'su'}){
    my $filtered_on_reads = read_support_filter($tumour_read_support, $filter_flags{'su'}, \@filter_reasons);
    @filter_reasons = @{ $filtered_on_reads };
  }

  if (not $filter_flags{'g'}){

    if (@samples > 1){ # In case there are no control samples...

      # for precise variants:
      if ($info_block !~ /IMPRECISE;/){

        # We want to be strict, so include all controls used for genotyping (and sum read support)
        for my $normal (@normals){
          $sample_info{$id}{$normal}{'PE'} eq '.' ? $sample_info{$id}{$normal}{'PE'} = '0' : $all_c_PE += $sample_info{$id}{$normal}{'PE'};
          $sample_info{$id}{$normal}{'SR'} eq '.' ? $sample_info{$id}{$normal}{'SR'} = '0' : $all_c_SR += $sample_info{$id}{$normal}{'SR'};
        }

        my $all_control_read_support = ( $all_c_PE + $all_c_SR );

        # Filter if there are more than 1 control reads
        if ( $all_control_read_support > 1 ){
          push @filter_reasons, "precise var with read support in a normal sample=" . $all_control_read_support;
        }
      }

      # for imprecise variants:
      if ($info_block =~ /IMPRECISE;/){

        # Get read support for direct control only
        $c_PE =  $sample_info{$id}{$control}{'PE'};
        $c_SR =  $sample_info{$id}{$control}{'SR'};

        my $direct_control_read_support = ( $c_PE + $c_SR );
        $pc_direct_control_read_support = $direct_control_read_support + 0.001;

        # Filter if # tumour reads supporting var is less than 5 * control reads
        # Or if there are more than 2 control reads
        if ( $pc_tumour_read_support/$pc_direct_control_read_support < 5 ){
          push @filter_reasons, 'imprecise var with less than 5 * more tum reads than control=' . $direct_control_read_support;
        }
        if ( $direct_control_read_support > 1 ){
          push @filter_reasons, 'imprecise var with control read support=' . $direct_control_read_support;
        }

      }
    }

  }

  # # if running in germline mode, require
  # if ($filter_flags{'g'}){
  #
  # }

  ######################
  # Read depth filters #
  ######################

  if ( exists $sample_info{$id}{$tumour}{'DP'} ){

    my $t_DP =  $sample_info{$id}{$tumour}{'DP'};

    if (@samples > 1){ # In case there are no control samples...

      my $c_DP =  $sample_info{$id}{$control}{'DP'};

      # Flag if either control or tumour has depth < 10 at site

      if ( exists $filter_flags{'dp'} and $c_DP <= $filter_flags{'dp'} ){
        push @filter_reasons, 'control_depth<' . $filter_flags{'dp'} . '=' . $c_DP;
      }
    }

    if ( exists $filter_flags{'dp'} and $t_DP <= $filter_flags{'dp'} ){
      push @filter_reasons, 'tumour_depth<' . $filter_flags{'dp'} . '=' . $t_DP;
    }


    # Subtract control reads from tumour reads
    # If this number of SU is less than 10% of tumour read_depth then filter
    if (not $filter_flags{'g'}){
      if ( exists $filter_flags{'rdr'} and ( $tumour_read_support - $pc_direct_control_read_support ) / ( $t_DP + 0.01 ) < $filter_flags{'rdr'} ){ # Maybe this is too harsh...
        push @filter_reasons, 'tumour_reads/tumour_depth<' . ($filter_flags{'rdr'}*100) . "%" . '=' . $tumour_read_support . "/" . $t_DP;
      }
    }

  }

  ##################
  # Quality filter #
  ##################

  if ( exists $sample_info{$id}{$tumour}{'SQ'} and exists $filter_flags{'sq'} ){
    $sample_info{$id}{$tumour}{'SQ'} = 0 if $sample_info{$id}{$tumour}{'SQ'} eq '.';

    if ( $sample_info{$id}{$tumour}{'SQ'} <= $filter_flags{'sq'} ){
      push @filter_reasons, "SQ<" . $filter_flags{'sq'} . '=' . $sample_info{$id}{$tumour}{'SQ'};
    }
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

  return ($SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, \@filter_reasons);
}


sub delly {
  my ($id, $info_block, $start, $SV_type, $filters, $filter_flags, $tumour_name, $sample_ref) = @_;

  my @filter_reasons = @{ $filters };

  my %filter_flags = %{ $filter_flags };

  my %sample_info = % { $sample_ref };

  my $dv = $sample_info{$id}{$tumour_name}{'DV'};
  my $dr = $sample_info{$id}{$tumour_name}{'DR'};

  my $ab = $dv/($dv+$dr);

  my ($stop) = $info_block =~ /;END=(.*?);/;

  my ($SV_length) = ($stop - $start);

  my ($t_SR, $t_PE) = (0,0);

    if ($info_block =~ /;SR=(\d+);/){
      $t_SR = $1;
    }

    if ($info_block =~ /;PE=(\d+);/){
      $t_PE = $1;
    }

  my $tumour_read_support = ( $t_PE + $t_SR );

  if (exists $filter_flags{'su'}){
    my $filtered_on_reads = read_support_filter($tumour_read_support, $filter_flags{'su'}, \@filter_reasons);
    @filter_reasons = @{$filtered_on_reads};
  }

  # if ($start > $stop){
  #   my $old_start = $start;
  #   my $old_stop = $stop;
  #   $start = $old_stop;
  #   $stop = $old_start;
  # }

  # my ($chr2) = 0;

  # if ($SV_type eq 'TRA'){
    my ($chr2) = $info_block =~ /CHR2=(.*?);/;
  # }

  return ($SV_length, $chr2, $stop, $t_SR, $t_PE, $ab, \@filter_reasons );
}


sub summarise_variants {
  my ( $SVs, $filter_switch, $chromosome ) = @_;

  my ($dels, $dups, $trans, $invs, $filtered) = (0,0,0,0,0);
  my ($tds, $CNVs, $ins) = (0,0,0);

  my ( $query_region, $query_start, $query_stop );

  my $specified_region = 0;

    if ( $chromosome ){

      if ( $chromosome =~ /:/ ){

        ($chromosome, $query_region) = split(/:/, $chromosome);

        if ( $query_region !~ /-/ ){
          die "Error parsing the specified region.\nPlease specify chromosome regions using the folloing format:\tchrom:start-stop\n";
        }

        ($query_start, $query_stop) = split(/-/, $query_region);
        $query_start =~ s/,//g;
        $query_stop =~ s/,//g;

        say "Limiting search to SVs within region '$chromosome:$query_start-$query_stop'";

        $specified_region = 1;
      }
      else {
        say "Limiting search to SVs on chromosome '$chromosome'";
      }

    }

  my %support_by_chrom;

  my $read_support;

  my %filtered_sv;

  for (keys %{ $SVs } ){

    my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $samples ) = @{ $SVs->{$_} };

    if ( $chromosome ){
      next if $chr ne $chromosome;
    }

    if ( $specified_region ){

      if ( ( $start < $query_start or $start > $query_stop ) and ( $stop < $query_stop or $stop > $query_stop ) ){
         next;
       }
    }

     # Need this re-assignment for novobreak - should be harmless for lumpy and delly
    $id = $_;

    my @filter_reasons = @{ $filters };

    foreach (@filter_reasons){
      my ($reason) = $_ =~ /(.+)=/;
      $filtered_sv{$reason}++ ;
    }

    if ( scalar @filter_reasons > 0 ){
      $filtered++;
      next if $filter_switch;
    }

    $read_support = ( $SR + $PE );
    $support_by_chrom{$id} = [ $read_support, $sv_type, $chr, $SV_length, $start ];

    $dels++ if $sv_type eq 'DEL';
    $dups++ if $sv_type eq 'DUP';
    $trans++ if $sv_type eq 'BND' or $sv_type eq 'TRA';
    $invs++ if $sv_type eq 'INV';
    $tds++ if $sv_type eq 'DUP:TANDEM';
    $CNVs++ if $sv_type eq 'CNV';
    $ins++ if $sv_type eq 'INS';

  }
  print "\n";

  if ($filter_switch){
    say "Running in filter mode: $filtered calls filtered out:";
    say " - $_: $filtered_sv{$_}" for sort {$filtered_sv{$b} <=> $filtered_sv{$a} } keys %filtered_sv;
    print "\n";
  }

  say "$dels deletions";
  say "$dups duplications";
  say "$trans translocations";
  say "$invs inversions";
  say "$tds tandem duplications" if $tds > 0;
  say "$CNVs CNV regions" if $CNVs > 0;
  say "$ins inserted sequences at BP" if $ins > 0;

  my $top_count = 0;
  my %connected_bps;

  print "\nTop SVs by read count:\n";
  for ( sort { $support_by_chrom{$b}[0] <=> $support_by_chrom{$a}[0] } keys %support_by_chrom ){
    my $bp_id = $_;

    if ($bp_id =~ /_/){
      ($bp_id) = $bp_id =~ /(.+)?_/;
    }

    # flatten connected bps into 1 id for summary
    next if $connected_bps{$bp_id}++;
    $top_count++;

    print join("\n",
    "ID: $_",
    "TYPE: $support_by_chrom{$_}[1]",
    "CHROM: $support_by_chrom{$_}[2]",
    "START: $support_by_chrom{$_}[4]",
    "READS: $support_by_chrom{$_}[0]",
    "LENGTH: $support_by_chrom{$_}[3]\n") . "\n";

    last if $top_count >= 5;
  }
}


sub get_variant {

  my ($id_lookup, $SVs, $info, $filter_flag) = @_;

  if (not $info->{$id_lookup}){
    say "Couldn't find any variant with ID: '$id_lookup' in file. Abort";
    exit;
  }

  my (@format)     = @{ $info->{$id_lookup}->[0]};
  my (%format_long)   = @{ $info->{$id_lookup}->[1]};
  my (%info_long)    = @{ $info->{$id_lookup}->[2]};

  my (@tumour_parts)   = @{ $info->{$id_lookup}->[3]};
  my (@normal_parts)   = @{ $info->{$id_lookup}->[4]};

  my (%information)  = @{ $info->{$id_lookup}->[5]};
  my (%sample_info)  = @{ $info->{$id_lookup}->[6]};

  my ($chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $samples ) = @{ $SVs->{$id_lookup} };

  my @filter_reasons = @{ $filters };

  my @samples = @{ $samples };

  # Should change so that it will only print filter reasons if user specifies them

  if (scalar @filter_reasons > 0 ){
  say "\n______________________________________________";
  say "Variant '$id_lookup' will be filtered for the following reasons:";
  say "* $_" foreach @filter_reasons;
  say "______________________________________________\n";
  }

  printf "%-10s %-s\n",        "ID:",     $id_lookup;
  printf "%-10s %-s\n",       "TYPE:",   $sv_type;
  $chr2 ? printf "%-10s %-s\n",    "CHROM1:",   $chr : printf "%-10s %-s\n",  "CHROM:",  $chr;
  printf "%-10s %-s\n",       "CHROM2:",   $chr2 if $chr2;
  printf "%-10s %-s\n",       "START:",   $start;
  printf "%-10s %-s\n",       "STOP:",    $stop;
  ($chr2 and ($chr2 ne $chr) ) ? printf "%-10s %-s\n",   "IGV:",     "$chr:$start" : printf "%-10s %-s\n", "IGV:", "$chr:$start-$stop";
  printf "%-10s %-s\n",       "LENGTH:",   $SV_length;
  printf "%-10s %-s\n",       "PE:",      $PE;
  printf "%-10s %-s\n",       "SR:",    $SR;
  printf "%-10s %-s\n",       "QUAL:",     $quality_score;
  printf "%-10s %-s\n",       "FILT:",     $filt;
  printf "%-10s %-s\n",       "REF:",     $ref;
  printf "%-10s %-s\n",       "ALT:",     $alt;

  say "__________________________________________________________________________________________________________________";
  printf "%-20s",         "INFO";
  printf "%-20s",         $_ for @samples;
  printf "%-s\n",         "EXPLAINER";
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
    # turn off warnings for badly formatted novobreak vcf
    no warnings;
    printf "%-20s %-20s %-s\n", $_, $information{$id_lookup}{$_}, $info_long{$_};
  }
  say "____________________________________________________________________________________";

}


sub dump_variants {
  my ( $SVs, $info, $filter_flag, $chromosome, $type ) = @_;
  my ( $query_region, $query_start, $query_stop );
  my $specified_region = 0;

  if ( $chromosome ){
    if ( $chromosome =~ /:/ ){
      ($chromosome, $query_region) = split(/:/, $chromosome);
      if ( $query_region !~ /-/ ){
        die "Error parsing the specified region.\nPlease specify chromosome regions using the folloing format:\tchrom:start-stop\n";
      }

      ($query_start, $query_stop) = split(/-/, $query_region);
      $query_start =~ s/,//g;
      $query_stop =~ s/,//g;

      say "Limiting search to SVs within region '$chromosome:$query_start-$query_stop'";

      $specified_region = 1;
    }
    else {
      say "Limiting search to SVs on chromosome '$chromosome'";
    }

  }

  say "Running in filter mode - not displaying filtered calls" if $filter_flag;
  say "\nEnter any key to start cycling through calls or enter 'q' to exit";

  for ( sort { @{ $SVs->{$a}}[0] cmp @{ $SVs->{$b}}[0] or
        @{ $SVs->{$a}}[1] <=> @{ $SVs->{$b}}[1]
      }  keys %{ $SVs } ){

    my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $samples ) = @{ $SVs->{$_} };

    $id = $_;

    if ( $chromosome ){
      next if $chr ne $chromosome;
    }

    if ( $specified_region ){
      if ( ( $start < $query_start or $start > $query_stop ) and ( $stop < $query_stop or $stop > $query_stop ) ){
        next;
      }
    }

    my (@format)        = @{ $info->{$_}->[0] };
    my (%format_long)   = @{ $info->{$_}->[1] };
    my (%info_long)     = @{ $info->{$_}->[2] };

    my (@tumour_parts)  = @{ $info->{$_}->[3] };
    my (@normal_parts)  = @{ $info->{$_}->[4] };

    my (%information)   = @{ $info->{$_}->[5] };
    my (%sample_info)   = @{ $info->{$_}->[6] };

    my @filter_reasons  = @{ $filters };

    my @samples         = @{ $samples };

    if ( scalar @filter_reasons > 0 ){
      next if $filter_flag;
    }

    my $next_line = <>;
    say "Displaying info for variant '$id'. Enter any key to go to the next variant or type 'q' to exit\n";

    if ( $next_line ){
      chomp($next_line);
      exit if $next_line eq 'q';

      if (scalar @filter_reasons > 0 ){
        say "______________________________________________";
        say "Variant '$id' will be filtered for the following reasons:";
        say "* $_" foreach @filter_reasons;
        say "______________________________________________\n";
      }
      if ($type ne 'snp'){
        printf "%-10s %-s\n",        "ID:",         $id;
        printf "%-10s %-s\n",        "TYPE:",       $sv_type;
        $chr2 ? printf "%-10s %-s\n",    "CHROM1:", $chr : printf "%-10s %-s\n",  "CHROM:",  $chr;
        printf "%-10s %-s\n",       "CHROM2:",      $chr2 if $chr2;
        printf "%-10s %-s\n",       "START:",       $start;
        printf "%-10s %-s\n",       "STOP:",        $stop;
        ($chr2 and ($chr2 ne $chr) ) ? printf "%-10s %-s\n",   "IGV:",     "$chr:$start" : printf "%-10s %-s\n", "IGV:", "$chr:$start-$stop";
        printf "%-10s %-s\n",       "LENGTH:",      $SV_length;
        printf "%-10s %-s\n",       "PE:",          $PE;
        printf "%-10s %-s\n",       "SR:",          $SR;
        printf "%-10s %-s\n",       "QUAL:",        $quality_score;
        printf "%-10s %-s\n",       "FILT:",        $filt;
        printf "%-10s %-s\n",       "REF:",         $ref;
        printf "%-10s %-s\n",       "ALT:",         $alt;
      }

      elsif ($type eq 'snp'){
        printf "%-10s %-s\n",    "CHROM:",  $chr;
        printf "%-10s %-s\n",    "POS:",    $start;
        printf "%-10s %-s\n",    "IGV:",    "$chr:$start";
        printf "%-10s %-s\n",    "FILT:",   $filt;
        printf "%-10s %-s\n",    "REF:",    $ref;
        printf "%-10s %-s\n",    "ALT:",    $alt;
        printf "%-10s %-s\n",    "MUT:",    "$ref>$alt";
      }

        say "__________________________________________________________________________________________________________________";
        printf "%-20s",         "INFO";
        printf "%-20s",         $_ for @samples;
        printf "%-s\n",         "EXPLAINER";
        say "__________________________________________________________________________________________________________________";

      foreach my $format_block (@format){
        printf "%-20s",       $format_block;
        foreach (@samples){
          printf "%-20s",     $sample_info{$id}{$_}{$format_block};
        }
        printf "%-s",         $format_long{$format_block};
        print "\n";
      }

      say "____________________________________________________________________________________";
      printf "%-20s %-20s %-s\n", "INFO", "VALUE", "EXPLAINER";
      say "____________________________________________________________________________________";

      for (sort keys %{$information{$id}}){
        # turn off warnings for badly formatted novobreak vcf
        no warnings;
        printf "%-20s %-20s %-s\n", $_, $information{$id}{$_}, $info_long{$_};
      }
      say "____________________________________________________________________________________";

    }
  }
}


sub print_variants {

  my ( $SVs, $filtered_SVs, $name, $output_dir, $germline ) = @_;
  my $out;

  if ($germline){
    open $out, '>', $output_dir . $name . ".germline_filtered.vcf" or die $!;
    say "Writing output to " . "'$output_dir" . $name . ".germline_filtered.vcf'";
  }
  else{
    open $out, '>', $output_dir . $name . ".filtered.vcf" or die $!;
    say "Writing output to " . "'$output_dir" . $name . ".filtered.vcf'";
  }

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


sub write_summary {
  my ( $SVs, $name, $summary_out, $type, $germline ) = @_;

  my $info_file;

  if ($germline){
    open $info_file, '>', $summary_out . $name . ".germline_filtered.summary.txt" or die $!;
    say "Writing useful info to " . "'$summary_out" . $name . ".germline_filtered.summary.txt'";
  }
  else{
    open $info_file, '>', $summary_out . $name . ".filtered.summary.txt" or die $!;
    say "Writing useful info to " . "'$summary_out" . $name . ".filtered.summary.txt'";
  }

  $type = "lumpy" if $type eq 'l';
  $type = "delly" if $type eq 'd';
  $type = "novobreak" if $type eq 'n';

  my %connected_bps;

  print $info_file join("\t", "source", "type", "chromosome1", "bp1", "chromosome2", "bp2", "split reads", "pe reads", "id", "length(Kb)", "position", "consensus|type", "microhomology", "configuration", "allele_frequency", "mechanism|log2(cnv)") . "\n";

  for ( sort { @{ $SVs->{$a}}[0] cmp @{ $SVs->{$b}}[0] or
        @{ $SVs->{$a}}[1] <=> @{ $SVs->{$b}}[1]
      }  keys %{ $SVs } ){
    my ( $chr, $start, $id, $ref, $alt, $quality_score, $filt, $info_block, $format_block, $tumour_info_block, $normal_info_block, $sv_type, $SV_length, $stop, $chr2, $SR, $PE, $ab, $filters, $samples ) = @{ $SVs->{$_} };

    if (scalar @{$filters} == 0){

      my $bp_id = $_;
      if ($bp_id =~ /_/){
        ($bp_id) = $bp_id =~ /(.+)?_/;
      }
      # flatten connected bps into 1 id for summary
      next if $connected_bps{$bp_id}++;

      my ($length_in_kb) = sprintf("%.1f", abs($SV_length)/1000);

      $ab = sprintf("%.2f", $ab) unless $type eq 'novobreak' or $ab eq '.';

      # Don't include DELS < 1kb with split read support == 0
      if ( ( $sv_type eq "DEL" and $length_in_kb < 1 ) and $SR == 0 ){
        say "Ommiting SV '$_' from '$name\.filtered.summary.txt' as $sv_type with length: $length_in_kb and split read support of $SR";
        next;
      }

      my ($consensus, $mh_length, $ct, $rdr, $rde );

      # Consensus seq
      if ($info_block =~ /CONSENSUS=(.*?);/){
         $consensus = $1;
      }
      else{
        $consensus = "-";
      }

      # Read depth ratio (delly)
      if ($info_block =~ /RDRATIO=(\d+\.?\d*)/){
        $rdr = log($1)/log(2);
        $rdr = sprintf("%.2f", $rdr)
      }
      else{
        $rdr = '-';
      }

      # # Read depth evidence (lumpy)
      # if ($info_block =~ /BD=(\d+)/){
      #   $rde = $1;
      # }
      # else {
      #   $rde = '-';
      # }

      # Microhology length (delly)
      if ($info_block =~ /HOMLEN=(\d+);/){
        $mh_length = $1;
      }
      else{
        $mh_length = "-";
      }

      # Configuration
      if ($info_block =~ /CT=(.*?);/){
        ($ct) = $1;
      }
      elsif ($alt =~ /\[|\]/) {
        $ct = $alt;
      }
      else {
        $ct = "-";
      }

      if ( $chr2 and ($chr2 ne $chr) ){
        print $info_file join("\t", $type, $sv_type, $chr, $start, $chr2, $stop, $SR, $PE, $_, $length_in_kb, "$chr:$start $chr2:$stop", $consensus, $mh_length, $ct, $ab, $rdr ) . "\n";
      }
      else {
        print $info_file join("\t", $type, $sv_type, $chr, $start, $chr, $stop, $SR, $PE, $_, $length_in_kb, "$chr:$start-$stop", $consensus, $mh_length, $ct, $ab, $rdr) . "\n";
      }

    }
  }
}


sub region_exclude_filter {
  my ( $chr1, $bp1, $chr2, $bp2, $exclude_regions, $filter_reasons ) = @_;

  my @filter_reasons = @{ $filter_reasons };

  my @bed = @{ $exclude_regions };

  foreach(@bed){
    chomp;
    my ($chromosome, $start, $stop) = split;

    if ( $bp1 >= $start and $bp1 <= $stop ) {
      next unless $chromosome eq $chr1;
      push @filter_reasons, 'bp1_in_unmappable_region=' . "$chromosome:$bp1\_in:" . $start . '-' . $stop;
    }
    if ( $bp2 >= $start and $bp2 <= $stop ) {
      next unless $chromosome eq $chr2;
      push @filter_reasons, 'bp2_in_unmappable_region=' . "$chromosome:$bp2\_in:" . $start . '-' . $stop;
    }
  }

  return (\@filter_reasons);
}


sub read_support_filter {
  my ($tumour_read_support, $read_support_flag, $filter_reasons ) = @_;
  my @filter_reasons = @{ $filter_reasons };

  # Filter if tum reads below specified threshold [default=4]
  if ( $tumour_read_support < $read_support_flag ){
    push @filter_reasons, 'tumour_reads<' . $read_support_flag . '=' . $tumour_read_support;
  }
  return(\@filter_reasons);
}


# Only for Drosophila so far...
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
