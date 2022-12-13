=head1 LICENSE

Copyright [1999-2014]

=cut

=head1 CONTACT

Lu.Liu@RoswellPark.org

=cut

=head1 NAME

Annotation - a script using three annotations programs to implement MNV annotaion.

Version 01

by Lu Liu (Lu.Liu@RoswellPark.org)

=cut

# This code would group the MNVs based on their codon position and give the corresponding
# annotations based on three different choices of annotation program: Annovar, SnpEff, Vep.

package Annotation;
use strict;
use POSIX;
use Bio::DB::Sam;
use IPC::System::Simple qw(system capture);
use Data::Dumper;
use Bio::DB::GFF;

our @input;
our %aoo            = ();
our %aono           = ();
our %aoi            = ();
our %ao_nonexon     = ();
our %aos            = ();
our $input_filename = "";
our $temp_filename  = "";
our %strand_table   = ();
our %cluster        = ();
our %clusters       = ();
our $option         = ();
our %mapping        = ();
our %mapping_vep    = ();

#-------------------------------------------------------------#

sub classify_variant {
    my $in = shift;
    if ( $in =~ /^[0-9A-Z_]+\.\d+\.([ATGCN]+)\.([ATGCN]+)$/i )
    {    #no dash: SNV/MNV/Indel
        if ( length($1) == length($2) ) {
            length($1) == 1 ? return 'SNV' : return 'MNV';
        }
        else {
            return 'Indel';
        }
    }
    elsif ( $in =~ /^[0-9A-Z_]+\.\d+\.[ATGCN-]+\.[ATGCN-]+$/i )
    {    #with dashL Indel
        return 'Indel';
    }
    else {
        return 'Junk';
    }
}

#-------------------------------------------------------------#

sub pre_annovar_process {

    my $fa_ref      = $_[1]{'assembly'};
    my $fa_infile   = Bio::DB::Sam::Fai->open($fa_ref);
    my $output_name = $temp_filename;
    open my $of, ">", $output_name or die("could not open file . $!");

    my ( $l_id, $new_l_id, $chr, $p_s, $p_e, $s_s, $s_e, $s_i, $d );
    my ( %original, @reformated );
    ( $l_id, $new_l_id ) = ( 1, 1 );
    my $hash_ref = $_[0];
    foreach my $line ( keys %$hash_ref ) {
        $original{$line} = $l_id;

        #parse the input file and build hash table
        $line =~ s/[\r\n]+$//g;
        my @columns = split( /,/, $line );
        my %table;

        for my $col ( 0 .. @columns - 1 ) {
            my $category = classify_variant( $columns[$col] );
            if ( $category eq 'Indel' ) {
                die("Indel found. $columns[$col]");
            }
            else {
                my @fields = split( /\./, $columns[$col] );
                push(
                    @{ ${ $table{ $fields[0] } }{ $fields[1] } },
                    @fields[ 2 .. 3 ]
                );
            }
        }

#process the hash table: extract sequence if necessary and then print accordingly
        foreach my $chr ( sort { $a cmp $b } keys %table ) {
            my ( $s_s, $s_e, $p_prev, $p_curr, $gap, $p_s, $p_e );
            my @sorted_pos = sort { $a <=> $b } keys %{ $table{$chr} };
            $p_s = $sorted_pos[0];
            $p_e = $sorted_pos[0];
            $s_s = @{ ${ $table{$chr} }{$p_s} }[0];
            $s_e = @{ ${ $table{$chr} }{$p_s} }[1];
            ( $p_prev, $p_curr ) = $p_s;
            for my $pos ( 1 .. @sorted_pos - 1 ) {
                $s_i    = "";
                $p_curr = $sorted_pos[$pos];
                $gap    = $p_curr - $p_prev;

                if ( $gap == 1 ) {
                    $s_s = $s_s . $s_i . @{ ${ $table{$chr} }{$p_curr} }[0];
                    $s_e = $s_e . $s_i . @{ ${ $table{$chr} }{$p_curr} }[1];
                    $p_e++;
                }
                elsif ( $gap == 2 ) {
                    $s_i =
                      $fa_infile->fetch(
                        $chr . ":" . ( $p_prev + 1 ) . "-" . ( $p_curr - 1 ) );
                    $s_s = $s_s . $s_i . @{ ${ $table{$chr} }{$p_curr} }[0];
                    $s_e = $s_e . $s_i . @{ ${ $table{$chr} }{$p_curr} }[1];
                    $p_e += 2;
                }
                else {
                    if ( $s_s eq "-" ) {
                        --$p_s;
                        --$p_e;
                        if ( $option ne "vep" ) {
                            print $of
                              "$chr\t$p_s\t$p_e\t$s_s\t$s_e\t+\t$l_id\n";
                        }
                        else {
                            print $of
                              "$chr\t$p_s\t$l_id\t$s_s\t$s_e\t+\t$l_id\n";
                        }
                        @{ $mapping_vep{$l_id}{$chr}{$p_s}{$p_e} } =
                          ( $s_s, $s_e )
                          ;    # keep a hash table for the generated input file
                        $mapping{"line$new_l_id"}{"mapto"} = $l_id;
                        $mapping{"line$new_l_id"}{"variant"} =
                          "$chr\t$p_s\t$p_e\t$s_s\t$s_e\t+\t$l_id"
                          ;    # keep a hash table for the generated input file
                    }
                    else {
                        if ( $option ne "vep" ) {
                            print $of
                              "$chr\t$p_s\t$p_e\t$s_s\t$s_e\t+\t$l_id\n";
                        }
                        else {
                            print $of
                              "$chr\t$p_s\t$l_id\t$s_s\t$s_e\t+\t$l_id\n";
                        }
                        @{ $mapping_vep{$l_id}{$chr}{$p_s}{$p_e} } =
                          ( $s_s, $s_e )
                          ;    # keep a hash table for the generated input file
                        $mapping{"line$new_l_id"}{"mapto"} = $l_id;
                        $mapping{"line$new_l_id"}{"variant"} =
                          "$chr\t$p_s\t$p_e\t$s_s\t$s_e\t+\t$l_id"
                          ;    # keep a hash table for the generated input file
                    }
                    $s_s = @{ ${ $table{$chr} }{$p_curr} }[0];
                    $s_e = @{ ${ $table{$chr} }{$p_curr} }[1];
                    $p_s = $p_curr;
                    $p_e = $p_curr;
                    $new_l_id++;
                }
                $p_prev = $p_curr;
            }
            if ( $option ne "vep" ) {
                print $of "$chr\t$p_s\t$p_e\t$s_s\t$s_e\t+\t$l_id\n";
            }
            else {
                print $of "$chr\t$p_s\t$l_id\t$s_s\t$s_e\t+\t$l_id\n";
            }
            @{ $mapping_vep{$l_id}{$chr}{$p_s}{$p_e} } =
              ( $s_s, $s_e );   # keep a hash table for the generated input file
            $mapping{"line$new_l_id"}{"mapto"} = $l_id;
            $mapping{"line$new_l_id"}{"variant"} =
              "$chr\t$p_s\t$p_e\t$s_s\t$s_e\t+\t$l_id"
              ;                 # keep a hash table for the generated input file
            $new_l_id++;
        }
        $l_id++;
    }
    close $of;
    return ( \%original );
}

#-------------------------------------------------------------#

#parse information from exonic_variant_function file and save to table
sub obtain_gene_exonic_variant {
    my %table;
    open( INFILE, $input_filename );
    while ( my $line = <INFILE> ) {
        $line =~ s/[\r\n]+$//g;
        my @columns  = split( /,/, $line );
        my $size_trs = $#columns - 1;
        my @temp     = split( /\t/, $columns[0] );
        my $lineid   = $temp[0];
        for my $i ( 0 .. $size_trs ) {
            $columns[$i] =~ m/([-_\w]*):(\w*):(\w*):c.(\w*)/;
            my ( $g_name, $t_name, $exon, $mut ) = ( $1, $2, $3, $4 );
            $exon =~ m/exon(\d*)/;
            $table{$lineid}{$t_name} =
              { 'gene' => $g_name, 'exon' => $1, 'mut' => $mut };
        }
        $table{"input"} = $columns[-1];
    }
    return \%table;
}

#-------------------------------------------------------------#

sub parse_predicted_annovar_output_Original {

	%aoo = ();
	%aono = ();
	%aoi = ();
    $Annotation::input_filename =
      $Annotation::temp_filename . ".exonic_variant_function";
    my $trs_gene_info = Annotation::obtain_gene_exonic_variant();

    $Annotation::input_filename = "predict_annovar_short.txt";
    open( INFILE, $input_filename );
    while ( my $line = <INFILE> ) {
        $line =~ s/[\r\n\t]+$//g;
        my @columns = split( / +/, $line );

        # get the original variations information
        my ( $inf_2, $mRNA_name, $ref_aa, $mut_aa ) =
          ( $columns[0], $columns[1], $columns[8], $columns[10] );
        $columns[5] =~ s/  / /g;
        $columns[5] =~ m/(\d*)-(\d*)/;
        my ( $aa_ids, $aa_ide ) = ( $1, $2 );
        ${$trs_gene_info}{$inf_2}{$mRNA_name}{"mut"} =~ m/(\d*)_(\d*)\w*/;
        my ( $mrna_s, $mrna_e ) = ( $1, $2 );
        my @input = split( /\t/, $mapping{$inf_2}{"variant"} );
        my ( $var_s, $var_e ) = ( $input[1], $input[2] );
        my $var_len = $var_e - $var_s;
        my $chr     = $input[0];
        my @in_s_s  = split( "", $input[3] );
        my @in_s_e  = split( "", $input[4] );
        my $inf_1   = $input[6];

        # get mutation from each transcript

        my $gene_name = ${$trs_gene_info}{$inf_2}{$mRNA_name}{'gene'};
        my $strand    = $strand_table{$gene_name}{$mRNA_name}{'strand'};
        my @ref_aa_l  = split( "", $ref_aa );
        my @mut_aa_l  = split( "", $mut_aa );
        my %aa_mu_T   = ();
        foreach ( $aa_ids .. $aa_ide ) {
            $aa_mu_T{$_} =
              [ $ref_aa_l[ $_ - $aa_ids ], $mut_aa_l[ $_ - $aa_ids ] ];
        }

        #create hash table
        my @aoa      = ();
        my $mrn_len  = $mrna_e - $mrna_s;
        my $diff_len = $var_len - $mrn_len;
        my $loop_len = ();

        if ( $diff_len >= 0 ) {
            $loop_len = $mrn_len;
        }
        else {
            $loop_len = $var_len;
        }

        # find the corresponding exon range
        # evaluate which side of the variant is located in intron
        # if the left side then
        # else if the right side
        my $exon_id = ${$trs_gene_info}{$inf_2}{$mRNA_name}{'exon'};
        my ( $exon_s, $exon_e ) =
          @{ $strand_table{$gene_name}{$mRNA_name}{'exon_list'}[ $exon_id - 1 ]
          };

        my $curr_mrn = 0;
        for my $n ( 0 .. $loop_len ) {
            if ( $strand eq '-' ) {
                $curr_mrn = $mrna_e - $n;
            }
            else {
                $curr_mrn = $mrna_s + $n;
            }
            my $curr_codon = ceil( $curr_mrn / 3.0 );

            my $shift_n = 0;
            if ( !in_range( $exon_s, $exon_e, $var_s ) ) {
                $shift_n = $diff_len + $n;
            }
            else {
                $shift_n = $n;
            }

            if ( $in_s_s[$shift_n] ne $in_s_e[$shift_n] ) {
                push @aoa,
                  [
                    $var_s + $shift_n, $in_s_s[$shift_n], $in_s_e[$shift_n],
                    $curr_mrn,         $curr_codon
                  ];
            }
        }
        my $last_ind = $aoa[-1][0];
        push @aoa,
          [ $last_ind + 1, $in_s_s[$mrn_len], $in_s_e[$mrn_len], $curr_mrn,
            -1 ];

        my $temp_key1 = $inf_1;

        for my $d ( 0 .. $diff_len - 1 ) {
            my $temp_key2 = "";
            my $shift_d   = 0;
            if ( !in_range( $exon_s, $exon_e, $var_s ) ) {
#    print "line 328\n"; #04032015
#    print Dumper $exon_s; #04032015
#    print Dumper $exon_e; #04032015
#    print Dumper $var_s; #04032015
                $shift_d = $d;
            }
            elsif ( !in_range( $exon_s, $exon_e, $var_e ) ) {
#    print "line 335\n"; #04032015
#    print Dumper $exon_s; #04032015
#    print Dumper $exon_e; #04032015
#    print Dumper $var_s; #04032015
                $shift_d = $d + $mrn_len + 1;
            }
            my $curr_coor = $var_s + $shift_d;
            $temp_key2 = "$chr.$curr_coor.$in_s_s[$shift_d].$in_s_e[$shift_d]";
            if ( $in_s_s[$shift_d] ne $in_s_e[$shift_d] ) {
                $aoi{$temp_key1}{$temp_key2} = 1;
            }
        }
#    print "line 347\n"; #04032015
#    print Dumper \@aoa; #04032015
#    print Dumper $#aoa; #04032015
        my $gid = 0;
        for my $m ( 0 .. $#aoa - 1 ) {
            my ( $pre, $nxt ) = ( $m - 1, $m + 1 );
            my $temp_key1 = $inf_1;
            my $temp_key2 = "$chr.$aoa[$m][0].$aoa[$m][1].$aoa[$m][2]";
            if ( $aoa[$m][-1] != $aoa[$pre][-1] ) {
                $gid += 1;
            }
            if ( ( $aoa[$m][-1] - $aoa[$pre][-1] ) *
                ( $aoa[$m][-1] - $aoa[$nxt][-1] ) == 0 )
            {

				if (exists $aono{$temp_key1}{$temp_key2}) {
					delete $aono{$temp_key1}{$temp_key2}
				}	
                $aoo{$temp_key1}{$temp_key2}{'functions'}
                  { $gene_name . ':' . $mRNA_name } =
                    @{ $aa_mu_T{ $aoa[$m][-1] } }[0]
                  . $aoa[$m][-1]
                  . @{ $aa_mu_T{ $aoa[$m][-1] } }[1];
                $aoo{$temp_key1}{$temp_key2}{'cluster'}{$gid.$aoa[$m][-1]} = 1; #04132015
                $cluster{$temp_key1}{$gid.$aoa[$m][-1]}{$temp_key2} = 1; #04132015

            }
            elsif ( not exists $aoo{$temp_key1}{$temp_key2} ) { #09/05
                $aono{$temp_key1}{$temp_key2}{ $gene_name . ':' . $mRNA_name }
                  = @{ $aa_mu_T{ $aoa[$m][-1] } }[0]
                  . $aoa[$m][-1]
                  . @{ $aa_mu_T{ $aoa[$m][-1] } }[1];
            }
        }
    }
    print "\n";
#    print "line 372\n"; #04032015
#    print Dumper \%aoi; #04032015
#    print Dumper \%cluster; #04032015
    return ( \%aoo, \%aono, \%aoi, \%cluster );
}

#-------------------------------------------------------------#

sub in_range {
    my $rs      = $_[0];
    my $re      = $_[1];
    my $nu      = $_[2];
    my $boolean = 0;
    if ( $rs <= $nu and $nu <= $re ) {
        $boolean = 1;
    }
    return $boolean;
}

#-------------------------------------------------------------#

sub process_input_variant_function {
    %ao_nonexon = ();
    open( INFILE, $input_filename );
    while ( my $line = <INFILE> ) {

        $line =~ s/[\r\n\s]+$//g;
        my @columns = split( /\t/, $line );
        my @temp =
          ( $columns[8], "$columns[2].$columns[3].$columns[5].$columns[6]" );
        $temp[0] =~ s/\s//g;
        $temp[1] =~ s/\s//g;
        if ( $columns[0] ne "exonic" ) {
            my $m_l = $columns[4] - $columns[3];
            foreach my $i ( 0 .. $m_l ) {
                $columns[2] =~ s/\s//g;
                $columns[5] =~ s/\s//g;
                $columns[6] =~ s/\s//g;
                my $s_nuc = substr( $columns[5], $i, 1 );
                my $e_nuc = substr( $columns[6], $i, 1 );
                my $c_pos = $i + $columns[3];
                $temp[1] = "$columns[2].$c_pos.$s_nuc.$e_nuc";
                $ao_nonexon{ $temp[0] }{ $temp[1] } = 1;
            }
        }
    }
    return ( \%ao_nonexon );
}

#-------------------------------------------------------------#

# this subroutine reads the snp annotation from
# the .exonic_variant_function file

sub process_exonic_variant_function_snp {
    my %aos = ();
    open( INFILE, $input_filename );
    my ( $num_line_total, $num_line_snp ) = (0) x 2;
    my $next_call = 0;
    while ( my $line = <INFILE> ) {
        $num_line_total += 1;
        $line =~ tr/[ \r\n]//ds;
        my @columns = split( /\t/, $line );
        if ( $columns[4] == $columns[5] && $columns[6] ne "-" ) {
            $num_line_snp += 1;
            my @function_ann = split( /\,/, $columns[2] );
            my %function_list = ();
            foreach my $i ( 0 .. $#function_ann ) {
                my @curr_ann = split( /\:/, $function_ann[$i] );
                my @amino_mu = split( /\./, $curr_ann[-1] );
                if ( defined $curr_ann[1] ) {
                    $function_list{"$curr_ann[0]:$curr_ann[1]"} = $amino_mu[-1];
                }
                else {
                    $function_list{"$curr_ann[0]"} = $amino_mu[-1];
                }
            }
            my $curr_mut = "$columns[3].$columns[4].$columns[6].$columns[7]";
            $aos{ $columns[-1] }{$curr_mut} = \%function_list;
        }
    }
    if ( $num_line_total > $num_line_snp ) {
        $next_call = 1;
    }
    return ( \%aos, $next_call );
}

#-----------------------------------------------------------#

sub find_overlapping_SNVs {
    pre_check_required_files( $_[1] );
    %cluster = ();
    my @chars = ( "A" .. "Z", "a" .. "z" );
    $Annotation::temp_filename = "input";
    $Annotation::temp_filename .= $chars[ rand @chars ] for 1 .. 15;
    $Annotation::temp_filename .= ".vcf";
    if ( defined $_[1] ) { $option = $_[1]{'method'} }
    my $original = Annotation::pre_annovar_process( $_[0], $_[1] );

    my ( $outa, $outb, $outi, $outc, $oute ) = ();
    if ( $option eq 'annovar' ) {
        Annotation::readRefGene( $_[1] );
        ( $outa, $outb, $outi, $outc, $oute ) = runAnnovar( $_[1] );
    }
    elsif ( $option eq 'snpeff' ) {
        Annotation::readEnsemblGTF( $_[1] );
        ( $outa, $outb, $outi, $outc ) = runSnpeff( $_[1] );
    }
    elsif ( $option eq 'vep' ) {
        Annotation::readEnsemblGTF( $_[1] );
        ( $outa, $outb, $outi, $outc ) = runVep( $_[1] );
    }
    else {
        die("option $option not implemented");
    }

#    foreach my $x ( keys %$outc ) {
#        foreach my $y ( keys ${$outc}{$x} ) {
#            ${$outb}{$x}{$y} = 1;
#        }
#    }
    foreach my $x ( keys %$oute ) {
        foreach my $y ( keys %{$oute->{$x}} ) {
            ${$outb}{$x}{$y} = ${$oute}{$x}{$y};
        }
    }

    # the following is to print out according to the desired data structure.
    my %output = ();

    foreach my $x ( keys %$original ) {
        my $line_id   = ${$original}{$x};
        my $curr_hash = {};
        push @{ ${$curr_hash}{'Overlapped'} }, ();
        ${$curr_hash}{'Non-exonic'} = ();
        ${$curr_hash}{'NonOverlapped'} = ();
        if ( exists ${$outb}{$line_id} ) {
            foreach my $variant ( keys %{$outb->{$line_id}} ) {
                ${$curr_hash}{'NonOverlapped'}{$variant} =
                  ${$outb}{$line_id}{$variant};
            }
        }
        if ( exists ${$outc}{$line_id} ) {
            foreach my $variant ( keys %{${$outc}->{$line_id}} ) {
                ${$curr_hash}{'Non-exonic'}{$variant} =
                  ${$outc}{$line_id}{$variant};
            }
        }
        if ( exists ${$outi}{$line_id} ) {
            foreach my $variant ( keys %{${$outi}->{$line_id}} ) {
                ${$curr_hash}{'Non-exonic'}{$variant} =
                  ${$outi}{$line_id}{$variant};
            }
        }
        if ( exists ${$outa}{$line_id} ) {
            my @overlap_arr = ();
            foreach my $c ( keys %{ $cluster{$line_id} } ) {
                my $curr_c = {};
                foreach my $v ( keys %{ $cluster{$line_id}{$c} } ) {
                    $curr_c->{'mutations'}->{$v} = 1;
                    $curr_c->{'functions'} = {};
                    foreach
                      my $f ( keys %{ ${$outa}{$line_id}{$v}{'functions'} } )
                    {
                        $curr_c->{'functions'}->{$f} =
                          ${$outa}{$line_id}{$v}{'functions'}{$f};
                    }
                }
                push @overlap_arr, $curr_c;
            }
            ${$curr_hash}{'Overlapped'} = \@overlap_arr;
        }

        $output{$x} = $curr_hash;
    }

    return \%output;
}

#-------------------------------------------------------------#

sub readRefGene {
    %strand_table = ();
    my $refGene = $_[0]{'annovar_refgene'}
      or die(
"Please provide the full path to gene reference file: hg19_refGene.txt\n"
      );
    open( INFILE, $refGene );

    while ( my $line = <INFILE> ) {
        $line =~ s/[\r\n]+$//g;
        my @columns   = split( /\t/, $line );
        my @exonstart = split( /,/,  $columns[9] );
        my @exonend   = split( /,/,  $columns[10] );
        my ( $gene, $mRna, $strand, $num_exon ) =
          ( $columns[12], $columns[1], $columns[3], $columns[8] );

        $strand_table{$gene}{$mRna}{'strand'} = $strand;

        my @arr_exon_coor = ();
        for my $e ( 0 .. $#exonstart ) {
            if ( $strand eq "+" ) {
                $arr_exon_coor[$e] = [ $exonstart[$e], $exonend[$e] ];
            }
            else {
                my $ec = $#exonstart - $e;
                $arr_exon_coor[$e] = [ $exonstart[$ec], $exonend[$ec] ];
            }
        }
        $strand_table{$gene}{$mRna}{'exon_list'} = \@arr_exon_coor;
    }
    close INFILE;
}

#-------------------------------------------------------------#

sub readEnsemblGTF {
    %strand_table = ();
    my $refGene = $_[0]{'ensembl_Homo_sapiens_GRCh37_75_gtf'}
      or die(
"Please provide the full path to the Ensembl file: ensembl_Homo_sapiens.GRCh37.75.gtf"
      );
    open( INFILE, $refGene );

    my @last_exon_posi = ();
    while ( my $line = <INFILE> ) {
        $line =~ s/[\r\n]+$//g;
        my @columns = split( /\t/, $line );
        my ( $chr, $feature, $start, $end, $strand ) =
          ( $columns[0], $columns[2], $columns[3], $columns[4], $columns[6] );
        my ( $gene, $mRna, $exon_rank ) = ();

        $columns[-1] =~ m/^.*transcript_id "([^"]*)".*$/;
        $mRna = $1;
        $columns[-1] =~ m/^.*gene_name "([^"]*)".*$/;
        $gene = $1;
        $columns[-1] =~ m/^.*exon_number "([^"]*)".*$/;
        $exon_rank = $1;
        if ( defined $mRna && defined $gene ) {
            if ( $feature eq "transcript" ) {
                $strand_table{$gene}{$mRna}{'strand'} = $strand;
                @{ $strand_table{$gene}{$mRna}{'position'} } = ( $start, $end );
            }
            elsif ( $feature eq "exon" ) {
                @last_exon_posi = ( $start, $end );
                if ( defined $exon_rank ) {
                    $strand_table{$gene}{$mRna}{'exon_list'}{$start}{'end_pos'}
                      = $end;
                    $strand_table{$gene}{$mRna}{'exon_list'}{$start}{'rank'} =
                      $exon_rank;
                }
            }
            else {
                @{ $strand_table{$gene}{$mRna}{'exon_list'}
                      { $last_exon_posi[0] }{'details'}{$feature} } =
                  ( $start, $end );
            }
        }
    }
    close INFILE;
}

#-----------------------------------------------------------#

sub runVep {
    my $script   = $_[0]{'vep_path_to_variant_effect_predictor'};
    my $VCF_BASE = $Annotation::temp_filename;
    system "perl $script --offline --cache --format vcf -i $VCF_BASE --force_overwrite --everything -o $VCF_BASE.out";
    system "grep -v \"#\" $VCF_BASE.out > $VCF_BASE.noheader.out";
    $Annotation::input_filename = "$VCF_BASE.noheader.out";
    my ( $outa, $outb, $outi, $outc ) = Annotation::parse_vep_output();
    system "rm -f $temp_filename*";
    return ( $outa, $outb, $outi, $outc );
}

#-----------------------------------------------------------#

sub parse_vep_output {
    ( %aoo, %aono, %aoi, %ao_nonexon, %cluster ) = ();
    %cluster = ();
    open( INFILE, $input_filename );
	my $firstline = <INFILE>;
 	close INFILE;
 	my @firstline_columns = split(/\t/, $firstline); 
 	my $prevline_inf_1 = $firstline_columns[0]; 
 	my ($prevline_chr, $prevline_var_p) = split(/:/, $firstline_columns[1]);
 	my ($prevline_var_s, $prevline_var_e) = split(/-/, $prevline_var_p);
 	$prevline_var_e = (defined $prevline_var_e)? $prevline_var_e : $prevline_var_s;
 	open(INFILE, $input_filename);
    while ( my $line = <INFILE> ) {
        $line =~ tr/[\r\n]//ds;
        my @columns = split( /\t/, $line );
        my $inf_1 = $columns[0];
        my ( $chr,   $var_p ) = split( /:/, $columns[1] );
        my ( $var_s, $var_e ) = split( /-/, $var_p );
        $var_e = ( defined $var_e ) ? $var_e : $var_s;
        my $var_len = $var_e - $var_s;
        my ( $temp_s, $temp_e ) =
          @{ $mapping_vep{$inf_1}{$chr}{$var_s}{$var_e} };
        my @in_s_s = split( "", $temp_s );
        my @in_s_e = split( "", $temp_e );
        my ( $gene_id, $mRNA_name ) = ( $columns[3], $columns[4] );
        $columns[-1] =~ m/[^STRAND]*STRAND=([-\d]*);.*/;
        my $strand = $1;
        $columns[-1] =~ m/[^SYMBOL]*SYMBOL=(\w*);.*/;
        my $gene_name = $1;
        my ( $cds_posi, $amin_posi, $amin_record ) =
          ( $columns[8], $columns[9], $columns[10] );

        if ( $amin_record ne "-" ) {
            my ( $mrna_s, $mrna_e ) = split( /-/, $cds_posi );
            $mrna_e = ( defined $mrna_e ) ? $mrna_e : $mrna_s;
            my ( $aa_ids, $aa_ide ) = split( /-/, $amin_posi );
            $aa_ide = ( defined $aa_ide ) ? $aa_ide : $aa_ids;
            my ( $aa_ref, $aa_mut ) = split( /\//, $amin_record );
            my @ref_aa_l = split( "", $aa_ref );
 			my @mut_aa_l = ();
 			if (defined $aa_mut) {
 				@mut_aa_l = split("", $aa_mut);
 			} else {
 				@mut_aa_l = @ref_aa_l;
 			}
            my %aa_mu_T  = ();
            foreach ( $aa_ids .. $aa_ide ) {

                if ( defined $mut_aa_l[ $_ - $aa_ids ] ) {
                    $aa_mu_T{$_} =
                      [ $ref_aa_l[ $_ - $aa_ids ], $mut_aa_l[ $_ - $aa_ids ] ];
                }
                else {
                    $aa_mu_T{$_} = [ $ref_aa_l[ $_ - $aa_ids ], "-" ];
                }
            }

            #create hash table
            my @aoa      = ();
            my $mrn_len  = $mrna_e - $mrna_s;
            my $diff_len = $var_len - $mrn_len;
            my $loop_len = ();

            if ( $diff_len >= 0 ) {
                $loop_len = $mrn_len;
            }
            else {
                $loop_len = $var_len;
            }
            my ( $exon_s, $exon_e ) = ();
            my $currexon_e = ();
            if ( defined $strand_table{$gene_name}{$mRNA_name}{'exon_list'} ) {
                foreach my $currexon_s ( sort { $a <=> $b }
                    keys %{ $strand_table{$gene_name}{$mRNA_name}{'exon_list'} }
                  )
                {
                    $currexon_e =
                      $strand_table{$gene_name}{$mRNA_name}{'exon_list'}
                      {$currexon_s}{'end_pos'};
                    if ( in_range( $currexon_s, $currexon_e, $var_s ) ) {
                        ( $exon_s, $exon_e ) = ( $currexon_s, $currexon_e );
                        if (
                            defined $strand_table{$gene_name}{$mRNA_name}
                            {'exon_list'}{$currexon_s}{'details'} )
                        {
                            ( $exon_s, $exon_e ) =
                              @{ $strand_table{$gene_name}{$mRNA_name}
                                  {'exon_list'}{$currexon_s}{'details'}{'CDS'}
                              };
                        }
                    }
                }
            }
            my $curr_mrn = 0;
            for my $n ( 0 .. $loop_len ) {
                if ( $strand eq '-' ) {
                    $curr_mrn = $mrna_e - $n;
                }
                else {
                    $curr_mrn = $mrna_s + $n;
                }
                my $curr_codon = ceil( $curr_mrn / 3.0 );

                my $shift_n = 0;
                if ( !in_range( $exon_s, $exon_e, $var_s ) ) {
                    $shift_n = $diff_len + $n;
                }
                else {
                    $shift_n = $n;
                }

                if ( $in_s_s[$shift_n] ne $in_s_e[$shift_n] ) {
                    push @aoa,
                      [
                        $var_s + $shift_n, $in_s_s[$shift_n],
                        $in_s_e[$shift_n], $curr_mrn,
                        $curr_codon
                      ];
                }
            }
            my $last_ind = $aoa[-1][0];
            push @aoa,
              [
                $last_ind + 1,     $in_s_s[$mrn_len],
                $in_s_e[$mrn_len], $curr_mrn,
                -1
              ];
            my $temp_key1 = $inf_1;

            for my $d ( 0 .. $diff_len - 1 ) {
                my $temp_key2 = "";
                my $shift_d   = 0;
                if ( !in_range( $exon_s, $exon_e, $var_s ) ) {
                    $shift_d = $d;
                }
                elsif ( !in_range( $exon_s, $exon_e, $var_e ) ) {
                    $shift_d = $d + $mrn_len + 1;
                }
                my $curr_coor = $var_s + $shift_d;
                $temp_key2 =
                  "$chr.$curr_coor.$in_s_s[$shift_d].$in_s_e[$shift_d]";
                if ( $in_s_s[$shift_d] ne $in_s_e[$shift_d] ) {
                    $aoi{$temp_key1}{$temp_key2} = 1;
                }
            }
            my $gid = 0;
            for my $m ( 0 .. $#aoa - 1 ) {
                my ( $pre, $nxt ) = ( $m - 1, $m + 1 );
                my $temp_key1 = $inf_1;
                my $temp_key2 = "$chr.$aoa[$m][0].$aoa[$m][1].$aoa[$m][2]";

                if ( $aoa[$m][-1] != $aoa[$pre][-1] ) {
                    $gid += 1;
                }
                if ( ( $aoa[$m][-1] - $aoa[$pre][-1] ) *
                    ( $aoa[$m][-1] - $aoa[$nxt][-1] ) == 0 )
                {
                    if ( defined $gene_name && defined $mRNA_name ) {
						if (exists $aono{$temp_key1}{$temp_key2}) {
							delete $aono{$temp_key1}{$temp_key2}
						}	
                        $aoo{$temp_key1}{$temp_key2}{'functions'}
                          { $gene_name . ':' . $mRNA_name } =
                            @{ $aa_mu_T{ $aoa[$m][-1] } }[0]
                          . $aoa[$m][-1]
                          . @{ $aa_mu_T{ $aoa[$m][-1] } }[1];
                        $aoo{$temp_key1}{$temp_key2}{'cluster'}{$gid} = 1;
                        $cluster{$temp_key1}{$gid}{$temp_key2} = 1;
                    }
                }
                else {
#                    if ( defined $gene_name && defined $mRNA_name ) {
                    if ( defined $gene_name && defined $mRNA_name && ( not exists $aoo{$temp_key1}{$temp_key2} ) ) { #09/05
                        $aono{$temp_key1}{$temp_key2}
                          { $gene_name . ':' . $mRNA_name } =
                            @{ $aa_mu_T{ $aoa[$m][-1] } }[0]
                          . $aoa[$m][-1]
                          . @{ $aa_mu_T{ $aoa[$m][-1] } }[1];
                    }
                }
            }
        }
 		if ($inf_1 != $prevline_inf_1 || $chr != $prevline_chr || $var_s != $prevline_var_s || $var_e != $prevline_var_e || eof) {
 			my $prevline_var_len = $prevline_var_e - $prevline_var_s;
 			my ($prevline_temp_s, $prevline_temp_e) = @{$mapping_vep{$prevline_inf_1}{$prevline_chr}{$prevline_var_s}{$prevline_var_e}};
 			my @prevline_in_s_s = split("", $prevline_temp_s);
 			my @prevline_in_s_e = split("", $prevline_temp_e);
 			if (eof) {
 				$prevline_var_len = $var_e - $var_s;
 				($prevline_temp_s, $prevline_temp_e) = @{$mapping_vep{$inf_1}{$chr}{$var_s}{$var_e}};
 				@prevline_in_s_s = split("", $temp_s);
 				@prevline_in_s_e = split("", $temp_e);
 			}
 			for my $i (0..$prevline_var_len) {
 				my $curr_coor = $prevline_var_s + $i;
 				my $temp_key2 = "$prevline_chr.$curr_coor.$prevline_in_s_s[$i].$prevline_in_s_e[$i]";
 				if (exists $aoo{$prevline_inf_1}{$temp_key2}) {
 				} elsif (exists $aono{$prevline_inf_1}{$temp_key2}) {
 				} elsif (exists $aoi{$prevline_inf_1}{$temp_key2}) {
 				} else {
 					$ao_nonexon{$prevline_inf_1}{$temp_key2} = 1;
 				}	
 			}
 			$prevline_inf_1 = $inf_1;
 			$prevline_chr = $chr;
 			$prevline_var_s = $var_s;
 			$prevline_var_e = $var_e;
 		}	
 	}
  	close INFILE;
  	return (\%aoo, \%aono, \%aoi, \%ao_nonexon); #added on 09/04/2014
}

#-----------------------------------------------------------#

sub runSnpeff {
    system
"cat $Annotation::temp_filename | vcf-sort > ${Annotation::temp_filename}_sorted.vcf ";
    my $INPUT_VCF = $Annotation::temp_filename . "_sorted.vcf";
    my $REFERENCE = "GRCh37.75";
    my $MEM       = "4G"
      ;    # Amount of memory to use (make sure there is enough physical memory)
    my $SNPEFF = "java -Xmx$MEM -jar " . $_[0]{'snpeff_path_to_snpEff_jar'}
      or die("Please provide the path to the snpEff.jar");
    my $VCF_BASE = $Annotation::temp_filename . "_sorted";
    system "$SNPEFF -v -hgvs $REFERENCE  $INPUT_VCF > $VCF_BASE.hgvs.eff.vcf ";
    system
      "grep -v \"##\" $VCF_BASE.hgvs.eff.vcf > $VCF_BASE.noheader.hgvs.eff.vcf";
    $Annotation::input_filename = "$VCF_BASE.noheader.hgvs.eff.vcf";
    my ( $outa, $outb, $outi, $outc ) = Annotation::parse_snpeff_output();
    system "rm -f $temp_filename*";
    system "rm -f snpEff_*";
    return ( $outa, $outb, $outi, $outc );
}

#-----------------------------------------------------------#

sub parse_snpeff_output {
    ( %aoo, %aono, %aoi, %ao_nonexon, %cluster ) = ();
    open( INFILE, $input_filename );
    while ( my $line = <INFILE> ) {
        $line =~ tr/[ \r\n]//ds;
        my @columns = split( /\t/, $line );
        my $inf_1   = $columns[6];
        my $chr     = $columns[0];
        my ( $var_s, $var_e ) = ( $columns[1], $columns[2] );
#        my $var_len = $var_e - $var_s;
        my $var_len = $var_e - $var_s + 1; ##04_01_2015
        my @in_s_s  = split( "", $columns[3] );
        my @in_s_e  = split( "", $columns[4] );
        $columns[-1] =~ tr/EFF=//;
        my @effs = split( /\,/, $columns[-1] );
        my $size_effs = $#effs;
		my $count_exon = 0;

        for my $i ( 0 .. $size_effs ) {

            # get info of each transcript
            my @trs_info = split( /\|/, $effs[$i] );
            if ( $trs_info[3] ne "" && index( $trs_info[3], "p." ) != -1 ) {
				$count_exon += 1;
                my $gene_name = $trs_info[5];
                my $mRNA_name = $trs_info[8];
                $trs_info[3] =~ m/p.(\D*)(\d*)(\D*)\/c.(\d*)(\w*)>(\w*)/;
                my ( $ref_aa, $aa_ids, $mut_aa, $mrna_s, $mrna_ref, $mrna_mut )
                  = ( $1, $2, $3, $4, $5, $6 );
                my $mrna_e   = $mrna_s + length($mrna_ref) - 1;
                my $strand   = $strand_table{$gene_name}{$mRNA_name}{'strand'};
                my @ref_aa_l = ( $ref_aa =~ m/.../g );
                my @mut_aa_l = ( $mut_aa =~ m/.../g );
                my $aa_size  = $#ref_aa_l;
                my $aa_ide   = $aa_ids + $aa_size;
                my %aa_mu_T  = ();

print "line 920\n"; #04_01_2015
print Dumper $ref_aa; #04_01_2015
print Dumper \@ref_aa_l; #04_01_2015
print Dumper $aa_ids; #04_01_2015
print Dumper $aa_ide; #04_01_2015
                foreach ( $aa_ids .. $aa_ide ) {
                    if ( defined $mut_aa_l[ $_ - $aa_ids ] ) {
                        $aa_mu_T{$_} = [
                            $ref_aa_l[ $_ - $aa_ids ],
                            $mut_aa_l[ $_ - $aa_ids ]
                        ];
                    }
                    else {
                        $aa_mu_T{$_} = [ $ref_aa_l[ $_ - $aa_ids ], "-" ];
                    }
                }

                #create hash table
                my @aoa      = ();
                my $mrn_len  = length($mrna_ref);
                my $diff_len = $var_len - $mrn_len;

				print "line 934\n";  #04_01_2015
				print Dumper $var_len;  #04_01_2015
				print Dumper $mrn_len;  #04_01_2015
                my $loop_len = ();

                if ( $diff_len >= 0 ) {
                    $loop_len = $mrn_len;
                }
                else {
                    $loop_len = $var_len;
                }

                # find the corresponding exon range
                # evaluate which side of the variant is located in intron
                # if the left side then
                # else if the right side
                my $exon_id = $trs_info[9];
                my ( $exon_s, $exon_e ) = ();
                my $currexon_e = ();
                if (
                    defined $strand_table{$gene_name}{$mRNA_name}{'exon_list'} )
                {
                    foreach my $currexon_s ( sort { $a <=> $b }
                        keys
                        %{ $strand_table{$gene_name}{$mRNA_name}{'exon_list'} }
                      )
                    {
                        $currexon_e =
                          $strand_table{$gene_name}{$mRNA_name}{'exon_list'}
                          {$currexon_s}{'end_pos'};
print "line 965\n"; #04_04_2015
print Dumper $currexon_s;
print Dumper $currexon_e;
print Dumper $var_s;
#                        if ( in_range( $currexon_s, $currexon_e, $var_s ) ) { #04_01_2015
                        if ( in_range( $currexon_s, $currexon_e, $var_s ) or in_range( $currexon_s, $currexon_e, $var_e) ) { #04_01_2015
print "line 977 in_range\n"; #04_04_2015
                            ( $exon_s, $exon_e ) = ( $currexon_s, $currexon_e );
                            if (
                                defined $strand_table{$gene_name}{$mRNA_name}
                                {'exon_list'}{$currexon_s}{'details'} )
                            {
                                ( $exon_s, $exon_e ) =
                                  @{ $strand_table{$gene_name}{$mRNA_name}
                                      {'exon_list'}{$currexon_s}{'details'}
                                      {'CDS'} };
                            }
                        }
                    }
                }

print "line 985\n"; #04_04_2015
print Dumper $exon_s;
print Dumper $exon_e;
                my $curr_mrn = 0;
print "line 997\n"; #04_04_2015
print Dumper $loop_len;
print Dumper \@in_s_s;
print Dumper \@in_s_e;
#                for my $n ( 0 .. $loop_len ) { #04_01_2015
                for my $n ( 0 .. $loop_len-1 ) {	#04_01_2015
                    if ( $strand eq '-' ) {
                        $curr_mrn = $mrna_e - $n;
                    }
                    else {
                        $curr_mrn = $mrna_s + $n;
                    }
                    my $curr_codon = ceil( $curr_mrn / 3.0 );

                    my $shift_n = 0;
print "line 1011\n"; #04_04_2015
print Dumper $exon_s;
print Dumper $exon_e;
                    if ( !in_range( $exon_s, $exon_e, $var_s ) ) {
                        $shift_n = $diff_len + $n;
                    }
                    else {
                        $shift_n = $n;
                    }

                    if ( $in_s_s[$shift_n] ne $in_s_e[$shift_n] ) {
                        push @aoa,
                          [
                            $var_s + $shift_n, $in_s_s[$shift_n],
                            $in_s_e[$shift_n], $curr_mrn,
                            $curr_codon
                          ];
                    }
                }
                my $last_ind = $aoa[-1][0];
                push @aoa,
                  [
#                    $last_ind + 1,     $in_s_s[$mrn_len],
#                    $in_s_e[$mrn_len], $curr_mrn,
                    $last_ind + 1,     "O",
                    "O", $curr_mrn,
                    -1
                  ];
                my $temp_key1 = $inf_1;

                for my $d ( 0 .. $diff_len - 1 ) {
                    my $temp_key2 = "";
                    my $shift_d   = 0;
print "line 1025\n"; #04_04_2015
print Dumper $exon_s;
print Dumper $exon_e;
                    if ( !in_range( $exon_s, $exon_e, $var_s ) ) {
                        $shift_d = $d;
                    }
                    elsif ( !in_range( $exon_s, $exon_e, $var_e ) ) {
print "line 1032\n"; #04_04_2015
print Dumper $exon_s;
print Dumper $exon_e;
                        $shift_d = $d + $mrn_len + 1;
                    }
                    my $curr_coor = $var_s + $shift_d;
                    $temp_key2 =
                      "$chr.$curr_coor.$in_s_s[$shift_d].$in_s_e[$shift_d]";
                    if ( $in_s_s[$shift_d] ne $in_s_e[$shift_d] ) {
                        $aoi{$temp_key1}{$temp_key2} = 1;
                    }
                }
                my $gid = 0;
                for my $m ( 0 .. $#aoa - 1 ) {
                    my ( $pre, $nxt ) = ( $m - 1, $m + 1 );
                    my $temp_key1 = $inf_1;
                    my $temp_key2 = "$chr.$aoa[$m][0].$aoa[$m][1].$aoa[$m][2]";

                    if ( $aoa[$m][-1] != $aoa[$pre][-1] ) {
                        $gid += 1;
                    }
                    if ( ( $aoa[$m][-1] - $aoa[$pre][-1] ) *
                        ( $aoa[$m][-1] - $aoa[$nxt][-1] ) == 0 )
                    {
                        if ( defined $gene_name && defined $mRNA_name ) {
							if (exists $aono{$temp_key1}{$temp_key2}) {
								delete $aono{$temp_key1}{$temp_key2}
							}	
                            $aoo{$temp_key1}{$temp_key2}{'functions'}
                              { $gene_name . ':' . $mRNA_name } =
                                @{ $aa_mu_T{ $aoa[$m][-1] } }[0]
                              . $aoa[$m][-1]
                              . @{ $aa_mu_T{ $aoa[$m][-1] } }[1];
                            $aoo{$temp_key1}{$temp_key2}{'cluster'}{$gid} = 1;
                            $cluster{$temp_key1}{$gid}{$temp_key2} = 1;
                        }
                    }
                    else {
                        if ( defined $gene_name && defined $mRNA_name && ( not exists $aoo{$temp_key1}{$temp_key2} ) ) { #09/05
print "line 1085\n"; #04_01_2015
print Dumper \%aa_mu_T; #04_01_2015
print Dumper $m; #04_01_2015
print Dumper \@aoa; #04_01_2015
print Dumper $aoa[$m]; #04_01_2015
print Dumper $aoa[$m][-1]; #04_01_2015
print Dumper $aa_mu_T{$aoa[$m][-1]}; #04_01_2015
print Dumper $temp_key1; #04_01_2015
print Dumper $temp_key2; #04_01_2015
							if (defined $aa_mu_T{$aoa[$m][-1]} ) {	#04_01_2015
	                            $aono{$temp_key1}{$temp_key2}
	                              { $gene_name . ':' . $mRNA_name } =
	                                @{ $aa_mu_T{ $aoa[$m][-1] } }[0]
	                              . $aoa[$m][-1]
	                              . @{ $aa_mu_T{ $aoa[$m][-1] } }[1];
							}
							else {#04_01_2015
								$ao_nonexon{$temp_key1}{$temp_key2} = 1;#04_01_2015
							}#04_01_2015
                        }
                    }
                }
            }
        }
 		if ($count_exon == 0) {
# 			for my $i (0..$var_len) { #04_01_2015	
 			for my $i (0..$var_len-1) {#04_01_2015	 	
 				my $temp_key2 = ""; 
 				my $curr_coor = $var_s + $i; 
 				$temp_key2 = "$chr.$curr_coor.$in_s_s[$i].$in_s_e[$i]"; 
 				if ($in_s_s[$i] ne $in_s_e[$i]) { 					
 					$ao_nonexon{$inf_1}{$temp_key2} = 1; 		
 				} 											
 			} 											
 		} 											
	}
    return ( \%aoo, \%aono, \%aoi, \%ao_nonexon );
}

#-----------------------------------------------------------#

sub runAnnovar {
    do {

        local @ARGV;
        my $annovar_db_dir = `dirname $_[0]{'annovar_refgene'}`;
        $annovar_db_dir =~ tr/[\r\n]//ds;
        @ARGV = (
            "-buildver", "hg19", $Annotation::temp_filename, $annovar_db_dir,
            "--exonicsplicing"
        );
        my $temp = capture( $^X, $_[0]{'annovar_annotate_variation'}, @ARGV )
          ;  #or die ("Please provide the full path to annovar_variation.pl\n");
    };

    $Annotation::input_filename =
      $Annotation::temp_filename . ".variant_function";
    my $outc = Annotation::process_input_variant_function();

    $Annotation::input_filename =
      $Annotation::temp_filename . ".exonic_variant_function";
    my ( $oute, $next_call ) =
      Annotation::process_exonic_variant_function_snp();

    my ( $outa, $outb, $outi, $clusters ) = ();
    if ( -z "$temp_filename.exonic_variant_function" ) {
        print "File $temp_filename.exonic_variant_function is empty";
    }
    elsif ($next_call) {
        my @out_2;
        do {
            local @ARGV;
            @ARGV = (
                "$temp_filename.exonic_variant_function",
                $_[0]{'annovar_refgene'},
                $_[0]{'annovar_refmrna'}
              )
              or
              die("Please provide the full path to refgene and refmrna files.");
            @out_2 = capture( $^X, $_[0]{'annovar_coding_change'}, @ARGV ); #04032015              or die("Please provide the full path to annovar script: coding_change.pl");
            open my $of_1, ">", "predict_annovar.txt"
              or die("could not open file. $!");
            print $of_1 @out_2;
            close $of_1;
            system
"grep \">\" predict_annovar.txt | grep -v \"WILDTYPE\" | sed -e 's/>//g' > predict_annovar_short.txt";
        };
        my $sizef = -s "predict_annovar.txt";
        if ( $sizef < 48 ) {
            print
              "File predict.txt is empty. Only single SNV in this sample. \n";
        }
        else {
            ( $outa, $outb, $outi, $clusters ) =
              Annotation::parse_predicted_annovar_output_Original();
        }
    }

#    system "rm -f $temp_filename*";
#    system "rm -f predict*.txt";
    return ( $outa, $outb, $outi, $outc, $oute );
}

#-----------------------------------------------------------#

sub pre_check_required_files {
    if ( defined $_[0]{'method'} ) {
        if ( $_[0]{'method'} eq "annovar" ) {
            my $refGene = $_[0]{'annovar_refgene'}
              or die(
"Please provide the full path to gene reference file: hg19_refGene.txt\n"
              );
            my $annovar_script_1 = $_[0]{'annovar_annotate_variation'}
              or die("Please provide the full path to annovar_variation.pl\n");
            my $refmrna = $_[0]{'annovar_refmrna'}
              or die("Please provide the full path to hg19_refGeneMrna.fa\n");
            my $annovar_script_2 = $_[0]{'annovar_coding_change'}
              or die(
"Please provide the full path to annovar script: coding_change.pl\n"
              );

        }
        elsif ( $_[0]{'method'} eq "snpeff" ) {
            my $vcfexit = `which vcf-sort`
              or die("Please install vcf-sort shell command.\n");
            my $refGene = $_[0]{'ensembl_Homo_sapiens_GRCh37_75_gtf'}
              or die(
"Please provide the full path to the Ensembl file: ensembl_Homo_sapiens.GRCh37.75.gtf\n"
              );
            my $SNPEFF = $_[0]{'snpeff_path_to_snpEff_jar'}
              or die("Please provide the path to the snpEff.jar\n");

        }
        elsif ( $_[0]{'method'} eq "vep" ) {
            my $refGene = $_[0]{'ensembl_Homo_sapiens_GRCh37_75_gtf'}
              or die(
"Please provide the full path to the Ensembl file: ensembl_Homo_sapiens.GRCh37.75.gtf\n"
              );
            my $script = $_[0]{'vep_path_to_variant_effect_predictor'}
              or die(
"Please provide the full path to the Vep script: variant_effect_    predictor.pl\n"
              );
        }
        else {
            die("Unknown annotation methods and not implemented.\n");
        }
    }
    else {
        die(
"Please provide with an annotation method: annovar, snpeff, or vep.\n"
        );
    }
}

#-----------------------------------------------------------#

sub output_QA {
	
	my $output_hash_ref = shift;
	foreach my $variants (keys %$output_hash_ref) {

		my  (%input_hash, %output_hash) = ();
		my @columns = split(/,/, $variants);
		my $size_variant = $#columns +1;
		if (defined ${$output_hash_ref}{$variants}{'Non-exonic'}) {
			foreach my $item (keys % {$output_hash_ref->{$variants}->{'Non-exonic'}} ) {
#				$output_hash{$item} = 1;
				(exists $output_hash{$item})? ($output_hash{$item} += 1) : ($output_hash{$item} = 1);
			}
		}
		if (defined ${$output_hash_ref}{$variants}{'NonOverlapped'}) {
			foreach my $item (keys % {$output_hash_ref->{$variants}->{'NonOverlapped'}} ) {
#				$output_hash{$item} = 1;
				(exists $output_hash{$item})? ($output_hash{$item} += 1) : ($output_hash{$item} = 1);
			}
		}
		if (defined ${$output_hash_ref}{$variants}{'Overlapped'}) {
			my $size_cluster = @{$output_hash_ref->{$variants}->{'Overlapped'}};
			foreach my $i (0..$size_cluster-1) {
#				foreach my $item (keys % {@{$output_hash_ref->{$variants}->{'Overlapped'}}[$i]}->{'mutations'}) {
				foreach my $item (keys %{$output_hash_ref->{$variants}{'Overlapped'}[$i]{'mutations'}}) {
#					$output_hash{$item} = 1;
					(exists $output_hash{$item})? ($output_hash{$item} += 1) : ($output_hash{$item} = 1);
				}
			}
		}
		foreach my $i (0..$#columns) {
			$input_hash{$columns[$i]} = 1;
		}

		if (%input_hash ne %output_hash) {
			print "line1233--WARNING: input and output don't match\n";
			print Dumper %input_hash;
			print Dumper %output_hash;
		} else {
			my %cmp = map { $_ => 1 } keys %input_hash;
			for my $key (keys %output_hash) {
				print "within check loop\n";
				print Dumper $key;
				print Dumper $output_hash{$key};
				print Dumper $cmp{$key};
				last unless exists $cmp{$key};
				delete $cmp{$key};
			}
			if (%cmp) {
				print "line1245--WARNING: input and output don't match\n";
				print Dumper \%cmp; #04032015
			}
		}
	}
	return 1;
}

#-----------------------------------------------------------#
1;
