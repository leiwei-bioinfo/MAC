package General;
use strict;
use Carp;
use FileHandle;
use Storable;
use Storable qw(dclone);
use Data::Dumper;

sub open_bam {
    my %ref_db;
    my ( $bam, $fasta ) = @_;
    my $sam = Bio::DB::Sam->new(
        -bam   => $bam,
        -fasta => $fasta
    );
    return $sam;
}

sub mutation_format_to_type {
    my ($in) = @_;
    if ( $in =~ /^[0-9A-Z_]+\.\d+\.([ATGCN]+)\.([ATGCN]+)$/i ) {
        if ( length($1) == length($2) ) {
            length($1) == 1 ? return 'SNV' : return 'MNV';
        }
        else {
            return 'Indel';
        }
    }
    elsif ( $in =~ /^[0-9A-Z_]+\.\d+\.[ATGCN-]+\.[ATGCN-]+$/i ) {
        return 'Indel';
    }
    else {
        croak "cannot figure out type for mutation id $in";
    }
}

sub get_allele_count {
    my ( $sam, $chr, $posi, $min_base_qual) = @_;
    croak "FATAL: incorrect input bam filehandle $sam ref ref($sam)"
      unless ref($sam) eq 'Bio::DB::Sam';
    croak "FATAL: incorrect position $posi" unless $posi =~ /^\d+$/;
    my @alignments = $sam->get_features_by_location(
        -seqid => $chr,
        -start => $posi,
        -end   => $posi
    );
    my ( %h, %h2 );
    for my $read (@alignments) {
        my $read_id    = $read->query->seq_id;
        my $self_start = $read->start;
        my $cigar      = $read->cigar_str;
        my $flag       = $read->flag;
        next if ( $flag & 0x0400 );
        my $read_strand   = ( $flag & 0x0010 ) ? 1 : 0;
        my $self_sequence = $read->query->dna;
        my @qual          = $read->qscore;
        my $match_qual    = $read->qual;
        my $attributes    = $read->aux;
        my ( $bases_match_ref, $base_shift_ref ) =
          parse_cigar( $self_start, $cigar, $self_sequence, \@qual );
        my %bases = %{$bases_match_ref};
        next unless defined( $bases{$posi}{'base'} );
        my $t = {
            'cigar'    => $cigar,
            'flag'     => $flag,
            'strand'   => $read_strand,
            'start'    => $self_start,
            'seq'      => $self_sequence,
            'qual_ref' => \@qual,
            'mapQ'     => $match_qual,
        };
        push( @{ $h{'broad'}{$read_id}{ $bases{$posi}{'base'} } }, $t );

        if ( $bases{$posi}{'qual'} >= $min_base_qual ) {
            push( @{ $h{'hq'}{$read_id}{ $bases{$posi}{'base'} } }, $t );
        }
    }
    for my $quality ( keys %h ) {
        for my $read_id ( keys %{ $h{$quality} } ) {
            next if ( scalar keys %{ $h{$quality}{$read_id} } > 1 );
            my $base = [ keys %{ $h{$quality}{$read_id} } ]->[0];
            $h2{$quality}{$base}{$read_id} = $h{$quality}{$read_id}{$base};
        }
    }
    return ( \%h2 );
}

sub parse_cigar {
    my ( $posi_start, $cigar_ori, $string_seq, $qual_info ) = @_;
    croak "FATAL: position contains non numbers $posi_start"
      unless $posi_start =~ /^\d+$/;
    my $cigar = $cigar_ori;
    my ( %positions, %base_shift );
    my ( $ref_posi, $string_posi ) = ( $posi_start, 0 );
    while ( $cigar =~ m/^(\d+)([MIDNSHPX=])(.*)$/ ) {
        my ( $number, $letter ) = ( $1, $2 );
        $cigar = $3;
        if ( $letter eq 'D' || $letter eq 'N' ) {
            $base_shift{$letter}{$ref_posi} = $number;
        }
        for ( my $i = 1 ; $i <= $number ; $i++ ) {
            if ( $letter eq 'M' || $letter eq '=' || $letter eq 'X' ) {
                $positions{$ref_posi}{'base'} =
                  substr( $string_seq, $string_posi, 1 );
                if ( ref($qual_info) eq "ARRAY" ) {
                    $positions{$ref_posi}{'qual'} = $qual_info->[$string_posi];
                }
                elsif ( $qual_info =~ /[!-~]+|\*/ ) {
                    $positions{$ref_posi}{'qual'} =
                      ord( substr( $qual_info, $string_posi, 1 ) ) - 33;
                }
                else {
                    croak
                      "FATAL: unrecognizable quality information $qual_info";
                }
                $ref_posi    += 1;
                $string_posi += 1;
            }
            elsif ( $letter eq 'I' ) {
                $base_shift{'I'}{$ref_posi}{'base'} .=
                  substr( $string_seq, $string_posi, 1 );
                my $qual_t =
                  ref($qual_info) eq "ARRAY"
                  ? $qual_info->[$string_posi]
                  : ord( substr( $qual_info, $string_posi, 1 ) ) - 33;
                push( @{ $base_shift{'I'}{$ref_posi}{'qual'} }, $qual_t );
                $string_posi += 1;
            }
            elsif ( $letter eq 'D' ) {
                $ref_posi += 1;
            }
            elsif ( $letter eq 'N' ) {
                $ref_posi += 1;
            }
            elsif ( $letter eq 'S' ) {
                $string_posi += 1;
            }
            elsif ( $letter eq 'H' ) {
            }
            elsif ( $letter eq 'P' ) {
            }
        }
    }
    my $bases_match_ref = \%positions;
    my $base_shift_ref  = \%base_shift;
    return ( $bases_match_ref, $base_shift_ref );
}

sub deep_defined {
    my ( $ref, @keys ) = @_;
    unless (@keys) {
        warn "deep_defined: no keys";
        return;
    }
    foreach my $key (@keys) {
        if ( ref $ref eq 'HASH' ) {
            return unless defined( $ref->{$key} );
            $ref = $ref->{$key};
            next;
        }
        if ( ref $ref eq 'ARRAY' ) {
            return unless 0 <= $key && $key < @{$ref};
            return unless defined( $ref->[$key] );
            $ref = $ref->[$key];
            next;
        }
        return;
    }
    return 1;
}

sub QA_input_output_one_to_one {
    my ( $r_in, $r_out ) = @_;
    croak "input hash reference is not an array ref"
      unless ref($r_in) eq 'HASH';
    croak "output hash reference is not an array ref"
      unless ref($r_out) eq 'HASH';
    for my $i ( keys %{$r_in} ) {
        croak "FATAL: input and output are not one-to-one, input uniq: $i"
          unless defined( $r_out->{$i} );
    }
    for my $i ( keys %{$r_out} ) {
        croak "FATAL: input and output are not one-to-one, output uniq: $i"
          unless defined( $r_in->{$i} );
    }
}

sub w_write {
    my ( $out_ref, $outfile, $option ) = @_;
    $option = 0 unless $option;
    croak "not an array ref, cannot write" unless ref($out_ref) eq 'ARRAY';
    croak "you must specify an output file to write" unless $outfile;
    my $bk_file = $outfile . "_bk_" . RandomString();
    my $md5_old = '';
    croak "bk_file already exist, unlikely with random string $bk_file"
      if -e $bk_file;
    if ( -e $outfile ) {
        croak "output file already exist, cannot use option 0" if $option == 0;
        if ( $option == 3 || $option == 4 ) {
            $md5_old = md5_file($outfile);
            rename $outfile, $bk_file;
        }
    }
    my $fh_out;
    if ( $option == 0 || $option == 1 || $option == 3 || $option == 4 ) {
        $fh_out = FileHandle->new("> $outfile")
          or croak "cannot write to output file $outfile";
    }
    elsif ( $option == 2 ) {
        $fh_out = FileHandle->new(">> $outfile")
          or croak "cannot write to output file $outfile";
    }
    else {
        croak "unrecognized option $option";
    }
    for ( @{$out_ref} ) {
        chomp;
        print $fh_out "$_\n";
    }
    if ( -e $bk_file ) {
        my $md5_new = md5_file($outfile);
        if ( $md5_old eq $md5_new ) {
            if ( $option == 3 ) {
                print
"output opt #3 - keep older output, remove newer outfile with the same md5\n";
                unlink($outfile);
                rename $bk_file, $outfile;
            }
            elsif ( $option == 4 ) {
                print
"output opt #4 - keep newer output, remove older outfile with the same md5\n";
                unlink($bk_file);
            }
            else {
                croak "FATAL: Strange error, option is not 3 or 4, -$option-";
            }
        }
    }
    $fh_out->close;
}

sub RandomString {
    my ($l) = @_;
    $l = 10 unless $l;
    return join( "", map { ( "a" .. "z" )[ rand 26 ] } ( 1 .. $l ) );
}

sub md5_file {
    my ($file) = @_;
    croak "file -$file- not exist" unless $file && -e $file;
    my $md5 = `md5sum $file|cut -f1 -d' '`;
    chomp $md5;
    return $md5;
}

1;
