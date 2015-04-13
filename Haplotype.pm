package Haplotype;

use strict;
use Carp;
use FileHandle;
use Storable;
use Storable qw(dclone);
use Data::Dumper;
use FindBin;
use lib "$FindBin::Bin/";
use General;

sub identify_haplotypes_by_mutation_and_reads {
    my ( $r_mut2read, $opts ) = @_;
    croak "FATAL: you must provide hash ref r_mut2read $r_mut2read"
      unless ref($r_mut2read) eq 'HASH';
    my $r_read2mut;
    for my $sample ( keys %{$r_mut2read} ) {
        for my $id_m ( keys %{ $r_mut2read->{$sample} } ) {
            for my $read_id ( keys %{ $r_mut2read->{$sample}{$id_m} } ) {
                $r_read2mut->{$sample}{$read_id}{$id_m}++;
            }
        }
    }
    my $r_cluster;
    for my $sample ( sort keys %{$r_mut2read} ) {
        my $r_reads_used;
        while ( keys %{ $r_mut2read->{$sample} } ) {
            my $id_m1 = ( keys %{ $r_mut2read->{$sample} } )[0];
            my $t     = {};
            $t->{$id_m1} = $r_mut2read->{$sample}{$id_m1};
            delete( $r_mut2read->{$sample}{$id_m1} );
            my %reads_new = map { $_ => 1 } keys %{ $t->{$id_m1} };
            my %reads_used;
            while ( keys %reads_new ) {
                my $read_id = ( keys %reads_new )[0];
                $reads_used{$read_id}++;
                delete( $reads_new{$read_id} );
                for my $id_m2 ( keys %{ $r_read2mut->{$sample}{$read_id} } ) {
                    if ( General::deep_defined( $r_mut2read, $sample, $id_m2 ) )
                    {
                        $t->{$id_m2} = $r_mut2read->{$sample}{$id_m2};
                        delete( $r_mut2read->{$sample}{$id_m2} );
                        for my $read_id2 ( keys %{ $t->{$id_m2} } ) {
                            $reads_new{$read_id2}++
                              unless defined( $reads_used{$read_id2} );
                        }
                    }
                }
            }
            my $ta;
            if ( defined( $opts->{'max_allowed_adjacent_distance'} ) ) {
                $ta =
                  separate_SNV_Indel_cluster_by_max_allowed_adjacent_distance(
                    $t, $opts->{'max_allowed_adjacent_distance'} );
            }
            else {
                push( @{$ta}, $t );
            }
            for my $t1 ( @{$ta} ) {
                next unless scalar keys %{$t1} > 1;
                my $t2;
                my @cluster_all_ids = sort keys %{$t1};
                my $cluster_joint_id = join( ",", @cluster_all_ids );
                for my $read_id3 ( keys %reads_used ) {
                    my @joint_status;
                    for my $id_m3 (@cluster_all_ids) {
                        my $status;
                        if ( General::deep_defined( $t1, $id_m3, $read_id3 ) ) {
                            my $status_raw =
                              scalar keys %{ $t1->{$id_m3}{$read_id3} } == 1
                              ? ( keys %{ $t1->{$id_m3}{$read_id3} } )[0]
                              : croak
"FATAL: multiple statuses found for id_m3 $id_m3 read_id3 $read_id3";
                            if ( $status_raw eq 'spt' ) {
                                $status = 1;
                            }
                            elsif ( $status_raw eq 'agt' ) {
                                $status = 0;
                            }
                            elsif ( $status_raw eq 'unknown' ) {
                                $status = '-';
                            }
                            else {
                                croak
"FATAL: unknown status_raw $status_raw for id_m3 $id_m3 read_id3 $read_id3";
                            }
                        }
                        else {
                            $status = '-';
                        }
                        push( @joint_status, $status );
                    }
                    my $joint_status = join( "", @joint_status );
					next if $joint_status=~/^-+$/;
                    $t2->{'joint_status'}{$joint_status}{$read_id3}++;
                }
                $t2->{'ids'} = \@cluster_all_ids;
                $r_cluster->{$sample}{$cluster_joint_id} = $t2;
            }
        }
    }
    return $r_cluster;
}

sub separate_SNV_Indel_cluster_by_max_allowed_adjacent_distance {
    my ( $r_in, $max_allowed_adjacent_distance ) = @_;
    croak "FATAL: input r must be a hash ref" unless ref($r_in) eq 'HASH';
    croak
"FATAL: max_allowed_adjacent_distance must be a non-negative number -$max_allowed_adjacent_distance-"
      unless defined($max_allowed_adjacent_distance)
      and $max_allowed_adjacent_distance >= 0;
    my $r1 = dclone($r_in);
    my $r1_by_chr;
    for my $id ( keys %{$r1} ) {
        my ( $chr, $posi, $ref, $mut ) = split /\./, $id;
        croak "FATAL: incorrect id for SNV/Indels -$id-"
          unless $posi =~ /^\d+$/;
        $r1_by_chr->{$chr}{$id} = $r1->{$id};
        delete( $r1->{$id} );
    }
    my $ra_out;
    my $r1_by_chr_posi;
    for my $chr ( sort keys %{$r1_by_chr} ) {
        my @ids =
          sort { ( split /\./, $a )[1] <=> ( split /\./, $b )[1] }
          keys %{ $r1_by_chr->{$chr} };
        my %t;
        my $id_last;
        for my $id (@ids) {
            if ( defined($id_last)
                and ( split /\./, $id )[1] - ( split /\./, $id_last )[1] >
                $max_allowed_adjacent_distance )
            {
                push( @{$ra_out}, dclone( \%t ) );
                %t = ();
            }
            $t{$id} = $r1_by_chr->{$chr}{$id};
            delete( $r1_by_chr->{$chr}{$id} );
            $id_last = $id;
        }
        push( @{$ra_out}, \%t );
    }
    return $ra_out;
}

sub extract_mutations_from_haplotype_cluster {
    my ( $r_in, $r_mut ) = @_;
    croak "FATAL: must provide input mutation hash r_mut"
      unless ref($r_mut) eq 'HASH';
    for my $key (qw/ids joint_status/) {
        croak "FATAL: cannot find expected key $key in input hash"
          unless General::deep_defined( $r_in, $key );
    }
    croak "FATAL: cannot find expected array ref"
      unless ref( $r_in->{ids} ) eq 'ARRAY';
    my $r_out;
    my $r_mut_order;
    for my $i ( 0 .. $#{ $r_in->{'ids'} } ) {
        my $mut = $r_in->{'ids'}[$i];
        if ( defined( $r_mut_order->{$mut} ) ) {
            croak
"FATAL: ordered mutation ids must be unique but containing duplicate mut $mut";
        }
        else {
            $r_mut_order->{$mut} = $i;
        }
    }
    my @mut_indexes_selected;
    for ( keys %{$r_mut} ) {
        croak "FATAL: r_mut contains mutations that is not in r_in $_"
          unless defined( $r_mut_order->{$_} );
        my $index = $r_mut_order->{$_};
        push( @mut_indexes_selected, $index );
    }
    my @mut_indexes_selected_sorted = sort { $a <=> $b } @mut_indexes_selected;
    my @cluster_ids_selected;
    for my $i (@mut_indexes_selected_sorted) {
        my $mut = $r_in->{'ids'}[$i];
        push( @cluster_ids_selected, $mut );
    }
    for my $joint_status ( keys %{ $r_in->{'joint_status'} } ) {
        my @statuses = split //, $joint_status;
        croak
"FATAL: joint_status $joint_status mutation N not equal to the number of ids"
          unless scalar @statuses == scalar @{ $r_in->{'ids'} };
        my @statuses_selected =
          map { $statuses[$_] } @mut_indexes_selected_sorted;
        my $joint_status_selected = join( "", @statuses_selected );
        for my $read_id ( keys %{ $r_in->{'joint_status'}{$joint_status} } ) {
            $r_out->{'joint_status'}{$joint_status_selected}{$read_id}++;
        }
    }
    $r_out->{'ids'} = \@cluster_ids_selected;
    return $r_out;
}

sub extract_present_mutations_from_haplotypes {
    my ( $combined_mutations, $haplotype ) = @_;
    croak "FATAL: you must provide combined_mutations $combined_mutations"
      unless defined($combined_mutations);
    croak
"FATAL: input haplotype $haplotype must only contain 1 or 0 but not unknown "
      unless $haplotype =~ /^[01]+$/;
    my @indexes   = split //,   $haplotype;
    my @mutations = split /\,/, $combined_mutations;
    print Dumper \@mutations
      and print Dumper \@indexes
      and croak
"FATAL: input haplotype $haplotype must have the same length as the mutation array"
      unless $#indexes == $#mutations;
    my @out;
    for my $i ( 0 .. $#indexes ) {
        my $status = $indexes[$i];
        my $mut    = $mutations[$i];
        push( @out, $mut ) if $status == 1;
    }
    my $combined_mutations_present = join( ",", @out );
    return $combined_mutations_present;
}

1;
