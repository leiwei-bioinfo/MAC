#!/usr/bin/env perl

use strict;
use warnings;
use Bio::DB::Sam;
use Carp;
use Getopt::Long;
use FileHandle;
use Data::Dumper;
use Storable;
use Storable qw(dclone);
use FindBin;
use lib "$FindBin::Bin/";
use Pod::Usage;
use IPC::System::Simple qw(system capture);
use General;
use Haplotype;
use Annotation;

my ( $infile, $outprefix, $bam, $reference, $o_print_incomplete_haplotype,
    $max_allowed_adjacent_distance, $min_base_qual,
    $printlog, $annotator );
my (
    $annovar_refgene, $annovar_annotate_variation,
    $annovar_refmrna, $annovar_coding_change
);
my ( $ensembl_Homo_sapiens_GRCh37_75_gtf, $snpeff_path_to_snpEff_jar );
my ($vep_path_to_variant_effect_predictor);
my ( $man, $usage, $help );
GetOptions(

    #minimum requirements
    'i=s'   => \$infile,
    'bam=s' => \$bam,
    'r=s'   => \$reference,

    #optional arguments
    'annotator=s'                  => \$annotator,
    'annovar_annotate_variation=s' => \$annovar_annotate_variation,
    'annovar_coding_change=s'      => \$annovar_coding_change,
    'annovar_refgene=s'            => \$annovar_refgene,
    'annovar_refmrna=s'            => \$annovar_refmrna,
    'ensembl_Homo_sapiens_GRCh37_75_gtf=s' =>
      \$ensembl_Homo_sapiens_GRCh37_75_gtf,
    'snpeff_path_to_snpEff_jar=s' => \$snpeff_path_to_snpEff_jar,
    'vep_path_to_variant_effect_predictor=s' =>
      \$vep_path_to_variant_effect_predictor,
    'man'                             => \$man,
    'usage'                           => \$usage,
    'help'                            => \$help,
    'o=s'                             => \$outprefix,
    'print_incomplete_haplotype!'     => \$o_print_incomplete_haplotype,
    'max_allowed_adjacent_distance=i' => \$max_allowed_adjacent_distance,
	'min_base_qual=i'				  => \$min_base_qual,
    'printlog!'                       => \$printlog,
);

pod2usage( -verbose => 1 ) if ($help);
pod2usage( -verbose => 2 ) if ( $man or $usage );

croak "FATAL: you must provide infile -$infile-" unless $infile && -e $infile;
croak "FATAL: you must provide bam -$bam-"       unless $bam    && -e $bam;
croak "FATAL: you must provide genome reference file -$reference-"
  unless $reference && -e $reference;
$outprefix = $infile . "_" unless defined($outprefix);
#$max_allowed_adjacent_distance = 2
#  unless defined($max_allowed_adjacent_distance);
$min_base_qual=20 unless defined($min_base_qual);

my $opts;
$opts->{'max_allowed_adjacent_distance'} = $max_allowed_adjacent_distance;

my @log;

#s1: load inputs
print "Load input SNVs, open bam...\n";
my $fh_bam = General::open_bam( $bam, $reference );
my $r_mut;
my $fh_in = FileHandle->new($infile);
while (<$fh_in>) {
    chomp;
    croak "FATAL: input must be SNV, please check $_"
      unless General::mutation_format_to_type($_) eq 'SNV';
    $r_mut->{SAMPLE}{$_}++;
}
$fh_in->close;

#log
my $n = scalar keys %{$r_mut};
push( @log, "s1 total mutations:\t$n\n" );

#step2: add raw reads
print "Add raw reads for each SNV...\n";
my $r_mut2read;
for my $id ( keys %{ $r_mut->{SAMPLE} } ) {
    my ( $chr, $posi, $ref, $mut ) = split /\./, $id;
    my $t = General::get_allele_count( $fh_bam, $chr, $posi, $min_base_qual);
    my $t1;
    for my $base ( keys %{ $t->{'hq'} } ) {
        for my $read_id ( keys %{ $t->{'hq'}{$base} } ) {
            for my $t2 ( @{ $t->{'hq'}{$base}{$read_id} } ) {
                my $flag = $t2->{'flag'};
                $t1->{$read_id}{$base}{$flag}++;
            }
        }
    }
    for my $read_id ( keys %{$t1} ) {
        my $combined_status;
        if ( scalar keys %{ $t1->{$read_id} } > 1 ) {
            $combined_status = 'unknown';
        }
        else {
            my $combined_base = ( keys %{ $t1->{$read_id} } )[0];
            $combined_status = $combined_base eq $mut ? 'spt' : 'agt';
        }
        for my $base ( keys %{ $t1->{$read_id} } ) {
            my $status = $base eq $mut ? 'spt' : 'agt';
            for my $flag ( keys %{ $t1->{$read_id}{$base} } ) {
                $r_mut2read->{SAMPLE}{$id}{$read_id}{$combined_status}{$flag}
                  {$status}++;
            }
        }
    }
}

#step3: Look for Block of Mutations
print "Look for Block of Mutations...\n";
my $t_r_mut_all;
for my $sample ( keys %{$r_mut2read} ) {
    for my $id_m ( keys %{ $r_mut2read->{$sample} } ) {
        $t_r_mut_all->{$id_m}++;
    }
}

my $r_cluster_all =
  Haplotype::identify_haplotypes_by_mutation_and_reads( $r_mut2read, $opts );

#output raw haplotype counts
my @header_raw =
  qw/BM_mutations Haplotype_index Joint_status N_reads Annotation/;
my @out_raw = keys %{$r_cluster_all} ? ( join( "\t", @header_raw ) ) : ();
for my $sample ( sort keys %{$r_cluster_all} ) {
    for my $block_joint_id_all ( sort keys %{ $r_cluster_all->{$sample} } ) {
        my $index_haplotype = 1;
        my $N_haplotypes    = scalar keys
          %{ $r_cluster_all->{$sample}{$block_joint_id_all}{'joint_status'} };
        for my $haplotype (
            sort { $b cmp $a }
            keys
            %{ $r_cluster_all->{$sample}{$block_joint_id_all}{'joint_status'} }
          )
        {
            my $function;
            if ( $haplotype =~ /-/ ) {
                $function = 'Incomplete';
            }
            elsif ( $haplotype =~ /^0+$/ ) {
                $function = 'Wt';
            }
            else {
                $function = 'Mutant';
            }
            my $n_reads =
              scalar keys
              %{ $r_cluster_all->{$sample}{$block_joint_id_all}{'joint_status'}
                  {$haplotype} };
            my $out = join(
                "\t",
                (
                    $block_joint_id_all, $index_haplotype . "/" . $N_haplotypes,
                    $haplotype, $n_reads, $function
                )
            );
            $index_haplotype++;
            push( @out_raw, $out );
        }
    }
}

#log
my @log_cluster_mut_N;
my @log_cluster_mut_full;
my $total_mut_N = 0;
for my $sample ( keys %{$r_cluster_all} ) {
    if ( General::deep_defined( $r_cluster_all, $sample ) ) {
        for my $cluster_joint_id ( keys %{ $r_cluster_all->{$sample} } ) {
            my $n =
              scalar @{ $r_cluster_all->{$sample}{$cluster_joint_id}{'ids'} };
            push( @log_cluster_mut_N, $n );
            my @joint_status;
            for my $joint_status (
                sort { $b cmp $a } keys %{
                    $r_cluster_all->{$sample}{$cluster_joint_id}{'joint_status'}
                }
              )
            {
                my $read_n =
                  scalar keys %{ $r_cluster_all->{$sample}{$cluster_joint_id}
                      {'joint_status'}{$joint_status} };
                push( @joint_status, $joint_status . "_" . $read_n );
            }
            push( @log_cluster_mut_full,
                    $sample . "\t"
                  . $n . "\t"
                  . $cluster_joint_id . "\t"
                  . join( "|", @joint_status ) );
            $total_mut_N += $n;
        }
    }
}
my $N_cluster = scalar @log_cluster_mut_N;
push( @log, "total clusters $N_cluster\ttotal mutations $total_mut_N" );
General::w_write( \@log_cluster_mut_N, 'log_cluster_mut_N.txt', 1 )
  if $printlog;
General::w_write( \@log_cluster_mut_full, 'log_cluster_mut_full.txt', 1 )
  if $printlog;

unless ( defined($annotator) ) {
    my $outfile2 = $outprefix . "raw_haplotype_counts.txt";
    print
"No annotator specified, output all identified raw MNVs without annotation...\nOutput:\t$outfile2\n";
    General::w_write( \@out_raw, $outfile2, 1 );
    exit 0;
}

#step4: for each mutation, find out which ones share the same codon
my $opts2;
$opts2->{'assembly'} = $reference;
$opts2->{'method'}   = $annotator;
if ( $annotator eq 'annovar' ) {
    $opts2->{'annovar_annotate_variation'} = $annovar_annotate_variation;
    $opts2->{'annovar_coding_change'}      = $annovar_coding_change;
    $opts2->{'annovar_refgene'}            = $annovar_refgene;
    $opts2->{'annovar_refmrna'}            = $annovar_refmrna;
}
elsif ( $annotator eq 'snpeff' ) {
    $opts2->{'snpeff_path_to_snpEff_jar'} = $snpeff_path_to_snpEff_jar;
    $opts2->{'ensembl_Homo_sapiens_GRCh37_75_gtf'} =
      $ensembl_Homo_sapiens_GRCh37_75_gtf;
}
elsif ( $annotator eq 'vep' ) {
    $opts2->{'vep_path_to_variant_effect_predictor'} =
      $vep_path_to_variant_effect_predictor;
    $opts2->{'ensembl_Homo_sapiens_GRCh37_75_gtf'} =
      $ensembl_Homo_sapiens_GRCh37_75_gtf;
}
else {
    croak "unknown annotator $annotator";
}

print "Find codon info for each Block of Mutation...\n";
my $r_all_combined_ids =
  { map { $_ => 1 }
      ( map { keys %{ $r_cluster_all->{$_} } } keys %{$r_cluster_all} ) };
print "\n\tstart:\trunning $annotator to find BMC\n";
my $r_all_combined_ids_codons =
  Annotation::find_overlapping_SNVs( $r_all_combined_ids, $opts2 );
Annotation::output_QA($r_all_combined_ids_codons);
print "\tfinished:\trunning $annotator to find BMC\n\n";
General::QA_input_output_one_to_one( $r_all_combined_ids,
    $r_all_combined_ids_codons );
my $t_r_QA_1;

for my $cluster_id ( keys %{$r_all_combined_ids_codons} ) {
    if ( keys %{ $r_all_combined_ids_codons->{$cluster_id}{'NonOverlapped'} } )
    {
        for my $mut_id (
            keys %{ $r_all_combined_ids_codons->{$cluster_id}{'NonOverlapped'} }
          )
        {
            if ( General::deep_defined( $t_r_QA_1, $cluster_id, $mut_id ) ) {
                print "t_r_QA_1->cluster_id\n";
                print Dumper $t_r_QA_1->{$cluster_id};
                print "r_all_combined_ids_codons->cluster_id\n";
                print Dumper $r_all_combined_ids_codons->{$cluster_id};
                croak
"FATAL: mut_id $mut_id appear more than once in cluster_id $cluster_id";
            }
            else {
                $t_r_QA_1->{$cluster_id}{$mut_id}++;
            }
        }
    }
    if (
        scalar @{ $r_all_combined_ids_codons->{$cluster_id}{'Overlapped'} } >
        0 )
    {
        for
          my $t ( @{ $r_all_combined_ids_codons->{$cluster_id}{'Overlapped'} } )
        {
            for my $mut_id ( keys %{ $t->{'mutations'} } ) {
                if ( General::deep_defined( $t_r_QA_1, $cluster_id, $mut_id ) )
                {
                    print "t_r_QA_1->cluster_id\n";
                    print Dumper $t_r_QA_1->{$cluster_id};
                    print "r_all_combined_ids_codons->cluster_id\n";
                    print Dumper $r_all_combined_ids_codons->{$cluster_id};
                    croak
"FATAL: mut_id $mut_id appear more than once in cluster_id $cluster_id";
                }
                else {
                    $t_r_QA_1->{$cluster_id}{$mut_id}++;
                }
            }
        }
    }
}

#step5: extract the overlapping(codons) in $r_all_combined_ids_codons from original $r_cluster_all
print
"Prepare haplotypes for each Block of Mutation inside Codon for annotator $annotator\n";
my $r_cluster_codon;
my $r_combined_codon_haplotype_ids = {};
for my $sample ( keys %{$r_cluster_all} ) {
    for my $cluster_joint_id ( keys %{ $r_cluster_all->{$sample} } ) {
        croak
"FATAL: cannot find cluster_joint_id $cluster_joint_id in r_all_combined_ids_codons"
          unless General::deep_defined( $r_all_combined_ids_codons,
            $cluster_joint_id );
        if (
            scalar
            @{ $r_all_combined_ids_codons->{$cluster_joint_id}{'Overlapped'} }
            > 0 )
        {
            for my $rh_codon (
                @{
                    $r_all_combined_ids_codons->{$cluster_joint_id}
                      {'Overlapped'}
                }
              )
            {
                croak
                  "FATAL: cannot find hashkey word of mutations in input hash"
                  unless General::deep_defined( $rh_codon, 'mutations' );
                my $t_out = Haplotype::extract_mutations_from_haplotype_cluster(
                    $r_cluster_all->{$sample}{$cluster_joint_id},
                    $rh_codon->{'mutations'} );
                my $cluster_joint_id_codon =
                  defined( $t_out->{'ids'} )
                  ? join( ",", @{ $t_out->{'ids'} } )
                  : croak
"FATAL: cannot find the key ids in output of Haplotype::extract_mutations_from_haplotype_cluster";
                $r_cluster_codon->{$sample}{$cluster_joint_id_codon} = $t_out;
                for my $joint_status_codon_haplotype (
                    keys %{ $t_out->{'joint_status'} } )
                {
                    next
                      if $joint_status_codon_haplotype =~ /-/
                      or $joint_status_codon_haplotype =~ /^0+$/;
                    my $combined_mutations_present =
                      Haplotype::extract_present_mutations_from_haplotypes(
                        $cluster_joint_id_codon,
                        $joint_status_codon_haplotype
                      );
                    $r_combined_codon_haplotype_ids
                      ->{$combined_mutations_present} = 1;
                }
            }
        }
    }
}
warn "No MNV overlapping codon identified in this data set\n"
  unless keys %{$r_combined_codon_haplotype_ids};

#log
my @log_cluster_mut_N_s6;
my $total_mut_N_s6 = 0;
for my $sample ( keys %{$r_cluster_codon} ) {
    if ( General::deep_defined( $r_cluster_codon, $sample ) ) {
        for my $cluster_joint_id ( keys %{ $r_cluster_codon->{$sample} } ) {
            my $n =
              scalar @{ $r_cluster_codon->{$sample}{$cluster_joint_id}{'ids'} };
            push( @log_cluster_mut_N_s6, $n );
            $total_mut_N_s6 += $n;
        }
    }
}
my $N_cluster_s6 = @log_cluster_mut_N_s6 ? scalar @log_cluster_mut_N_s6 : 0;
push( @log,
    "total clusters_s6 $N_cluster_s6\ttotal mutations_s6 $total_mut_N_s6" );

#step6: rerun Annovar for the codon haplotypes
print "Predict haplotype amino acid change using $annotator\n";
print "\n\tstart:\trunning $annotator to annotate haplotypes\n";
my $r_combined_codon_haplotype_ids_anno =
  Annotation::find_overlapping_SNVs( $r_combined_codon_haplotype_ids, $opts2 );
Annotation::output_QA($r_combined_codon_haplotype_ids_anno);
print "\tfinished:\trunning $annotator to annotate haplotypes\n\n";
General::QA_input_output_one_to_one( $r_combined_codon_haplotype_ids,
    $r_combined_codon_haplotype_ids_anno );
unless ($o_print_incomplete_haplotype) {
    print "clean incomplete haplotypes\n";
    for my $sample ( keys %{$r_cluster_codon} ) {
        for my $cluster_joint_id_codon ( keys %{ $r_cluster_codon->{$sample} } )
        {
            for my $joint_status_codon_haplotype (
                keys %{
                    $r_cluster_codon->{$sample}{$cluster_joint_id_codon}
                      {'joint_status'}
                }
              )
            {
                delete( $r_cluster_codon->{$sample}{$cluster_joint_id_codon}
                      {'joint_status'}{$joint_status_codon_haplotype} )
                  if $joint_status_codon_haplotype =~ /-/;
            }
        }
    }
}

#step7: output
print "Output...\n";
my @header = qw/BMC_mutations Haplotype_index Joint_status N_reads Annotation/;
my @out =
  keys %{$r_combined_codon_haplotype_ids} ? ( join( "\t", @header ) ) : ();
for my $sample ( sort keys %{$r_cluster_codon} ) {
    for
      my $cluster_joint_id_codon ( sort keys %{ $r_cluster_codon->{$sample} } )
    {
        my $index_haplotype = 1;
        my $N_haplotypes =
          scalar keys %{ $r_cluster_codon->{$sample}{$cluster_joint_id_codon}
              {'joint_status'} };
        for my $joint_status_codon_haplotype (
            sort { $b cmp $a } keys %{
                $r_cluster_codon->{$sample}{$cluster_joint_id_codon}
                  {'joint_status'}
            }
          )
        {
            my $function;
            if ( $joint_status_codon_haplotype =~ /-/ ) {
                $function = 'Unknown';
            }
            elsif ( $joint_status_codon_haplotype =~ /^0+$/ ) {
                $function = 'Wt';
            }
            else {
                my $combined_mutations_present =
                  Haplotype::extract_present_mutations_from_haplotypes(
                    $cluster_joint_id_codon, $joint_status_codon_haplotype );
                my $t;
                if (
                    scalar @{
                        $r_combined_codon_haplotype_ids_anno
                          ->{$combined_mutations_present}{'Overlapped'}
                    } == 1
                  )
                {
                    for my $mrna (
                        sort keys %{
                            $r_combined_codon_haplotype_ids_anno
                              ->{$combined_mutations_present}{'Overlapped'}[0]
                              {'functions'}
                        }
                      )
                    {
                        my $aachange =
                          $r_combined_codon_haplotype_ids_anno
                          ->{$combined_mutations_present}{'Overlapped'}[0]
                          {'functions'}{$mrna};
                        push( @{ $t->{$mrna} }, $aachange );
                    }
                }
                elsif (
                    scalar @{
                        $r_combined_codon_haplotype_ids_anno
                          ->{$combined_mutations_present}{'Overlapped'}
                    } > 1
                  )
                {
                    print Dumper $r_combined_codon_haplotype_ids_anno;
                    croak
"FATAL: pre-selected codon_haplotype did not yield unique output in r_combined_codon_haplotype_ids_anno for combined_mutations_present -$combined_mutations_present-";
                }
                if (
                    General::deep_defined(
                        $r_combined_codon_haplotype_ids_anno,
                        $combined_mutations_present,
                        'NonOverlapped'
                    )
                  )
                {
                    for my $sinlge_mutation (
                        sort keys %{
                            $r_combined_codon_haplotype_ids_anno
                              ->{$combined_mutations_present}{'NonOverlapped'}
                        }
                      )
                    {
                        for my $mrna (
                            sort keys %{
                                $r_combined_codon_haplotype_ids_anno
                                  ->{$combined_mutations_present}
                                  {'NonOverlapped'}{$sinlge_mutation}
                            }
                          )
                        {
                            my $aachange =
                              $r_combined_codon_haplotype_ids_anno
                              ->{$combined_mutations_present}{'NonOverlapped'}
                              {$sinlge_mutation}{$mrna};
                            push( @{ $t->{$mrna} }, $aachange );
                        }
                    }
                }
                my @t;
                for my $mrna ( sort keys %{$t} ) {
                    my $aachange_joint = join( "|", @{ $t->{$mrna} } );
                    push( @t, $mrna . ":" . $aachange_joint );
                }
                $function = @t ? join( ",", @t ) : 'NonCoding';
            }
            my $n_reads =
              scalar
              keys %{ $r_cluster_codon->{$sample}{$cluster_joint_id_codon}
                  {'joint_status'}{$joint_status_codon_haplotype} };
            my $out = join(
                "\t",
                (
                    $cluster_joint_id_codon,
                    $index_haplotype . "/" . $N_haplotypes,
                    $joint_status_codon_haplotype,
                    $n_reads,
                    $function
                )
            );
            $index_haplotype++;
            push( @out, $out );
        }
    }
}

my $outfile1 = $outprefix . "updated_haplotype_anno_by_" . $annotator . ".txt";
General::w_write( \@out, $outfile1, 1 );
print "Output file with annotation: $outfile1\n";
my $logfile = $outprefix . "log.txt";
General::w_write( \@log, $logfile, 1 ) if $printlog;

=head1 NAME

	MAC_v1.0.pl - Multi-nucleotide Variation Annotation Corrector


=head1 VERSION

	This documentation refers to MAC_v1.0.pl version 1.0.


=head1 USAGE

	#run without annotation:
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta

	#run with one of the three precompiled annotation programs
	#Annovar
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta -annotator annovar -annovar_annotate_variation FULLPATH/annotate_variation.pl -annovar_coding_change FULLPATH/coding_change.pl -annovar_refgene FULLPATH/hg19_refGene.txt -annovar_refmrna FULLPATH/hg19_refGeneMrna.fa
	
	#SnpEff
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta -annotator snpeff -snpeff_path_to_snpEff FULLPATH/snpEff.jar -ensembl_Homo_sapiens_GRCh37_75_gtf FULLPATH/Homo_sapiens.GRCh37.75.gtf

	#VEP
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta -annotator vep -vep_path_to_variant_effect_predictor FULLPATH/variant_effect_predictor.pl -ensembl_Homo_sapiens_GRCh37_75_gtf FULLPATH/Homo_sapiens.GRCh37.75.gtf

=head1 MINIMUM REQUIRED ARGUMENTS

	-i: input SNV file, each line contains one SNV in the format of "chr.posi.ref.mut", e.x. "17.7577084.T.A"
	-bam: input BAM file
	-r: reference genome sequence file in fasta format, must be the same one as the input BAM file

=head1 OPTIONS

	#Annotation related options
		#First select annotator
				-annotator						"annovar", "snpeff" or "vep"

		#Then provide full paths to required files and programs for working with each annotator
			#Annovar
				-annovar_annotate_variation		FULLPATH/annotate_variation.pl					(source: Annovar Package)
				-annovar_coding_change 			FULLPATH/coding_change.pl						(source: Annovar Package)
				-annovar_refgene 				FULLPATH/hg19_refGene.txt 						(source: Annovar Package)
				-annovar_refmrna 				FULLPATH/hg19_refGeneMrna.fa					(source: Annovar Package)

			#SnpEff
				-snpeff_path_to_snpEff			FULLPATH/snpEff.jar								(source: SnpEff Package)
				-ensembl_Homo_sapiens_GRCh37_75_gtf	FULLPATH/Homo_sapiens.GRCh37.75.gtf			(source: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)

			#VEP
				-vep_path_to_variant_effect_predictor	FULLPATH/variant_effect_predictor.pl	(source: VEP Package)
				-ensembl_Homo_sapiens_GRCh37_75_gtf	FULLPATH/Homo_sapiens.GRCh37.75.gtf			(source: ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz)
				
	#Other options
				-o:output prefix, default is to use input SNV file name as prefix.
				-max_allowed_adjacent_distance: max allowed distance between two adjacent SNVs in one Block of Mutation
				--print_incomplete_haplotype: select this option to keep haplotypes containing unknown status in annotation
				--printlog: print number counts at each step for debugging purpose
			    --man: Print the man page.
			    --usage: Print usage information.

=head1 DESCRIPTION

	This is a program to identify multiple nucleotide variation (MNV) from a list of 
	user-provided SNVs and the matching BAM file. If an annotation program is selected,
	the program will further determine which MNVs may contain multiple bases within the
	same codon, and re-annotate such MNV on a haplotype basis.


=head1 DEPENDENCIES

	The program depend on several packages:
	1.  Perl v5 (tested with Perl v5.16.1)
    2.  Bio::DB::Sam (http://search.cpan.org/~lds/Bio-SamTools/), which also requires SamTools library (http://sourceforge.net/projects/samtools/files/).
    3.  IPC::System::Simple (http://search.cpan.org/~pjf/IPC-System-Simple-1.25/lib/IPC/System/Simple.pm)
	
	To run with an annotator, additional packages will be needed:
	1. to run with ANNOVAR
		* ANNOVAR package
	2. to run with SnpEff
		* SnpEff package
		* ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
	3. to run with VEP
		* VEP package
		* ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
		* VCFTOOLs

=head1 AUTHORS

	Lei Wei (Lei.Wei@roswellpark.org)
	Lu Liu (Lu.Liu@RoswellPark.org)

=head1 LICENCE AND COPYRIGHT

Copyright (c) 2014 by Roswell Park Cancer Institue.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 2 of the License, or
(at your option) any later version. 

This program is distributed in the hope that it will be useful, 
but WITHOUT ANY WARRANTY; without even the implied warranty of 
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
