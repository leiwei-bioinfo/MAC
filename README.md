#########################################################
Multi-nucleotide Variation Annotation Corrector (MAC)	
													
CONTACT: Lei.Wei@roswellpark.org						
														
VERSION: 1.2											
VERSION DATE: 03/27/2014								
													
Copyright (c) 2014 by Roswell Park Cancer Institue.		
#########################################################

This document provides information on how to run MAC to screen through a list 
of SNVs called by any existing SNV-based variant caller to identify and fix 
incorrect amino acid prediction caused by multi-nucleotide variations (MNV).

------------
REQUIREMENTS
------------
MAC is implemented in Perl and tested with Perl v5.16.1, Bio::DB::Sam 1.36, 
Linux version 2.6.18. MAC can be run under any *nix OS. Before running MAC, 
please make sure the following packages are properly installed:

	1. Perl v5 (tested with Perl v5.16.1). Please ensure your running environment 
		points to the right Perl.
	2. Bio::DB::Sam (http://search.cpan.org/~lds/Bio-SamTools/), which also 
		requires SamTools library (http://sourceforge.net/projects/samtools/files/).
	3. IPC::System::Simple 
		(http://search.cpan.org/~pjf/IPC-System-Simple-1.25/lib/IPC/System/Simple.pm)

	To run with existing annotators, additional packages/files will be needed 
	depending on which annotator(s) will be used:

	4a. to run with Annovar
		Annovar:
			http://www.openbioinformatics.org/annovar/annovar_download_form.php

	4b. to run with SnpEff
		SnpEff:
			http://snpeff.sourceforge.net/download.html#download
		Ensembl gene annotation GTF file:   
			ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
		* Please note that Java is required when running SnpEff.

	4c. to run with VEP
		Vep:
			https://github.com/Ensembl/ensembl-tools/archive/release/75.zip
		VCFTOOLs
			http://vcftools.sourceforge.net/
		Ensembl gene annotation GTF file:
			ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz

	The users can install any one or all of the above 3 programs.

---------------------
PREPARING INPUT FILES
---------------------
MAC requires 3 input files:
	1.	a list of SNVs, each line contains one SNV in the format of 
		"chr.posi.ref.mut",	e.x. "17.7577084.T.A". Please make sure that the 
		nomenclature of chromosome name matches to the BAM and reference files,
		 such as with or without "chr". Currently only SNVs supported, any 
		Indels or MNVs need to be removed from the input file.
	2.	the BAM file from which these SNVs were detected. The BAM file must be 
		sorted and indexed.
	3.	the reference genome sequence file in fasta format, e.x. hg19.fasta. 
		Must be the same as the reference of above BAM file.

-------------------
RNNING THE PIPELINE
-------------------
MAC can be run either with or without a selected annotation program. For user's
 convenience, we precompiled our program to work with three popular annotators: 
ANNOVAR, SnpEff and VEP (Cingolani, et al., 2012; McLaren, et al., 2010; Wang, et al., 2010).
If the user prefers another program, simply not specify any annotator, the 
program will output identified Block of Mutations(BM), defined as a group of 
mutations where every mutation contains at least one read or mate pair 
(regardless of mutation status) that is shared with at least one other mutation
 in that group, in the format of haplotypes with read counts. The users can 
then annotate any contained MNVs using their own selection of annotator. Before
 running MAC, please make sure all required perl modules are contained in @INC.
 One simple way is to put all .pm and .pl scripts in the same directory and run
 from there. The additional program packages for Annovar/SnpEff/VEP can be 
installed in separate location, with the full path of required programs/files 
provided when running MAC.

Output file name prefix: MAC use input SNV file name as prefix by default. User 
can change it with "-o" option. 

MAC's 4 features:
	1. run without annotation (output without information about codons/amino 
acid prediction, user can then do annotation using their own selection of 
annotator):
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta

	#Output format when running MAC without annotation:
	Output file name: ${prefix}_raw_haplotype_counts.txt
	Example:
		BM_mutations	Haplotype_index	Joint_status	N_reads	Annotation
		17.7577084.T.A,17.7577085.C.G	1/3	11	17	Mutant
		17.7577084.T.A,17.7577085.C.G	2/3	00	18	Wt
		17.7577084.T.A,17.7577085.C.G	3/3	0-	1	Incomplete
	Each row contains one haplotype of one BM (Block of Mutations).
	Description of columns:
		BM_mutations: all SNVs contained in that BM. In the above example, 
			this BM contains 2 SNVs: 17.7577084.T.A and 17.7577085.C.G. Please 
			note these SNVs may or may not share any protein codon, such 
			information must be obtained by using one of the annotators.
		Haplotype_index: in example "1/3", contains the index number of 
			current haplotype (1) and the number of total haplotypes (3)
		Joint_status: the joint statuses of individual SNVs in current 
			haplotype. The total number of digits in "Joint_status" should 
			equal to the total number of SNVs in "BM_mutations", with each 
			digit corresponding to one SNV, in the same order. There are 3 
			possible statuses: 1, 0, and dash (-), to represent mutant, 
			non-mutant, and unknown, respectively. For example, a typical 
			dinucleotide SNV will be reported as 11, while two consecutive 
			SNVs on different reads will be reported as “10” or “01”.
		N_reads: total number of unique read pairs supporting the current 
			haplotype
		Annotation: when run without an annotator, MAC only reports Mutant, Mt,
			 or Incomplete (any haplotype containing unknown "-" status)


	#run with one of the three precompiled annotation programs	
	2. Annovar
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta -annotator 
		annovar -annovar_annotate_variation FULLPATH/annotate_variation.pl 
		-annovar_coding_change FULLPATH/coding_change.pl -annovar_refgene 
		FULLPATH/hg19_refGene.txt -annovar_refmrna FULLPATH/hg19_refGeneMrna.fa
	Please note the following programs/files are from ANNOVAR package:
		*	FULLPATH/annotate_variation.pl (MAC tested with "Revision: 527")
		*	FULLPATH/coding_change.pl (MAC tested with "Revision 
			545fca3432e51d34037d6055b34f10f4eb1b72fc")
		*	FULLPATH/hg19_refGene.txt
		*	FULLPATH/hg19_refGeneMrna.fa

	3. SnpEff
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta -annotator 
		snpeff -snpeff_path_to_snpEff FULLPATH/snpEff.jar 
		-ensembl_Homo_sapiens_GRCh37_75_gtf FULLPATH/Homo_sapiens.GRCh37.75.gtf
	Please note the following program is from SnpEff package:
		*	FULLPATH/snpEff.jar (MAC tested with "v3.6c")
	Please note the following file is from Ensembl website:
		*	FULLPATH/Homo_sapiens.GRCh37.75.gtf

	4. VEP
	./MAC_v1.0.pl -i input_SNVs.txt -bam sample.bam -r hg19.fasta -annotator 
		vep -vep_path_to_variant_effect_predictor 
		FULLPATH/variant_effect_predictor.pl -ensembl_Homo_sapiens_GRCh37_75_gtf
		 FULLPATH/Homo_sapiens.GRCh37.75.gtf
	Please note the following program is from VEP package:
	* FULLPATH/variant_effect_predictor.pl (MAC tested with "Version 75")
	Please note the following file is from Ensembl website:
	* FULLPATH/Homo_sapiens.GRCh37.75.gtf

	#Output format when running MAC with a selected annotator
	Output file name: ${prefix}_updated_haplotype_anno_by_${annotator}.txt, 
		where ${annotator} will be one of the 3 selected ones: annovar/snpeff/vep
	Example (using SnpEff as annotator):
		BMC_mutations	Haplotype_index	Joint_status	N_reads	Annotation
		17.7577084.T.A,17.7577085.C.G	1/2	11	17	TP53:ENST00000269305:Glu285Leu,TP53:ENST00000359597:Glu285Leu,TP53:ENST00000420246:Glu285Leu,TP53:ENST00000445888:Glu285Leu,TP53:ENST00000455263:Glu285Leu,TP53:ENST00000509690:Glu153Leu
		17.7577084.T.A,17.7577085.C.G	2/2	00	18	Wt
	The description of columns are the same as the above section of 
	"running MAC without annotation". The only difference is the column 
	"Annotation" now contains amino acid prediction, in the format of 
	"gene:mRNA transcript:amino acid change" (e.x. TP53:ENST00000269305:Glu285Leu).
	 When multiple transcripts are available, all transcripts will be reported 
	and deliminated by comma.
	*Please note that when running with an annotator, by default "Incomplete" 
	haplotypes ("0-" in the above example) are not reported. The user can use 
	the option of "--print_incomplete_haplotype" to include those.

Additional options
	-o:output prefix, default is to use input SNV file name as prefix.
	--print_incomplete_haplotype: select this option to keep haplotypes 
		containing unknown status in annotation
	--max_allowed_adjacent_distance: max allowed distance between two adjacent 
		SNVs in one Block of Mutation: for MNV annotation, use 2; for SNP phasing,
		use 1000 (and do not select any annotator). If set at 1, only continous SNVs
		are considered as BM.
	--printlog: print number counts at each step for debugging purpose
	--h: print brief help information
	--man: Print the man page
	--usage: Print usage information.

------------
EXAMPLE DATA
------------

An example data set is included for testing purpose under "./example/".

Input:
"input_SNVs.txt"	a testing input SNV file containing two consecutive SNVs
"sample.bam"		a mini BAM file containing NGS reads overlapping with the site
					 of	the input SNVs
"sample.bam.bai"	the index file for the above mini BAM file
* The reference genome for this test data set is not included but can be 
downloaded from:
	http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

Output:
When successfully finished, the expected output files are included under
"./example/output/":

	#when running MAC without annotation
	"input_SNVs.txt_raw_haplotype_counts.txt"
	
	#when running MAC with ANNOVAR
	"input_SNVs.txt_updated_haplotype_anno_by_annovar.txt"

	#when running MAC with SnpEff
	"input_SNVs.txt_updated_haplotype_anno_by_snpeff.txt"

	#when running MAC with VEP
	"input_SNVs.txt_updated_haplotype_anno_by_vep.txt"

----------
REFERENCES
----------
Cingolani, P., et al. (2012) A program for annotating and predicting the 
	effects of single nucleotide polymorphisms, SnpEff: SNPs in the genome of 
	Drosophila melanogaster strain w1118; iso-2; iso-3, Fly, 6, 80-92.
McLaren, W., et al. (2010) Deriving the consequences of genomic variants with 
	the Ensembl API and SNP Effect Predictor, Bioinformatics, 26, 2069-2070.
Wang, K., Li, M. and Hakonarson, H. (2010) ANNOVAR: functional annotation of 
	genetic variants from high-throughput sequencing data, Nucleic Acids Res, 38, e164.
