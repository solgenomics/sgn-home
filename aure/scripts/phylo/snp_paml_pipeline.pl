#!/usr/bin/perl

=head1 NAME

 snp_paml_pipeline.pl
 Pipeline to analyze Non-synonymous/Synonymous ratio from an assembly file in maf format (version.0.1).

=cut

=head1 SYPNOSIS

 snp_paml_pipeline.pl [-h] -i <maf_file> -c <control_file> [-s <step>] [-d <input_dir>] 
                           [-o <output_basename>][-X][-D]

=head1 EXAMPLE

 snp_paml_pipeline.pl [-h] -i coffea.maf -c snpsel.ctr

=head2 I<Flags:>

=over


=item -i

B<maf_file>                     assembly file in maf format format (mandatory)

=item -o

B<output_basename>              output basename (input_file.output by default)

=item -c 

B<control_file>                 control file with the program executables and the its options (mandatory)

=item -s

B<step>                         step to start the pipeline (1 by default)

=item -d 

B<input_dir>                    input dir with the files to start in an step different from 1

=item -X

B<create_control_file>          create the control file with the name maf_to_phylo.ctr

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

 This script is a pipeline that calculate omega (dN/dS ratio) for the proteins
 predicted in an unigene assembly.
 
 It use maf file produced by MIRA assembler, but can be used from different 
 steps.

 It have the following steps:

    0- Parse controler file
       Load input files (maf or step files)
       * Input: snp_selection.ctr
       * Output: none
 
    1- Extract the contigs and the reads from maf file
       * Input: .maf
       * Output: .fasta (per single contig)
                 .fasta (multifasta for all the reads, extracted from the contig)
                 .tab (read mapping data)

    2- Calculate predicted cds and proteins for contigs and reads, 
       using EstScan executable.
       * Input: .fasta (contig and multifasta for reads)
       * Output: .fasta (cds with X for unknowns)
                 .fasta (protein with X for unknows)

    3- Change X for N in cds and remove ';' from id predicted by estscan.
       Filter sequences based in protein size 
       or the in absent of the first met.
       * Input: .fasta (cds with X)
       * Output: .fasta (cds)

    4- Align proteins using Muscle executable
       * Input: .fasta (protein for contig + reads)
                .fasta (only reads)
       * Output: .clw (protein contig + reads)
                 .clw (protein reads).

    5- Calculate the coverage of the reads protein alignment using
       contig protein as reference.
       * Input: .clw (protein contig + reads).
       * Output: .txt (result)

    6- Create paml and codon alignment format from reads aligment
       * Input: .clw (protein reads)
       * Output: .paml (cds reads)
                 .codon (protein + cds reads)

    7- Generate PAML controler files, run PAML and get omega (dN/dS ratio)
       * Input: .paml (cds reads)
       * Output: .txt (ym file)

    8- Parse PAML ym file
       * Input: .txt (single ym file)
       * Output: .tab (results for all the assembly)
                 .txt (summary of analysis)

 This pipeline can start from any step supplying the input files.
 In that case the controler need to have the output_format parameter
 set for the previous step. 

 To create a control file: 
    snp_selection_pipeline.pl -X

 

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 contig_phyloanalysis.pl

=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use List::Util 'shuffle';
use Getopt::Std;
use Cwd;
use Path::Class;

use Bio::Seq;
use Bio::Seq::Quality;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Align::DNAStatistics;
use Math::BigFloat;

our ($opt_i, $opt_o, $opt_c, $opt_s, $opt_d, $opt_X, $opt_h);
getopts("i:o:c:s:d:Xh");
if (!$opt_i && !$opt_o && !$opt_c && !$opt_s && !$opt_d && !$opt_X && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}
elsif ($opt_X) {
    create_ctr();
}

## Get the arguments and check the mandatory ones

my $maf = $opt_i;
my $baseout = $opt_o 
    || $opt_i . '_output';
my $ctr = $opt_c
    || die("MANDATORY ARGUMENT ERROR: -c <control_file> was not supplied.\n");
my $step = $opt_s
    || 1;
my $input_dir = $opt_d;

print STDERR "\n\n========================";
print STDERR "\nSTARTING PIPELINE...\n";
print STDERR "========================\n\n";

## Check if exists outfile, outtree, infile and intree
## In can give problems to programs from Phylip

my @default_filename = ('rub', 'lnf', 'rst', 'rates', '2YN.dN', '2YN.dS', '2YN.t', 'yn00.ctl');

foreach my $def_file (@default_filename) {
    if (-e $def_file) {
	print STDERR "\nWARNING: $def_file exists into the working dir. It could be a problem with PAML programs.\n\n";
    }
}


## Define alway the input_file array
## This array always will be overwrite by the result of a step, but before do it
## the script will store the list into a hash as %files = ( step => { input or output => \@list_of_files ) 

my (%master_files, %count_files);

if ($step > 1) {

    print STDERR "START POINT MOVE TO STEP:$step\n";

    unless (defined $input_dir) {
	die("ARGUMENT ERROR: -s <step> option only can be used if it is specified a dir of input arguments.\n");
    }
    else {

	## It will take the files from the input dir.

	print STDERR "\n0) DATA LOADING\n";

	opendir(my $idir_fh, $input_dir);
	
	print STDERR "\nLoading files into the input file variable.\n";

	my @files = readdir($idir_fh);

	if ($step == 2) {
	    print STDERR "\nMESSAGE: Only files with the extensions .contig.fasta or .reads.fasta will be loaded.\n";
	    foreach my $file (@files) {
		unless ($file =~ m/^\.+$/) {
		    my $basename = basename($file);
		    if ($basename =~ m/(.+).contig\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'assembly_contig' => $file };
		    }
		    elsif ($basename =~ m/(.+).reads\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'assembly_reads' => $file };
		    }		
		}
	    }
	}
	elsif ($step == 3) {
	    print STDERR "\nMESSAGE: Only files with the extensions .cds/prot_contig.fasta or .cds/prot_reads.fasta will be loaded.\n";
	    foreach my $file (@files) {
		unless ($file =~ m/^\.+$/) {
		    my $basename = basename($file);
		    if ($basename =~ m/(.+).cds.contig\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'cds_contig' => $file };
		    }
		    elsif ($basename =~ m/(.+).cds.reads\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'cds_reads' => $file };
		    }
		    elsif ($basename =~ m/(.+).prot.+contig\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'protein_contig' => $file };
		    }
		    elsif ($basename =~ m/(.+).prot.+reads\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'protein_reads' => $file };
		    }
		}
	    }
	}
	elsif ($step =~ m/^[4|5|6]$/ ) {
	    print STDERR "\nMESSAGE: Only files with the extensions .cds/prot_contig.fasta OR.cds/prot_reads.fasta AND .aln/.clw will be loaded.\n";
	    foreach my $file (@files) {
		unless ($file =~ m/^\.+$/) {
		    my $basename = basename($file);
		    if ($basename =~ m/(.+).cds.contig\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'formated_cds_contig' => $file };
		    }
		    elsif ($basename =~ m/(.+).cds.reads\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'formated_cds_reads' => $file };
		    }
		    elsif ($basename =~ m/(.+).prot.+contig\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'formated_protein_contig' => $file };
		    }
		    elsif ($basename =~ m/(.+).prot.+reads\.[fasta|fas|seq]$/) {
			$master_files{$1} = { 'formated_protein_reads' => $file };
		    }
		    elsif ($basename =~ m/(.+)contig\.[aln|clw]$/) {
			$master_files{$1} = { 'contig_protein_alignment' => $file };
		    }
		    elsif ($basename =~ m/(.+)\.[aln|clw]$/) {
			$master_files{$1} = { 'protein_alignment' => $file };
		    }
		}
	    }
	}
	elsif ($step =~ m/^7$/ ) {
	    print STDERR "\nMESSAGE: Only files with the extensions .paml will be loaded.\n";
	    foreach my $file (@files) {
		unless ($file =~ m/^\.+$/) {
		    my $basename = basename($file);
		    if ($basename =~ m/(.+)\.paml$/) {
			$master_files{$1} = { 'cds_alignment_paml' => $file };
		    }
		}
	    }
	}
	else {
	    die("ARGUMENT ERROR: -s <step> option only can be used with values from 1 to 7.\n\n");
	}

	my $f_count = scalar(@files);

	if ($f_count == 0) {
	    die("\nINPUT ERROR: None file was found in the input dir: $input_dir\n");
	}
	print STDERR "\n$f_count files have been loaded to start the pipeline in the step $step.\n";
    }
}
else {
    ## Check if exists input file.

    unless (defined $maf) {
	die("ARGUMENT ERROR: -i <input_maf_file> was not supplied. It must to be supplied for step = 1.\n");
    }
}



## STEP 0
## First, parse the control file
    
print STDERR "\n0) PARSING CONTROL FILE.\n\n";

my @ctr_arguments_href = parse_ctr_file($ctr);

my $curr_path = getcwd();

print STDERR "\n\tDone...control file -c $ctr has been parsed.\n";

## Define the control variables 

my %assem_exec = %{$ctr_arguments_href[0]};
my %cds_exec = %{$ctr_arguments_href[1]};
my %prot_filter = %{$ctr_arguments_href[2]};
my %align_exec = %{$ctr_arguments_href[3]};
my %prot_cov_filter = %{$ctr_arguments_href[4]};
my %pal2nal_exec = %{$ctr_arguments_href[5]};
my %paml_exec = %{$ctr_arguments_href[6]};

my %strain_list = ();


## Check strains and extract from the maf file if extract tag is used

if (defined $assem_exec{'strain_data'}) {
    
    %strain_list = %{$assem_exec{'strain_data'}};

    if (defined $strain_list{'{extract'} ) {
	
	%strain_list = extract_strain_from_maf($maf, \%strain_list);
    }	
    $assem_exec{'strain_data'} = \%strain_list;
}
else {
    print STDERR "\n\tMSG: None strain constrain was defined. ";
    print STDERR "\n\tIt will take as many sequences as possible for the best overlaping region.\n";
}


############################################
## STEP 1 (ALIGNMENT EXTRACTION)          ##
##  Input type  => a single maf file      ##
##  Output type => a list of fasta files  ##
############################################

if ($step == 1) {

    ## Create the preprocessing_folder

    print STDERR "\n1) EXTRACTING ALIGMENTS FROM MAF FILE.\n";
    
    print STDERR "\n  1.1) CREATING PREPROCESSING FOLDER.\n";
    my $curr_path = getcwd();
    my $preprocessing_dir = $curr_path . '/dSdN_1_align_extract';

    mkdir($preprocessing_dir);
    print "\n\tDone ($preprocessing_dir).\n";


    ## Open the input maf file and parse

    print STDERR "\n  1.2) PARSING MAF FILE.\n";

    my ($fasta_files_href, $contig_href, $read_href) = extract_align_from_maf($maf, \%assem_exec, $preprocessing_dir);

    %count_files = file_class_counting($fasta_files_href);
    print STDERR "\n\n\tDone... $count_files{'assembly_reads'} reads and $count_files{'assembly_contigs'} contigs files have been created in the dir:$preprocessing_dir.\n\n";

    %master_files = %{$fasta_files_href};

    ## Finally add a new step to the step var

    $step++;
}


############################################################
## STEP 2 (CDS prediction)                                ##
##  Input type  => a list of fasta files                  ##
##  Output type => a list of fasta files                  ##
############################################################

if ($step == 2) {

    ## The program will be runned over the fasta file writting the output in the aligment folder

    print STDERR "\n2) CDS PREDICTION PROCESS.\n";

    if (defined $cds_exec{'cds_exec'}) {
    
	## It will realign the sequences using an specific program, 
	## so in that case it will create a bunch of new alignments.

	print STDERR "\n  2.1) CREATING CDS PREDICTION FOLDER.\n";
	my $cds_dir = $curr_path . '/dSdN_2_cds_prediction';
    
	mkdir($cds_dir);
	print "\n\tDone ($cds_dir).\n";    

	print STDERR "\n  2.2) RUNNING CDS PREDICTION SEQUENCES.\n";

	my %out_files = run_cds_prediction(\%cds_exec, \%master_files, $cds_dir);

	%count_files = file_class_counting(\%out_files);
	my $total = $count_files{'cds_contig'} + $count_files{'cds_reads'} + $count_files{'protein_contig'} + $count_files{'protein_reads'};
	print STDERR "\n\n\tDone... $total files \n";
	print STDERR "\t\t(cds_contigs=$count_files{'cds_contig'} | cds_reads=$count_files{'cds_reads'} | protein_contigs=$count_files{'protein_contig'} | ";
	print STDERR "protein_reads=$count_files{'protein_reads'})\n\t\thave been created in the dir:$cds_dir.\n";
		
	%master_files = %out_files;	
    }
    else {
	print STDERR "\n  2.1) CREATING CDS PREDICTION FOLDER.\n";
	print STDERR "\n\tNone cds prediction process have been asigned. Skipping step.\n";
	
	print STDERR "\n  2.2) RUNNING CDS PREDICTION SEQUENCES.\n";
	print STDERR "\n\tNone cds prediction process have been asigned. Skipping step\n\n";
    }

    ## Finally add a new step to the step var
    $step++;
}

############################################################
## STEP 3 (fasta seq reformat)                            ##
##  Input type  => a list of fasta files                  ##
##  Output type => a list of fasta files                  ##
############################################################

if ($step == 3) {

    ## The program will be runned over the fasta file writting the output in the aligment folder

    print STDERR "\n3) FASTA SEQUENCE REFORMAT AND FILTERING.\n";
    
    ## It will realign the sequences using an specific program, 
    ## so in that case it will create a bunch of new alignments.

    print STDERR "\n  3.1) CREATING FASTA REFORMAT FOLDER.\n";
    my $ref_dir = $curr_path . '/dSdN_3_fasta_reformat';
    
    mkdir($ref_dir);
    print "\n\tDone ($ref_dir).\n";    
    
    print STDERR "\n  3.2) RUNNING REFORMAT OF THE SEQUENCES AND FILTERING.\n";

    my %out_files = run_seqformating(\%prot_filter, \%master_files, $ref_dir);
    
    %count_files = file_class_counting(\%out_files);
    my $total2 = $count_files{'formated_cds_contig'} + $count_files{'formated_cds_reads'} + $count_files{'formated_protein_contig'} + $count_files{'formated_protein_reads'};
    print STDERR "\n\tDone... $total2\n\t\t(formated_cds_contigs=$count_files{'formated_cds_contig'} | formated_cds_reads=$count_files{'formated_cds_reads'} | ";
    print STDERR "formated_protein_contig=$count_files{'formated_protein_contig'} | formated_protein_reads=$count_files{'formated_protein_reads'})\n\t";
    print STDERR "have been created in the dir:$ref_dir.\n";
    
    %master_files = %out_files;
	
    ## Finally add a new step to the step var
    $step++;
}


############################################################
## STEP 4 (protein aligment)                              ##
##  Input type  => a list of fasta files                  ##
##  Output type => a list of clw files                    ##
############################################################

if ($step == 4) {

    ## The program will be runned over the fasta file writting the output in the aligment folder

    print STDERR "\n4) PROTEIN ALIGNMENT.\n";
    
    ## It will realign the sequences using an specific program, 
    ## so in that case it will create a bunch of new alignments.

    print STDERR "\n  4.1) CREATING PROTEIN ALIGNMENT FOLDER.\n";
    my $align_dir = $curr_path . '/dSdN_4_protein_alignment';
    
    mkdir($align_dir);
    print "\n\tDone ($align_dir).\n";    
    
    print STDERR "\n  4.2) RUNNING ALIGNMENT FOR PROTEIN SEQUENCES.\n";

    my %out_files = run_alignment(\%align_exec, \%master_files, $align_dir);
    
    %count_files = file_class_counting(\%out_files);
    print STDERR "\n\tDone... $count_files{'protein_alignment'} reads alignment and $count_files{'contig_protein_alignment'} global alignment files have been created";
    print STDERR " in the dir:$align_dir.\n";
    
    %master_files = %out_files;
	
    ## Finally add a new step to the step var
    $step++;
}


############################################################
## STEP 5 (protein coverage)                              ##
##  Input type  => a list of clw files                    ##
##  Output type => a list of clw files and txt list       ##
############################################################

if ($step == 5) {

    ## The program will be runned over the fasta file writting the output in the aligment folder

    print STDERR "\n5) PROTEIN COVERAGE CALCULATION.\n";
    
    ## It will calculate the coverage with a reference sequence, usually
    ## the contig

    print STDERR "\n  5.1) CREATING PROTEIN COVERAGE FOLDER.\n";
    my $cov_dir = $curr_path . '/dSdN_5_protein_coverage';
    
    mkdir($cov_dir);
    print "\n\tDone ($cov_dir).\n";    
    
    print STDERR "\n  5.2) RUNNING COVERAGE FOR PROTEIN ALIGNMENTS.\n";

    my %out_files = run_coverage(\%prot_cov_filter, \%master_files, $cov_dir);
    
    %count_files = file_class_counting(\%out_files);
    my $prot_coverage_files =  $count_files{'protein_coverage'} || 0;
    print STDERR "\n\tDone... $prot_coverage_files files have been created";
    print STDERR " in the dir:$cov_dir.\n";
    
    %master_files = %out_files;
	
    ## Finally add a new step to the step var
    $step++;
}

############################################################
## STEP 6 (converting to paml format)                     ##
##  Input type  => a list of clw files                    ##
##  Output type => a list of paml files                   ##
############################################################

if ($step == 6) {

    ## The program will be runned over the fasta file writting the output in the aligment folder

    print STDERR "\n6) PAML FORMAT.\n";
    
    ## It will take the alignment file and converting to paml format using the script pal2nal.pl

    print STDERR "\n  6.1) CREATING PAML FORMAT FOLDER.\n";
    my $pal2nal_dir = $curr_path . '/dSdN_6_paml_converting';
    
    mkdir($pal2nal_dir);
    print "\n\tDone ($pal2nal_dir).\n";    
    
    print STDERR "\n  6.2) RUNNING PAML FORMATING SCRIPT.\n";

    my %out_files = run_pal2nal(\%pal2nal_exec, \%master_files, $pal2nal_dir);
    
    %count_files = file_class_counting(\%out_files);
    my $paml_input_files =  $count_files{'cds_alignment_paml'} || 0;
    print STDERR "\n\tDone... $paml_input_files files have been created";
    print STDERR " in the dir:$pal2nal_dir.\n";
    
    %master_files = %out_files;
	
    ## Finally add a new step to the step var
    $step++;
}

############################################################
## STEP 7 (running paml program)                          ##
##  Input type  => a list of paml files                   ##
##  Output type => a list of paml result files            ##
############################################################

if ($step == 7) {

    ## The program will be runned over the fasta file writting the output in the aligment folder

    print STDERR "\n7) PAML PROGRAM.\n";
    
    ## It will take paml program and it will run paml a paml program

    print STDERR "\n  7.1) CREATING PAML PROGRAM FOLDER.\n";
    my $paml_dir = $curr_path . '/dSdN_7_paml_running';
    
    mkdir($paml_dir);
    print "\n\tDone ($paml_dir).\n";    
    
    print STDERR "\n  7.2) RUNNING PAML PROGRAM.\n";

    my %out_files = run_paml(\%paml_exec, \%master_files, $paml_dir);
    
    %count_files = file_class_counting(\%out_files);
    my $paml_output_files =  $count_files{'paml_result'} || 0;
    print STDERR "\n\tDone... $paml_output_files files have been created";
    print STDERR " in the dir:$paml_dir.\n";
    
    %master_files = %out_files;
	
    ## Finally add a new step to the step var
    $step++;
}

############################################################
## STEP 8 (parsing paml result)                           ##
##  Input type  => a list of paml result                  ##
##  Output type => a tab file with paml result            ##
############################################################

if ($step == 8) {

    print STDERR "\n8) PARSING PAML RESULTS.\n";
    
    ## It will parse paml result to get the omega values

    print STDERR "\n  8.1) CREATING PAML PARSING FOLDER.\n";
    my $paml_r_dir = $curr_path . '/dSdN_8_paml_parsing';
    
    mkdir($paml_r_dir);
    print "\n\tDone ($paml_r_dir).\n";    
    
    print STDERR "\n  8.2) PARSING PAML RESULT.\n";

    my ($result, $distribution) = run_paml_parsing(\%master_files, $paml_r_dir);
    
    print STDERR "\n\tDone... $result AND $distribution files have been created\n";

	
    ## Finally add a new step to the step var
    $step++;
}



print STDERR "\n\n";    

=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script is a pipeline that calculate omega (dN/dS ratio) for the proteins
     predicted in an unigene assembly.
 
     It use maf file produced by MIRA assembler, but can be used from different 
     steps.

     It have the following steps:

     0- Parse controler file
         Load input files (maf or step files)
         * Input: snp_selection.ctr
         * Output: none
 
     1- Extract the contigs and the reads from maf file
         * Input: .maf
         * Output: .fasta (per single contig)
                   .fasta (multifasta for all the reads, extracted from the contig)
                   .tab (read mapping data)

     2- Calculate predicted cds and proteins for contigs and reads, 
        using EstScan executable.
         * Input: .fasta (contig and multifasta for reads)
         * Output: .fasta (cds with X for unknowns)
                   .fasta (protein with X for unknows)

     3- Change X for N in cds and remove ';' from id predicted by estscan.
        Filter sequences based in protein size or the in absent of the first met.
         * Input: .fasta (cds with X)
         * Output: .fasta (cds)

     4- Align proteins using Muscle executable. If a sequence reference file an
        extra alignment will be done with these references, if not an alignment
	with the contig sequence will be used as reference
         * Input: .fasta (protein for contig + reads)
                  .fasta (only reads)
         * Output: .clw (protein contig + reads)
                   .clw (protein reads).

     5- Calculate the coverage of the reads protein alignment using
        contig protein as reference.
         * Input: .clw (protein contig + reads).
         * Output: .txt (result)

     6- Create paml and codon alignment format from reads aligment
         * Input: .clw (protein reads)
         * Output: .paml (cds reads)
                   .codon (protein + cds reads)

     7- Generate PAML controler files, run PAML and get omega (dN/dS ratio)
         * Input: .paml (cds reads)
         * Output: .txt (ym file)

     8- Parse PAML ym file
         * Input: .txt (single ym file)
         * Output: .tab (results for all the assembly)
                   .txt (summary of analysis)

     This pipeline can start from any step supplying the input files.
     In that case the controler need to have the output_format parameter
     set for the previous step. 

     To create a control file: 
     snp_selection_pipeline.pl -X
      
         
    Usage:

      snp_paml_pipeline.pl [-h] -i <maf_file> -c <control_file> 
                                [-s <step>] [-d <input_dir>] 
                                [-o <output_basename>][-X]
      
    Example:

      snp_paml_pipeline.pl [-h] -i coffea.maf -c snpsel.ctr
       
    Flags:

      -i <maf_file>               assembly file in maf format format (mandatory)
      -o <output_basename>        output basename (input_file.output by default)
      -c <control_file>           control file (mandatory)
      -s <step>                   step to start the pipeline (1 by default)
      -d <input_dir>              input dir to start in an step > 1
      -X <create_control_file>    create the control file: maf_to_phylo.ctr
      -h <help>                   print the help
      

EOF
exit (1);
}


=head2 create_ctr

  Usage: create_ctr
  Desc: create_ctr file
  Ret: none
  Args: filename
  Side_Effects: exit of the script
  Example: if ($opt_X) {
               create_ctr();
           }

=cut

sub create_ctr {    
    my $ctr_filename = shift || 'dsdn_snp.ctr';
 
    open my $ctr_fh, '>', $ctr_filename 
	|| die("OPEN FILE ERROR: $ctr_filename file can not be openned (system error: $!),\n");

    
    print $ctr_fh <<EOF;


############################################################
## This pipeline is composed by eight steps:
##
## 0- Parse this controler file
##       Load input files (maf or step files)
##       * Input: snp_selection.ctr
##       * Output: none
## 
##    1- Extract the contigs and the reads from maf file
##       * Input: .maf
##       * Output: .fasta (per single contig)
##                 .fasta (multifasta for all the reads, 
##                  extracted from the contig)
##                 .tab (read mapping data)
##
##    2- Calculate predicted cds and proteins for contigs 
##       and reads, using EstScan executable.
##       * Input: .fasta (contig and multifasta for reads)
##       * Output: .fasta (cds with X for unknowns)
##                 .fasta (protein with X for unknows)
##
##    3- Change X for N in cds and remove ';' from id 
##       predicted by estscan. Size sequence filtering or 
##       based the in absent of the first met.
##       * Input: .fasta (cds with X)
##       * Output: .fasta (cds)
##
##    4- Align proteins using Muscle executable. If a sequence 
##       reference file an extra alignment will be done with 
##       these references, if not an alignment with the contig 
##       sequence will be used as reference
##       * Input: .fasta (protein for contig + reads)
##                .fasta (only reads)
##       * Output: .clw (protein contig + reads)
##                 .clw (protein reads).
##
##    5- Calculate the coverage of the reads protein alignment 
##       using contig protein as reference.
##       * Input: .clw (protein contig + reads).
##       * Output: .txt (result)
##
##    6- Create paml and codon alignment format from reads 
##       aligment
##       * Input: .clw (protein reads)
##       * Output: .paml (cds reads)
##                 .codon (protein + cds reads)
##
##    7- Generate PAML controler files, run PAML and get omega 
##       (dN/dS ratio)
##       * Input: .paml (cds reads)
##       * Output: .txt (ym file)
##
##    8- Parse PAML ym file
##       * Input: .txt (single ym file)
##       * Output: .tab (results for all the assembly)
##                 .txt (summary of analysis)
## 
############################################################

############################################################
### 1) ASSEMBLY FILE DATA EXTRACTION                     ###
############################################################

## Note: The script parse maf (mira assembly format) files
##       To use other formats like caf or ace, you can 
##       change the format using pharp2caf script
##       (http://www.sanger.ac.uk/resources/software/caf/)
##       or convert_project script from mira 
##       (http://sourceforge.net/projects/mira-assembler/)

<STRAIN_LIST>="list of strains sepparated by comma"
## Requeriment: one per strain in the file by default 
## Example:     <ALIGMENT_EXEC>="str1:1,str2:1,str3:2"
## Note:        List of strains with the sequence number
##              as cutoff to consider an alignment
##              If none strains are detailed it will take
##              all the sequences.
##              If '{extract:n}' is used, it will parse the
##              different strains from the maf file and use
##              as n sequences per strain.

<STRAIN_SELECTION>="method to select the different strain"
## Requeriment: none by default
## Example:     <STRAIN_SELECTION>="Distance:str1"
## Note:        method used to select the sequences in the
##              matching region of the alignment.
##              None   => parsing order
##              Random => use random function so every time
##                        that the script is run, it will 
##                        take different set of sequences 
##                        with the specified requeriments
##              Distance => It will calculate the distance 
##                          between sequences and it will 
##                          rename them according kimura 
##                          distance to the seed

<MIN_ALIGNMENT_LENGTH>="minimum alignment length"
## Requeriment: none, 100 by default
## Example:     <MIN_ALIGNMENT_LENGTH>="50"
## Note:        All the aligments with length 
##              < this value will be discarted 


############################################################
### 2) PREDICT CDS                                       ###
############################################################

## Note: It use EstScan to predict cds for nucleotide 
##       sequences
    
<CDS_EXEC>="Predict cds tool executable path" 
## Requeriment: None, estscan by default 
## Example:     <CDS_EXEC>="~/programs/estscan"
## Note:        Path to access to executable

<CDS_MATRIX>="Output format for the alignment" 
## Requeriment: none by default 
## Example:     <CDS_MATRIX>="tomato.dat"
## Note:        Score matrix file


############################################################
### 3) CHANGE X for N and REMOVE ';' FROM IDs            ###
############################################################

<PROT_SIZE_FILTER>="Protein size filtering"
## Requeriment: None
## Example: <PROT_SIZE_FILTER>="33"
## Note: It will remove all the sequences with size less 
##       than this.

<FIRST_MET_FILTER>="Yes/No | Y/N | 1/0"
## Requeriment: None, No by default
## Example: <FIRST_MET_FILTER>="Y"
## Note: It will remove all the sequences that have not a 
##       first met


############################################################
### 4) ALIGN PROTEINS USING MUSCLE OR CLUSTAW            ###
############################################################

<ALIGN_EXEC>="Align tool executable"
## Requeriment: muscle by defalut
## Example:     <ALIGN_EXEC>="muscle"
## Note:        Alignment tool. The output should be clustaw
##              strict format. Values allowed muscle or clustalw
##
## muscle always will run with the following arguments:
## -clwstrict
## -stable
## clustalw always will run with the following arguments:
## -OUTORDER=INPUT
##

<ADD_ALIGN_ARGS>="list of arguments to add to alignment"
## Example: <ADD_ALIGN_ARGS>="-diags" 
## Note: To add more options to the arguments supplied to
##       alignment tool. It still will use clustalw output
##       and stable/outorder=input to be able to be used by
##       pal2nal.pl script

<BLAST_EXEC>="blast executable" 
## Example: <BLAST_EXEC>="blastall" 
## Note: Blast executable to do the comparisons with the 
##       reference

<REFERENCE_PROTEIN_FILE>="A file with a set of ref. proteins" 
## Example: <REFENCE_PROTEIN_FILE>="swisprot" 
## Note: It will use a referebce protein dataset. It will do
##       one tblastx per contig sequence with a single match
##       and it will take the protein to do an extra alignment
##       with the reads 

############################################################
### 5) PROTEIN COVERAGE CALCULATION                      ###
############################################################

<PROT_COVERAGE_FILTER>="Protein coverage filtering"
## Requeriment: None
## Example: <PROT_COVERAGE_FILTER>="75"
## Note: It will remove all the alignments where the 
## sequence coverage will be less than this value

<PROT_IDENTITY_FILTER>="Protein identity filtering"
## Requeriment: None
## Example: <PROT_IDENTITY_FILTER>="75"
## Note: It will remove all the alignments where the 
## sequence indentity of the global alignment will be less 
## than this value


############################################################
### 6) CREATE PAML AND CODON PROTEIN ALIGNMENT           ###
############################################################

      
<PAL2NAL_EXEC>="executable for pal2nal.pl script"
## Requeriment: pal2nal.pl by default
## Example:     <PAL2NAL_EXEC>="pal2nal.pl"
## Note:        It always will be run two times, one with
##              the argument -output=paml and other with
##              -output=codon

<ADD_PAL2NAL_ARGS>="list of arguments to add to pal2nal"
## Example: <ADD_PAL2NAL_ARGS>="-nogap" 
## Note: To add more options to the arguments supplied to
##       pal2nal.pl script.


############################################################
### 7) CREATE PAML CRL FILES AND RUN PAML                ###
############################################################

      
<PAML_EXEC>="PAML program path"
## Requeriment: yn00 by default
## Example:     <PAML_EXEC>="yn00"
## Note:        

<ADD_PAML_ARGS>="Add arguments to ctr file used by PAML" 
## Requeriment: none by default 
## Example:     <CONSTREE_IFORMAT>="icode = 1"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. 


############################################################
### 8) PARSE PAML OUTPUT                                 ###
############################################################

## None variable is requested at this step  

EOF

print STDERR "\nControl file have been created ($ctr_filename)\n\n";

exit (1);

}

=head2 parse_ctr_file

  Usage: my @ctr_arguments_href = parse_ctr_file($ctr_file);

  Desc: parse control file

  Ret: An array of hash references:
       $ctr_args[0] = $assembly_opts_href
       $ctr_args[1] = $predict_cds_args_href
       $ctr_args[2] = $protein_length_filter_href
       $ctr_args[3] = $align_args_href
       $ctr_args[4] = $protein_cov_filter_href
       $ctr_args[5] = $pal2nal_args_href
       $ctr_args[6] = $paml_arg_href

  Args: $ctr_file, a scalar with the filename

  Side_Effects: Add default values
                Die if a mandatory argument is not specified

  Example: my @ctr_arguments_href = parse_ctr_file($ctr_file);

=cut

sub parse_ctr_file {    
    
    my $ctrfile = shift 
	|| die("ARGUMENT ERROR: None control file name was supplied to parse_ctr_file function.\n");

    ## Define the control variables with the default values

    my %assembly_args =    ( 'min_alignment' => 100          );
    my %cds_predict_args = ( 'cds_exec'      => 'estscan'    );
    my %prot_size_filter = ();
    my %align_args =       ( 'align_exec'    => 'muscle'     );
    my %prot_cov_filter =  ();
    my %pal2nal_args =     ( 'pal2nal_exec'  => 'pal2nal.pl' );
    my %paml_args =        ( 'paml_exec'     => 'yn00'       );

    ## Open the file and parse it (use autodie)

    open( my $ctr_fh, '<', $ctrfile);

    while(<$ctr_fh>) {
	chomp($_);

	## Ignore the lines that start with #

	unless ($_ =~ m/^\s*#/) {

	    ## ASSEMBLY ARGUMENTS

	    if ($_ =~ m/<STRAIN_LIST>="(.+)"/) {
		my @strains = split(/,/, $1);
		foreach my $strain (@strains) {
		    my @str_data = split(/:/, $strain);

		    $assembly_args{'strain_data'}->{$str_data[0]} = $str_data[1];
		}
	    }
	    if ($_ =~ m/<STRAIN_SELECTION>="(.+)"/) {
		my $strain_select = $1;

		## check that it have the right format for distance

		if ($strain_select =~ m/distance:(.+)/i) {
		    my @par = split(',', $1);
		    foreach my $par (@par) {
			unless ($par =~ m/.+?>.+?\*\d+/) {
			    die("PARAMETER ERROR: The parents detailed for distance selection need to have format StrainParent>StrainChildren*Digit\n");
			}
		    }
		}

		$assembly_args{'strain_selection'} = $strain_select;
	    }
	    
	    ## Overwrite the default value with the new one

	    if ($_ =~ m/<MIN_ALIGNMENT_LENGTH>="(.+)"/) {
		$assembly_args{'min_alignment'} = $1;
	    }
	    
	    ## CDS PREDICTION ARGUMENTS
	     	    
	    if ($_ =~ m/<CDS_EXEC>="(.+)"/) {
		$cds_predict_args{'cds_exec'} = $1; 
	    }
	    if ($_ =~ m/<CDS_MATRIX>="(.+)"/) {
		$cds_predict_args{'cds_matrix'} = $1; 
	    }

	    ## PROTEIN FORMAT AND FILTER
	    
	    if ($_ =~ m/<PROT_SIZE_FILTER>="(.+)"/) {
		$prot_size_filter{'prot_size_filter'} = $1; 
	    }
	    if ($_ =~ m/<FIRST_MET_FILTER>="(.+)"/) {
		$prot_size_filter{'first_met_filter'} = $1; 
	    }

	    ## ALIGNMENT ARGUMENTS
	    	    
	    if ($_ =~ m/<ALIGN_EXEC>="(.+)"/) {
		$align_args{'align_exec'} = $1; 
	    }
	    if ($_ =~ m/<ADD_ALIGN_ARGS>="(.+)"/) {
		$align_args{'more_align_args'} = $1; 
	    }
	     if ($_ =~ m/<BLAST_EXEC>="(.+)"/) {
		$align_args{'blast_exec'} = $1; 
	    }
	    if ($_ =~ m/<REFERENCE_PROTEIN_FILE>="(.+)"/) {
		$align_args{'ref_protein_file'} = $1; 
	    }

	    ## PROTEIN COVERAGE ARGUMENTS

	    if ($_ =~ m/<PROT_COVERAGE_FILTER>="(.+)"/) {
		$prot_cov_filter{'prot_cov_filter'} = $1; 
	    }
	    if ($_ =~ m/<PROT_IDENTITY_FILTER>="(.+)"/) {
		$prot_cov_filter{'prot_ident_filter'} = $1; 
	    }

	    ## PAL2NAL ARGUMENTS
	    
	    if ($_ =~ m/<PAL2NAL_EXEC>="(.+)"/) {
		$pal2nal_args{'pal2nal_exec'} = $1; 
	    }
	    if ($_ =~ m/<ADD_PAL2NAL_ARGS>="(.+)"/) {
		$pal2nal_args{'more_pal2nal_args'} = $1; 
	    }
	    
	    ## PAML ARGUMENTS

	    if ($_ =~ m/<PAML_EXEC>="(.+)"/) {
		$paml_args{'paml_exec'} = $1; 
	    }
	    if ($_ =~ m/<ADD_PAML_ARGS>="(.+)"/) {
		$paml_args{'more_paml_args'} = $1; 
	    }
	}
    }

    return (\%assembly_args, \%cds_predict_args, \%prot_size_filter, \%align_args, \%prot_cov_filter, \%pal2nal_args, \%paml_args);
}


=head2 extract_strain_from_maf

  Usage: my %strain_list = extract_strain_from_maf($maf, \%strain_list);

  Desc: parse the maf (mira alignment file) getting only the strains

  Ret: A hash with keys=strain and value=extract value from argument

  Args: $maf_file, maf filename
        \%strain_list, a hash reference with the list of strains (extract tag)

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %strain_list = extract_strain_from_maf($maf, \%strain_list);

=cut

sub extract_strain_from_maf {
    my $maf = shift ||
	die("ARGUMENT ERROR: maf filename was not supplied to the extract_strain_from_maf function.\n");
    my $strain_href = shift ||
	die("ARGUMENT ERROR: strain argument hash ref. was not supplied to the extract_strain_from_maf function.\n");

    my %strain_list = %{$strain_href};

    my $str_n = $strain_list{'{extract'};

    ## remove the { symbol from the parsing
    $str_n =~ s/\}//;
    
    unless ($str_n =~ m/\d+/) {
	die("ARGUMENT ERROR: Only integers can be used with {extract:$str_n} tag.\n");	    
    }

    print STDERR "\n\tMSG: Extract tag was supplied into the control file.\n\tIt will take strains from the maf file...\n";
    open (my $maf0_fh, '<', $maf);
    
    while(<$maf0_fh>) {
	chomp($_);
	if ($_ =~ m/^SN\s+(.+?)$/) { 
	    unless($strain_list{$1}) {
		$strain_list{$1} = $str_n;
	    }
	}
    }
	
    ## Delete the old tag
    delete $strain_list{'{extract'};

    $str_n = scalar(keys %strain_list);
    if ($str_n == 0) {
	die("ERROR: There are not strain data in the maf file.\n\n");
    }
    else {
	my $strain_list_p = join(',', keys %strain_list);
	print STDERR "\t$str_n strains have been obtained from the maf file ($strain_list_p).\n";
    }
    
    return %strain_list;
}

=head2 extract_align_from_maf

  Usage: my ($align_files_aref, $contig_href, $read_href) = extract_align_from_maf($maf_file, \%align_args, $preprocessing_dir);

  Desc: parse the maf (mira alignment file) getting only the tags that 
        needs and create the alignment files for the strains detailed into
        the comand file.

  Ret: $output_files_href, a hash reference with keys => contig_id
                                                 value => hash reference with:
                                                          keys => type (assembly_contig/assembly_reads)
                                                          value => filename
       $contig_href, a hash ref with: keys   => contig_id, 
                                      values => hash ref with:
                                                keys  => type
                                                value => data
       $read_href, a hash ref with:   keys   => read_id, 
                                      values => hash ref with:
                                                keys  => type
                                                value => data

  Args: $maf_file, maf filename
        \%align_args, a hash reference with the extract alignment args 

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my @align_hrefs = extract_align_from_maf($maf_file, \%align_args, $preprocessing_dir);

=cut

sub extract_align_from_maf {
    my $maf = shift ||
	die("ARGUMENT ERROR: maf filename was not supplied to the extract_align_from_maf function.\n");
    my $align_href = shift ||
	die("ARGUMENT ERROR: aligment argument hash ref. was not supplied to the extract_align_from_maf function.\n");
    my $preprocessing_dir = shift ||
	die("ARGUMENT ERROR: aligment preprocessing dir was not supplied to the extract_align_from_maf function.\n");
    

    ## Get the data from the alignment hash ref

    my %strain_list;

    if (defined $assem_exec{'strain_data'}) {
	%strain_list = %{$assem_exec{'strain_data'}};
    }

    my $min_align_l = $assem_exec{'min_alignment'};



    open( my $maf_fh, '<', $maf);                 ## (use autodie)

    ## Define the catch variables

    ## It only will need catch the ids and the sequences

    my ($contig_id, $contig_seq_pad, $contig_seq_unpad, $read_id, $read_seq);
    my ($l, $n_co, $n_rd) = (0, 0, 0);

    ## It will create one file per contig with the name of the 
    ## contig and the sequences of it. It will put into the
    ## fasta_files array

    my (%contig, %read, %start, %end);

    my %output_files;
    my $tag;
    my $seqio;

    while(<$maf_fh>) {
	chomp($_);
	my $line = $_;
	$l++;   

	## Maf file works with tags and tabs. To catch a data type only needs to catch the tag
	## and the data after the tag

	## All the reads in the contig finish with //
	## Now it will calculate the alignments

	if ($line =~ m/^\/\//) {

	    my $contig_length = length($contig_seq_pad);

	    ## Scan the contig positions from contig_start = 1 to contig_end = X
	    ## Use the start and the end with the position

	    my @reads = @{$contig{$contig_id}->{'read'}};

	    my (%starts, %ends);

	    my ($n, $k) = (0, 0);
	    my %p_count;

	    ## Create the hash to count nt for strains

	    my %strain_c;
	    foreach my $str (keys %strain_list) {
		$strain_c{$str} = 0;
	    }

	    my $acc_count;
	    my %curr_rds;
	    my $rds_analyzed = '';
	    my $e = 0;
	    my %elegible = ();
	

	    while ($n < $contig_length) {
		$n++;

		my @st_reads = (); 
		my @en_reads = ();
	    
		## Count of totals and strains

		if (exists $start{$n}) {
		    @st_reads = @{$start{$n}};
		    foreach my $rd (@st_reads) {
			my $st = $read{$rd}->{'strain'};
			$strain_c{'total'}++;
			$strain_c{$st}++;
			$curr_rds{$rd} = 1;
		    }

		    ## Check if it elegible, if it is store.
		    ## If strain_data is not defined all the sequences will be elegible from 
		    ## the strain constraint

		    my $is_elegible = 1;

		    if (defined $assem_exec{'strain_data'}) { 
			foreach my $check_str (keys %strain_list) {
			    if (defined $strain_c{$check_str}) {
				unless ($strain_c{$check_str} >= $strain_list{$check_str}) {
				    $is_elegible = 0;
				}
			    }
			}
		    }

		    ## While the curr_reads are the same it should not change.
		    ## When change the curr_reads should open or close elegible zone

		    my @curr_reads = sort keys %curr_rds;
		    my $curr_rds = join(',', @curr_reads);
				
		    if ($is_elegible == 1 && $rds_analyzed ne $curr_rds) {
		    
			## It will be elegible for this sequences enabled.	    
			$e++;
			$elegible{$e}->{'start'} = $n;
			$elegible{$e}->{'reads'} = \@curr_reads;
			$rds_analyzed = $curr_rds;
		    }			
		}
       
		if (exists $end{$n}) {
		    @en_reads = @{$end{$n}};		
		    foreach my $rd (@en_reads) {
			my $st = $read{$rd}->{'strain'};
			$strain_c{'total'}--;
			$strain_c{$st}--;
			delete $curr_rds{$rd};
			
			## It will add the end to the elegible zone when one of the current reads will be delete 
   
			foreach my $e (keys %elegible) {
			    my @reads = @{$elegible{$e}->{'reads'}};
			    foreach my $reads (@reads) {
				if ($reads eq $rd && !$elegible{$e}->{'end'}) {  ## and it still is opened
				    
				    ## It should close the elegible
				    $elegible{$e}->{'end'} = $n;
				    $elegible{$e}->{'length'} = $elegible{$e}->{'end'} - $elegible{$e}->{'start'};
				}
			    }
			}
		    }
		}
	    
		## Cluster the counts by groups

		foreach my $strain (keys %strain_c) {
		    my $count = $strain_c{$strain};

		    if (exists $p_count{$count}) {
			if (exists $p_count{$count}->{$strain}) {
			    push @{$p_count{$count}->{$strain}}, $n;
			}
			else {
			    $p_count{$count}->{$strain} = [$n];
			}
		    }
		    else {
			$p_count{$count}->{$strain} = [$n];
		    }
		}	    
	    }
    
	    ## Now it should have a list of elegible zones. 
	    ## It will take the bigger... with at least two sequences

	    my %selected_ovlp = ();
	    my $overlap = 0;
	    foreach my $ovlps (keys %elegible) {
		my $ovlp_length = $elegible{$ovlps}->{'length'};
		my $ovlp_n = scalar(@{$elegible{$ovlps}->{'reads'}});

		if ($overlap < $ovlp_length && $ovlp_n > 1) {
		    $overlap = $ovlp_length;
		    %selected_ovlp = %{$elegible{$ovlps}};
		}
	    }

	    ## Once the bigger overlap is chossen it will get each trimmed sequence to print into a file
	    
	    ## it can happens that do not exists sequences for a concrete contig

	    ## Before create the file also it will check that the aligment have the right length
	
		
	    if (defined $selected_ovlp{'reads'}) {

		my $align_l = $selected_ovlp{'end'} - $selected_ovlp{'start'};

		if ($align_l >= $min_align_l) {	
		
		    my $filename = $preprocessing_dir . '/' . $contig_id . '.reads.fasta';	

		    $seqio = Bio::SeqIO->new( -file => ">$filename", -format => "fasta" );

		    my $filename2 = $preprocessing_dir . '/' . $contig_id . '.overlaps.tab';
		    open(my $out1_fh, '>', $filename2);

		    print $out1_fh "#co_id\t#co_length\t#ovlp_length\t#ovlp_start\t#ovlp_end\t";
		    print $out1_fh "#rd_id\t#short_rd_id\t#rd_st\t#rd_en\t#snp_count\t#snp_list\n";

		    my @rds = @{$selected_ovlp{'reads'}};

		    my %curr_request_status = ();
		    if (defined $assem_exec{'strain_data'}) {
			%curr_request_status = %strain_list;
		    }


		    ## Now it will take the selection function
		    
		    my @selected_reads = @{$selected_ovlp{'reads'}};
		    if (defined $assem_exec{'strain_selection'}) {
			if ($assem_exec{'strain_selection'} =~ m/random/i) {
			    @selected_reads = shuffle(@{$selected_ovlp{'reads'}});
			}
			elsif ($assem_exec{'strain_selection'} =~ m/distance:(.+)/i) {
			    my @parent_ids = split(',', $1);

			    if (defined $1) {
				@selected_reads = order_by_distance($selected_ovlp{'reads'}, \%read, \@parent_ids);
			    }
			}
		    }

		    foreach my $sel_rd_id (@selected_reads) {
			
			## It can happens that there are more than one sequence from the same strain
			## (With the request variable at least it have selected that number, but it can have more)

			$k++;
			my $select_rd_strain = $read{$sel_rd_id}->{'strain'};
			my $selected_by_strain = 1;

			if (defined $curr_request_status{$select_rd_strain}) {
			    if ($curr_request_status{$select_rd_strain} > 0) {				
				$read{$sel_rd_id}->{'short_id'} = $select_rd_strain . '-'. $curr_request_status{$select_rd_strain};
				$curr_request_status{$select_rd_strain}--;
			    }
			    else {
				$selected_by_strain = 0;
			    }
			}
			else {
			    if (defined $assem_exec{'strain_data'}) {
				$selected_by_strain = 0;
			    }
			    $read{$sel_rd_id}->{'short_id'} = $select_rd_strain . '-'. $k;
			}
			    
			if ($selected_by_strain == 1) {
			    ## To write the sequences it should replace pad symbol (*) for a dash (-) 

			    my $seq = $read{$sel_rd_id}->{'pad_seq'};
			    $seq =~ s/\*/-/g;

			    my $rd_obj = Bio::Seq->new( -id => $read{$sel_rd_id}->{'short_id'}, -seq => $seq);
			
			    if ($read{$sel_rd_id}->{'sense'} == 0) {
				$rd_obj = $rd_obj->revcom();
			    }
	
			    my $ovlp_l = $selected_ovlp{'end'} - $selected_ovlp{'start'}; 
			    my $st = $selected_ovlp{'start'} - $read{$sel_rd_id}->{'co_start'} + $read{$sel_rd_id}->{'rd_start'};
			    my $en = $selected_ovlp{'end'} - $read{$sel_rd_id}->{'co_start'} + $read{$sel_rd_id}->{'rd_start'};
	
			    my $trim_seqobj = $rd_obj->trunc($st, $en);
		    		  
			    $seqio->write_seq($trim_seqobj);

			    ## Also it is interesting to have a list of the sequences of each contig with the overlapping region

			    my $contig_l = length($contig{$contig_id}->{'pad_seq'});
			    my $contig_s = $contig{$contig_id}->{'pad_seq'};

			    my @snps;
			    if (defined $contig{$contig_id}->{'snp_position'}) {
				my @snp_positions = @{$contig{$contig_id}->{'snp_position'}};
				my @snp_types = @{$contig{$contig_id}->{'snp_type'}};
				
				## It will add the snps of the overlapping region

				my $c = 0;
				foreach my $snp_p (@snp_positions) {
				    if ($snp_p >= $selected_ovlp{'start'} && $snp_p <= $selected_ovlp{'end'}) {
					my $contig_snp = substr($contig_s, $snp_p-1, 1);
					
					my $rd_p = $snp_p - $read{$sel_rd_id}->{'co_start'} + $read{$sel_rd_id}->{'rd_start'};
					my $read_snp = $rd_obj->subseq($rd_p, $rd_p);
					
					my $snp_tag = $snp_p . '-[' . $contig_snp . '=>' . $read_snp . ']';

					push @snps, $snp_tag;
				    }
				    $c++;
				}
				my $snp_n = scalar(@snps);
				my $snp_list = join(';', @snps);
				
				print $out1_fh "$contig_id\t$contig_l\t$ovlp_l\t$selected_ovlp{'start'}\t$selected_ovlp{'end'}\t";
				print $out1_fh "$sel_rd_id\t$read{$sel_rd_id}->{'short_id'}\t$st\t$en\t$snp_n\t$snp_list\n";
			    }		    			
			}       	    
		    }
		    $output_files{$contig_id} = { 'assembly_reads' => $filename };

		    my $contig_filename = $preprocessing_dir . '/' . $contig_id . '.contig.fasta';	
		    my $contig_seqio = Bio::SeqIO->new( -file => ">$contig_filename", -format => "fasta" );
		    my $contig_seq = Bio::Seq->new( -id => $contig_id, -seq => $contig_seq_unpad);
		    $contig_seqio->write_seq($contig_seq);
		    $output_files{$contig_id}->{'assembly_contig'} = $contig_filename;
		}
		else {
		    print STDERR "\n\t\tPROCESS MSG I: $contig_id alignment ($align_l) do not pass min alignment length ($min_align_l)."; 
		    print STDERR "SKIPPING...\n";
		}	
	    }
	    else {
		print STDERR "\n\t\tPROCESS MSG II: $contig_id alignment have not enought reads to pass strain requeriment. ";
		print STDERR "SKIPPING...\n";  
	    }
	}
	
	## Catch the contig_id and create a new file to store the sequences

	if ($line =~ m/^CO\s+(.+?)$/) { 
	    $contig_id = $1;
		
	    $contig{$contig_id} = {};

	    $n_co++;
	    %start = ();
	    %end = ();
	}

	if ($line =~ m/^CS\t(.+?)$/) { 
	    $contig_seq_pad = $1;
	    $contig{$contig_id}->{'pad_seq'} = $contig_seq_pad;
	    $contig_seq_unpad = $contig_seq_pad;
	    $contig_seq_unpad =~ s/\*//g;
	    $contig{$contig_id}->{'unpad_seq'} = $contig_seq_unpad;
	}

	if ($line =~ m/^LC\t(.+?)$/) {
    
	    $contig{$contig_id}->{'contig_length'} = $1;
	}


	## Catch the SNPs

	if ($line =~ m/^CT\s+(S\w+)\s+(\d+)\s+(\d+)\s+(\w+)/) { 
	    my $snp_type = $1;
	    my $snp_position = $2;
	    
	    if (exists $contig{$contig_id}->{'snp_type'}) {
		push @{$contig{$contig_id}->{'snp_type'}}, $snp_type;	
	    }
	    else {
		$contig{$contig_id}->{'snp_type'} = [$snp_type]
	    }
	    
	    if (exists $contig{$contig_id}->{'snp_position'}) {
		push @{$contig{$contig_id}->{'snp_position'}}, $snp_position;	
	    }
	    else {
		$contig{$contig_id}->{'snp_position'} = [$snp_position]
	    }
	}

	## Catch the read_id

	if ($line =~ m/^RD\s+(.+?)$/) { 
	    $read_id = $1;
	    
	    $read{$read_id} = {};

	    if (exists $contig{$contig_id}->{'read'}) {
		push @{$contig{$contig_id}->{'read'}}, $read_id;	
	    }
	    else {
		$contig{$contig_id}->{'read'} = [$read_id]
	    }
	    
	    $n_rd++;
	}

	if ($line =~ m/^RS\s+(.+?)$/) { 
	    $read_seq = $1;
	    $read{$read_id}->{'pad_seq'} = $read_seq;
	    my $read_seq_unpad = $read_seq;
	    $read_seq_unpad =~ s/\*//g;    
	    $read{$read_id}->{'unpad_seq'} = $read_seq_unpad;
	}

	## Strain tag 

	if ($line =~ m/^SN\s+(.+?)$/) { 
	    $read{$read_id}->{'strain'} = $1;
	}

	## The assembly coordinates will be defined for the AT tag.
	## AT contig_start contig_end read_start read_end

	if ($line =~ m/^AT\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	    if ($1 <= $2) {
		$read{$read_id}->{'co_start'} = $1;
		$read{$read_id}->{'co_end'} = $2;
		$read{$read_id}->{'sense'} = 1;
	    }
	    else {
		$read{$read_id}->{'co_start'} = $2;
		$read{$read_id}->{'co_end'} = $1;
		$read{$read_id}->{'sense'} = 0;
	    }
	
	    $read{$read_id}->{'rd_start'} = $3;
	    $read{$read_id}->{'rd_end'} = $4;
	    $read{$read_id}->{'overlap_length'} = $4 - $3;
	    
	    if (exists $start{$read{$read_id}->{'co_start'}}) {
		push @{$start{$read{$read_id}->{'co_start'}}}, $read_id;
	    }
	    else {
		$start{$read{$read_id}->{'co_start'}} = [$read_id];
	    }
	    
	    if (exists $end{$read{$read_id}->{'co_end'}}) {
		push @{$end{$read{$read_id}->{'co_end'}}}, $read_id;
	    }
	    else {
		$end{$read{$read_id}->{'co_end'}} = [$read_id];
	    }
	}
	
	print STDERR "\tParsing line=$l for file=$maf. ($n_co contigs and $n_rd reads have been parsed).      \r";
    }

    return (\%output_files, \%contig, \%read);
}

=head2 order_by_distance

  Usage: my @order_seq = order_by_distance($ids_aref, $seq_href, $seed_id);

  Desc: Run some blast comparissons between sequences and order them based in the distance
        to the seed_id

  Ret: An array with the sequences ordered

  Args: $ids_aref, an array reference with the sequence ids
        $seq_href, a hash reference with key=id and hash=hash ref seq data
        $parent_ids_aref, array reference for parent selection  

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my @order_seq = order_by_distance($ids_aref, $seq_href, $parents_aref);

=cut

sub order_by_distance {
    my $ids_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None ids_aref array reference have been supplied to order_by_distance function.\n");
    
    my $reads_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None read hash reference have been supplied to order_by_distance function.\n");

    my $parent_ids_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None seed_id have been supplied to order_by_distance function.\n");

    my @ordered_ids = ();
    
    my %selection;
    my @parent_ids = ();
    foreach my $parent (@{$parent_ids_aref}) {
	my @par = split('>', $parent);
	my @chi = split('\*', $par[1]);
	$selection{$par[0]} = { $chi[0] => $chi[1] };
    }

    ## First, put the sequences into seq objects

    my @reads = ();
    my @parents = ();
    my %children = ();

    foreach my $id (@{$ids_aref}) {

	my $strain = $reads_href->{$id}->{'strain'};
	if (defined $selection{$strain}) {
	    push @parents, $id;
	}
	else {
	    $children{$id} = 1;
	}
	my $rd_obj = Bio::Seq->new( -id => $id, -seq => $reads_href->{$id}->{'unpad_seq'});
	push @reads, $rd_obj;
    }

    ## Run the alignment

    my @params = ('quiet' => 1, 'matrix' => 'BLOSUM');
    my $factory = Bio::Tools::Run::Alignment::Clustalw->new(@params);
    my $align = $factory->align(\@reads);
      
    ## Calculate distances

    my $stats = Bio::Align::DNAStatistics->new();
    my $jcmatrix = $stats->distance( -align => $align, 
				     -method => 'Kimura');
    my %selected_ids;

    ## Now it will get distances
    foreach my $parent_id (@parents) {
	
	push @ordered_ids, $parent_id;

	my $pstrain =  $reads_href->{$parent_id}->{'strain'};

	my %selec = %{$selection{$pstrain}};
	
	foreach my $chi (keys %selec) {
	    
	    my $pass = $selec{$chi};
	    while ( $pass > 0) {

		if (scalar(keys %children) > 0) {

		    my %dist = ();
		    foreach my $children_id (keys %children) {
			unless (exists $selected_ids{$children_id}) {
			    my $distance_value = $jcmatrix->get_entry($parent_id, $children_id);
			    $dist{$distance_value} = $children_id;
			}
		    }

		    ## Now order
		    my @ord_children = sort {$a <=> $b} keys %dist;

		    if (scalar(@ord_children) > 0) {
		    
			push @ordered_ids, $dist{$ord_children[0]};
			$selected_ids{$dist{$ord_children[0]}} = 1;

			## And delete for the hash
			delete $children{$ord_children[0]};
		    }
		}
		$pass--;
	    }
	   
	}
    }

    return @ordered_ids;
}

=head2 run_cds_prediction

  Usage: my %out_files = run_cds_prediction(\%cds_predict_args, \%fasta_files, $out_dirname);

  Desc: Run the cds prediction program over the two set of fasta files, input read files and input contig
        files 

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%align_args, a hash reference with keys=argument_type and value=value
        \%fasta_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_cds_prediction(\%align_args, \%input_files, $out_dirname);

=cut

sub run_cds_prediction {
    my $cds_predict_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None cds_prediction_argument hash reference have been supplied to run_cds_prediction function.\n");
    
    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to run_cds_prediction function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_cds_prediction function.\n");

    my %input_files = %{$input_files_href};

    my %outfiles;

    ## First, create a log file

    my $log_filename = $output_dir . '/snpsel_pipeline.step2.log';
    open my $log, '>', $log_filename;


    my $n = 0;
    my $N = scalar(keys %input_files);
    print STDERR "\n";

    foreach my $contig_id (sort keys %input_files) {
	$n++;

	## First get the basename

	my %files = %{$input_files{$contig_id}};
	
	if (defined $files{'assembly_contig'} && defined $files{'assembly_reads'}) {

	    my $contig_basename = basename($files{'assembly_contig'});
	    my $reads_basename = basename($files{'assembly_reads'});

	    ## Create the output_basename

	    my $contig_cds = $output_dir . '/' . $contig_basename;
	    $contig_cds =~ s/\.fasta/\.estscan_cds\.fasta/;
	    my $reads_cds = $output_dir . '/' . $reads_basename;
	    $reads_cds =~ s/\.fasta/\.estscan_cds\.fasta/;
	    my $contig_proteins = $output_dir . '/' . $contig_basename;
	    $contig_proteins =~ s/\.fasta/\.estscan_protein\.fasta/;
	    my $reads_proteins = $output_dir . '/' . $reads_basename;
	    $reads_proteins =~ s/\.fasta/\.estscan_protein\.fasta/;
	    
	    ## Prepare the command line as executable + options + input/output

	    print STDERR "\tRunning Estscan for contig $contig_id ($n/$N)\r";

	    my $command_contig = $cds_predict_args_href->{'cds_exec'} 
	                        . ' -M ' . $cds_predict_args_href->{'cds_matrix'}
	                        . ' -o ' . $contig_cds
                                . ' -t ' . $contig_proteins
			        . ' -n ' . $files{'assembly_contig'};
			     
	    print $log "RUNNING ESTSCAN over contig=$contig_id:\n\t$command_contig\n";
	    system($command_contig);

	    my $command_reads = $cds_predict_args_href->{'cds_exec'} 
	                       . ' -M ' . $cds_predict_args_href->{'cds_matrix'}
	                       . ' -o ' . $reads_cds
                               . ' -t ' . $reads_proteins
			       . ' -n ' . $files{'assembly_reads'};
			     
	    print $log "RUNNING ESTSCAN over reads from $contig_id:\n\t$command_reads\n\n";
	    system($command_reads);

	    ## Store the new files

	    $input_files{$contig_id}->{'cds_contig'} = $contig_cds;
	    $input_files{$contig_id}->{'cds_reads'} = $reads_cds;
	    $input_files{$contig_id}->{'protein_contig'} = $contig_proteins;	
	    $input_files{$contig_id}->{'protein_reads'} = $reads_proteins;
	}
	else {
	    print STDERR "WARNING: contig=$contig_id has not associated any assembly_contig or assembly_reads file.\n";
	}
    }
    print STDERR "\n\n";

    return %input_files;
}

=head2 file_class_counting

  Usage: my %count_files = file_class_counting(\%master_files);

  Desc: Count the number of files in the master_file hash

  Ret: %count_files, a hash with keys=type and value=count

  Args: \%master_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %count_files = file_class_counting(\%master_files);

=cut

sub file_class_counting { 
    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to file_class_counting function.\n");

    my %input_files = %{$input_files_href};

    my %count;

    foreach my $id (keys %input_files) {

	if (defined $input_files{$id}) {
	    my %files = %{$input_files{$id}};
	    
	    foreach my $class (keys %files) {
		if (defined $count{$class}) {
		    $count{$class}++;
		}
		else {
		    $count{$class} = 1;
		}	
	    }
	}
    }
    
    return %count; 
}



=head2 run_seqformating

  Usage: my %out_files = run_seqformating(\%prot_filter_args, \%master_files, $out_dirname);

  Desc: Replace X for N in the cds sequences and remove ';' from the
        id in proteins and cds sequences

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%prot_filter, a hash reference with keys=argument and value=filter value
        \%master_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_seqformating(\%prot_filter, \%master_files, $out_dirname);

=cut

sub run_seqformating { 
    my $prot_filter_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None protein_filter_argument hash reference have been supplied to run_seqformating function.\n");

    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to run_seqformating function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_seqformating function.\n");


    my %prot_filter = %{$prot_filter_args_href};
    my %input_files = %{$input_files_href};

    ## Define filter

    my ($prot_size, $cds_size);
    if (defined $prot_filter{'prot_size_filter'}) {
	$prot_size = $prot_filter{'prot_size_filter'};
	$cds_size = $prot_size * 3;
    }

    my %outfiles;

    ## First, create a log file

    my $log_filename = $output_dir . '/snpsel_pipeline.step3.log';
    open my $log, '>', $log_filename;


    my $n = 0;
    my $N = scalar(keys %input_files);
    print STDERR "\n";

    foreach my $contig_id (sort keys %input_files) {
	$n++;

	## First get the basename

	my %files = %{$input_files{$contig_id}};

	foreach my $type (keys %files) {

	    print STDERR "\tReformating [$n/$N] contig $contig_id (file=$files{$type})\r";

	    ## Take only the cds_ and proteins_ files

	    
	    if ($type =~ m/^[cds|protein]/) {
		    
		## To catch files with proteins or cds

		if (-s $files{$type}) {

		    my $seqin_io = Bio::SeqIO->new( -file => "$files{$type}", -format => "fasta" );
				
		    ## It will count the number of sequences and how many sequences pass the filter.
		    ## If both numbers are the same it will create the output file and add the sequences
		
		    my $seq_n = 0;
		    my $nomet = 0;
		    my @seqs;

		    while( my $seq = $seqin_io->next_seq() ) {
			$seq_n++;

			my $sequence = $seq->seq();
			my $id = $seq->id();
			my $seqtype = $seq->alphabet();
			if ($seqtype eq 'dna') {
			    $sequence =~ s/X/N/gi;
			}
			$id =~ s/;//;
			$seq->id($id);
			$seq->seq($sequence);
			my $seqlength = length($sequence);
			
			## It will filter the seq by size, if none size_protein
			## is defined, it will skip this part

			if (defined $prot_size) {
			    if ($seqtype eq 'dna') {
				if ($seqlength > $cds_size) {			
				    push @seqs, $seq;
				}
				else {
				    print $log "\ncontig_id=$contig_id with type=$type has sequence length = $seqlength (< $cds_size). SKIPPING SEQUENCE.\n";
				}
			    }
			    else {
				if ($seqlength > $prot_size) {
				    if (defined $prot_filter{'first_met_filter'} && $prot_filter{'first_met_filter'} =~ m/^Y|1/i) {
					if ($sequence =~ m/^M/i) {
					    push @seqs, $seq;
					}
					else {
					    $nomet = 1;
					    print $log "\ncontig_id=$contig_id with type=$type has not first metionine. SKIPPING SEQUENCE.\n";
					}
				    }
				    else {
					push @seqs, $seq;
				    }
				}
				else {
				    print $log "\ncontig_id=$contig_id with type=$type has sequence length = $seqlength (< $prot_size). SKIPPING SEQUENCE.\n";
				}
			    }
			}
			else {
			    if (defined $prot_filter{'first_met_filter'} && $prot_filter{'first_met_filter'} =~ m/^Y|1/i) {
				if ($sequence =~ m/^M/i) {
				    push @seqs, $seq;
				}
				else {
				    $nomet = 1;
				    print $log "\ncontig_id=$contig_id with type=$type has not first metionine. SKIPPING SEQUENCE.\n";
				}
			    }
			    else {
				push @seqs, $seq;			
			    }
			}
		    }
		    if (scalar(@seqs) == $seq_n) {
		    
			## Create the file
			
			my $new_type = 'formated_'. $type;
			
			my $basename = basename($files{$type});
			$basename =~ s/protein\./formated_protein\./;
			$basename =~ s/cds\./formated_cds\./;
			$files{$new_type} = $output_dir . '/' . $basename;
			
			my $seqout_io = Bio::SeqIO->new( -file => ">$files{$new_type}", -format => "fasta" );
			foreach my $seq_filter (@seqs) {
			    $seqout_io->write_seq($seq_filter);
			}
		    }
		    else {
			print STDERR "\n\tMESSAGE: sequence contig_id=$contig_id and type=$type have predicted protein size < $prot_size or cds_size < $cds_size. SKIPPING SEQUENCE.\n";
			if ($nomet == 1) {
			    print STDERR "\n\tMESSAGE: contig_id=$contig_id with type=$type has not first metionine. SKIPPING SEQUENCE.\n";
			}
		    }
		

		    $input_files{$contig_id} = \%files;
		}
		else {
		    print STDERR "\n\tMESSAGE: Sequence file:$files{$type} have not any sequence. SKIPPING FILE.\n";
		}
	    }
	}
    }
    print STDERR "\n\n";

    return %input_files;
}

=head2 run_alignment

  Usage: my %out_files = run_alignment(\%align_args, \%master_files, $out_dirname);

  Desc: Run the aligment program over the fasta sequences

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%align_args, a hash reference with keys=argument and value=value
        \%master_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_alignment(\%align_args, \%master_files, $out_dirname);

=cut

sub run_alignment { 
    my $align_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None alignment_argument hash reference have been supplied to run_alignment function.\n");

    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to run_alignment function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_alignment function.\n");

    my %input_files = %{$input_files_href};

    my %outfiles;

    ## First, create a log file

    my $log_filename = $output_dir . '/snpsel_pipeline.step4.log';
    open my $log, '>', $log_filename;


    my $n = 0;
    my $N = scalar(keys %input_files);
    print STDERR "\n";

    foreach my $contig_id (sort keys %input_files) {
	$n++;

	## First get the basename

	my %files = %{$input_files{$contig_id}};

	## Calculate the blast if a reference file is defined

	my %match;

	$files{'protein_ref'} = $files{'formated_protein_contig'};

	if (defined $align_args_href->{'blast_exec'} &&  $align_args_href->{'ref_protein_file'} ) {
	    my $blast_output = $output_dir . '/' . $contig_id . '.blastx.ref.m8';
	    my $blast_command = "$align_args_href->{'blast_exec'} -p blastx -i $files{'assembly_contig'} -d $align_args_href->{'ref_protein_file'} -o $blast_output -m 8";

	    print $log "RUNNING BLAST over contig=$contig_id:\n\t$blast_command\n";
	    system($blast_command);

	    ## Now it will get the results
	    open my $blast_fh, '<', $blast_output;
	    while(<$blast_fh>) {
		chomp($_);
		my @data = split(/\t/, $_);

		if ($data[0] eq $contig_id && !$match{$data[0]}) {
		    $match{$data[0]} = $match{$data[1]};
		}
	    }

	    ## To get the sequence it will use fastacmd

	    if (defined $match{$contig_id}) {
		my $ref_sequence_file =  $output_dir . '/' . $contig_id . '.ref_protein.fasta';
		my $extract_command = "fastacmd -d $align_args_href->{'ref_protein_file'} -s $match{$contig_id} -o $ref_sequence_file";
		
		print $log "RUNNING FASTACMD over contig=$contig_id:\n\t$extract_command\n";
		system($blast_command);
	    
		$files{'protein_ref'} = $ref_sequence_file;
	    }
	}

	## It will run the script over the protein_reads types

	if (defined $files{'formated_protein_reads'} && defined $files{'formated_protein_contig'}) {

	    my $basename = basename($files{'formated_protein_reads'});
	    $basename =~ s/\.fasta$/\.clw/;
	    my $contig_basename = basename($files{'formated_protein_contig'});
	    my $global_basename = $contig_basename;
	    $contig_basename =~ s/\.fasta$/\.clw/;
	    $global_basename =~ s/contig/global/;
	
	    $files{'formated_protein_global'} = $output_dir . '/'. $global_basename; 
	
	    ## It always cat with protein ref, but it can be the contig or a reference protein depending of the options used

	    my $cat_command = "cat $files{'protein_ref'} $files{'formated_protein_reads'} > $files{'formated_protein_global'} ";
	    system($cat_command);
	    print $log "\n RUNNING CAT: $cat_command\n\n";

	    $files{'protein_alignment'} = $output_dir . '/' . $basename;
	    $files{'contig_protein_alignment'} = $output_dir . '/' . $contig_basename;

	    my ($command, $contig_command);

	    if ($align_args_href->{'align_exec'} =~ m/muscle/) {
		$command = $align_args_href->{'align_exec'} 
		       . ' -in ' . $files{'formated_protein_reads'}
	               . ' -out ' . $files{'protein_alignment'}
	               . ' -clwstrict -stable -quiet'
	               . ' -loga ' . $log_filename;
		$contig_command = $align_args_href->{'align_exec'} 
		       . ' -in ' . $files{'formated_protein_global'}
		       . ' -out ' . $files{'contig_protein_alignment'}
		       . ' -clwstrict -stable -quiet'
		       . ' -loga ' . $log_filename;
	    
	    }
	    elsif ($align_args_href->{'align_exec'} =~ m/clustalw/) {
		$command = $align_args_href->{'align_exec'} 
	               . ' -INFILE=' . $files{'formated_protein_reads'}
	               . ' -OUTFILE=' . $files{'protein_alignment'}
	               . ' -OUTORDER=INPUT';
		$contig_command = $align_args_href->{'align_exec'} 
	               . ' -INFILE=' . $files{'formated_protein_global'}
	               . ' -OUTFILE=' . $files{'contig_protein_alignment'}
	               . ' -OUTORDER=INPUT';
	    }
	    else {
		print STDERR "WARNING: alignment executable is not muscle or clustaw. SKIPPING alignment.\n";
	    }
	
	    if (defined $align_args_href->{'add_align_args'}) {
		$command .= ' ' . $align_args_href->{'add_align_args'};
		$contig_command .= ' ' . $align_args_href->{'add_align_args'};
	    }

	    print $log "\n RUNNING $command\n\n";
	    print STDERR "\tRunning Protein Alignment for contig $contig_id ($n/$N)\r";

	    system($command);
	    system($contig_command);

	    $input_files{$contig_id} = \%files;

	}
    }
    print STDERR "\n\n";

    return %input_files;
}


=head2 run_coverage

  Usage: my %out_files = run_coverage(\%coverage_args, \%fasta_files, $out_dirname);

  Desc: Run the coverage calculation over the global alignments

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%coverage_args, a hash reference with keys=argument_type and value=value
        \%fasta_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_coverage(\%align_args, \%input_files, $out_dirname);

=cut

sub run_coverage {
    my $coverage_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None coverage_argument hash reference have been supplied to run_coverage function.\n");
    
    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to run_coverage function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_coverage function.\n");

    my %input_files = %{$input_files_href};
    my %cov_args = %{$coverage_args_href};

    ## It will take each alignment 

    my $coverage_file =  $output_dir . '/snpsel_pipeline.step5.coverage.tab';
    open my $cov_fh, '>', $coverage_file;
    print $cov_fh "#alignment\t#id_list(ref+elem)\t#alignment_length\t#n_sequences\t#identity\t#coverage\t#qualification\n";

    foreach my $id (sort keys %input_files) {
	
	my %file = %{$input_files{$id}};

	if (defined $file{'contig_protein_alignment'}) {
	    
	    my $align_io = Bio::AlignIO->new(-file   => $file{'contig_protein_alignment'} , -format => 'clustalw');

	    print STDERR "\tRunning alignment coverage filtering over alignment=$file{'contig_protein_alignment'} \r";

	    while (my $align = $align_io->next_aln() ) {

		## Take the alignment data

		my $align_length = $align->length();
		my $n_seq = $align->num_sequences();
		my $ident = $align->percentage_identity();
		my $ident_obj = Math::BigFloat->new($ident);
		my $f_ident = $ident_obj->bfround(-3);

		my $consensus = $align->consensus_string(50);
		$consensus =~ s/\?//g;
		my $consensus_length = length($consensus);
		my $coverage = ($consensus_length / $align_length)*100;
		my $coverage_obj = Math::BigFloat->new($coverage);
		my $f_coverage = $coverage_obj->bfround(-3);
		
		my @ids;
		foreach my $seq ($align->each_seq) {
		    my $seq_id = $seq->id();
		    push @ids, $seq_id;
		}

		my $id_list = join(',', @ids);

		## check if pass or not the filter
		my $qualif = 'PASS';

		if (defined $cov_args{'prot_cov_filter'} ) {
		    if ($coverage <  $cov_args{'prot_cov_filter'}) {
			$qualif = 'NO PASS COVERAGE FILTER';
		    }
		}
			
		if (defined $cov_args{'prot_ident_filter'} ) {
		    if ($ident <  $cov_args{'prot_ident_filter'}) {
			if ($qualif ne 'PASS') {
			    $qualif .= '+ NO PASS IDENTITY FILTER';
			}
			else {
			    $qualif = 'NO PASS IDENTITY FILTER';
			}
		    }
		}
		
		## Finally copy as new files into the new folder files that have pass the filter

		if ($qualif eq 'PASS') {
		    my $basename = basename($file{'protein_alignment'});
		    my $filename = $output_dir . '/' . $basename;
		    $file{'protein_coverage'} = $filename;
		    
		    my $cp_cmd = "cp $file{'protein_alignment'} $filename";
		    system($cp_cmd);
		}
		else {
		    print STDERR "\n\talignment=$file{'contig_protein_alignment'} have not pass coverage filtering. SKIPPING FILE.\n";
		}

		my $infile = basename($file{'contig_protein_alignment'});
		print $cov_fh "$infile\t$id_list\t$align_length\t$n_seq\t$f_ident\t$f_coverage\t$qualif\n";
	    }
	}
	$input_files{$id} = \%file;
    }
    print STDERR "\n\n";
    return %input_files;
}


=head2 run_pal2nal

  Usage: my %out_files = run_pal2nal(\%pal2nal_args, \%master_files, $out_dirname);

  Desc: Run the pal2nal script over the alignments to get files with paml format

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%pal2nal_args, a hash reference with keys=argument_type and value=value
        \%master_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_pal2nal(\%pal2nal_args, \%input_files, $out_dirname);

=cut

sub run_pal2nal {
    my $pal2nal_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None pal2nal_argument hash reference have been supplied to run_pal2nal function.\n");
    
    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to run_pal2nal function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_pal2nal function.\n");

    my %input_files = %{$input_files_href};
    my %pal2nal_args = %{$pal2nal_args_href};

    ## Create a log file to redirect the commands used 
    
    my $log_filename = $output_dir . '/snpsel_pipeline.step6.log';
    open my $log, '>', $log_filename;

    ## Pal2nal script check

    unless (defined $pal2nal_args{'pal2nal_exec'}) {
	die("CONTROL FILE ERROR: pal2nal_executable is not defined.\n");
    }

    ## Pal2nal needs protein alignment file and a cds file. Protein files must be
    ## the translation of the cds file. X are not permited in cds file and the
    ## sequence id must to be in the same order in both files

    foreach my $id (sort keys %input_files) {
	
	my %file = %{$input_files{$id}};

	## pal2nal will use two files: 'protein_coverage' from the last step
	## filtered by coverage and identity and a cds file of the sequences 
	## used to predict the proteins: 'formated_cds_reads'

	if (defined $file{'protein_coverage'} && $file{'formated_cds_reads'}) {
	    
	    ## Output should be:

	    my $basename = basename($file{'protein_coverage'});
	    $basename =~ s/\.clw/\.paml/;
	    my $outfile = $output_dir . '/' . $basename;

	    ## prepare the command

	    my $pal2nal_cmd = "$pal2nal_args{'pal2nal_exec'} $file{'protein_coverage'} $file{'formated_cds_reads'} -output paml";

	    if (defined $pal2nal_args{'more_pal2nal_args'}) {
		$pal2nal_cmd .= " $pal2nal_args{'more_pal2nal_args'}";
	    }

	    $pal2nal_cmd .= "> $outfile ";

	    print $log "RUNNING $pal2nal_cmd\n";
	    system($pal2nal_cmd);

	    ## Add to the file hash

	    $file{'cds_alignment_paml'} = $outfile;

	    ## Also it will run the same command to get codon file

	    my $pal2nal_codon_cmd = "$pal2nal_args{'pal2nal_exec'} $file{'protein_coverage'} $file{'formated_cds_reads'} -output codon";

	    if (defined $pal2nal_args{'more_pal2nal_args'}) {
		$pal2nal_codon_cmd .= " $pal2nal_args{'more_pal2nal_args'}";
	    }

	    my $outfile2 = $outfile;
	    $outfile2 =~ s/\.paml/\.codon/;

	    $pal2nal_codon_cmd .= "> $outfile2 ";

	    print $log "RUNNING CODON $pal2nal_codon_cmd\n";
	    system($pal2nal_codon_cmd);

	    $file{'cds_alignment_codon'} = $outfile2;
	 
	    $input_files{$id} = \%file;
	}
    }
    print STDERR "\n\n";
    return %input_files;
}

=head2 run_paml

  Usage: my %out_files = run_paml(\%paml_args, \%master_files, $out_dirname);

  Desc: Create the ctr files, run the paml program over the paml files and
        move the output to another name before rerun the program with another file.

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%paml_args, a hash reference with keys=argument_type and value=value
        \%master_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_paml(\%paml_args, \%input_files, $out_dirname);

=cut

sub run_paml {
    my $paml_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None pal2nal_argument hash reference have been supplied to run_paml function.\n");
    
    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to run_paml function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_paml function.\n");

    my %input_files = %{$input_files_href};
    my %paml_args = %{$paml_args_href};

    ## Create a log file to redirect the commands used 
    
    my $log_filename = $output_dir . '/snpsel_pipeline.step7.log';
    open my $log, '>', $log_filename;

    ## Paml script check

    unless (defined $paml_args{'paml_exec'}) {
	die("CONTROL FILE ERROR: paml_executable is not defined.\n");
    }

    ## Paml works with a control (ctr) file that supply all the info
    ## about how should be run the program. The ctr file alway have the
    ## same name, executable + . + ctr
    ## for example ym00 will have the ctr file ym00.ctr
    
    my $ctr_file = basename($paml_args{'paml_exec'}) . '.ctl';

    ## After the creation of the ctr file and it is use, it will move to 
    ## ctr folder

    my $ctr_folder = $output_dir . '/ctr_files';
    mkdir($ctr_folder);
    my $paml_folder = $output_dir . '/paml_files';
    mkdir($paml_folder);

    foreach my $id (sort keys %input_files) {
	
	my %file = %{$input_files{$id}};

	if (defined $file{'cds_alignment_paml'}) {

	    ## First create a ctr file.
	    ## For now it will deal with yn00 ctr files but 
	    ## in some point we should implement more methods 
	    ## like codeml

	    if ($paml_args{'paml_exec'} =~ m/yn00/) {

		print STDERR "RUNNING $paml_args{'paml_exec'} with $ctr_file over the file: $file{'cds_alignment_paml'}\n";

		my $outfile = $output_dir . '/' . basename($file{'cds_alignment_paml'});
		$outfile =~ s/\.paml/\.yn00_output\.txt/;

		open my $ctr_fh, '>', $ctr_file;
		print $ctr_fh "seqfile = $file{'cds_alignment_paml'}\n";
		print $ctr_fh "outfile = $outfile\n";
		
		my %default_args = ( 'verbose' => 0, 'icode' => 0, 'weighting' => 0, 'commonf3x4' => 0);
		if (defined $paml_args{'more_paml_args'} ) {
		    
		    my @args = split(';', $paml_args{'more_paml_args'});
		    foreach my $arg (@args) {
			if ($arg =~ m/(\w+)\s*=\s*(\d+)/) {
			    if (defined $default_args{$1}) {
				$default_args{$1} = $2;
			    }
			    else {
				die("CONTROL FILE ERROR: <ADD_PAML_ARGS> = $arg is not a valid control file line for yn00 paml program.\n");
			    }
			}
		    }
		}

		## Now complete the ctr file

		foreach my $key (keys %default_args) {
		    print $ctr_fh "$key = $default_args{$key} \n";
		}
		close $ctr_fh;

		## The ctr file should be complete now, so run the paml program

		system($paml_args{'paml_exec'});
		print $log "RUNNING PAML PROGRAM: $paml_args{'paml_exec'}\n\tCONTROL FILE: $ctr_file\n\tOUTPUT_FILE: $outfile\n\n";

		## It should have produce a result file, so now it will load the filename into the file hash

		$file{'paml_result'} = $outfile;

		## Finally it will change the name of the control file and the paml files

		my $new_ctr_file = $outfile;
		$new_ctr_file =~ s/\.yn00_output\.txt/\.yn00_ctr/;
		my $mv_cmd = "mv $ctr_file $new_ctr_file";
		system($mv_cmd);
		$file{'ctr_file'} = $new_ctr_file;
		print $log "RUNNING MOVING CMD: $mv_cmd\n";

		my @paml_files = ('2YN.dN', 'rub', '2YN.dS', 'rst', '2YN.t', 'rst1');
		foreach my $paml_file (@paml_files) {
		    my $new_name = $paml_folder . '/' . basename($outfile);
		    $new_name =~ s/\.yn00_output\.txt/\.yn00_$paml_file/;
		    my $mv2_cmd = "mv $paml_file $new_name";
		    system($mv2_cmd);
		    print $log "RUNNING MOVING 2 CMD: $mv2_cmd\n";
		    $file{'paml_'.$paml_file} = $new_name;
		}
	    }

	 
	    $input_files{$id} = \%file;
	}
    }
    print STDERR "\n\n";
    return %input_files;
}

=head2 run_paml_parsing

  Usage: my ($result_filename, $distrib) = run_paml_parsing(\%master_files, $out_dirname);

  Desc: Parse the paml output files (yn00) to get omega (dS/dN) values

  Ret: An array with the names of the result filenames

  Args: \%master_files, a hash reference with keys=contig_id and value=hash ref
                       with key=type and value=filename
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my ($result_filename, $distrib) = run_paml_parsing(\%input_files, $out_dirname);

=cut

sub run_paml_parsing {
    my $input_files_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file hash reference have been supplied to run_paml function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_paml function.\n");

    my %input_files = %{$input_files_href};

    ## Create a log file to redirect the commands used 
    
    my $log_filename = $output_dir . '/snpsel_pipeline.step8.log';
    open my $log, '>', $log_filename;

    ## This will take each file, parse them, and put into a results file.

    my $result_filename = $output_dir . '/snpsel_result.tab';
    open my $result_fh, '>', $result_filename;

    ## Also it will take pairs to calculate the distribution of omega values based into
    ## amount of relations

    my %pairs;

    my $distrib_filename = $output_dir . '/snpsel_distribution.tab';
    open my $distr_fh, '>', $distrib_filename;

    ## Finally another file that it will produce will be aq columnar table with -f1 contig
    ## and -f2 to X omega values for different comparissons

    my %contigs;
    my $contig_filename = $output_dir . '/snpsel_contig.tab';
    open my $contig_fh, '>', $contig_filename;


    foreach my $id (sort keys %input_files) {
	
	my %file = %{$input_files{$id}};

	my %seqid;

	if (defined $file{'paml_result'} && $file{'paml_2YN.t'}) {

	    open my $paml_2yn_fh, '<', $file{'paml_2YN.t'};

	    my $n = 1;
	    while(<$paml_2yn_fh>) {
		chomp($_);
		      
		## First, get the sequences id in order
		## using the file 2YN.t (more easy to parse)

 		if ($_ =~ m/^(\w.+?)\s+/) {
		    $seqid{$n} = $1;
		    $n++;
		}
	    }
	    close($paml_2yn_fh);

	    ## Now it will take the data from the result file 
	    
	    open my $paml_fh, '<', $file{'paml_result'};

	    my $match_region = 0;
	    while(<$paml_fh>) {
		chomp($_);

		if ($_ =~ m/\(B\)\s+Yang\s+&\s+Nielsen\s+\(2000\)\s+method/) {
		    $match_region = 1;
		}

		if ($match_region == 1 && $_ =~ m/^\s+\w+/) {
		    my @data = split(/\s+/, $_);
		    my $seq1 = $data[1];
		    my $seq2 = $data[2];
		    my $S = $data[3];
		    my $N = $data[4];
		    my $t = $data[5];
		    my $kappa = $data[6];
		    my $omega = $data[7];
		    my $dN = $data[8] . ' ' . $data[9] . ' ' . $data[10];
		    my $dS = $data[11] . ' ' . $data[12] . ' ' . $data[13];

		    my $selection;
		    if ($omega == 99.0000) {
		     	$selection = "NA (SAME CDS SEQUENCE)";
		    }
		    elsif ($omega < 99.0000 && $omega > 1.5) {
		     	$selection = "STRONG POSITIVE SELECTION";
		    } 
		    elsif ($omega < 1.5 && $omega > 1.2) {
		     	$selection = "POSITIVE SELECTION";
		    }
		    elsif ($omega < 1.2 && $omega > 1) {
		     	$selection = "CLOSE NEUTRALITY";
		    }
		    elsif ($omega == 1) {
		     	$selection = "NEUTRALITY";
		    }
		    elsif ($omega < 1 && $omega > 0.8) {
		     	$selection = "CLOSE NEUTRALITY";
		    }
		    elsif ($omega < 0.8 && $omega > 0.5) {
		     	$selection = "NEGATIVE SELECTION";
		    }
		    elsif ($omega < 0.5 && $omega > 0) {
		     	$selection = "STRONG NEGATIVE SELECTION";
		    }
		    elsif ($omega == 0) {
		     	$selection = "NA (SAME PROTEIN, DIFFERENT CDS SEQUENCE)";
		    }
		    
		    print $result_fh "$id\t$seqid{$seq1}\t$seqid{$seq2}\t$omega\t$dN\t$dS\t$selection **\n";

		    my @pairs = ($seqid{$seq1}, $seqid{$seq2});
		    my $seqpair = join(',', sort @pairs);

		    ## take the data for the contigs

		    if (exists $contigs{$id}) {
			$contigs{$id}->{$seqpair} = $omega;
		    }
		    else {
			$contigs{$id} = { $seqpair => $omega};
		    }

		    ## Now it will calculate the pair for the distribution

		    my @omega_fractions = (99.0, 2.0, 1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0.0);
		    foreach my $fract (@omega_fractions) {
			if ($omega <= $fract) {
			    if (exists $pairs{$seqpair}) {
				
				if (exists $pairs{$seqpair}->{$fract}) {
				    $pairs{$seqpair}->{$fract}++;
				}
				else {
				    $pairs{$seqpair}->{$fract} = 1;
				}
			    }
			    else {
				$pairs{$seqpair} = { $fract => 1 };
			    }
			}
		    }
		}
		
		if ($_ =~ m/\(C\)\s+LWL85,\s+LPB93\s+&\s+LWLm\s+methods/) {
		    $match_region = 0;
		}
	    }
	    close($paml_fh);
	}
    }

    ## Now it will print the distributions

    foreach my $pair (sort keys %pairs) {
	my $p = 0;
	my @fraction = sort keys %{$pairs{$pair}};
	foreach my $fraction (@fraction) {
	    my $value;
	    if ($p > 0) {
		$value = $pairs{$pair}->{$fraction} - $pairs{$pair}->{$fraction[$p-1]};
	    }
	    else {
		$value = $pairs{$pair}->{$fraction};
	    }
	    $p++;
	    print $distr_fh "$pair\t$fraction\t$value\n";
	}
	print $distr_fh "\n";
    }

    ## And print also the contig file

    ## first print headers

    print $contig_fh "#ID";

    foreach my $omega_seqpair (sort keys %pairs) {
	print $contig_fh "\t#$omega_seqpair";
    }
    print $contig_fh "\n";
    

    foreach my $contig_id (sort keys %contigs) {
	print $contig_fh "$contig_id";
	my %omega_pairs = %{$contigs{$contig_id}};

	foreach my $omega_seqpair (sort keys %pairs) {
	    print $contig_fh "\t";
	    if (defined $omega_pairs{$omega_seqpair}) {
		print $contig_fh "$omega_pairs{$omega_seqpair}";
	    }
	}
	print $contig_fh "\n";
    }


    print STDERR "\n\n";
    return ($result_filename, $distrib_filename, $contig_filename);
}


#####
1; ##
#####
