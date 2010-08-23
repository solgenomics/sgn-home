#!/usr/bin/perl

=head1 NAME

 phylotypes_pipeline.pl
 Pipeline to create phylogenetic tree from maf assembly files (version.0.1).

=cut

=head1 SYPNOSIS

 phylotypes_pipeline.pl [-h] -i <maf_file> -c <control_file> [-s <step>] [-d <input_dir>] 
                          [-o <output_basename>] [-f <input_format>] [-X][-D]

=head1 EXAMPLE

 contigs_phyloanalysis.pl [-h] -i coffea.maf -c maf_to_phylo.default.ctr

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

=item -f

B<input_format>                 input format for files, when the pipeline start in a step different from 1
                                (more than one format can be specified, it will match the formats, the first
				format must to be the input format, secondary formats will be used for other tasks)

=item -X

B<create_control_file>          create the control file with the name maf_to_phylo.ctr

=item -D

B<create_default_control_file>  create the control file with the name maf_to_phylo.ctr

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

 This script is a pipeline that calculate the phylotree for contigs and
 the reads of these contigs from a MIRA assembly. It will create at least
 file per contig.

 It have the following steps:

    0- Parse controler file
       Load input files (maf or step files)
 
    1- Extract the contigs and the reads from maf file

    2- Rerun alignments.
       (Programs like clustalw or muscle will be used as executables)

    3- Apply bootstrapping or delete-half-jackknifing process to the 
       alignments.
       (Programs like bootseq will be used)

    4- Calculate distance matrix based in the alignments
       (Programs like dnadist or protdist will be used as executables)

    5- Calculate the tree
       (Programs like neighbor will be used as executables)
    
    6- Modify the tree (for example reroot it)
       (Programs like retree could be used)

    7- Calculate the consensus tree and the bootstrap values
       (Programs like consense will be used as executables)

    8- Parse the results (trees), filter and calculate the
       tree types.

 This pipeline can start from any step supplying the input files.
 In that case the controler need to have the output_format parameter
 set for the previous step. 

 To create a control file: 
    contigs_phyloanalysis.pl -X

 To create a default control file
    contigs_phyloanalysis.pl -D

 This default control file will use the following programs and arguments.

    1) Extract maf file with 1 seq/strain and all the strains 
       and minimum overlapping length of 50 bp.

    2) Skip the realignment

    3) Bootstrap with seqboot program and 1000 replicates.
       Bootstrap cutoff value for final tree parsing: 500

    4) Distance matrix calculated with dnadist and Kimura algorithm.
       
    5) Tree calculated with neighbor

    6) Retree the result to root the trees using the middle point
       with retree program

    7) Calculate the consensus tree and use the following tags for
       the branch length: 0.001=>SIM (similar) and 1.0 => DIF
       (different)
 

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
use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Align::DNAStatistics;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;
use Math::BigFloat;
use Statistics::Basic qw(:all);

our ($opt_i, $opt_o, $opt_c, $opt_s, $opt_d, $opt_f, $opt_X, $opt_D, $opt_h);
getopts("i:o:c:s:d:f:XDh");
if (!$opt_i && !$opt_o && !$opt_c && !$opt_s && !$opt_d && !$opt_f && !$opt_X && !$opt_D && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}
elsif ($opt_X) {
    create_ctr();
}
elsif ($opt_D) {
    create_default_ctr();
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
my $curr_format = $opt_f;

print STDERR "\n\n========================";
print STDERR "\nSTARTING PIPELINE...\n";
print STDERR "========================\n\n";


## Check if exists outfile, outtree, infile and intree
## In can give problems to programs from Phylip

my @default_filename = ('infile', 'intree', 'outfile', 'outtree');

foreach my $def_file (@default_filename) {
    if (-e $def_file) {
	print STDERR "\nWARNING: $def_file exists into the working dir. It could be a problem with phylip programs.\n\n";
    }
}

my %master_files;
my @input_files;

if ($step > 1) {

    print STDERR "START POINT MOVE TO STEP:$step\n";

    unless (defined $input_dir) {
	die("ARGUMENT ERROR: -s <step> option only can be used if it is specified a dir of input arguments.\n");
    }
    else {

	unless (defined $curr_format) {
	    die("ARGUMENT ERROR: -f <input_format> was not supplied. It must to be supplied when step > 1.\n");
	}
	## It will take the files from the input dir.

	print STDERR "\n0) DATA LOADING\n";

	opendir(my $idir_fh, $input_dir);

	my @current_formats = split(/,/, $curr_format);

	print STDERR "\nInput format: $curr_format.\n";
	
	print STDERR "\nLoading files into the input file variable.\n";

	my @files = readdir($idir_fh);

	my $prev_step = $step - 1;
	my $step_tag = 'step'.$prev_step;
	$master_files{$step_tag} = {};

	foreach my $file (@files) {
	    unless ($file =~ m/^\.+$/) {
		my $f = 0;
		foreach my $format (@current_formats) {
		    $f++;
		    if ($file =~ m/$format$/) {
			if (exists $master_files{$step_tag}->{$format}) {
			    push @{$master_files{$step_tag}->{$format}}, $input_dir.'/'.$file;
			}
			else {
			    $master_files{$step_tag}->{$format} = [$input_dir.'/'.$file];
			}

			## Also it will use as input file the first format
			if ($f == 1) {
			    push @input_files, $input_dir.'/'.$file;
			}
		    }
		}		
	    }
	}

	## Now it will overwrite the input format for the first format, keeping secondary formats for
	## other functions (it will check the only the first input format for change format purposes)

	$curr_format = $current_formats[0];

	my $f_count = scalar(@input_files);

	if ($f_count == 0) {
	    die("\nINPUT ERROR: None file was found in the input dir: $input_dir\n");
	}
	print STDERR "\n$f_count files have been loaded to start the pipeline in the step $step.\n";
    }
}
else {
    ## Check if exists input file.

    unless (defined $maf) {
	die("ARGUMENT ERROR: -i <input_maf_file> was not supplied. It must to be supplied for steps = 1.\n");
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
my %align_exec = %{$ctr_arguments_href[1]};
my %boot_exec = %{$ctr_arguments_href[2]};
my %dist_exec = %{$ctr_arguments_href[3]};
my %tree_exec = %{$ctr_arguments_href[4]};
my %modftree_exec = %{$ctr_arguments_href[5]};
my %constree_exec = %{$ctr_arguments_href[6]};
my %parsetree_exec = %{$ctr_arguments_href[7]};

my %strain_list = ();


unless ($step > 1) {    
    $curr_format = 'maf';
}

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


## Define alway the input_file array
## This array always will be overwrite by the result of a step, but before do it
## the script will store the list into a hash as %files = ( step => { input or output => \@list_of_files ) 



my @ref_input_files;



############################################
## STEP 1 (ALIGNMENT EXTRACTION)          ##
##  Input type  => a single maf file      ##
##  Output type => a list of fasta files  ##
############################################

if ($step == 1) {

    ## Add input to master file
 
    $master_files{'step1'} = { 'input' => [$maf] };
    $master_files{'step1'}->{'in_format'} = 'maf';

    ## Create the preprocessing_folder

    print STDERR "\n1) EXTRACTING ALIGMENTS FROM MAF FILE.\n";
    
    print STDERR "\n  1.1) CREATING PREPROCESSING FOLDER.\n";
    my $curr_path = getcwd();
    my $preprocessing_dir = $curr_path . '/pmaf_1_align_extract';

    mkdir($preprocessing_dir);
    print "\n\tDone ($preprocessing_dir).\n";


    ## Open the input maf file and parse

    print STDERR "\n  1.2) PARSING MAF FILE.\n";

    my ($fasta_files_href, $contig_href, $read_href) = extract_align_from_maf($maf, \%assem_exec, $preprocessing_dir);

    my @fasta_files;
    my $n_fasta_files = 0;
    my %fasta_files = %{$fasta_files_href};
    foreach my $contig_id (sort keys %fasta_files) {
	my %filetypes = %{$fasta_files{$contig_id}};
	push @fasta_files, $filetypes{'assembly_reads'};
	
	foreach my $type (keys %filetypes) {
	    $n_fasta_files++;
	}
    }

    print STDERR "\n\n\tDone... $n_fasta_files have been created in the dir:$preprocessing_dir.\n\n";

    ## Add output to master file
 
    $master_files{'step1'}->{'output'} = \@fasta_files;
    $master_files{'step1'}->{'out_format'} = 'fasta';
    $curr_format = 'fasta';

    @input_files = @fasta_files;

    ## Finally add a new step to the step var

    $step++;
}


############################################################
## STEP 2 (REALIGNMENT)                                   ##
##  Input type  => a list of fasta files                  ##
##  Output type => a list of input_next_step_format files ##
############################################################

if ($step == 2) {

    ## Add input to master file
 
    $master_files{'step2'} = { 'input' => \@input_files };
    $master_files{'step2'}->{'in_format'} = 'fasta';

    ## The program will be runned over the fasta file writting the output in the aligment folder

    print STDERR "\n2) REALIGNMENT PROCESS.\n";

    if (defined $align_exec{'executable'}) {
    
	## It will realign the sequences using an specific program, 
	## so in that case it will create a bunch of new alignments.

	print STDERR "\n  2.1) CREATING REALIGNMENT FOLDER.\n";
	my $alignment_dir = $curr_path . '/pmaf_2_realignment';
    
	mkdir($alignment_dir);
	print "\n\tDone ($alignment_dir).\n";    

	print STDERR "\n  2.2) REALIGNING SEQUENCES.\n";

	my %out_files = run_realignment(\%align_exec, \@input_files, $alignment_dir);

	my @realign_files = @{$out_files{'output'}};
	my $n_realign_files = scalar(@realign_files);
	print STDERR "\n\n\tDone... $n_realign_files have been created in the dir:$alignment_dir.\n";
	
	## Add output to master file
 
	$master_files{'step2'}->{'output'} = \@realign_files;
	$master_files{'step2'}->{'out_format'} = $align_exec{'output_format'};
	$curr_format = $align_exec{'output_format'};

	@input_files = @realign_files;		
    }
    else {
	print STDERR "\n  2.1) CREATING REALIGNMENT FOLDER.\n";
	print STDERR "\n\tNone realigment process have been asigned. Skipping step.\n";
	
	print STDERR "\n  2.2) REALIGNING SEQUENCES.\n";
	print STDERR "\n\tNone realigment process have been asigned. Using the aligment data from maf file.\n\n";
    }

    ## Finally add a new step to the step var
    $step++;
}


########################################
## INTERSTEP: Changing Align format   ##
##  From: $align_exec{output_format}  ##
##  To:   $boots_exec{input_format}   ##
########################################

if (defined $boot_exec{'executable'}) {

    if ($step == 3 && $curr_format ne $boot_exec{'input_format'}) {
    
	my @newformat_files = change_format(\@input_files, $curr_format, $boot_exec{'input_format'}, $step);
	
	@input_files = @newformat_files;
	$curr_format = $boot_exec{'input_format'};
    }
}

## It will get the input files without bootstraping to use 
## as reference.

@ref_input_files = @input_files;


############################################################
## STEP 3 (BOOTSTRAP)                                     ##
##  Input type  => $boots_exec{input_format}              ##
##  Output type => $boots_exec{output_format}             ##
############################################################

## Now it will incluide the option of do bootstraping
## It will need to create the command files

if ($step == 3) {

    ## Add input to master file
 
    $master_files{'step3'} = { 'input' => \@input_files };
    $master_files{'step3'}->{'in_format'} = $boot_exec{'input_format'};

    ## Now it will have different options to run the bootstrap depending of the
    ## program
    
    print STDERR "\n\n3) RUNNING BOOTSTRAP CALCULATION.\n";

    if (defined $boot_exec{'executable'}) {

	print STDERR "\n  3.1) CREATING BOOTSTRAP FOLDER.\n";
	my $bootstrap_dir = $curr_path . '/pmaf_3_bootstrap';
    
	mkdir($bootstrap_dir);
	print "\n\tDone ($bootstrap_dir).\n";    

	print STDERR "\n  3.2) BOOTSTRAP ALIGNMENT.\n";

	my %out_files = run_bootstrap(\%boot_exec, \@input_files, $bootstrap_dir);

	my @bootstrap_files = @{$out_files{'outfile'}};
	my $n_bootstrap_files = scalar(@bootstrap_files);
	print STDERR "\n\n\tDone... $n_bootstrap_files have been created in the dir:$bootstrap_dir.\n";
	
	## Add output to master file
 
	$master_files{'step3'}->{'output'} = \@bootstrap_files;
	$master_files{'step3'}->{'out_format'} = $boot_exec{'output_format'};
	$curr_format = $boot_exec{'output_format'};

	@input_files = @bootstrap_files;
    }
    else {
	print STDERR "\n  3.1) CREATING BOOTSTRAP FOLDER.\n";
	print STDERR "\n\tNone bootstrap process have been asigned. Skipping step.\n";
	
	print STDERR "\n  3.2) BOOTSTRAP ALIGNMENT.\n";
	print STDERR "\n\tNone bootstrap process have been asigned. Skipping step.\n\n";
    }

    ## Finally add a new step to the step var	
    $step++;	
}

########################################
## INTERSTEP: Changing Dist format    ##
##  From: $boots_exec{output_format}  ##
##  To:   $dist_exec{input_format}    ##
########################################

if (defined $dist_exec{'executable'}) {
    if ($step == 4 && $curr_format ne $dist_exec{'input_format'}) {
    
	my @newformat_files = change_format(\@input_files, $curr_format, $dist_exec{'input_format'}, $step);

	@input_files = @newformat_files;
	$curr_format = $dist_exec{'input_format'};
    }
}

############################################################
## STEP 4 (SEQUENCE DISTANCE MATRIX)                      ##
##  Input type  => $dist_exec{input_format}               ##
##  Output type => $dist_exec{output_format}              ##
############################################################


## To run the phylip tools it will need to create command files 
## generally with the name of the file and the options into other 
## lines, and finally Y
## for example:
## input_file.phyl
## D Kimura
## Y

if ($step == 4) {

    ## Add input to master file
 
    $master_files{'step4'} = { 'input' => \@input_files };
    $master_files{'step4'}->{'in_format'} = $dist_exec{'input_format'};

    print STDERR "\n\n4) RUNNING MATRIX DISTANCE CALCULATION.\n";

    if (defined $dist_exec{'executable'}) {
    
	print STDERR "\n  4.1) CREATING DISTANCE MATRIX FOLDER.\n";
	my $matrix_dir = $curr_path . '/pmaf_4_dist_matrix';    
	mkdir($matrix_dir);
	print "\n\tDone ($matrix_dir).\n";   
	
	print STDERR "\n  4.2) CREATING DISTANCE MATRIX.\n";

	my %out_files = run_distmatrix(\%dist_exec, \@input_files, $matrix_dir);

	my @distmatrix_files = @{$out_files{'outfile'}};

	my $n_distmatrix_files = scalar(@distmatrix_files);
	print STDERR "\n\n\tDone... $n_distmatrix_files have been created in the dir:$matrix_dir.\n";
	
	## Add output to master file
 
	$master_files{'step4'}->{'output'} = \@distmatrix_files;
	$master_files{'step4'}->{'out_format'} = $dist_exec{'output_format'};
	$curr_format = $dist_exec{'output_format'};

	@input_files = @distmatrix_files;
    }
    else {
	print STDERR "\n  4.1) CREATING DISTANCE MATRIX FOLDER.\n";
	print STDERR "\n\tNone distance matrix calculation process have been asigned. Skipping step.\n";
	
	print STDERR "\n  4.2) CREATING DISTANCE MATRIX.\n";
	print STDERR "\n\tNone distance matrix calculation process have been asigned. Skipping step.\n\n";
    }
    ## Finally add a new step to the step var	
    $step++;
}

########################################
## INTERSTEP: Changing tree format    ##
##  From: $dist_exec{output_format}   ##
##  To:   $tree_exec{input_format}    ##
########################################

if (defined $tree_exec{'executable'}) {
    if ($step == 5 && $curr_format ne $tree_exec{'input_format'}) {
    
	my @newformat_files = change_format(\@input_files, $curr_format, $tree_exec{'input_format'}, $step);

	@input_files = @newformat_files;
	$curr_format = $dist_exec{'input_format'};
    }
}


############################################################
## STEP 5 (TREE)                                          ##
##  Input type  => $tree_exec{input_format}               ##
##  Output type => $tree_exec{output_format}              ##
############################################################
	
if ($step == 5) {

    ## Add input to master file
 
    $master_files{'step5'} = { 'input' => \@input_files };
    $master_files{'step5'}->{'in_format'} = $tree_exec{'input_format'};

    print STDERR "\n\n5) RUNNING TREE CALCULATION.\n";

    if (defined $tree_exec{'executable'}) {
    
	print STDERR "\n  5.1) CREATING TREE FOLDER.\n";
	my $tree_dir = $curr_path . '/pmaf_5_tree';    
	mkdir($tree_dir);
	print "\n\tDone ($tree_dir).\n";   
	
	print STDERR "\n  5.2) CALCULATING TREES.\n";

	my %out_files = run_treecalc(\%tree_exec, \@input_files, $tree_dir);

	my @tree_files = @{$out_files{'outtree'}};
	my $n_tree_files = scalar(@tree_files);
	print STDERR "\n\n\tDone... $n_tree_files have been created in the dir:$tree_dir.\n";
	
	## Add output to master file
 
	$master_files{'step5'}->{'output'} = \@tree_files;
	$master_files{'step5'}->{'outfile'} = $out_files{'outfile'}; 
	$master_files{'step5'}->{'out_format'} = $tree_exec{'output_format'};
	$curr_format = $tree_exec{'output_format'};

	@input_files = @tree_files;	
    }
    else {
	print STDERR "\n  5.1) CREATING TREE FOLDER.\n";
	print STDERR "\n\tNone tree calculation process have been asigned. Skipping step.\n";
	
	print STDERR "\n  5.2) CALCULATING TREES.\n";
	print STDERR "\n\tNone tree calculation process have been asigned. Skipping step.\n\n";
    }
    ## Finally add a new step to the step var	
    $step++;
}




#########################################
## INTERSTEP: Changing tree format     ##
##  From: $tree_exec{output_format}    ##
##  To:   $modftree_exec{input_format} ##
#########################################

if (defined $modftree_exec{'executable'}) {
    if ($step == 6 && $curr_format ne $modftree_exec{'input_format'}) {
    
	my @newformat_files = change_format(\@input_files, $curr_format, $modftree_exec{'input_format'}, $step);

	@input_files = @newformat_files;
	$curr_format = $modftree_exec{'input_format'};
    }
}


############################################################
## STEP 6 (RETREE)                                        ##
##  Input type  => $retree_exec{input_format}             ##
##  Output type => $retree_exec{output_format}            ##
############################################################

if ($step == 6) {

    ## Add input to master file
 
    $master_files{'step6'} = { 'input' => \@input_files };
    $master_files{'step6'}->{'in_format'} = $modftree_exec{'input_format'};

    print STDERR "\n\n6) RUNNING TREE MODIFICATION.\n";

    if (defined $tree_exec{'executable'}) {
    
	print STDERR "\n  6.1) CREATING TREE MODIFICATION FOLDER.\n";
	my $treemodf_dir = $curr_path . '/pmaf_6_tree_modif';    
	mkdir($treemodf_dir);
	print "\n\tDone ($treemodf_dir).\n";   
	
	print STDERR "\n  6.2) MODIFICATING TREES.\n";

	my %out_files = run_treemodf(\%modftree_exec, \@input_files, $treemodf_dir);

	my @treemodf_files = @{$out_files{'outtree'}};
	my $n_treemodf_files = scalar(@treemodf_files);
	print STDERR "\n\n\tDone... $n_treemodf_files have been created in the dir:$treemodf_dir.\n";
	
	## Add output to master file
 
	$master_files{'step6'}->{'output'} = \@treemodf_files;
	$master_files{'step6'}->{'out_format'} = $modftree_exec{'output_format'};
	$curr_format = $modftree_exec{'output_format'};

	@input_files = @treemodf_files;

    }
    else {
	print STDERR "\n  6.1) CREATING TREE MODIFICATION FOLDER.\n";
	print STDERR "\n\tNone tree calculation process have been asigned. Skipping step.\n";
	
	print STDERR "\n  6.2) MODIFICATING TREES.\n";
	print STDERR "\n\tNone tree calculation process have been asigned. Skipping step.\n\n";
    }

    ## Finally add a new step to the step var
    $step++;
}


#########################################
## INTERSTEP: Changing tree format     ##
##  From: $retree_exec{output_format}  ##
##  To:   $constree_exec{input_format} ##
#########################################

if (defined $constree_exec{'executable'}) {
    if ($step == 7 && $curr_format ne $constree_exec{'input_format'}) {
    
	my @newformat_files = change_format(\@input_files, $curr_format, $constree_exec{'input_format'}, $step);

	@input_files = @newformat_files;
	$curr_format = $constree_exec{'input_format'};
    }
}


############################################################
## STEP 7 (CONSENSUS TREE)                                ##
##  Input type  => $constree_exec{input_format}           ##
##  Output type => $constree_exec{output_format}          ##
############################################################

my @consensus_tree = ();

if ($step == 7) {

    ## Add input to master file
 
    $master_files{'step7'} = { 'input' => \@input_files };
    $master_files{'step7'}->{'in_format'} = $constree_exec{'input_format'};

    print STDERR "\n\n7) RUNNING CONSENSUS TREE FOR BOOTSTRAP METHODS.\n";

    if (defined $tree_exec{'executable'}) {
    
	print STDERR "\n  7.1) CREATING CONSENSUS TREE FOLDER.\n";
	my $treeconsensus_dir = $curr_path . '/pmaf_7_constree';    
	mkdir($treeconsensus_dir);
	print "\n\tDone ($treeconsensus_dir).\n";   
	
	print STDERR "\n  7.2) CALCULATING CONSENSUS TREES.\n";

	my %out_files = run_consensus_tree(\%constree_exec, \@input_files, $treeconsensus_dir);

	my @constree_files = @{$out_files{'outtree'}};
	my @ref_trees = @{$out_files{'reftree'}}; 
	
	my $n_constree_files = scalar(@constree_files);
	print STDERR "\n\n\tDone... $n_constree_files have been created in the dir:$treeconsensus_dir.\n";
	
	## Add output to master file
 
	$master_files{'step7'}->{'output'} = \@constree_files;
	$master_files{'step7'}->{'out_format'} = $constree_exec{'output_format'};
	$curr_format = $constree_exec{'output_format'};

	## In this case oit will not overwrite the input files with these consensus trees
	## because it they haven't distances, they are just bootstrap values

	@consensus_tree = @constree_files;
	@input_files = @ref_trees;
    }
    else {
	print STDERR "\n  7.1) CREATING CONSENSUS TREE FOLDER.\n";
	print STDERR "\n\tNone consensus tree process have been asigned. Skipping step.\n";
	
	print STDERR "\n  7.2) CALCULATING CONSENSUS TREES.\n";
	print STDERR "\n\tNone consensus tree process have been asigned. Skipping step.\n\n";
    }
    ## Finally add a new step to the step var
    $step++;
}





############################################################
## STEP 8 (PARSE TREES)                                   ##
##                                                        ##
##                                                        ##
############################################################


## Finally it will parse the trees and cluster them into categories according the followings parameters:
## 1) Shape.
## 2) Distance between members.
## So following the number of strains it will calculate the number of different possibilities, for example
## for 4 members you can find (without consider the distances) for rooted trees
##
## i)       +-E    ((E1:int1,E2:int2)):intA,(E3:int3,E4:int4):intB);     
##        +-1            
##        | +-E    Comparison between E strains      
##       -2        Comparison between int
##        | +-E
##        +-1  
##          +-E
##                   
##
## ii)      +-E    (((E1:int1,E2:int2)):intA,E3:int3):intB,E4:int4); 
##        +-1
##        | +-E
##      +-2
##     -3 +---E
##      |
##      +-----E
##
##
##  It will open the tree files, replace the variables by strains and ==0 or !=0 and compare how many different types are
##
##  Using the bootstraping, it will take the first tree (it should be the original without bootstrapping, if the option 1 yes, 
##  is used during the boostraping
##




my %tree_types;
my %co_tree;

print STDERR "\n8) PARSING TREES\n\n";

print STDERR "\n  8.1) CREATING RESULTS FOLDER.\n";
my $results_dir = $curr_path . '/pmaf_8_results';    
mkdir($results_dir);
print STDERR "\n\tDone ($results_dir).\n";

my $tree_c = 0;

if ($step == 8) {
    
    my $x = 0;
    
    print STDERR "\n  8.2) PARSING BOOTSTRAPS AND ANALYZING TREES.\n";

    foreach my $tree_filename (@input_files) {
	
	my $treeio = Bio::TreeIO->new(-format => $tree_exec{'output_format'}, -file => $tree_filename);
	
	my $bootstrap_pass = 1;
	my $bootstrap_line;

	if (scalar(@consensus_tree) > 0) {
	    my $boots_treeio = Bio::TreeIO->new( -format => $constree_exec{'output_format'}, 
                                                 -file   => $consensus_tree[$x],
		                                 -internal_node_id => 'bootstrap' );
	    
	    while (my $boot_tree = $boots_treeio->next_tree() ) {

		print STDERR "\tAnalyzing bootstrap values for $consensus_tree[$x]\r";

		my @nodes = $boot_tree->get_nodes();
		    
		foreach my $node (@nodes) {
			
		    ## The bootstrap value will be added to the tree as branch value
		    my $node_bl = $node->branch_length();
		    
		    if (defined $node_bl && $node_bl <= $boot_exec{'cutoff'}) {
			$bootstrap_pass = 0;
			$bootstrap_line = $boot_tree->as_text('newick');
			chomp($bootstrap_line);
		    }
		}
	    } 
	}
	    
	my $t = 0;
	
	## Now it will analyze the results

	if ($bootstrap_pass == 1) {
	    while( my $tree = $treeio->next_tree ) {
	    		

		## Store the tree asociated to a contig
		
		my $tree_string = $tree->as_text('newick');		
	   
		chomp($tree_string);	

		## store only the contig name
		my $contig_id = $tree_filename;
		if ($tree_filename =~ m/ref_(.+?)\.\w+/) {
		    $contig_id = $1;
		}

		$co_tree{$contig_id} = $tree_string;

		## It will change a couple of things inside the nodes.	

		my @tr_nodes = $tree->get_nodes();

		foreach my $tr_node (@tr_nodes) {
		    
		    ## 1) Modify the id, removing the integer (to get only the strain)

		    my $id = $tr_node->id();
		    $id =~ s/-\d+$//;
		    $tr_node->id($id);

		    ## 2) Change the branch by tags according the value.

		    my $branch = $tr_node->branch_length();
		    my $match = 0;
		    
		    my %branch_tags = %{$parsetree_exec{'branch_tags'}};

		    if (defined $branch) {
			foreach my $branch_value ( sort {$a <=> $b} keys %branch_tags ) {
			    my $branch_tag = $branch_tags{$branch_value};

			    my $branch_v_obj = Math::BigFloat->new($branch_value);
			    my $branch_obj = Math::BigFloat->new($branch);

			    my $test1 = $branch_v_obj->bcmp($branch_obj);
			    my $test2 = $branch_obj->bcmp($branch_v_obj);

			    if ($branch_obj->bcmp($branch_v_obj) <= 0 && $match == 0) {

				## The branch values will be order asc. It only will replace the first match
				## for example, if there are 3 tags 0.0=>EQ, 0.1=>SI, 1.0=>DI
				## and it is checking the following values 0.8, 0.3, 0.001 and 0.0
				## 0.8 and 0.3 will pass only the last tag (DI)
				## 0.001 will pass SI and DI tags but SI is before DI, so it will tag as SI
				## 0.0 will pass EQ, SI and DI, but the first is EQ.

				$tr_node->branch_length($branch_tag);
				$match = 1;
			    }		    
			}
		    }
		}

		## Now it need order the nodes
		
		my $tree_as_string = $tree->as_text('newick');		
	   
		chomp($tree_as_string);	
	
		my @data = split(/,\(/, $tree_as_string);
		my @ord_b_data;
		foreach my $b_data (@data) {
		    my @bb_data = split(/\)/, $b_data);
		    my @ord_bb_data;
		    foreach my $bb_data (@bb_data) {
			my @bbb_data = split(/\(/, $bb_data);
			my @ord_bbb_data;
			foreach my $bbb_data (@bbb_data) {
			    my @bbbb_data = split(/,/, $bbb_data);						    			
			    my $ord_data = join(',', sort @bbbb_data);		  
			    push @ord_bbb_data, $ord_data;
			}
			my $ord_bb_data = join('(', @ord_bbb_data);
			push @ord_bb_data, $ord_bb_data;
		    }
		    my $ord_b_data = join(')', @ord_bb_data);
		    push @ord_b_data, $ord_b_data;
		}
		$tree_as_string = join(',(', @ord_b_data);

		
		## store only the contig name
		my $recontig_id = $tree_filename;
		if ($tree_filename =~ m/ref_(.+?)\.\w+/) {
		    $recontig_id = $1;
		}

		## Add another tree to the tree count

		$tree_c++;

		if (exists $tree_types{$tree_as_string}) {
		    push @{$tree_types{$tree_as_string}}, $recontig_id;
		}
		else {
		    $tree_types{$tree_as_string} = [$recontig_id];
		}
	    }
	}
	else {
	    my $b_tree_name = basename($tree_filename);
	    print STDERR "\n\t$b_tree_name have not passed bootstrat filter (value:$boot_exec{'cutoff'})\n\t(tree:$bootstrap_line)\n";
	}   
	$x++;
    }
    
    ## Finally print the tree types.

    my $treetypes_file = $results_dir .'/tree_analysis_result.tab';
    open(my $ofh, '>', $treetypes_file);

    my $treesum_file = $results_dir .'/tree_analysis_summary.tab';
    open(my $sufh, '>', $treesum_file);

    
    my $tree_types_n = scalar(keys %tree_types);

    print STDERR "\n\nThere are $tree_types_n tree types ($tree_c trees):\n";

    my $t = 0;
   
    foreach my $treetype (sort keys %tree_types) {

	my $cat = "type_" . $t;

	my @contigs = @{$tree_types{$treetype}};
	my $tree__n = scalar(@contigs);
	my $f_cat = sprintf '%10s', $cat; 
	my $f_tree_n = sprintf '%5s', $tree__n;

	## To calculate the porcentage

	my $perc = ($tree__n/$tree_c)*100;
	my $perc_obj = Math::BigFloat->new($perc);
	my $f_perc = sprintf '%10s', $perc_obj->bfround(-3);

	my $f_treetype = sprintf '%10s', $treetype;
	print STDERR "$f_cat\t$f_tree_n\t$f_perc%\t$f_treetype\n";
	print $sufh "$f_cat\t$f_tree_n\t$f_perc%\t$f_treetype\n";

	foreach my $co (@contigs) {
	    print $ofh "$co\t$cat\t$co_tree{$co}\n";
	}
	
	$t++;
    }

    ## One last function, if dnaml was used as tree_calc{'executable'} it will take
    ## the ln likehood values

    if (defined  $tree_exec{'executable'} &&  $tree_exec{'executable'} =~ m/dnaml/) {
	extract_ln_likehood($master_files{'step5'}->{'outfile'}, $results_dir);
    }

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

      This script is a pipeline that calculate the phylotree for contigs and
    the reads of these contigs from a MIRA assembly. It will create at least
    file per contig.

    It have the following steps:

      0- Parse controler file
         Load input files (maf or step files)
 
      1- Extract the contigs and the reads from maf file

      2- Rerun alignments.
         (Programs like clustalw or muscle will be used as executables)

      3- Apply bootstrapping or delete-half-jackknifing process to the 
         alignments.
         (Programs like bootseq will be used)

      4- Calculate distance matrix based in the alignments
         (Programs like dnadist or protdist will be used as executables)

      5- Calculate the tree
         (Programs like neighbor will be used as executables)
    
      6- Modify the tree (for example reroot it)
         (Programs like retree could be used)

      7- Calculate the consensus tree and the bootstrap values
         (Programs like consense will be used as executables)

      8- Parse the results (trees), filter and calculate the
         tree types.

    This pipeline can start from any step supplying the input files.
    In that case the controler need to have the output_format parameter
    set for the previous step. 

    To create a control file: 
      contig_phyloanalysis.pl -X

    To create a default control file
      contig_phyloanalysis.pl -D

    This default control file will use the following programs and arguments.

      1) Extract maf file with 1 seq/strain and all the strains 
         and minimum overlapping length of 100 bp.

      2) Skip the realignment

      3) Bootstrap with seqboot program and 1000 replicates.
         Bootstrap cutoff value for final tree parsing: 500

      4) Distance matrix calculated with dnadist and Kimura algorithm.
       
      5) Tree calculated with neighbor

      6) Retree the result to root the trees using the middle point
         with retree program

      7) Calculate the consensus tree and use the following tags for
         the branch length: 0.001=>SIM (similar) and 1.0 => DIF
         (different)      
         
    Usage:

      phylotypes_pipeline.pl [-h] -i <maf_file> -c <control_file> 
                              [-s <step>] [-d <input_dir>] [-f <input_format>] 
                              [-o <output_basename>] [-X] [-D]    

    Example:

      phylotypes_pipeline.pl [-h] -i coffea.maf -c maf_to_phylo.default.ctr
       
    Flags:

      -i <maf_file>               assembly file in maf format format (mandatory)
      -o <output_basename>        output basename (input_file.output by default)
      -c <control_file>           control file (mandatory)
      -s <step>                   step to start the pipeline (1 by default)
      -d <input_dir>              input dir to start in an step > 1
      -f <input_format>           input format for files, when starting step > 1
                                  (more than one format can be specified, it will
                                   match the formats, the first format must to be
                                   the input format, secondary formats will be used 
                                   for other tasks)
      -X <create_control_file>    create the control file: maf_to_phylo.ctr
      -D <default_control_file>   create control file: maf_to_phylo.default.ctr
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
    my $ctr_filename = shift || 'maf_to_phylo.ctr';
 
    open my $ctr_fh, '>', $ctr_filename 
	|| die("OPEN FILE ERROR: $ctr_filename file can not be openned (system error: $!),\n");

    
    print $ctr_fh <<EOF;


############################################################
## This pipeline is composed by four steps:
##
## 1) Extract the sequence data from the assembly file
##
## 2) Realign the sequences from the assembly file
##    (if is needed)
##
## 3) Create the similarity matrix 
##    (if is needed)
##
## 4) Compute the tree
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
### 2) SEQUENCE REALIGNMENT                              ###
############################################################

## Note: Many realignment tools can be used, like muscle or
##       clustaw. There are five variables that the pipeline 
##       needs to run the sequence realignment. 
##       * executable: path for the sequence alignment 
##                     program
##       * output_format: output format produced by the 
##                        aligment program
##       * options: used by the program
##       * input argument: 
##       * output_argument:

    
<REALIGNMENT_EXEC>="Aligment tool executable path" 
## Requeriment: none by default 
## Example:     <REALIGMENT_EXEC>="~/programs/muscle"
## Note:        Path to access to executable
##              To skip this step, leave undef


<REALIGNMENT_OFORMAT>="Output format for the alignment" 
## Requeriment: none by default 
## Example:     <REALIGMENT_OFORMAT>="clustalw"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. The alignment tool 
##              usually produce a concrete format, so the 
##              pipeline will change this format using 
##              bioperl. The options are: bl2seq, clustalw, 
##              emboss, fasta, maf, mase, mega, meme, 
##              metafasta, msf, nexus, pfam, phylip, po, 
##              prodom, psi, selex, stockholm, XMFA and arp. 

<REALIGNMENT_OPT>="List of options except input and output"
## Requeriment: none
## Example:     <ALIGMENT_OPT>="-clw -maxhours 12"
## Note:        Input and output tags should be detailed as 
##              <ALIGNMENT_IN> and <ALIGNMENT_OUT>
 
<REALIGNMENT_IN>="Input argument"
## Requeriment: none by default
## Example:     <ALIGMENT_IN>="-in"
## Note:        For STDIN leave empty
##              For command file use "{cmd}"

<REALIGNMENT_OUT>="Output argument"
## Requeriment: none by default
## Example:     <ALIGMENT_OUT>="-out"
## Note:        For STDOUT leave empty
##              For predefined output file by the program
##              use '{predefined=name}'


############################################################
### 3) BOOTSTRAP SAMPLE CALCULATION                      ###
############################################################

<BOOTSTRAP_EXEC>="Bootstrap program path"
## Requeriment: dnadist by default
## Example:     <BOOTSTRAP_EXEC>="/home/aure/seqboot"
## Note:        

<BOOTSTRAP_IFORMAT>="Input alignment format used" 
## Requeriment: none by default 
## Example:     <BOOTSTRAP_IFORMAT>="phylip"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. The alignment tool
##              usually produce a concrete format, so the 
##              pipeline will change this format using 
##              bioperl. The options are: bl2seq, clustalw, 
##              emboss, fasta, maf, mase, mega, meme, 
##              metafasta, msf, nexus, pfam, phylip, po, 
##              prodom, psi, selex, stockholm, XMFA and arp. 

<BOOTSTRAP_OFORMAT>="Output alignment format used" 
## Requeriment: none by default 
## Example:     <BOOTSTRAP_OFORMAT>="phylip"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. The alignment tool
##              usually produce a concrete format, so the 
##              pipeline will change this format using 
##              bioperl. The options are: bl2seq, clustalw, 
##              emboss, fasta, maf, mase, mega, meme, 
##              metafasta, msf, nexus, pfam, phylip, po, 
##              prodom, psi, selex, stockholm, XMFA and arp. 

<BOOTSTRAP_OPT>="Bootstrap program options"
## Requeriment: Optional**
## Example:     <BOOTSTRAP_OPT>="1 Yes"
## Note:        Separate the options by commas

<BOOTSTRAP_REP>="Bootstrap replicates"
## Requeriment: Optional**
## Example:     <BOOTSTRAP_REP>="1000"
## Note:        Overwrite the R option in the seqboot program
##              It is 100 by default
##              Also it will replace the M options in other
##              executables

    
<BOOTSTRAP_IN>="Input argument"
## Requeriment: Optional
## Example:     <BOOTSTRAP_IN>="{cmd}"
## Note:        For STDIN leave empty
##              For command file use "{cmd}"

<BOOTSTRAP_OUT>="Output argument"
## Requeriment: Optional
## Example:     <BOOTSTRAP_OUT>='{predefined=outfile}'
## Note:        For STDOUT leave empty
##              For predefined output file by the program
##              use '{predefined=name}'

<BOOTSTRAP_CUTOFF>="Cutoff value to discard a tree"
## Requeriment: Optional
## Example:     <BOOTSTRAP_OUT>='{predefined=outfile}'
## Note:        For STDOUT leave empty
##              For predefined output file by the program
##              use '{predefined=name}'


############################################################
### 4) MATRIX DISTANCE CALCULATION                       ###
############################################################

<DISTANCE_EXEC>="Distance matrix program path"
## Requeriment: dnadist by default
## Example:     <DISTANCE_EXEC>="/home/aure/dnadist"
## Note:        For phylogenetic analysis programs that 
##              doesn't need distance matrix calculation
##              it should be undef.

<DISTANCE_IFORMAT>="Input alignment format used" 
## Requeriment: none by default 
## Example:     <DISTANCE_IFORMAT>="phylip"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. The alignment tool
##              usually produce a concrete format, so the 
##              pipeline will change this format using 
##              bioperl. The options are: bl2seq, clustalw, 
##              emboss, fasta, maf, mase, mega, meme, 
##              metafasta, msf, nexus, pfam, phylip, po, 
##              prodom, psi, selex, stockholm, XMFA and arp. 

<DISTANCE_OPT>="Distance matrix program path"
## Requeriment: Optional**
## Example:     <DISTANCE_OPT>="D Kimura,2 No"
## Note:        Separate the options by commas
    
<DISTANCE_IN>="Input argument"
## Requeriment: Optional
## Example:     <DISTANCE_IN>="{cmd}"
## Note:        For STDIN leave empty
##              For command file use "{cmd}"

<DISTANCE_OUT>="Output argument"
## Requeriment: Optional
## Example:     <DISTANCE_OUT>='{predefined=outfile}'
## Note:        For STDOUT leave empty
##              For predefined output file by the program
##              use '{predefined=name}'


############################################################
### 5) TREE CALCULATION                                  ###
############################################################

      
<PHYLOTREE_EXEC>="Phylogenetic analysis program path"
## Requeriment: neighbor by default
## Example:     <PHYLOTREE_EXEC>="/home/aure/neighbor"
## Note:        

<PHYLOTREE_IFORMAT>="Input alignment format used" 
## Requeriment: none by default 
## Example:     <PHYLOTREE_IFORMAT>="phylip"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. The alignment tool
##              usually produce a concrete format, so the 
##              pipeline will change this format using 
##              bioperl. The options are: bl2seq, clustalw, 
##              emboss, fasta, maf, mase, mega, meme, 
##              metafasta, msf, nexus, pfam, phylip, po, 
##              prodom, psi, selex, stockholm, XMFA and arp. 

<PHYLOTREE_OFORMAT>="Output tree format used" 
## Requeriment: none by default 
## Example:     <PHYLOTREE_OFORMAT>="newick"
## Note:        Phylogenetic tools like phylip poduce a 
##              specific tree format but other programs produce
##              others like nexus.
  
<PHYLOTREE_OPT>="Phylogenetic analysis program path"
## Requeriment: Optional
## Example:     <PHYLOTREE_OPT>="M Yes,1000,1001"
## Note:        When the input is an alignment distance tool
##              should be empty
##              For M multiple option in neighbor or other
##              phylip programs, it will need after the switch
##              the multiple amount number and the seed (odd)
    
<PHYLOTREE_IN>="Input argument"
## Requeriment: Optional
## Example:     <PHYLOTREE_IN>="{cmd}"
## Note:        For STDIN leave empty
##              For command file use "{cmd}"

<PHYLOTREE_OUT>="Output argument"
## Requeriment: Optional
## Example:     <PHYLOTREE_OUT>="{predefined=outfile}"
## Note:        For STDOUT leave empty
##              For predefined output file by the program
##              use '{predefined=name}'
  

############################################################
### 6) TREE MODIFICATION                                 ###
############################################################

      
<RETREE_EXEC>="Phylogenetic analysis program path"
## Requeriment: neighbor by default
## Example:     <RETREE_EXEC>="/home/aure/retree"
## Note:        

<RETREE_IFORMAT>="Input alignment format used" 
## Requeriment: none by default 
## Example:     <RETREE_IFORMAT>="newick"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. 

<RETREE_OFORMAT>="Output tree format used" 
## Requeriment: none by default 
## Example:     <RETREE_OFORMAT>="newick"
## Note:        Phylogenetic tools like phylip poduce a 
##              specific tree format but other programs produce
##              others like nexus.
  
<RETREE_GEN_OPT>="General options for the retree program"
## Requeriment: Optional
## Example:     <RETREE_GEN_OPT>=""
## Note:        General options of the program
    
<RETREE_TREE_OPT>="Specific options for the retree program"
## Requeriment: Optional
## Example:     <RETREE_TREE_OPT>="M"
## Note:        Specific options for each tree during the
##              program running

<RETREE_IN>="Input argument"
## Requeriment: Optional
## Example:     <RETREE_IN>="{cmd}"
## Note:        For STDIN leave empty
##              For command file use "{cmd}"

<RETREE_OUT>="Output argument"
## Requeriment: Optional
## Example:     <RETREE_OUT>="{predefined=outtree}"
## Note:        For STDOUT leave empty
##              For predefined output file by the program
##              use '{predefined=name}'


############################################################
### 7) TREE CONSENSUS                                    ###
############################################################

      
<CONSTREE_EXEC>="Phylogenetic analysis program path"
## Requeriment: neighbor by default
## Example:     <CONSTREE_EXEC>="/home/aure/consensus"
## Note:        

<CONSTREE_IFORMAT>="Input alignment format used" 
## Requeriment: none by default 
## Example:     <CONSTREE_IFORMAT>="newick"
## Note:        Phylogenetic tools like phylip use a 
##              specific alignment format. 

<CONSTREE_OFORMAT>="Output tree format used" 
## Requeriment: none by default 
## Example:     <CONSTREE_OFORMAT>="newick"
## Note:        Phylogenetic tools like phylip poduce a 
##              specific tree format but other programs produce
##              others like nexus.
  
<CONSTREE_OPT>="General options for the consensus program"
## Requeriment: Optional
## Example:     <CONSTREE_OPT>="R"
## Note:        General options of the program
    

<CONSTREE_IN>="Input argument"
## Requeriment: Optional
## Example:     <CONSTREE_IN>="{cmd}"
## Note:        For STDIN leave empty
##              For command file use "{cmd}"

<CONSTREE_OUT>="Output argument"
## Requeriment: Optional
## Example:     <CONSTREE_OUT>="{predefined=outfile,outtree}"
## Note:        For STDOUT leave empty
##              For predefined output file by the program
##              use '{predefined=name}'


############################################################
### 8) PARSE TREE                                        ###
############################################################

<BRANCH_LENGTH_TAGS>="values follow by an arrow and tag"
## Requeriment: Optional (0.0=>EQ,1.0=>NE by default)
## Example:     <BRANCH_LENGTH_TAGS>="0.0=>EQ,1.0=>NE";
## Note:        The branch lengths will be replaced by 
##              the tags will less or equal to the value.

EOF

print STDERR "\nControl file have been created ($ctr_filename)\n\n";

exit (1);

}

=head2 parse_ctr_file

  Usage: my @ctr_arguments_href = parse_ctr_file($ctr_file);

  Desc: parse control file

  Ret: An array of hash references:
       $ctr_args[0] = $assembly_opts_href
       $ctr_args[1] = $realignment_args_href
       $ctr_args[2] = $boots_args_href
       $ctr_args[3] = $mat_distance_args_href
       $ctr_args[4] = $phylotree_arg_href
       $ctr_args[5] = $retree_args_href
       $ctr_args[6] = $constree_args_href
       $ctr_args[7] = $parsetree_args_href

  Args: $ctr_file, a scalar with the filename

  Side_Effects: Add default values
                Die if a mandatory argument is not specified

  Example: my @ctr_arguments_href = parse_ctr_file($ctr_file);

=cut

sub parse_ctr_file {    
    
    my $ctrfile = shift 
	|| die("ARGUMENT ERROR: None control file name was supplied to parse_ctr_file function.\n");

    ## Define the control variables with the default values

    my %assembly_args = ( 'min_alignment' => 100 );
    my %align_args = ();
    my %boot_args = ();
    my %dist_args = ();
    my %tree_args = ();
    my %retree_args = ();
    my %constree_args =();
    my %parsetree_args =();

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
	    
	    ## REALIGNMENT ARGUMENTS
	    	    
	    if ($_ =~ m/<REALIGNMENT_EXEC>="(.+)"/) {
		$align_args{'executable'} = $1; 
	    }
	    if ($_ =~ m/<REALIGNMENT_OFORMAT>="(.+)"/) {
		$align_args{'output_format'} = $1; 
	    }
	    if ($_ =~ m/<REALIGNMENT_OPT>="(.+)"/) {
		$align_args{'options'} = $1; 
	    }
	    if ($_ =~ m/<REALIGNMENT_IN>="(.+)"/) {
		$align_args{'input'} = $1; 
	    }
	    if ($_ =~ m/<REALIGNMENT_OUT>="(.+)"/) {
		$align_args{'output'} = $1; 
	    }

	    ## BOOTSTRAP ARGUMENTS BOOTSTRAP

	     if ($_ =~ m/<BOOTSTRAP_EXEC>="(.+)"/) {
		$boot_args{'executable'} = $1; 
	    }
	    if ($_ =~ m/<BOOTSTRAP_IFORMAT>="(.+)"/) {
		$boot_args{'input_format'} = $1;
	    }
	    if ($_ =~ m/<BOOTSTRAP_OFORMAT>="(.+)"/) {
		$boot_args{'output_format'} = $1;
	    }
	    if ($_ =~ m/<BOOTSTRAP_OPT>="(.+)"/) {
		$boot_args{'options'} = $1; 
	    }
	    if ($_ =~ m/<BOOTSTRAP_REP>="(.+)"/) {
		$boot_args{'boots_replicates'} = $1; 
	    }
	    if ($_ =~ m/<BOOTSTRAP_IN>="(.+)"/) {
		$boot_args{'input'} = $1; 
	    }
	    if ($_ =~ m/<BOOTSTRAP_OUT>="(.+)"/) {
		$boot_args{'output'} = $1; 
	    }
	     if ($_ =~ m/<BOOTSTRAP_CUTOFF>="(.+)"/) {
		$boot_args{'cutoff'} = $1; 
	    }

	    ## MATRIX DISTANCE ARGUMENTS

	    if ($_ =~ m/<DISTANCE_EXEC>="(.+)"/) {
		$dist_args{'executable'} = $1; 
	    }
	    if ($_ =~ m/<DISTANCE_IFORMAT>="(.+)"/) {
		$dist_args{'input_format'} = $1;
	    }
	    if ($_ =~ m/<DISTANCE_OFORMAT>="(.+)"/) {
		$dist_args{'output_format'} = $1;
	    }
	    if ($_ =~ m/<DISTANCE_OPT>="(.+)"/) {
		$dist_args{'options'} = $1; 
	    }
	    if ($_ =~ m/<DISTANCE_IN>="(.+)"/) {
		$dist_args{'input'} = $1; 
	    }
	    if ($_ =~ m/<DISTANCE_OUT>="(.+)"/) {
		$dist_args{'output'} = $1; 
	    }
	    

	    ## PHYLOTREE ARGUMENTS

	    if ($_ =~ m/<PHYLOTREE_EXEC>="(.+)"/) {
		$tree_args{'executable'} = $1; 
	    }	
	    if ($_ =~ m/<PHYLOTREE_IFORMAT>="(.+)"/) {
		$tree_args{'input_format'} = $1;
	    }
	    if ($_ =~ m/<PHYLOTREE_OFORMAT>="(.+)"/) {
		$tree_args{'output_format'} = $1;
	    }
	    if ($_ =~ m/<PHYLOTREE_OPT>="(.+)"/) {
		$tree_args{'options'} = $1; 
	    }		
	    if ($_ =~ m/<PHYLOTREE_IN>="(.+)"/) {
		$tree_args{'input'} = $1; 
	    }
	    if ($_ =~ m/<PHYLOTREE_OUT>="(.+)"/) {
		$tree_args{'output'} = $1; 
	    }

	    ## REETREE ARGUMENTS

	    if ($_ =~ m/<RETREE_EXEC>="(.+)"/) {
		$retree_args{'executable'} = $1;
	    }
	    if ($_ =~ m/<RETREE_IFORMAT>="(.+)"/) {
		$retree_args{'input_format'} = $1;
	    }
	    if ($_ =~ m/<RETREE_OFORMAT>="(.+)"/) {
		$retree_args{'output_format'} = $1;
	    }
	    if ($_ =~ m/<RETREE_GEN_OPT>="(.+)"/) {
		$retree_args{'general_options'} = $1;
	    }
	    if ($_ =~ m/<RETREE_TREE_OPT>="(.+)"/) {
		$retree_args{'retree_options'} = $1;
	    }
	    if ($_ =~ m/<RETREE_IN>="(.+)"/) {
		$retree_args{'input'} = $1;
	    }
	    if ($_ =~ m/<RETREE_OUT>="(.+)"/) {
		$retree_args{'output'} = $1;
	    }

	     ## CONSTREE ARGUMENTS

	    if ($_ =~ m/<CONSTREE_EXEC>="(.+)"/) {
		$constree_args{'executable'} = $1;
	    }
	    if ($_ =~ m/<CONSTREE_IFORMAT>="(.+)"/) {
		$constree_args{'input_format'} = $1;
	    }
	    if ($_ =~ m/<CONSTREE_OFORMAT>="(.+)"/) {
		$constree_args{'output_format'} = $1;
	    }
	    if ($_ =~ m/<CONSTREE_OPT>="(.+)"/) {
		$constree_args{'options'} = $1;
	    }
	    if ($_ =~ m/<CONSTREE_IN>="(.+)"/) {
		$constree_args{'input'} = $1;
	    }
	    if ($_ =~ m/<CONSTREE_OUT>="(.+)"/) {
		$constree_args{'output'} = $1;
	    }

	    ## PARSE TREE ARGS

	    if ($_ =~ m/<BRANCH_LENGTH_TAGS>="(.+)"/) {
		my @branch_tags_pair = split(/,/, $1);

		my %branch_tags;
		foreach my $tag (@branch_tags_pair) {
		    if ($tag =~ m/^(\d\.?\d*)=>(.+)$/) {
			$branch_tags{$1} = $2;
		    }
		}
		$parsetree_args{'branch_tags'} = \%branch_tags;
	    }
	}
    }

    ## Now it will overwrite the M and R option based in the number of replicates

    if (defined $boot_args{'boots_replicates'}) {
	my $bootsrep = $boot_args{'boots_replicates'};
	my $bootsrep_i = $bootsrep+1;

	## First overwrite bootstrap

	if (defined $boot_args{'options'}) {

	    if ($boot_args{'options'} =~ m/R\s+(\d+)/) {
		print STDERR "\n\tMESSAGE: R $1 option for bootstrap step will be replaced by R $bootsrep according replicates option.\n";
		$boot_args{'options'} =~ s/R\s+(\d+)/R $bootsrep/;
	    }
	    else {
		$boot_args{'options'} .= ",R $bootsrep";
		print STDERR "\n\tMESSAGE: R $bootsrep option for bootstrap step will be added according replicates option.\n";
	    }
	}
	else {
	    $boot_args{'options'} = "R $bootsrep";
	    print STDERR "\n\tMESSAGE: R $bootsrep option for bootstrap step will be added according replicates option.\n";
	}

	## Second overwrite distmatrix and add a new one (bootstrap + reference)

	if (defined $dist_args{'options'}) {

	    if ($dist_args{'options'} =~ m/M\s+Yes\s*,\s*D\s*,\s*(\d+)/) {
		print STDERR "\n\tMESSAGE: M $1 option for distance step will be replaced by M $bootsrep according replicates option.\n";
		$dist_args{'options'} =~ s/M\s+Yes\s*,\s*D\s*,\s*(\d+)/M Yes,D,$bootsrep_i/;
	    }
	    else {
		$dist_args{'options'} .= ",M Yes,D,$bootsrep_i";
		print STDERR "\n\tMESSAGE: R $bootsrep option for distance step will be added according replicates option.\n";
	    }
	}
	else {
	    $dist_args{'options'} = "M Yes,D,$bootsrep_i";
	    print STDERR "\n\tMESSAGE: R $bootsrep option for distance step will be added according replicates option.\n";
	}

	## Third overwrite tree calc and add a new one (bootstrap + reference)

	if (defined $tree_args{'options'}) {

	    if ($tree_args{'options'} =~ m/M\s+Yes\s*,\s*(\d+),*\s*\d*/i) {
		print STDERR "\n\tMESSAGE: M $1 option for tree step will be replaced by M $bootsrep according replicates option.\n";

		## It will use as seed value 1807.
		$tree_args{'options'} =~ s/M\s+Yes\s*,\s*(\d+),*\s*\d*/M Yes,$bootsrep_i,1807/;
	    }
	    else {
		if ($tree_args{'executable'} =~ m/dnaml/) {
		    $tree_args{'options'} .= ",M Yes,D,$bootsrep_i,1807,1";
		}
		else {
		    $tree_args{'options'} .= ",M Yes,$bootsrep_i,1807";
		}
		print STDERR "\n\tMESSAGE: R $bootsrep option for tree step will be added according replicates option.\n";
	    }
	}
	else {
	    if ($tree_args{'executable'} =~ m/dnaml/) {
		$tree_args{'options'} .= "M Yes,D,$bootsrep_i,1807,1";
	    }
	    else {
		$tree_args{'options'} = "M Yes,$bootsrep_i,1807";
	    }
	    print STDERR "\n\tMESSAGE: R $bootsrep option for tree step will be added according replicates option.\n";
	}
    }

    ## Default value for branch tags

    unless (defined $parsetree_args{'branch_tags'}) {
	$parsetree_args{'branch_tags'} = { '0.0' => 'EQ', '1.0' => 'NE' };
    }

    return (\%assembly_args, \%align_args, \%boot_args, \%dist_args, \%tree_args, \%retree_args, \%constree_args, \%parsetree_args);
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

  Usage: my ($align_files_href, $contig_href, $read_href) = extract_align_from_maf($maf_file, \%align_args, $preprocessing_dir);

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




=head2 run_realignment

  Usage: my %out_files = run_realignment(\%align_args, \@input_files, $out_dirname);

  Desc: Run the realignment programs over an array of fasta sequences.
        It return the sequences in the format spscified in the arguments

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%align_args, a hash reference with keys=argument_type and value=value
        \@input_file, an array of input filenames
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_realignment(\%align_args, \@input_files, $out_dirname);

=cut

sub run_realignment {
    my $align_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None alignment_argument hash reference have been supplied to run_realignment function.\n");
    
    my $input_files_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file array reference have been supplied to run_realignment function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_realignment function.\n");

    my @input_files = @{$input_files_aref};

    my %outfiles;

    foreach my $fasta (@input_files) {

	## First, redefine dirs and outnames

	my $in_basename = basename($fasta);

	## It will replace the fasta extension by the output format

	my $outformat = $align_args_href->{'output_format'} || '.dat';
	
	my $out_basename = $in_basename;
	$out_basename =~ s/\.\w+$/$out_basename/;	
	my $align_output = $output_dir . '/' . $out_basename;
    
	## Now, build the run command as executable + options + input + output

	my $command = $align_exec{'executable'} . " " . $align_exec{'options'};

	if (defined $align_exec{'input'}) {
	    $command .= ' ' . $align_exec{'input'} . ' ' . $fasta; 	
	}
	else {
	    $command .= ' ' . $fasta; 
	}	
	
	if (defined $align_exec{'output'} && $align_args_href->{'output'} !~ m/\{predefined/) {
	    $command .= ' ' . $align_args_href->{'output'} . ' ' . $align_output; 	
	}
	elsif (defined $align_exec{'output'}) {
	    ## If it is defined an output as predefined, it will not add nothing to the
	    ## command but it will change the names after
	}
	else {
	    $command .= ' >' . $align_output; 
	}
	
	## Print the command and run it as system. If something fail, with autodie should be printed a message

	print STDERR "\n\tRUNNING ALIGMENT: $command\n";
	system($command);		

	## I more than one output is produced, it will add each of them to the hash
	## They will be defined as {tag=list comma separated}

	if (defined $align_exec{'output'} && $align_exec{'output'} =~ m/\{predefined=(.+)\}/) {
	    my @list = split(/,/, $1);
	    
	    foreach my $filetype (@list) {
		## Change the name adding the filetype before the extension

		my $new_filename = $align_output;
		$new_filename =~ s/\.(\w+)$/$filetype\.$1/; 
		rename($filetype, $new_filename);
		
		if (exists $outfiles{$filetype}) {
		    push @{$outfiles{$filetype}}, $new_filename;
		}
		else {
		    $outfiles{$filetype} = [$new_filename];
		}
	    }
	}
	else {
	    ## By default it will store as filetype=output
	    
	    if (exists $outfiles{'output'}) {
		push @{$outfiles{'output'}}, $align_output;
	    }
	    else {
		$outfiles{'output'} = [$align_output];
	    }
	}
    }
    return %outfiles;
}

=head2 change_format

  Usage: my @newformat_files = change_format($input_aref, $curr_format, $next_format, $step);

  Desc: Change the format of an array of files from curr_format to next format
        Alignment formats supported are: 'fasta', 'mase', 'selex', 'clustalw', 
	  'msf', 'phylip', 'po', 'stockholm', 'XMFA' and 'metafasta'.   
        Tree formats supported are: 'netwick', 'nexus', 'nhx', 'lintree', 
          'treecluster', 'pag', 'tab', and 'svggraph'.

  Ret: @newformat_files, an array with the filenames of the new format files

  Args: $input_aref, an array reference with input filenames
        $curr_format, format of these input filenames
        $next_format, format of the output files
        $step, step of the pipeline to create a different dirnames

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my @newformat_files = change_format($input_aref, $curr_format, $next_format, $step);

=cut

sub change_format {

    my $input_files_aref = shift ||
	die("ARGUMENT ERROR: None input files array reference variable was supplied to change_format function.\n");

    my $curr_format  = shift ||
	die("ARGUMENT ERROR: None curr_format variable was supplied to change_format function.\n");

    my $next_format = shift ||
	die("ARGUMENT ERROR: None next_format variable was supplied to change_format function.\n");

    my $step = shift ||
	die("ARGUMENT ERROR: None step variable was supplied to change_format function.\n");

    my @input_files = @{$input_files_aref};

    my @output_files;

    ## Now it will change the alignment format to phylip (or the format specified in the next step)  
    ## and also it will create the phylip option files

    print STDERR "\n  *) REFORMATING FILES TO NEXT INPUT FORMAT.\n";

    ## Previous format will be fasta (extracted from maf file) except if the realignment
    ## produce a new one    

    print STDERR "\n    *.i) CREATING REFORMAT FOLDER.\n";
    my $format_dir = $curr_path . '/pmaf_changing_format_to_' . $step;
    
    mkdir($format_dir);
    print "\n\tDone ($format_dir).\n";    

    print STDERR "\n    *.ii) FORMATING ALIGMENT FILES.\n\n";


    ## Depending of the file type it will change format for aligments or trees

    my $bioperl_trees = { 
	                  'newick'     => 1, 
			  'nexus'       => 1, 
			  'nhx'         => 1, 
			  'lintree'     => 1, 
			  'treecluster' => 1, 
			  'pag'         => 1, 
			  'tab'         => 1, 
			  'svggraph'    => 1 
	                };

    my $bioperl_aligns = {
	                   'fasta'     => 2,
	                   'mase'      => 2, 
	                   'selex'     => 2, 
                           'clustalw'  => 2, 
	                   'msf'       => 2, 
	                   'phylip'    => 2, 
	                   'po'        => 2, 
	                   'stockholm' => 2, 
	                   'XMFA'      => 2, 
	                   'metafasta' => 2   
	                 };

    my $matrix_format = { 
	                  'matrix' => 3, 
			  'mat'    => 3 
                        };
   
    my @formats = ($bioperl_trees, $bioperl_aligns, $matrix_format);

    ## Now test if it is possible change the format based in current and next

    my ($curr_format_val, $next_format_val) = (0, 0);
    foreach my $formattype (@formats) {

	$curr_format_val += $formattype->{$curr_format} || 0;
	$next_format_val += $formattype->{$next_format} || 0;
    }

    if ($curr_format_val == 0) {
	die("FORMAT ERROR: From:$curr_format is not supported by this script.\n");
    }
    elsif ($next_format_val == 0) {
	die("FORMAT ERROR: To:$next_format is not supported by this script.\n");
    }
    elsif ($curr_format_val != $next_format_val) {
	die("FORMAT ERROR: current format($curr_format) and next format($next_format) have different types (aligment, matrix or tree)\n");
    }
    else {
	if ($curr_format_val == 2) {
	    
            ## It means that the formats are both alignments 

	    foreach my $align_file (@input_files) {
    
		my $out_basename = basename($align_file);	
		$out_basename =~ s/\.\w+$/\.$next_format/;	
		my $align_output = $format_dir . '/' . $out_basename;


		print STDERR "\tCHANGING FORMAT for file: $align_file   \r";	
	
		my $f_in = Bio::AlignIO->new( -file => "$align_file" ,
					      -format => $curr_format );
		
		my $f_out = Bio::AlignIO->new( -file => ">$align_output",
					       -format => $next_format );
	
		while ( my $align = $f_in->next_aln() ) { 
		    $f_out->write_aln($align); 
		}    
		push @output_files, $align_output;
	    }
	}
	elsif ($curr_format_val == 1) {

	    ## It means that both formats are trees

	    foreach my $tree_file (@input_files) {
    
		my $out_basename = basename($tree_file);	
		$out_basename =~ s/\.\w+$/\.$next_format/;	
		my $tree_output = $format_dir . '/' . $out_basename;


		print STDERR "\tCHANGING FORMAT for file: $tree_file   \r";	
	
		my $f_in = Bio::TreeIO->new( -file => "$tree_file" ,
					     -format => $curr_format );
		
		my $f_out = Bio::TreeIO->new( -file => ">$tree_output",
					      -format => $next_format );
	
		while ( my $tree = $f_in->next_tree() ) { 
		    $f_out->write_tree($tree); 
		}    
		push @output_files, $tree_output;
	    }
	}
	elsif ($curr_format_val == 3) {
	
	    ## The format is a matrix. This format can not be changed for now
	    
	    @output_files = @input_files;
	}
    }

    my $n_output_files = scalar(@output_files);
    print STDERR "\n\n\n\tDone... $n_output_files have been created in the dir:$format_dir.\n";

    return @output_files;
}




=head2 run_bootstrap

  Usage: my %out_files = run_bootstrap(\%boot_args, \@input_files, $out_dirname);

  Desc: Run the bootstrap programs over an array of fasta sequences.
        It return the sequences in the format specified in the arguments

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%boot_args, a hash reference with keys=argument_type and value=value
        \@input_file, an array of input filenames
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_bootstrap(\%align_args, \@input_files, $out_dirname);

=cut

sub run_bootstrap {
    my $boot_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None bootstrap_argument hash reference have been supplied to run_bootstrap function.\n");
    
    my $input_files_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file array reference have been supplied to run_bootstrap function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_bootstrap function.\n");

    my %outfiles;

    ## First, check if it need command files (for example using seqboot from phylip)

    my @input_files;
    my %ref_input_files;

    if ($boot_args_href->{'input'} =~ m/\{cmd\}/) {

	## It will create the command file, and run the bootstrap program using them

	my $cmd_dir = $output_dir . '/cmd';
	mkdir($cmd_dir); 

	## Now it will print them. Each cmd file will have the same name than 
	## the input file with the prefix 'cmd_' and the suffix 'txt'

	foreach my $infile (@{$input_files_aref}) {
	    
	    my $basename = basename($infile);
	    $basename =~ s/\.(\w+)$/\.txt/;
	    my $filename = $cmd_dir . '/cmd_' . $basename;

	    open(my $cmd_fh, '>', $filename);

	    print $cmd_fh "$infile\n";
	    
	    ## Now it will print the options (of there are any).
	    ## The bootstrap command (at least for seqboot from phylip)
	    ## needs one line per option AND result

	    if (defined $boot_args_href->{'options'}) {
		my @opts = split(/,/, $boot_args_href->{'options'});

		foreach my $opt (@opts) {
		    my @opt_data = split(/\s+/, $opt);
		    
		    foreach my $optdata (@opt_data) {
			print $cmd_fh "$optdata\n";
		    }
		}
	    }
	    
	    ## After the options it need confirm them with 'Y'

	    print $cmd_fh "Y\n";

	    ## And finally it will ask for a seed number.
	    ## It should be odd and between 1 and 32767
	    ## It is irrelevant, so as example, we can use
	    ## 1809 (Darwin Birth)
	
	    print $cmd_fh "3333\n";
	    
	    ## Finally it will add to infiles

	    push @input_files, $filename;
	    $ref_input_files{$filename} = $infile;
	}

	## To keep a record of the command files 

	$outfiles{'command'} = \@input_files;
    }
    else {

	## In it does not use command files
	@input_files = @{$input_files_aref};
    }

    ## Now the input files are set it will run one command per input file

    foreach my $in (@input_files) {
    
	my $command = $boot_args_href->{'executable'};

	## Now add the options

	if (defined $boot_args_href->{'options'} && $boot_args_href->{'input'} !~ m/\{cmd\}/) {
	    
	    my @options = split(/,/, $boot_args_href->{'options'});
	    my $option_string = join(' ', @options);
	    $command .= ' ' . $option_string;
	}

	## Add the input, depending of the input type 

	unless (defined $boot_args_href->{'input'}) { ## It will add as STDIN
	    
	    $command .= ' ' . $in
	}
	else {
	    if ($boot_args_href->{'input'} =~ m/\{cmd\}/) {
		$command .= ' <' . $in;
	    }
	    else {
		$command .= $boot_args_href->{'input'} . ' ' . $in;
	    }
	}

	## Create an output filename

	my $outformat = $boot_args_href->{'output_format'} || 'phylip';

	my $basename = basename($in);

	$basename =~ s/^cmd_//;                             ## Remove the prefix
	$basename =~ s/\.(\w+)$/_bootsp\.$outformat/;        ## Add outformat as suffix
	my $out_filename = $output_dir . '/' . $basename;
	
	## Add the output, depending of the output type

	if (defined $boot_args_href->{'output'} && $boot_args_href->{'output'} !~ m/\{predefined/) {
	    $command .= ' ' . $boot_args_href->{'output'} . ' ' . $out_filename; 	
	}
	elsif (defined $boot_args_href->{'output'}) {
	    ## If it is defined an output as predefined, it will not add nothing to the
	    ## command but it will change the names after... but it is interesing store 
	    ## the logs

	    $command .= ' >> ' . $output_dir . '/bootstrap_screenlog.txt';  
	}
	else {
	    $command .= ' >' . $out_filename; 
	}

	## Now the command is complete (executable + options + input + output)
	
	print STDERR "\tRUNNING: $command     \r";
	system($command);

	## Okay, it have been run. Now get the output for the predefined outputs

	## The bootstrap do not incluide the original filee so, it will add it using 
	## the following action

	if (defined $boot_args_href->{'output'} && $boot_args_href->{'output'} =~ m/\{predefined=(.+)\}/) {
	    my @list = split(/,/, $1);
	    
	    foreach my $filetype (@list) {

		## Change the name adding the filetype before the extension

		my $new_filename = $out_filename;
		$new_filename =~ s/\.(w+)$/$filetype\.$1/;
		
		## Join the input file before the bootstrap to the bootstrap
		## It will get before do the consensus

		system("cat $ref_input_files{$in} $filetype > $new_filename");
		
		## After cat, it will remove the file

		system("rm $filetype");

		if (exists $outfiles{$filetype}) {
		    push @{$outfiles{$filetype}}, $new_filename;
		}
		else {
		    $outfiles{$filetype} = [$new_filename];
		}
	    }
	}
	else {
	    ## By default it will store as filetype=output
	    
	    if (exists $outfiles{'output'}) {
		push @{$outfiles{'output'}}, $out_filename;
	    }
	    else {
		$outfiles{'output'} = [$out_filename];
	    }
	}
    }

    ## Finally, return the hash
    return %outfiles;
}


=head2 run_distmatrix

  Usage: my %out_files = run_distmatrix(\%dist_args, \@input_files, $out_dirname);

  Desc: Run the bootstrap programs over an array of fasta sequences.
        It return the sequences in the format specified in the arguments

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%dist_args, a hash reference with keys=argument_type and value=value
        \@input_file, an array of input filenames
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_distmatrix(\%dist_args, \@input_files, $out_dirname);

=cut

sub run_distmatrix {
    my $dist_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None dist_argument hash reference have been supplied to run_distmatrix function.\n");
    
    my $input_files_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file array reference have been supplied to run_distmatrix function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_distmatrix function.\n");

    my %outfiles;

    ## First, check if it need command files (for example using dnadist from phylip)

    my @input_files;

    if ($dist_args_href->{'input'} =~ m/\{cmd\}/) {

	## It will create the command file, and run the bootstrap program using them

	my $cmd_dir = $output_dir . '/cmd';
	mkdir($cmd_dir); 

	## Now it will print them. Each cmd file will have the same name than 
	## the input file with the prefix 'cmd_' and the suffix 'txt'

	foreach my $infile (@{$input_files_aref}) {
	    
	    my $basename = basename($infile);
	    $basename =~ s/_outfile//;
	    $basename =~ s/_outtree//;
	    $basename =~ s/\.(\w+)$/\.txt/;
	    my $filename = $cmd_dir . '/cmd_' . $basename;

	    open(my $cmd_fh, '>', $filename);

	    print $cmd_fh "$infile\n";
	    
	    ## Now it will print the options (of there are any).
	    ## The bootstrap command (at least for dnadist from phylip)
	    ## needs one line per option result pair

	    if (defined $dist_args_href->{'options'}) {
		my @opts = split(/,/, $dist_args_href->{'options'});

		foreach my $opt (@opts) {		    

		    print $cmd_fh "$opt\n";		    
		}
	    }
	    
	    ## After the options it need confirm them with 'Y'

	    print $cmd_fh "Y\n";
	   	    
	    ## Finally it will add to infiles

	    push @input_files, $filename;
	}

	## To keep a record of the command files 

	$outfiles{'command'} = \@input_files;
    }
    else {

	## In it does not use command files
	@input_files = @{$input_files_aref};
    }

    ## Now the input files are set it will run one command per input file

    foreach my $in (@input_files) {
    
	my $command = $dist_args_href->{'executable'};

	## Now add the options

	if (defined $dist_args_href->{'options'} && $dist_args_href->{'input'} !~ m/\{cmd\}/) {
	    
	    my @options = split(/,/, $dist_args_href->{'options'});
	    my $option_string = join(' ', @options);
	    $command .= ' ' . $option_string;
	}

	## Add the input, depending of the input type 

	unless (defined $dist_args_href->{'input'}) { ## It will add as STDIN
	    
	    $command .= ' ' . $in
	}
	else {
	    if ($dist_args_href->{'input'} =~ m/\{cmd\}/) {
		$command .= ' <' . $in;
	    }
	    else {
		$command .= $dist_args_href->{'input'} . ' ' . $in;
	    }
	}

	## Create an output filename

	my $outformat = $dist_args_href->{'output_format'} || 'matrix';

	my $basename = basename($in);

	$basename =~ s/_bootsp//;
	$basename =~ s/^cmd_//;                             ## Remove the prefix
	$basename =~ s/\.(w+)$/_distmatrix\.$outformat/;        ## Add outformat as suffix
	my $out_filename = $output_dir . '/' . $basename;
	
	## Add the output, depending of the output type

	if (defined $dist_args_href->{'output'} && $dist_args_href->{'output'} !~ m/\{predefined/) {
	    $command .= ' ' . $dist_args_href->{'output'} . ' ' . $out_filename; 	
	}
	elsif (defined $dist_args_href->{'output'}) {
	    ## If it is defined an output as predefined, it will not add nothing to the
	    ## command but it will change the names after... but it is interesing store 
	    ## the logs

	    $command .= ' >> ' . $output_dir . '/distance_matrix_screenlog.txt';  
	}
	else {
	    $command .= ' >' . $out_filename; 
	}

	## Now the command is complete (executable + options + input + output)
	
	print STDERR "\tRUNNING: $command\r";
	system($command);

	## Okay, it have been run. Now get the output for the predefined outputs

	if (defined $dist_args_href->{'output'} && $dist_args_href->{'output'} =~ m/\{predefined=(.+)\}/) {
	    my @list = split(/,/, $1);
	    
	    foreach my $filetype (@list) {

		## Change the name adding the filetype before the extension

		my $new_filename = $out_filename;
		$new_filename =~ s/\.(\w+)$/_$filetype\.$1/; 
		rename($filetype, $new_filename);
		
		if (exists $outfiles{$filetype}) {
		    push @{$outfiles{$filetype}}, $new_filename;
		}
		else {
		    $outfiles{$filetype} = [$new_filename];
		}
	    }
	}
	else {
	    ## By default it will store as filetype=output
	    
	    if (exists $outfiles{'output'}) {
		push @{$outfiles{'output'}}, $out_filename;
	    }
	    else {
		$outfiles{'output'} = [$out_filename];
	    }
	}
    }

    ## Finally, return the hash
    return %outfiles;
}


=head2 run_treecalc

  Usage: my %out_files = run_treecalc(\%dist_args, \@input_files, $out_dirname);

  Desc: Run the bootstrap programs over an array of fasta sequences.
        It return the sequences in the format specified in the arguments

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%tree_args, a hash reference with keys=argument_type and value=value
        \@input_file, an array of input filenames
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_treecalc(\%tree_args, \@input_files, $out_dirname);

=cut

sub run_treecalc {
    my $tree_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None tree_argument hash reference have been supplied to run_treecalc function.\n");
    
    my $input_files_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file array reference have been supplied to run_treecalc function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_treecalc function.\n");

    my %outfiles;

    ## First, check if it need command files (for example using neighbor from phylip)

    my @input_files;

    if ($tree_args_href->{'input'} =~ m/\{cmd\}/) {

	## It will create the command file, and run the bootstrap program using them

	my $cmd_dir = $output_dir . '/cmd';
	mkdir($cmd_dir); 

	## Now it will print them. Each cmd file will have the same name than 
	## the input file with the prefix 'cmd_' and the suffix 'txt'

	foreach my $infile (@{$input_files_aref}) {
	    
	    my $basename = basename($infile);
	    $basename =~ s/_outfile//;
	    $basename =~ s/_outtree//;
	    $basename =~ s/\.(\w+)$/\.txt/;
	    my $filename = $cmd_dir . '/cmd_' . $basename;

	    open(my $cmd_fh, '>', $filename);

	    print $cmd_fh "$infile\n";
	    
	    ## Now it will print the options (of there are any).
	    ## The bootstrap command (at least for dnadist from phylip)
	    ## needs one line per option result pair

	    if (defined $tree_args_href->{'options'}) {
		my @opts = split(/,/, $tree_args_href->{'options'});

		foreach my $opt (@opts) {		    
		    print $cmd_fh "$opt\n";		    
		}
	    }
	    
	    ## After the options it need confirm them with 'Y'

	    print $cmd_fh "Y\n";
	   	    
	    ## Finally it will add to infiles

	    push @input_files, $filename;
	}

	## To keep a record of the command files 

	$outfiles{'command'} = \@input_files;
    }
    else {

	## In it does not use command files
	@input_files = @{$input_files_aref};
    }

    ## Now the input files are set it will run one command per input file

    foreach my $in (@input_files) {
    
	my $command = $tree_args_href->{'executable'};

	## Now add the options

	if (defined $tree_args_href->{'options'} && $tree_args_href->{'input'} !~ m/\{cmd\}/) {
	    
	    my @options = split(/,/, $tree_args_href->{'options'});
	    my $option_string = join(' ', @options);
	    $command .= ' ' . $option_string;
	}

	## Add the input, depending of the input type 

	unless (defined $tree_args_href->{'input'}) { ## It will add as STDIN
	    
	    $command .= ' ' . $in
	}
	else {
	    if ($tree_args_href->{'input'} =~ m/\{cmd\}/) {
		$command .= ' <' . $in;
	    }
	    else {
		$command .= $tree_args_href->{'input'} . ' ' . $in;
	    }
	}

	## Create an output filename

	my $outformat = $tree_args_href->{'output_format'} || 'tree';

	my $basename = basename($in);

	$basename =~ s/^cmd_//;                             ## Remove the prefix
	$basename =~ s/_distmatrix//; 
	$basename =~ s/\.(\w+)$/_tree\.$outformat/;        ## Add outformat as suffix
	my $out_filename = $output_dir . '/' . $basename;
	
	## Add the output, depending of the output type

	if (defined $tree_args_href->{'output'} && $tree_args_href->{'output'} !~ m/\{predefined/) {
	    $command .= ' ' . $tree_args_href->{'output'} . ' ' . $out_filename; 	
	}
	elsif (defined $tree_args_href->{'output'}) {
	    ## If it is defined an output as predefined, it will not add nothing to the
	    ## command but it will change the names after... but it is interesing store 
	    ## the logs

	    $command .= ' >> ' . $output_dir . '/tree_screenlog.txt';  
	}
	else {
	    $command .= ' >' . $out_filename; 
	}

	## Now the command is complete (executable + options + input + output)
	
	print STDERR "\tRUNNING: $command\r";
	system($command);

	## Okay, it have been run. Now get the output for the predefined outputs

	if (defined $tree_args_href->{'output'} && $tree_args_href->{'output'} =~ m/\{predefined=(.+)\}/) {
	    my @list = split(/,/, $1);
	    
	    foreach my $filetype (@list) {

		## Change the name adding the filetype before the extension

		my $new_filename = $out_filename;
		$new_filename =~ s/\.(\w+)$/_$filetype\.$1/; 
		rename($filetype, $new_filename);
		
		if (exists $outfiles{$filetype}) {
		    push @{$outfiles{$filetype}}, $new_filename;
		}
		else {
		    $outfiles{$filetype} = [$new_filename];
		}
	    }
	}
	else {
	    ## By default it will store as filetype=output
	    
	    if (exists $outfiles{'output'}) {
		push @{$outfiles{'output'}}, $out_filename;
	    }
	    else {
		$outfiles{'output'} = [$out_filename];
	    }
	}
    }

    ## Finally, return the hash
    return %outfiles;
}


=head2 run_treemodf

  Usage: my %out_files = run_treemodf(\%treemodf_args, \@input_files, $out_dirname);

  Desc: Run the tree programs over an array of input files.
        Depending of the program, the input file can be alignments or distance matrix
        It return the sequences in the format specified in the arguments

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%treemodf_args, a hash reference with keys=argument_type and value=value
        \@input_file, an array of input filenames
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_treemodf(\%treemodf_args, \@input_files, $out_dirname);

=cut

sub run_treemodf {
    my $treemodf_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None tree modif argument hash reference have been supplied to run_treemodf function.\n");
    
    my $input_files_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file array reference have been supplied to run_treemodf function.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_treemodf function.\n");

    my %outfiles;

    ## First, check if it need command files (for example using retree from phylip)

    my @input_files;

    if ($treemodf_args_href->{'input'} =~ m/\{cmd\}/) {

	## It will create the command file, and run the retree program using them
	## In that case the retree infile will be compose by three diferent parts
	## 1) General options
	## 2) Tree modification options
	## 3) End the infile
	## 4) change the name to output name

	my $cmd_dir = $output_dir . '/cmd';
	mkdir($cmd_dir); 

	## Now it will print them. Each cmd file will have the same name than 
	## the input file with the prefix 'cmd_' and the suffix 'txt'

	## Check if exists intree, if it is the case print a warning because the 
	## program will read that

	if (-e 'intree') {
	    print STDERR "\n\nWARNING: Exists an intree file. Some phylip programs will read that ";
	    print STDERR "as input instead the tree detailed in the command file.\n";
	}

	foreach my $infile (@{$input_files_aref}) {
	    
	    ## Now it will count how many tree there are into the infile
	    ## using grep with system

	    my $tree_n = `grep -c ';' $infile`;
	    chomp($tree_n);

	    if ($tree_n > 0) {

		my $basename = basename($infile);
		$basename =~ s/_outfile//;
		$basename =~ s/_outtree//;
		$basename =~ s/\.(\w+)$/\.txt/;
		my $filename = $cmd_dir . '/cmd_' . $basename;
		
		open(my $cmd_fh, '>', $filename);

		## The command file should start with Y

		print $cmd_fh "Y\n";

		print $cmd_fh "$infile\n";

		## Now it will print the options (of there are any).
		## The bootstrap command (at least for dnadist from phylip)
		## needs one line per option result pair

		if (defined $treemodf_args_href->{'general_options'}) {
		    my @opts = split(/,/, $treemodf_args_href->{'general_options'});
		    
		    foreach my $opt (@opts) {		    
			print $cmd_fh "$opt\n";		    
		    }
		}
	    
		## After the general options it need confirm them with 'Y'

		print $cmd_fh "Y\n";	    

		## It will use the tree_options and after that + (tree_n - 1)
		## Also after the first tree, it will append the results to the old one
		## (outtree)(characters +,Y,R for first and +,Y,A,R for the rest)

		my $n = 1;
		
		if (defined $treemodf_args_href->{'retree_options'}) {
		    my @rt_opts = split(/,/, $treemodf_args_href->{'retree_options'});

		    foreach my $rt_opt (@rt_opts) {		    
			print $cmd_fh "$rt_opt\n";		    
		    }
		}

		print $cmd_fh "+\nY\nR\n";

		while($n < $tree_n) {

		    if (defined $treemodf_args_href->{'retree_options'}) {
			my @rt_opts = split(/,/, $treemodf_args_href->{'retree_options'});
			
			foreach my $rt_opt (@rt_opts) {		    
			    print $cmd_fh "$rt_opt\n";		    
			}
		    }
		    print $cmd_fh "+\nY\nA\nR\n";
		    $n++;
		}	    
	  
		## Finally it will need to add the last tree

		if (defined $treemodf_args_href->{'retree_options'}) {
		    my @rt_opts = split(/,/, $treemodf_args_href->{'retree_options'});
		    
		    foreach my $rt_opt (@rt_opts) {		    
			print $cmd_fh "$rt_opt\n";		    
		    }
		}

		## And close it

		print $cmd_fh "X\nY\nA\nR\n";
	    
		## Finally it will add to infiles

		push @input_files, $filename;
	    }	    
	    else {
		print STDERR "\n\tMSG: $infile have not any tree. Skipping file.\n";
	    }
	}

	## To keep a record of the command files 

	$outfiles{'command'} = \@input_files;
    }
    else {

	## In it does not use command files
	@input_files = @{$input_files_aref};
    }

    ## Now the input files are set it will run one command per input file

    foreach my $in (@input_files) {
    
	my $command = $treemodf_args_href->{'executable'};

	## Now add the options

	if (defined $treemodf_args_href->{'options'} && $treemodf_args_href->{'input'} !~ m/\{cmd\}/) {
	    
	    my @options = split(/,/, $treemodf_args_href->{'options'});
	    my $option_string = join(' ', @options);
	    $command .= ' ' . $option_string;
	}

	## Add the input, depending of the input type 

	unless (defined $treemodf_args_href->{'input'}) { ## It will add as STDIN
	    
	    $command .= ' ' . $in
	}
	else {
	    if ($treemodf_args_href->{'input'} =~ m/\{cmd\}/) {
		$command .= ' <' . $in;
	    }
	    else {
		$command .= $treemodf_args_href->{'input'} . ' ' . $in;
	    }
	}

	## Create an output filename

	my $outformat = $treemodf_args_href->{'output_format'} || 'tree';

	my $basename = basename($in);

	$basename =~ s/^cmd_//;                             ## Remove the prefix
	$basename =~ s/_tree//;
	$basename =~ s/\.(\w+)$/_retree\.$outformat/;        ## Add outformat as suffix
	my $out_filename = $output_dir . '/' . $basename;
	
	## Add the output, depending of the output type

	if (defined $treemodf_args_href->{'output'} && $treemodf_args_href->{'output'} !~ m/\{predefined/) {
	    $command .= ' ' . $treemodf_args_href->{'output'} . ' ' . $out_filename; 	
	}
	elsif (defined $treemodf_args_href->{'output'}) {
	    ## If it is defined an output as predefined, it will not add nothing to the
	    ## command but it will change the names after... but it is interesing store 
	    ## the logs

	    $command .= ' >> ' . $output_dir . '/retree_screenlog.txt';  
	}
	else {
	    $command .= ' >' . $out_filename; 
	}

	## Now the command is complete (executable + options + input + output)
	
	print STDERR "\tRUNNING: $command\r";
	system($command);

	## Okay, it have been run. Now get the output for the predefined outputs

	if (defined $treemodf_args_href->{'output'} && $treemodf_args_href->{'output'} =~ m/\{predefined=(.+)\}/) {
	    my @list = split(/,/, $1);
	    
	    foreach my $filetype (@list) {

		## Change the name adding the filetype before the extension

		my $new_filename = $out_filename;
		$new_filename =~ s/\.(\w+)$/_$filetype\.$1/; 
		rename($filetype, $new_filename);
		
		if (exists $outfiles{$filetype}) {
		    push @{$outfiles{$filetype}}, $new_filename;
		}
		else {
		    $outfiles{$filetype} = [$new_filename];
		}
	    }
	}
	else {
	    ## By default it will store as filetype=output
	    
	    if (exists $outfiles{'output'}) {
		push @{$outfiles{'output'}}, $out_filename;
	    }
	    else {
		$outfiles{'output'} = [$out_filename];
	    }
	}
    }

    ## Finally, return the hash
    return %outfiles;
}

=head2 run_consensus_tree

  Usage: my %out_files = run_consensus_tree(\%constree_args, \@input_files, $out_dirname);

  Desc: Run the tree programs over an array of input files.
        It return the sequences in the format specified in the arguments

  Ret: %out_files, a hash with keys=file_type and value=filename

  Args: \%constree_args, a hash reference with keys=argument_type and value=value
        \@input_file, an array of input filenames
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: my %out_files = run_consensus_tree(\%constree_args, \@input_files, $out_dirname);

=cut

sub run_consensus_tree {
    my $constree_args_href = shift
	|| die("FUNCTION ARGUMENT ERROR: None consensus tree argument hash reference have been supplied to run_consensus_tree func.\n");
    
    my $input_files_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file array reference have been supplied to run_consensus_tree func.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_consensus_tree func.\n");

    my %outfiles;

    ## First, check if it need command files (for example using retree from phylip)

    my @input_files;
    my @ref_files;
    my @boots_files;

    ## before run the consensus tree it need to remove the first tree. It is not a bootstrap
    ## it is the original 

    foreach my $input_file (@{$input_files_aref}) {
	my $basename_o = basename($input_file);
	$basename_o =~ s/_outfile//;
	$basename_o =~ s/_outtree//;
	my $boot_file = $output_dir . '/boots_' . $basename_o;
	my $ref_file = $output_dir . '/ref_' . $basename_o;
	
	open(my $input_fh, '<', $input_file);
	open(my $ref_fh, '>', $ref_file);
	open(my $cons_fh, '>', $boot_file);
	
	my $tree = 0;
	
	while(<$input_fh>) {	    
	    chomp($_);

	    if ($_ =~ m/\(/) { ## The begining of a tree
		$tree++;
		
		if ($tree == 1) {
		    print $ref_fh "$_\n";
		}
		else {
		    print $cons_fh "$_\n";
		}
	    }
	}
	push @ref_files, $ref_file;
	push @boots_files, $boot_file;    
    }

    $outfiles{'reftree'} = \@ref_files;

    

    if ($constree_args_href->{'input'} =~ m/\{cmd\}/) {

	## It will create the command file, and run the consensus program using them


	my $cmd_dir = $output_dir . '/cmd';
	mkdir($cmd_dir); 

	## Now it will print them. Each cmd file will have the same name than 
	## the input file with the prefix 'cmd_' and the suffix 'txt'

	foreach my $infile (@boots_files) {
	    
	    my $basename = basename($infile);
	    $basename =~ s/_outfile//;
	    $basename =~ s/_outtree//;
	    $basename =~ s/\.(\w+)$/\.txt/;
	    my $filename = $cmd_dir . '/cmd_' . $basename;

	    open(my $cmd_fh, '>', $filename);

	    print $cmd_fh "$infile\n";

	    ## Now it will print the options 

	    if (defined $constree_args_href->{'options'}) {
		my @opts = split(/,/, $constree_args_href->{'options'});

		foreach my $opt (@opts) {		    
		    print $cmd_fh "$opt\n";		    
		}
	    }
	    
	    ## After the general options it need confirm them with 'Y'

	    print $cmd_fh "Y\n";
	    
	    ## Finally it will add to infiles

	    push @input_files, $filename;
	}

	## To keep a record of the command files 

	$outfiles{'command'} = \@input_files;
    }
    else {

	## In it does not use command files
	@input_files = @{$input_files_aref};
    }

    ## Now the input files are set it will run one command per input file

    foreach my $in (@input_files) {
    
	my $command = $constree_args_href->{'executable'};

	## Now add the options

	if (defined $constree_args_href->{'options'} && $constree_args_href->{'input'} !~ m/\{cmd\}/) {
	    
	    my @options = split(/,/, $constree_args_href->{'options'});
	    my $option_string = join(' ', @options);
	    $command .= ' ' . $option_string;
	}

	## Add the input, depending of the input type 

	unless (defined $constree_args_href->{'input'}) { ## It will add as STDIN
	    
	    $command .= ' ' . $in
	}
	else {
	    if ($constree_args_href->{'input'} =~ m/\{cmd\}/) {
		$command .= ' <' . $in;
	    }
	    else {
		$command .= $constree_args_href->{'input'} . ' ' . $in;
	    }
	}

	## Create an output filename

	my $outformat = $constree_args_href->{'output_format'} || 'tree';

	my $basename = basename($in);

	$basename =~ s/^cmd_//;                             ## Remove the prefix
	$basename =~ s/_tree//;
	$basename =~ s/\.(\w+)$/_consensus\.$outformat/;        ## Add outformat as suffix
	my $out_filename = $output_dir . '/' . $basename;
	
	## Add the output, depending of the output type

	if (defined $constree_args_href->{'output'} && $constree_args_href->{'output'} !~ m/\{predefined/) {
	    $command .= ' ' . $constree_args_href->{'output'} . ' ' . $out_filename; 	
	}
	elsif (defined $constree_args_href->{'output'}) {
	    ## If it is defined an output as predefined, it will not add nothing to the
	    ## command but it will change the names after... but it is interesing store 
	    ## the logs

	    $command .= ' >>' . $output_dir . '/constree_screenlog.txt';  
	}
	else {
	    $command .= ' >' . $out_filename; 
	}

	## Now the command is complete (executable + options + input + output)
	
	print STDERR "\tRUNNING: $command\r";
	system($command);

	## Okay, it have been run. Now get the output for the predefined outputs

	if (defined $constree_args_href->{'output'} && $constree_args_href->{'output'} =~ m/\{predefined=(.+)\}/) {
	    my @list = split(/,/, $1);
	    
	    foreach my $filetype (@list) {

		## Change the name adding the filetype before the extension

		my $new_filename = $out_filename;
		$new_filename =~ s/\.(\w+)$/_$filetype\.$1/; 
		rename($filetype, $new_filename);
		
		if (exists $outfiles{$filetype}) {
		    push @{$outfiles{$filetype}}, $new_filename;
		}
		else {
		    $outfiles{$filetype} = [$new_filename];
		}
	    }
	}
	else {
	    ## By default it will store as filetype=output
	    
	    if (exists $outfiles{'output'}) {
		push @{$outfiles{'output'}}, $out_filename;
	    }
	    else {
		$outfiles{'output'} = [$out_filename];
	    }
	}
    }    
    ## Finally, return the hash
    return %outfiles;
}

=head2 extract_ln_likehood

  Usage: extract_ln_likehood(\@input_files, $out_dirname);

  Desc: Extract the ln likehood variable from outfile produced by dnaml

  Ret: nothing, it will print the results as a file

  Args: \@input_file, an array of input filenames
        $out_dirname, a scalar with output dirname

  Side_Effects: Die if something is wrong
                Print parsing status messages

  Example: extract_ln_likehood(\@input_files, $out_dirname);

=cut

sub extract_ln_likehood {

    my $input_files_aref = shift
	|| die("FUNCTION ARGUMENT ERROR: None input file array reference have been supplied to run_consensus_tree func.\n");

    my $output_dir = shift
	|| die("FUNCTION ARGUMENT ERROR: None output dirname have been supplied to run_consensus_tree func.\n");


    ## it will add the result to a file

    my $filename = $output_dir . '/contigs_ln_likehood.txt';
    open my $ofh, '>', $filename;

    foreach my $infile (@{$input_files_aref}) {
    
	my $ctg;
	my $basename = basename($infile);
	if ($basename =~ m/(.+?)\..+/) {
	    $ctg = $1;
	}
	else {
	    $ctg = $basename;
	}

	open my $infh, '<', $infile;
	
	my %first_ln;
	my @ln;

	while(<$infh>) {
	
	    chomp($_);

	    if ($_ =~ m/^Ln\s+Likelihood\s+=\s+(.+)$/) {
		my $ln_value = $1;
		unless (exists $first_ln{$ctg}) {
		    $first_ln{$ctg} = $ln_value;
		}
		push @ln, $ln_value;
	    }
	}

	my $median = median(@ln);
	my $stddv = stddev(@ln);
	
	if (defined $ctg && defined $first_ln{$ctg}) {
	    print $ofh "$ctg\t$first_ln{$ctg}\t$median\t$stddv\n";
	}

	close $infh;
    }
}


=head2 create_default_ctr

  Usage: create_ctr
  Desc: create_ctr file
  Ret: none
  Args: filename
  Side_Effects: exit of the script
  Example: if ($opt_D) {
               create_default_ctr();
           }

=cut

sub create_default_ctr {    
    my $ctr_filename = shift || 'maf_to_phylo.default.ctr';
 
    open my $ctr_fh, '>', $ctr_filename 
	|| die("OPEN FILE ERROR: $ctr_filename file can not be openned (system error: $!),\n");

    
    print $ctr_fh <<EOF;

#############################################################
## This is the default controller file.                    ##
## To print a complete file with instruction use -X option ##
#############################################################

<STRAIN_LIST>="{extract}"
<MIN_ALIGNMENT_LENGTH>="100"
    
<REALIGNMENT_EXEC>="" 
<REALIGNMENT_OFORMAT>="" 
<REALIGNMENT_OPT>=""
<REALIGNMENT_IN>=""
<REALIGNMENT_OUT>=""

<BOOTSTRAP_EXEC>="seqboot"
<BOOTSTRAP_IFORMAT>="phylip"
<BOOTSTRAP_OFORMAT>="phylip" 
<BOOTSTRAP_OPT>=""
<BOOTSTRAP_REP>="1000"
<BOOTSTRAP_IN>="{cmd}"
<BOOTSTRAP_OUT>="{predefined=outfile}"
<BOOTSTRAP_CUTOFF>="500"

<DISTANCE_EXEC>="dnadist"
<DISTANCE_IFORMAT>="phylip"
<DISTANCE_OFORMAT>="matrix"
<DISTANCE_OPT>="D Kimura,2 No"
<DISTANCE_IN>="{cmd}"
<DISTANCE_OUT>="{predefined=outfile}"

<PHYLOTREE_EXEC>="neighbor"
<PHYLOTREE_IFORMAT>="matrix"
<PHYLOTREE_OFORMAT>="newick"
<PHYLOTREE_OPT>=""
<PHYLOTREE_IN>="{cmd}"
<PHYLOTREE_OUT>="{predefined=outtree,outfile}"

<RETREE_EXEC>="retree" 
<RETREE_IFORMAT>="newick"
<RETREE_OFORMAT>="newick"
<RETREE_GEN_OPT>=""
<RETREE_TREE_OPT>="M"
<RETREE_IN>="{cmd}"
<RETREE_OUT>="{predefined=outtree}"

<CONSTREE_EXEC>="consense"
<CONSTREE_IFORMAT>="newick" 
<CONSTREE_OFORMAT>="newick"
<CONSTREE_OPT>="R"
<CONSTREE_IN>="{cmd}"
<CONSTREE_OUT>="{predefined=outfile,outtree}"

<BRANCH_LENGTH_TAGS>="0.001=>SIM,1.0=>DIF"

EOF

print STDERR "\nControl file have been created ($ctr_filename)\n\n";

exit (1);

}




#####
1; ##
#####
