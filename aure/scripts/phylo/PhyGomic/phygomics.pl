#!/usr/bin/perl

=head1 NAME

 phygomics.pl
 Script to run PhyGomics pipeline.

=cut

=head1 SYPNOSIS

 phygomics.pl [-h] -c config.file -i input.dir -o output.dir [-C] [-S]

=head1 EXAMPLE 

 phygomics.pl -c phygo1.txt -i analysis1 -o results1

 ## To print an ampty configuration file:

 phygomics.pl -C
                                  
=head2 I<Flags:>

=over


=item -c

B<configuration.file>          Configuration file printed using -C.

=item -i

B<input.dir>                   Input dir. with input datasets (mandatory)

=item -o

B<output.dir>                  Output dir. to create output datasets (mandatory)

=item -h

B<help>                        Print the help

=item -C

B<create_config.file>         Create an empty configuration file with the name:
                              phygomics.conf 
=item -S

B<print_pipeline_status>      Print as STDOUT the pipeline status for parsing
                              and executing methods.

=back

=cut

=head1 DESCRIPTION

   phygomics.pl is the scripts that executes the PhyGomics pipeline.

   It will run the pipeline in two three steps:
     1) Data loading [1 cycle]
        1.1) Cluster extraction from one source, assembly file or selfblast.
        1.2) Sequence members loading from a fasta file.
        1.3) Strains for members loading from a tabular file.

     2) Data processing [as many cycles as paths are described]
        2.1) Search homologous using blast                   [optional]
        2.2) Run alignments                                  [mandatory]
        2.3) Run distances                                   [mandatory for NJ]
        2.4) Prune members, it will rerun 2.2 and 2.3        [optional]
        2.5) Run trees                                       [mandatory]
        2.6) Tree rerooting                                  [optional]
        2.7) Run bootstrapping                               [optional]
        2.8) Run topoanalysis                                [mandatory]

     3) Data integration and analysis [1 cycle]
        3.1) Create graphs and tables.

    (see -C, configuration file for more details about data)


=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 phygomics.pl

=cut

use strict;
use warnings;
use autodie;

use Cwd;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;

use FindBin;
use lib "$FindBin::Bin/lib";

use PhyGeCluster;

our ($opt_c, $opt_i, $opt_o, $opt_h, $opt_C, $opt_S);
getopts("c:i:o:hCS");
if (!$opt_c && !$opt_i && !$opt_o && !$opt_C && !$opt_S) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
elsif ($opt_C) {
    print_config();
}



## Mandatory arguments (die if they are not supplied)

my $err_m = "MANDATORY ARGUMENT ERROR:";

my $control = $opt_c 
    || die("\n$err_m -c <control_file> argument was not supplied.\n\n");

my $indir = $opt_i
    || die("\n$err_m -i <input_dir> argument was not supplied.\n\n");

my $outdir = $opt_o
    || die("\n$err_m -o <out_dir> argument was not supplied.\n\n");


## Check files and dirs

unless (-s $control) {
    die("\n$err_m $control configuration file doesnt exist or is empty.\n\n");
}
unless (-d $indir) {
    die("\n$err_m $indir input dir. doesnt exist.\n\n");
}
unless (-d $outdir) {
    print STDOUT ("\nMESSAGE: Outdir doesnt exist. Creating outdir=$outdir\n");
    mkdir($outdir);
}

##############################################
#### PIPELINE EXECUTION  #####################
##############################################

############################################
## 0) Print start message and start to parse

print STDERR "\n** STARTING PhyGomic PIPELINE (" . date() . "):\n\n";


#############################################################################
## 0) Parse configuration file, check files and PhyGeCluster object creation.

print STDERR "0) PARSING CONFIGURATION FILE (" .  date() . "):\n";

my %conf = parse_config($control);
my $cv = scalar(keys %conf);

print STDERR "\tPARSED. $cv configuration arguments have been parsed.\n\n";

## Add indir and outdir to the conf.

$conf{indir} = $indir;
$conf{outdir} = $outdir;

filecheck(\%conf);

## Create an empty phygecluster object

my $phyg = PhyGeCluster->new();


###############################################
## 1.1) Add the cluster source file and parse it.

print STDERR "1.1) PARSING CLUSTER SOURCE FILE (" .  date() . "):\n\n";

my $sc = $conf{cl_sc} || '';

my %clusters;
if ($sc eq 'blast') {

    my $blast_args = { blastfile => $indir . '/' . $conf{cl_fn} };

    ## Add report_status if $opt_S

    if ($opt_S) {
	$blast_args->{report_status} = 1;
    }

    ## Add cluster values if they are defined

    my $cls_val = parse_clustervals($conf{cl_vl}); 
    if (defined $cls_val) {
	$blast_args->{clustervalues} = $cls_val;
    }

    ## Choose the parsing according 'cl_fast' switch

    if (exists $conf{cl_fast}) {
	%clusters = PhyGeCluster::fastparse_blastfile($blast_args);
    }
    else {
	%clusters = PhyGeCluster::parse_blastfile($blast_args);
    }
}
else {  ## By default it will .ace file
    
    my $ace_args = { acefile => $indir . '/' . $conf{cl_fn} };

    ## Add report_status if $opt_S

    if ($opt_S) {
	$ace_args->{report_status} = 1;
    }

    %clusters = PhyGeCluster::parse_acefile($ace_args);
}

## Set the cluster for PhyGeCluster

$phyg->set_clusters(\%clusters);

my $cl_v = scalar(keys %clusters);
print STDERR "\tPARSED. $cl_v clusters have been extracted.\n\n";


###################################################
## 1.2) Parse the strain file and add to phygecluster

my $seq_args = { sequencefile => $indir . '/' . $conf{seq_fn} };
if ($opt_S) {
    $seq_args->{report_status} = 1;
}

print STDERR "1.2) PARSING SEQUENCE MEMBER FILE (" .  date() . "):\n\n";
$phyg->load_seqfile($seq_args);

my $seqn = 0;
my %cls = %{$phyg->get_clusters()};
foreach my $clid (keys %cls) {
    my @members = $cls{$clid}->get_members();
    $seqn += scalar(@members);
}

print STDERR "\tPARSED. $seqn sequences have been extracted.\n\n";


###################################################
## 1.3) Parse the strain file and add to phygecluster

my $str_args = { strainfile => $indir . '/' . $conf{str_fn} };
if ($opt_S) {
    $str_args->{report_status} = 1;
}

print STDERR "1.3) PARSING STRAIN FILE (" .  date() . "):\n\n";
$phyg->load_strainfile($str_args);
my $str_n = scalar(keys %{$phyg->get_strains()});

print STDERR "\tPARSED. $str_n strains have been extracted.\n\n";


##############################################################################
## 2) Now the pipeline will run different paths depending of the configuration
##    file. Each path will clone phygecluster and run a sepparate analysis.
##    at the end, it will integrate all of them in a phygestat object.

my %paths = %{$conf{paths}};
my $pn = scalar(keys %paths);

print STDERR "\n** $pn PATHS has been defined (".date()."):\n\n";

my %phygecls = ();

foreach my $path_idx (sort {$a <=> $b} keys %paths) {
    my %pargs = %{$paths{$path_idx}};
    my $phname = $pargs{pa_name} || $path_idx; 

    print STDERR "\t$path_idx ==> INIT. PATH=$phname (" .  date() . "):\n\n";

    ## 2.0) Clone the PhyGeCluster

    print STDERR "\t\t2.0) CLONNING CLUSTER DATA (" .  date() . "):\n\n";

    my $paphyg = $phyg->clone();

    ## 2.1) Homologous search

    print STDERR "\t\t2.1) HOMOLOGOUS SEARCH (" .  date() . "):\n\n";


    if (defined $pargs{hs_dts} && defined $pargs{hs_str}) {
	
	my %hom_args = (
	    -blast  => [ -p => 'blastn', -d => $indir . '/' . $pargs{hs_dts} ],
	    -strain => $pargs{hs_dts},
	    );
	
	if (exists $pargs{hs_arg}) {
	    my @bargs = split(/;/, $pargs{hs_arg});
	    foreach my $barg (@args) {
		if ($bargs =~ m/\s*(-\w)\s+(.+?)\s*/) {
		    my $b_ar = $1;
		    my $b_vl = $2;
		    unless ($b_ar =~ m/-(o|i|p|d)/) {
			push @{$hom_args{-blast}}, $b_ar => $b_vl; 
		    }
		}
	    }
	}
	if (exists $pargs{hs_fil}) {
	    my @bfils = split(/;/, $pargs{hs_fil});
	    foreach my $bfil (@fils) {
		if ($bfil =~ m/\s*(.+?)\s*(.+?)\s*(.+?)\s*/) {
		    my ($arg, $ope, $val) = ($1, $2, $3);
		    if (exists $hom_args{-filter}) {
			$hom_args{-filter}->{$arg} = [$2, $3];
		    }
		    else {
			$hom_args{-filter} = { $arg = [$2, $3] };
		    }
		}
	    }
	}
	
	## And run homologous search
	
	$paphyg->homologous_search($hom_args_href);

    }
    else {
    }



    print STDERR "\t$path_idx ==> END. PATH=$phname (" .  date() . "):\n\n";
}














##############################################
#### METHODS #################################
##############################################

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
 
      phygomics.pl is the scripts that executes the PhyGomics pipeline.

      It will run the pipeline in two three steps:
     
      1) Data loading [1 cycle]
        1.1) Cluster extraction from one source, assembly file or selfblast.
        1.2) Sequence members loading from a fasta file.
        1.3) Strains for members loading from a tabular file.

      2) Data processing [as many cycles as paths are described]
        2.1) Search homologous using blast                   [optional]
        2.2) Run alignments                                  [mandatory]
        2.3) Run distances                                   [mandatory for NJ]
        2.4) Prune members, it will rerun 2.2 and 2.3        [optional]
        2.5) Run trees                                       [mandatory]
        2.6) Tree rerooting                                  [optional]
        2.7) Run bootstrapping                               [optional]
        2.8) Run topoanalysis                                [mandatory]

      3) Data integration and analysis [1 cycle]
        3.1) Create graphs and tables.

      (see -C, configuration file for more details about data)


    Usage:
  
     phygomics.pl [-h] -c config.file -i input.dir -o output.dir [-C]
                              
    Example:

     phygomics.pl -c phygo1.txt -i analysis1 -o results1

    Flags:
     
      -c <configuration.file>  Configuration file printed using -C.
      -i <input.dir>           Input dir. with input datasets (mandatory)
      -o <output.dir>          Output dir. to create output datasets (mandatory)
      -h <help>                Print the help
      -C <create_config.file>  Create an empty configuration file with the name:
                               phygomics.conf 
      -S <print_status>        Print as STDOUT the pipeline status for parsing
                               and executing methods.

EOF
exit (1);
}

=head2 date

  Usage: my $date = date()

  Desc: return a string with the current date run by unix command 'date'

  Ret: A string

  Args: none

  Side_Effects: none

  Example: print STDERR date()

=cut

sub date {
    my $date = `date`;
    chomp($date);
    return $date;
}


=head2 filecheck

  Usage: filecheck(\%conf);

  Desc: Check if exists all the files from the conf. hash. If not, die.

  Ret: None

  Args: \%conf, a hash reference with the configuration data.

  Side_Effects: Die if no conf. hashref. is supplied to filecheck().

  Example: filecheck(\%conf);

=cut

sub filecheck {
    my $confhref = shift ||
	die("ERROR: No config. hashref. was supplied to filecheck function.\n");

    if(ref($confhref) ne 'HASH') {
	die("ERROR: Config. hashref supplied to filecheck() isnt a hashref.\n");
    }

    my %reqfiles = (
	cl_fn  => 'cluster_file',
	seq_fn => 'membersequence_file',
	str_fn => 'memberstrain_file',
	);

    foreach my $req_k (keys %reqfiles) {
	unless (exists $confhref->{$req_k}) {
	    my $e = "ERROR: Required file:$reqfiles{$req_k} didnt was";
	    die("$e supplied to the configuration file.\n");
	}
    }

    my $indir = $confhref->{indir} || '';
    foreach my $confk (sort keys %{$confhref}) {
	if ($confk =~ m/_fn$/) {
	    my $file = $indir . '/' . $confhref->{$confk};
	    unless (-s $file) {
		die("ERROR: file=$file doesnt exist or is empty.\n");
	    }
	}
    }
}



##############################################
#### CONFIGURATION METHODS ###################
##############################################


=head2 print_config

  Usage: print_config()

  Desc: create and print a configuration file.

  Ret: none

  Args: none

  Side_Effects: none

  Example: print_config()

=cut

sub print_config {
    
    my $dir = getcwd;
    my $contrfile = $dir . '/phygomics.conf';
    open my $cfh, '>', $contrfile;

    my $info = q/## Configuration file start.
## This is a control file for PhyGomics pipeline:
## Fill the data between square brackets without break each line.

###############################################################################
## 1) DATA LOADING ############################################################
###############################################################################
##
## It can be based in two approaches: 
##   1) A sequence assembly. Each contig will be a cluster. 
##      Recomended for species close related where homologous genes can be
##      coassembled in the same contig.
##
##   2) A selfblast. Cluster will be created based in sequence match results.
##      Recomended for species where homologous genes cannot be coasssembled.

######################
## 1.1 ## MANDATORY ##
######################

<>CLUSTER_DATASOURCE:          []
##  
## Source type to extract the clusters
##
## format:  [{value}]
## example: [ace]
## values:  blast, ace

<>CLUSTER_FILENAME:            []  
##
## Name of the file to extract the clusters.
##
## format:  [{filename_without_dir}]
## example: [assembly_out.test.ace]


<>CLUSTER_VALUES:              []  
##
## Arguments to use to analyze the blast result file to get sequence clusters
##
## format:  [{parameter}{space}{operator}{space}{value}{semicolon}...]
## example: [score > 100; percent_identity > 90]
## parameters: evalue, expect, frac_identical, frac_conserved, gaps, 
##             hsp_length, num_conserved, num_identical, score, bits, 
##             percent_identity, 
## If fastblastparser is enabled, the valid values are:
##             query_id, subject_id, percent_identity, align_length, 
##             mismatches, gaps_openings, q_start, q_end, s_start, s_end, 
##             e_value, bit_score


<>FASTBLASTPARSER:             []
##
## Switch to use the fast blast parser (to enable add any character inside [])
##
## format: [1]
## example [1]


######################
## 1.2 ## MANDATORY ##
######################

<>MEMBERSEQ_FILENAME:          []
##
## Sequence file for the members of the blast or ace file in fasta format.
##
## format:  [{filename_without_dir}]
## example: [seq.test.fasta]


######################
## 1.3 ## MANDATORY ##
######################

<>MEMBERSTRAIN_FILENAME:       []
## File with 2 columns: -f1 member_id, and -f2 strain name.
##
## format:  [{filename_without_dir}]
## example: [strains.test.tab]


###############################################################################
## ANALYSIS PATHS #############################################################
###############################################################################
##
## The phygomics pipeline is designed to work with different analysis paths,
## so it could run different alignment methods, distance calculations and\/or
## tree analysis. To do that it defines analysis paths using an integer 
## between diamond brackets <>
##
## For example, analysis path 1 will be annotated as <1> 

<>PATH_NAME:                   []
##
## A name to define the path.
##
## format:  [{name}]
## example: ['ML']

#####################
## 2.1 ## OPTIONAL ##
#####################

<>HOMOLOGOUS_SEARCH_ARGUMENTS: []
##
## Blast arguments to use with the homologous search
## 
## format:  [{blast_argument}{space}{value}{space}{semicolon}...]
## example: [-e 1e-10; -a 2]
## arguments: all the blast arguments except: -i, -o, -p, -d

<>HOMOLOGOUS_SEARCH_DATASET:   []
##
## Dataset to use with the homologous search
##
## format:  [{blastdb_file_without_dir}]
## example: [blastref.test.fasta]

<>HOMOLOGOUS_SEARCH_STRAIN:    []
##
## Strain to define as homologous
##
## format:  [{strain}]
## example: [Sly]

<>HOMOLOGOUS_SEARCH_FILTER:    []
##
## To filter the blast output before parse and analyze to find the homologous
##
## format:  [{parameter}{space}{operator}{space}{value}{semicolon}...]
## example: [score > 100; percent_identity > 90]
## parameters: evalue, expect, frac_identical, frac_conserved, gaps, hsp_length,
##             num_conserved, num_identical, score, bits, percent_identity

######################
## 2.2 ## MANDATORY ##
######################

<>RUN_ALIGNMENT_PROGRAM:       []
##
## Name of the sequence alignment program to use. 
##
## format:  [{program}]
## example: [clustalw]
## programs: clustalw, kalign, MAFFT, muscle, tcoffee

<>RUN_ALIGNMENT_ARGUMENTS:     []
##
## Arguments to be used with the program. For more infomation try perldoc
## Bio::Tools::Run::Alignment.
##
## format:  [{parameter}{=}{value}{semicolon}...]
## example: [quiet = yes; matrix = BLOSUM]


###################################
## 2.3 ## MANDATORY FOR NJ TREES ##
###################################

<>RUN_DISTANCE_FUNCTION:       []
##
## Name of the distance method to use. Required for NJ trees and ignored for ML.
##
## format:  [{function}] 
## example: [Kimura]
## functions: JukesCantor, Uncorrected, F81, Kimura,  Tamura or TajimaNei.  


######################
## 2.4 ## OPTIONALS ##
######################

<>PRUNE_ALIGN_ARGUMENTS:       []
##
## Prune (remove) elements based in the alignment data.
##
## format:    [{argument}{space}{operator}{space}{value}{semicolon}...]
## example:   [ num_sequences < 4; length < 100 ]
## arguments: score, length, num_residues, num_sequences or percentage_identity

<>PRUNE_STRAINS_ARGUMENTS:     []
##
## Prune (remove) elements based in the strain data that dont have the 
## specified conditions
##
## format:   [composition => {strain}={value}{comma}...{semicolon}...]
##           [min_distance => {strain1}={strain2}{comma}...{semicolon}]
## example:  [composition => Sly=1,Nta=2,Nto=1,Nsy=1]
##           [min_distance => Nta=Nta,Nsy=Nta,Nto=Nta,Nsy=Nto,Nta=Sly]
## arguments: composition, min_distance, max_distance

<>PRUNE_OVERLAPS_ARGUMENTS:    []
##
## Prune (remove) elements based in the alignment overlap values. It takes the
## best overlap
##
## format:   [composition => {strain}={value}{comma}...{semicolon}...]
##           [random => {integer}]
             [trim   => {1}]
## example:  [composition => Sly=1,Nta=2,Nto=1,Nsy=1]
##           [trim => 1 ]
## arguments: composition, random, trim


######################
## 2.5 ## MANDATORY ##
######################

<>RUN_TREE_METHOD:             []
##
## Name of the tree method to run (NJ, UPGMA, ML)
##
## format:  [{method}]
## example: [ML]
## methods: NJ, UPGMA, ML

<>RUN_TREE_ARGUMENTS:          []
##
## List of arguments to use with the tree calculation
##
## format:  [{argument}{space}{=}{space}{value}{semicolon}...]
## example: [quiet = 1; outgroup_strain = Sly]
## arguments for NJ or UPGMA: lowtri (0 or 1), uptri (0 or 1), subrep (0 or 1)
##                            jumble (integer), quiet (0 or 1) or 
##                            outgroup_strain (strain).
## arguments for ML: phyml, dnaml and outgroup_strain (strain)


#####################
## 2.6 ## OPTIONAL ##
#####################

<>REROOT_TREE:                 []
##
## Reroot tree based the specified method 
##
## format:  [{method}{=}{value}]
## example: [midpoint=1]
## methods: midpoint (1), strainref (strain) and longest (1)


#####################
## 2.7 ## OPTIONAL ##
#####################

<>RUN_BOOTSTRAPPING:           []
##
## Number of replicates. It will use the same distance and tree method than
## the rest of the path
##
## format:  [{integer}]
## example: [1000]

<>FILTER_BOOTSTRAPPING:        []
##
## Minimum bootstrapping value to discard a tree from a topology analysis.
## Values are normalized to 100.
##
## format:  [{integer}]
## example: [1000]


######################
## 2.8 ## MANDATORY ##
######################

<>RUN_TOPOANALYSIS:            []
##
## Run the topology analysis over the trees set. 
##
## format:  [branch_cutoff{=>}{integer}{=}{integer}{comma}...]
## example: [branch_cutoff => 0.1=1]






## Configuration file end.
/;

    print $cfh $info;
    print STDOUT "\nCONTROL FILE: $contrfile has been created (".date().")\n\n";
    exit (1);
}


=head2 parse_config

  Usage: my %config = parse_config($control_file)

  Desc: parse configuration file and return a hash with the configuration 
        variables.

  Ret: %config, a hash with key = configuration_variable_name
                            val = configuration_data

  Args: $control_file, a file with the configuration variables

  Side_Effects: die if no control file is supplied

  Example: my %config = parse_config($control_file)

=cut

sub parse_config {
    my $file = shift ||
	die("ERROR: No control file was supplied to parse_config.\n");
    
    ## Define variables

    my %config = ( paths => {} );
    my $region = '';

    ## Define the regexp and the keys:

    my %match = (
	'CLUSTER_DATASOURCE'           => 'cl_sc',
	'CLUSTER_FILENAME'             => 'cl_fn',
	'CLUSTER_VALUES'               => 'cl_vl',
	'MEMBERSEQ_FILENAME'           => 'seq_fn',
	'MEMBERSTRAIN_FILENAME'        => 'str_fn',
	'FASTBLASTPARSER'              => 'cl_fast',
	'PATH_NAME'                    => 'pa_name',
	'HOMOLOGOUS_SEARCH_ARGUMENTS'  => 'hs_arg',
	'HOMOLOGOUS_SEARCH_DATASET'    => 'hs_dts',
	'HOMOLOGOUS_SEARCH_STRAIN'     => 'hs_str',
	'HOMOLOGOUS_SEARCH_FILTER'     => 'hs_fil',	
	'RUN_ALIGNMENT_PROGRAM'        => 'al_prg',
	'RUN_ALIGNMENT_ARGUMENTS'      => 'al_arg',
	'RUN_DISTANCE_FUNCTION'        => 'di_fun',
	'PRUNE_ALIGN_ARGUMENTS'        => 'pr_aln',
	'PRUNE_STRAINS_ARGUMENTS'      => 'pr_str',
	'PRUNE_OVERLAPS_ARGUMENTS'     => 'pr_ovl',
	'RUN_TREE_METHOD'              => 'tr_met',
	'RUN_TREE_ARGUMENTS'           => 'tr_arg',
	'REROOT_TREE'                  => 'tr_rer',
	'RUN_BOOTSTRAPPING'            => 'bo_run',
	'FILTER_BOOTSTRAPPING'         => 'bo_fil',
	'RUN_TOPOANALYSIS'             => 'tp_run',	
	);


    open my $cfh, '<', $file;
    while(<$cfh>) {
	chomp($_);

	my $line = $_;
	if ($_ =~ m/^(.+?)#/) {  ## Get the data between start and '#'
	    $line = $1;
	}
	
	unless ($line =~ m/^(#|)$/) {  ## Just ignore lines that start with #
                                       ## or are empty

	    if ($line =~ m/^<(\d*)>(.+?):\s+\[\s*(.+?)\s*\]/) {  ## ignore []
		my $path = $1;
		my $arg = $2;
		my $match_k = $match{$arg};
		my $match_v = $3;
		if (defined $match_k) {
		    if ($path =~ m/^\d+$/) {
			if (defined $config{paths}->{$path}) {
			    $config{paths}->{$path}->{$match_k} = $match_v;
			}
			else {
			    $config{paths}->{$path} = { $match_k => $match_v };
			}
		    }
		    else {
			$config{$match_k} = $match_v;
		    }
		}
		else {
		    print STDERR "\nWARNING: $arg isnt a valid conf. arg.\n";
		}
	    }
	}
    }

    return %config;
}

=head2 parse_clustervals

  Usage: my $cluster_val_href = parse_clustervals($cluster_val_line);

  Desc: parse cluster value line and return a hash reference with the
        right format to be used by parse_blastfile and fastparse_blastfile

  Ret: $cluster_val_href, a hash reference with key=blast_argument
                                                val=arrayref. with operator and 
                                                    value.

  Args: $cluster_val_line, a cluster line to be parsed.

  Side_Effects: Return undef if undef is used as argument.

  Example: my $cluster_val_href = parse_clustervals('evalue < 1e-100');

=cut

sub parse_clustervals {
    my $cl_line = shift;

    my $clval_href;

    if (defined $cl_line) {
	my %clval = ();
	my @sublines = split(/;/,  $cl_line);
	foreach my $subline (@sublines) {
	    if ($subline =~ m/^\s*(.+?)\s+(.+?)\s+(.+?)\s*$/) {
		$clval{$1} = [$2, $3];
	    }	
	}
	$clval_href = \%clval;
    }

    return $clval_href;
}


=head2 parse_homolog_args

  Usage: my $hom_val_href = parse_homolog_args(\%conf);

  Desc: parse homologous_search arguments and return a hash ref. usable by
        homologous_search function.

  Ret: $hom_val_href, a hash reference with homologous_search arguments

  Args: \%conf, a hash ref. with the configuration data.

  Side_Effects: None

  Example: my $hom_val_href = parse_homolog_args(\%conf);

=cut

sub parse_homolog_args {
    
}












####
1; #
####
