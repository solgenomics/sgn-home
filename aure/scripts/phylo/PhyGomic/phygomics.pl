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
if (!$opt_c && !$opt_i && !$opt_o && !$opt_h && !$opt_C && !$opt_S) {
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
## 1) Add the cluster source file and parse it.

print STDERR "1) PARSING CLUSTER SOURCE FILE (" .  date() . "):\n\n";

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
## 2) Parse the strain file and add to phygecluster

my $seq_args = { sequencefile => $indir . '/' . $conf{seq_fn} };
if ($opt_S) {
    $seq_args->{report_status} = 1;
}

print STDERR "2) PARSING SEQUENCE MEMBER FILE (" .  date() . "):\n\n";
$phyg->load_seqfile($seq_args);

my $seqn = 0;
my %cls = %{$phyg->get_clusters()};
foreach my $clid (keys %cls) {
    my @members = $cls{$clid}->get_members();
    $seqn += scalar(@members);
}

print STDERR "\tPARSED. $seqn sequences have been extracted.\n\n";


###################################################
## 3) Parse the strain file and add to phygecluster

my $str_args = { strainfile => $indir . '/' . $conf{str_fn} };
if ($opt_S) {
    $str_args->{report_status} = 1;
}

print STDERR "3) PARSING STRAIN FILE (" .  date() . "):\n\n";
$phyg->load_strainfile($str_args);
my $str_n = scalar(keys %{$phyg->get_strains()});

print STDERR "\tPARSED. $str_n strains have been extracted.\n\n";


##############################################################################
## 4) Now the pipeline will run different paths depending of the configuration
##    file. Each path will clone phygecluster and run a sepparate analysis.
##    at the end, it will integrate all of them in a phygestat object.


















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
## CLUSTERS CREATION ##########################################################
###############################################################################
##
## It can be based in two approaches: 
##   1) A sequence assembly. Each contig will be a cluster. 
##      Recomended for species close related where homologous genes can be
##      coassembled in the same contig.
##
##   2) A selfblast. Cluster will be created based in sequence match results.
##      Recomended for species where homologous genes cannot be coasssembled.


>CLUSTER_DATASOURCE:           []  
## Two possible sources, blast or ace files

>CLUSTER_FILENAME:             []  
## Name of the file to extract the clusters.

>CLUSTER_VALUES:               []  
## If cluster_datasource: [blast], a list separated by semicolons of arguments
## operators and values of different blast datatypes. Valid values are:
##   evalue, expect, frac_identical, frac_conserved, gaps, hsp_length, 
##   num_conserved, num_identical, score, bits, percent_identity, 
## If fastblastparser is enabled, the valid values are:
##   query_id, subject_id, percent_identity, align_length, mismatches, 
##   gaps_openings, q_start, q_end, s_start, s_end, e_value, bit_score
## Example: [ hsp_length > 100; percent_identity > 90, gaps < 5 ]

>FASTBLASTPARSER:              []
## Switch to use the fast blast parser (to enable add any character inside [])

>MEMBERSEQ_FILENAME:           []
## Sequence file for the members of the blast or ace file in fasta format.

>MEMBERSTRAIN_FILENAME:        []
## File with 2 columns: -f1 member_id, and -f2 strain name.



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

    my %config = ();
    my $region = '';

    ## Define the regexp and the keys:

    my %match = (
	'CLUSTER_DATASOURCE'           => 'cl_sc',
	'CLUSTER_FILENAME'             => 'cl_fn',
	'CLUSTER_VALUES'               => 'cl_vl',
	'MEMBERSEQ_FILENAME'           => 'seq_fn',
	'MEMBERSTRAIN_FILENAME'        => 'str_fn',
	'FASTBLASTPARSER'              => 'cl_fast',
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

	    if ($line =~ m/^>(.+?):\s+\[\s*(.+?)\s*\]/) {  ## ignore empty []
		my $match_k = $match{$1};                 
		if (defined $match_k) {
		    $config{$match_k} = $2;
		}
		else {
		    print STDERR "\nWARNING: $1 isnt a valid conf. argument.\n";
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















####
1; #
####
