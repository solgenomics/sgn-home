#!/usr/bin/perl

=head1 NAME

 phygomics.pl
 Script to run PhyGomics pipeline.

=cut

=head1 SYPNOSIS

 phygomics.pl [-h] -c config.file -i input.dir -o output.dir [-C] [-S] [-O]

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

=item -O

B<create_steps_files>         Create the files for each step (alignments, 
                              distances, trees, consensus and topologies)

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
        2.1) Search homologous using blast                [optional]
        2.2) Run alignments, 'clustalw' by default        [mandatory]*
        2.3) Run distances, 'kimura' by default           [mandatory for NJ]*
        2.4) Prune members, it will rerun 2.2 and 2.3     [optional]
        2.5) Run trees, 'ML' & 'phyml' by default         [mandatory]*
        2.6) Tree rerooting                               [optional]
        2.7) Run bootstrapping                            [optional]*
        2.8) Run topoanalysis                             [mandatory]*

     3) Data integration and analysis [1 cycle]
        3.1) Create graphs and tables.

    (see -C, configuration file for more details about data)
    (using -O argument will create output files for each of these steps)

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
use PhyGeConfig qw( write_conf_file read_conf_file );
use PhyGeTopo;
use PhyGeStats;
use PhyGePaml;

use Strain::AlleleIdentification qw/ identify_alleles /;

our ($opt_c, $opt_i, $opt_o, $opt_h, $opt_C, $opt_S, $opt_O);
getopts("c:i:o:hCSO");

if (!$opt_c && !$opt_i && !$opt_o && !$opt_C && !$opt_S && !$opt_O) {

    ## If no options are used, print help and exit

    print "There are n\'t any tags. Print help\n\n";
    help();
}
elsif ($opt_C) {

    ## If option C is used, print the configuration file and exit

    my $dir = getcwd;
    my $conffile = $dir . '/phygomics.conf';
    
    write_conf_file($conffile);

    print STDERR "\nEMPTY CONFIGURATION FILE: $conffile has been created ("
	. date() 
	. ")\n\n";
 
    exit (1);
}



## Mandatory arguments (die if they are not supplied)

my $err_m = "MANDATORY ARGUMENT ERROR:";

my $conffile = $opt_c 
    || die("\n$err_m -c <configuration_file> argument was not supplied.\n\n");

my $indir = $opt_i
    || die("\n$err_m -i <input_dir> argument was not supplied.\n\n");

my $outdir = $opt_o
    || die("\n$err_m -o <out_dir> argument was not supplied.\n\n");


## Check files and dirs

unless (-s $conffile) {
    die("\n$err_m $conffile configuration file doesnt exist or is empty.\n\n");
}
unless (-d $indir) {
    die("\n$err_m $indir input dir. doesnt exist.\n\n");
}


##############################################
#### PIPELINE EXECUTION  #####################
##############################################

############################################
## 0) Print start message and start to parse

print STDERR "\n\n";
print STDERR "*************************************************************";
print STDERR "\n** STARTING PhyGomic PIPELINE (" . date() . "):\n";
print STDERR "*************************************************************";
print STDERR "\n\n";

unless (-d $outdir) {
    print STDERR ("\nMESSAGE: Outdir doesnt exist.Creating outdir=$outdir\n\n");
    mkdir($outdir);
}

#############################################################################
## 0) Parse configuration file, check files and PhyGeCluster object creation.

print STDERR "0) READING CONFIGURATION FILE (" .  date() . "):\n\n";

my %conf = read_conf_file($conffile);

my $cv = scalar(keys %conf);
my %cpaths = %{$conf{path}};
foreach my $pi (keys %cpaths) {
    $cv += scalar(keys %{$cpaths{$pi}});
}
$cv--;

print STDERR "\tPARSED. $cv configuration arguments have been parsed.\n\n";

## Add indir and outdir to the conf.

$conf{indir} = $indir;
$conf{outdir} = $outdir;

filecheck(\%conf);


## Create an empty phygecluster object

my $phyg = PhyGeCluster->new();

 if ($opt_S) {
     $phyg->enable_reportstatus();
}


###############################################
## 1.1) Add the cluster source file and parse it.

print STDERR "1.1) PARSING CLUSTER SOURCE FILE (" .  date() . "):\n\n";

my $sc = $conf{general}->{input_filenames}->{source_filetype} || '';

my %clusters;
if ($sc eq 'blast') {

    ## Get the complete path for the blastfile

    my $blastfile = File::Spec->catfile($indir, 
					$conf{general}->{input_filenames}
					              ->{source}           );

    ## Check the cluster arguments

    my $fast = $conf{general}->{input_filenames}->{fastcluster_values};
    my $regular = $conf{general}->{input_filenames}->{cluster_values};

    if (defined $fast && defined $regular) {
	my $err = "ERROR: cluster_values and fastcluster_values arguments are";
	$err .= "used at the same time.\nPlease choose one of them and remove";
	$err .= "the other from the configuration file.\n\n";
	die($err);
    }
    elsif (defined $fast && $regular eq undef) {
    
	$fast->{blastformat} = 'blasttable';
	$fast->{blastfile} = $blastfile;
	
	if ($opt_S) {
	    $fast->{report_status} = 1;
	}

	%clusters = PhyGeCluster::fastparse_blastfile($fast);
    }
    elsif (defined $regular && $fast eq undef) {

	$regular->{blastformat} = 'blasttable';
	$regular->{blastfile} = $blastfile;
	
	if ($opt_S) {
	    $regular->{report_status} = 1;
	}

	%clusters = PhyGeCluster::parse_blastfile($regular);
    }
    else {
	die("ERROR: No cluster_value or fastcluster_values args. were used.\n");
    }
}
else {  ## By default it will .ace file
    
    my $acefile = File::Spec->catfile($indir, $conf{general}->{input_filenames}
					                    ->{source});
    
    my $ace_args = { acefile => $acefile };

    ## Add report_status if $opt_S

    if ($opt_S) {
	$ace_args->{report_status} = 1;
    }

    %clusters = PhyGeCluster::parse_acefile($ace_args);
}

## Set the cluster for PhyGeCluster

$phyg->set_clusters(\%clusters);

my $cl_v = scalar(keys %clusters);
print STDERR "\n\tPARSED. $cl_v clusters have been extracted.\n\n";

if ($opt_O) {

    ## Create a folder:
    my $clusterdir = $outdir . '/00_cluster_prepath';
    mkdir($clusterdir);

    ## Print files
    my $clbasename = $clusterdir . '/clusters.prepath';
    my %clfiles = $phyg->out_clusterfile({ rootname     => $clbasename,
					   distribution => 'single',
					 });
    
    my $clf_n = scalar(keys %clfiles);
    print STDERR "\tOPTION -O enabled: $clf_n prepath clusters files have been";
    print STDERR " created with basename:\n\t$clbasename\n\n";
}


###################################################
## 1.2) Parse the strain file and add to phygecluster

my $memberseq = File::Spec->catfile($indir, $conf{general}->{input_filenames}
                                                          ->{memberseq}); 

my $seq_args = { sequencefile => $memberseq };
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

print STDERR "\n\tPARSED. $seqn sequences have been extracted.\n\n";


###################################################
## 1.3) Parse the strain file and add to phygecluster

my $memberstr = File::Spec->catfile($indir, $conf{general}->{input_filenames}
                                                          ->{memberstrain}); 

my $str_args = { strainfile => $memberstr };
if ($opt_S) {
    $str_args->{report_status} = 1;
}

print STDERR "1.3) PARSING STRAIN FILE (" .  date() . "):\n\n";
$phyg->load_strainfile($str_args);
my $str_n = scalar(keys %{$phyg->get_strains()});

print STDERR "\n\tPARSED. $str_n strains have been extracted.\n\n";


##############################################################################
## 2) Now the pipeline will run different paths depending of the configuration
##    file. Each path will clone phygecluster and run a sepparate analysis.
##    at the end, it will integrate all of them in a phygestat object.

my %paths = %{$conf{path}};
my $pn = scalar(keys %paths);

print STDERR "\n** $pn PATHS has been defined (".date()."):\n\n";

my %topostats = ();

foreach my $path_idx (sort {$a <=> $b} keys %paths) {
    my %pargs = %{$paths{$path_idx}};
    my $phname = $pargs{pa_name} || $path_idx; 

    print STDERR "\t$path_idx ==> INIT. PATH=$phname (" .  date() . "):\n\n";

    my $pathoutdir = $outdir . '/path_' . $path_idx;
    if ($opt_O) {
	mkdir($pathoutdir);
    }


    ## 2.0) Clone the PhyGeCluster #######################################

    print STDERR "\t\t2.0) CLONING CLUSTER DATA (" .  date() . "):\n\n";

    my $paphyg = $phyg->clone();

    if (defined $opt_S) {
	print STDERR "\n\n";
    }

    if ($opt_O) {
	
	## Create a folder:
	my $alndir = $pathoutdir . '/00_alignments';
	mkdir($alndir);
	
	## Print files
	my $alnbase = $alndir . '/alignments';
	my %alnfiles = $paphyg->out_alignfile({ 
	    'rootname'     => $alnbase,
	    'distribution' => 'single',
	    'format'       => 'clustalw',
						 });
	
	my $aln_n = scalar(keys %alnfiles);
	print STDERR "\t\t\tOPTION -O enabled: ";
	print STDERR "$aln_n alignment files have been created";
	print STDERR " with basename:\n\t\t\t$alnbase\n\n";
    }

    ## 2.1) Homologous search ############################################

    print STDERR "\t\t2.1) HOMOLOGOUS SEARCH (" .  date() . "):\n\n";

    if (exists $conf{path}->{$path_idx}->{homologous_search}) {
	
	## Get the arguments and set the defaults
	
	my %homsearch_args = %{$conf{path}->{$path_idx}->{homologous_search}};
	$homsearch_args{indir} = $indir;

	my %homsearch = parse_homolog_args(\%homsearch_args);

	## And run homologous search
	
	$paphyg->homologous_search(\%homsearch);
	
	## Count homologous
	my $hom_c = 0;

	my %cls = %{$paphyg->get_clusters()};
	my %str = %{$paphyg->get_strains()};
	foreach my $clid (keys %cls) {
	    my @membs = $cls{$clid}->get_members();
	    foreach my $memb (@membs) {
		my $m_id = $memb->id();
		if (defined $m_id) {
		    if ($str{$m_id} eq $homsearch_args{homologous_strain}) {
			$hom_c++;
		    }
		}
	    }
	}

	print STDERR "\t\t\tDONE. $hom_c homologous have been assigned.\n\n";

	if ($opt_O) {

	    ## Create a folder:
	    my $homologdir = $pathoutdir . '/01_cluster_homologs';
	    mkdir($homologdir);
	    
	    ## Print files
	    my $hoclbase = $homologdir . '/clusters.homologs';
	    my %hoclfiles = $paphyg->out_clusterfile({ 
		rootname     => $hoclbase,
		distribution => 'single',
						     });

	    my $hof_n = scalar(keys %hoclfiles);
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "$hof_n clusters with homolog files have been created";
	    print STDERR " with basename:\n\t\t\t$hoclbase\n\n";
	}

    }
    else {
	
	print STDERR "\t\t\tNo HOMOLOGOUS SEARCH arguments were set:\n";
	print STDERR "\t\t\t\tHOMOLOGOUS_SEARCH_STRAIN\n";
	print STDERR "\t\t\t\tHOMOLOGOUS_SEARCH_DATASET\n";
	print STDERR "\t\t\t\tSKIP STEP.\n\n";
    }

    ## 2.2) Run alignments ################################################

    print STDERR "\t\t2.2) RUN ALIGNMENTS (" .  date() . "):\n\n";

    my $align_args;

    if (defined $conf{path}->{$path_idx}->{alignment}) {
	
	 $align_args = $conf{path}->{$path_idx}->{alignment};		
    }
    else {
	 print STDERR "\t\t\t USING DEFAULT ALIGNMENT (clustalw)\n\n";
	 $align_args = { program    => 'clustalw', 
			 parameters => { quiet => undef, matrix => 'BLOSUM'},
	 };
    }


    ## Run the alignment

    my @failed_align = $paphyg->run_alignments($align_args);

    my $alig_n = 0;
    my %a_cls = %{$paphyg->get_clusters()};
    foreach my $a_clid (keys %a_cls) {
	if (defined $a_cls{$a_clid}->alignment()) {
	    $alig_n++;
	}
    }
    
    print STDERR "\n\t\t\tDONE. $alig_n alignments have been added.\n\n";
    
    if ($opt_O) {
	
	## Create a folder:
	my $alndir = $pathoutdir . '/02_alignments';
	mkdir($alndir);
	
	## Print files
	my $alnbase = $alndir . '/alignments';
	my %alnfiles = $paphyg->out_alignfile({ 
	    'rootname'     => $alnbase,
	    'distribution' => 'single',
	    'format'       => 'clustalw',
						 });
	
	my $aln_n = scalar(keys %alnfiles);
	print STDERR "\t\t\tOPTION -O enabled: ";
	print STDERR "$aln_n alignment files have been created";
	print STDERR " with basename:\n\t\t\t$alnbase\n\n";

	if (scalar(@failed_align) > 0) {
	    open my $ffh_aln, '>', $alndir . '/failed_alignments.err';
	    foreach my $f_align (@failed_align) {
		print $ffh_aln "$f_align\n";
	    }
	    close($ffh_aln);
	}

	## Print alignments stats
	my $alnstatbase = $alndir . '/stats_aln.tab';
	my $stats = $paphyg->out_alignstats({basename => $alnstatbase});
	print STDERR "\t\t\tOPTION -O enabled: ";
	print STDERR "a stats alignment file have been created";
	print STDERR " with filename:\n\t\t\t$stats\n\n";
    }

    ## 2.3) Prune alignments ##############################################

    print STDERR "\t\t2.3) PRUNE BY ALIGNMENTS (" .  date() . "):\n\n";
    
    my $cls_n = scalar( keys %{$paphyg->get_clusters()});

    if (defined $conf{path}->{$path_idx}->{prune_alignments}) {

	my $prun_align_args = $conf{path}->{$path_idx}->{prune_alignments};

	my %rem_cls = $paphyg->prune_by_align($prun_align_args);
	my $rem_c = scalar(keys %rem_cls);

	print STDERR "\n\t\t\t\tDONE. $rem_c clusters have been removed.\n\n";
	
	print STDERR "\t\t\t\tRERUNING alignments.\n\n";
	my @failed_pr_aln = $paphyg->run_alignments($align_args);
 
        if ($opt_O) {
	
	    ## Create a folder:
	    my $pru_alndir = $pathoutdir . '/03_pruned_alignments';
	    mkdir($pru_alndir);
	
	    ## Print files
	    my $prubase = File::Spec->catfile($pru_alndir, 'pruned_alignments');
	    my %prunalnfiles = $paphyg->out_alignfile({ 
		'rootname'     => $prubase,
		'distribution' => 'single',
		'format'       => 'clustalw',
						  });
	    
	    my $prunaln_n = scalar(keys %prunalnfiles);
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "$prunaln_n alignment files have been created after ";
	    print STDERR "pruning by aligment ";
	    print STDERR " with basename:\n\t\t\t$prubase\n\n";
	    
	    ## Print failed alignments
	    if (scalar(@failed_pr_aln) > 0) {
		open my $ffh_aln, '>', $pru_alndir . '/failed_alignments.err';
		foreach my $f_align (@failed_pr_aln) {
		    print $ffh_aln "$f_align\n";
		}
		close($ffh_aln);
	    }
	    
	    ## Print alignments stats
	    my $pru_alnstat = $pru_alndir . '/stats_aln_after_pruning.tab';
	    my $pru_stats = $paphyg->out_alignstats({basename => $pru_alnstat});
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "a stats alignment file have been created";
	    print STDERR " with filename:\n\t\t\t$pru_stats\n\n";
	}

	$cls_n = scalar( keys %{$paphyg->get_clusters()});
	print STDERR "\t\t\t$cls_n clusters remain after prunning.\n\n";
    }
    else {
	print STDERR "\t\t\tNO PRUNE ALIGNMENTS ARGUMENT WAS USED.\n";
	print STDERR "\t\t\tSKIPPING PRUNE BY ALIGNMENTS.\n\n";
    }

    ## 2.4) Prune by overlaps ###############################################

    print STDERR "\t\t2.4) PRUNE BY OVERLAPINGS (" .  date() . "):\n\n";

    if (defined $conf{path}->{$path_idx}->{prune_overlaps}) {
    
	my %pruovl_args = %{$conf{path}->{$path_idx}->{prune_overlaps}};

	## Format the arguments

	if (exists $pruovl_args{composition}) {
	    my %comp = ();
	    if ($pruovl_args{composition} =~ m/(,|;)/) {
		my @comp_list = split(/$1/, $pruovl_args{composition});
		foreach my $comp_element (@comp_list) {
		    $comp_element =~ s/\s+//g;
		    my @comp_fact = split(/=/, $comp_element);
		    if (scalar(@comp_fact) != 2) {
			die("ERROR: wrong format: prune_overlaps>composition");
		    } 
		    $comp{$comp_fact[0]} = $comp_fact[1]; 
		}
	    }
	    $pruovl_args{composition} = \%comp;
	}
	if (exists $pruovl_args{filter}) {
	    my %filt = ();
	    my $sep = ',';
	    if ($pruovl_args{filter} =~ m/(,|;)/) {
		$sep = $1;
	    }
	    my @filt_list = split(/$sep/, $pruovl_args{filter});
	    foreach my $filt_element (@filt_list) {

		$filt_element =~ s/\s+//g;  ## Remove spaces

		my $sep2 = '=';
		if ($filt_element =~ m/(=>|->|=|>)/) {
		    $sep2 = $1;
		}
		my @filfact = split(/$sep2/, $filt_element);
		if (scalar(@filfact) != 2) {
		    my $line = join("$sep2", @filfact);
		    my $test = scalar(@filfact);
		    die("ERROR: wrong format: prune_overlaps>filter ($line)\n");
		} 
		$filt{$filfact[0]} = $filfact[1]; 
	    }
	    $pruovl_args{filter} = \%filt;
	}
	

	my ($rem_o_cls, $rem_o_mbs) = $paphyg->prune_by_overlaps(\%pruovl_args);
	my $rem_o_c = scalar(keys %{$rem_o_cls});
	my $rem_o_m = scalar(keys %{$rem_o_mbs});

	print STDERR "\n\t\t\t\tDONE. $rem_o_c clusters and $rem_o_m members ";
	print STDERR "have been removed.\n\n";

	my @failed_pr_aln = ();
	if ($pruovl_args{trim} == 0) {
	    print STDERR "\t\t\t\tRERUNING alignments.\n\n";
	    @failed_pr_aln = $paphyg->run_alignments($align_args);
	}

	if ($opt_O) {
	    
	    ## Create a folder:
	    my $pru_alndir = $pathoutdir . '/04_pruned_overlaps';
	    mkdir($pru_alndir);
	
	    ## Print files
	    my $prubase = File::Spec->catfile($pru_alndir, 'pruned_overlaps');
	    my %prunalnfiles = $paphyg->out_alignfile({ 
		'rootname'     => $prubase,
		'distribution' => 'single',
		'format'       => 'clustalw',
						  });
	    
	    my $prunaln_n = scalar(keys %prunalnfiles);
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "$prunaln_n alignment files have been created after ";
	    print STDERR "pruning by overlaps";
	    print STDERR " with basename:\n\t\t\t$prubase\n\n";
	    
	    ## Print failed alignments
	    if (scalar(@failed_pr_aln) > 0) {
		open my $ffh_aln, '>', $pru_alndir . '/failed_alignments.err';
		foreach my $f_align (@failed_pr_aln) {
		    print $ffh_aln "$f_align\n";
		}
		close($ffh_aln);
	    }
	    
	    ## Print alignments stats
	    my $pru_alnstat = $pru_alndir . '/stats_aln_after_pruning.tab';
	    my $pru_stats = $paphyg->out_alignstats({basename => $pru_alnstat});
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "a stats alignment file have been created";
	    print STDERR " with filename:\n\t\t\t$pru_stats\n\n";
	}

	$cls_n = scalar( keys %{$paphyg->get_clusters()});
	print STDERR "\t\t\t$cls_n clusters remain after prunning.\n\n";
    }
    else {
	print STDERR "\t\t\tNO PRUNE OVERLAP ARGUMENT WAS USED.\n";
	print STDERR "\t\t\tSKIPPING PRUNE BY OVERLAPS.\n\n";
    }

    ## 2.5) Run distances  ###############################################

    ## It will run distances only if tree>method = ML OR exists prune_distances

    print STDERR "\t\t2.5) RUN DISTANCES (" .  date() . "):\n\n";

    my $dist_args;

    ## Get the conditions to run distances

    my $run_distances = 0;

    if (defined $conf{path}->{$path_idx}->{prune_distances}) {
	$run_distances = 1;
    }
    elsif ($conf{path}->{$path_idx}->{tree}->{method} eq 'ML') {
	$run_distances = 1;
    }
    elsif (defined $conf{path}->{$path_idx}->{distance}) {
	$run_distances = 1;
    }

    if ($run_distances == 1) {
	if (defined $conf{path}->{$path_idx}->{distance}) {
	
	    $dist_args = $conf{path}->{$path_idx}->{distance};	
	}
	else {
	
	    print STDERR "\t\t\t USING DEFAULT DISTANCE (Kimura)\n\n";
	    $dist_args = { method => 'Kimura', quiet => 1};
	}

	## And now use the Run distances method

	my @failed_dist = $paphyg->run_distances($dist_args);

	my $dist_n = scalar( keys %{$paphyg->get_distances()});    
	print STDERR "\n\t\t\tDONE. $dist_n distance matrices have been ";
	print STDERR "added.\n\n";

	if ($opt_O) {

	    ## Create a folder:
	    my $disdir = $pathoutdir . '/05_distances';
	    mkdir($disdir);
	
	    ## Print files
	    my $disbase = $disdir . '/distances';
	    my %disfiles = $paphyg->out_distancefile({ 
		'rootname'     => $disbase,
		'distribution' => 'single',
		'format'       => 'phylip',
						     });
	    
	    my $dis_n = scalar(keys %disfiles);
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "$dis_n distance files have been created";
	    print STDERR " with basename:\n\t\t\t$disbase\n\n";
	
	    if (scalar(@failed_dist) > 0) {
		open my $ffh_dist, '>', $disdir . '/failed_distances.err';
		foreach my $f_dist (@failed_dist) {
		    print $ffh_dist "$f_dist\n";
		}
		close($ffh_dist);
	    }	
	}
    }
    else {
	print STDERR "\t\t\tNO DISTANCE ARGUMENTS (OR RELATED) WAS USED.\n";
	print STDERR "\t\t\tSKIPPING RUN DISTANCES.\n\n";
    }
	

    ## 2.6) Prune by distances ###########################################

    print STDERR "\t\t2.6) PRUNE BY DISTANCES (" .  date() . "):\n\n";

    if (defined $conf{path}->{$path_idx}->{prune_distances}) {
    
	my %prudis_args = %{$conf{path}->{$path_idx}->{prune_distances}};

	## Format the arguments (need to be refactored)

	if (exists $prudis_args{composition}) {
	    my %comp = ();
	    if ($prudis_args{composition} =~ m/(,|;)/) {
		my @comp_list = split(/$1/, $prudis_args{composition});
		foreach my $comp_element (@comp_list) {
		    $comp_element =~ s/\s+//g;
		    my @comp_fact = split(/=/, $comp_element);
		    if (scalar(@comp_fact) != 2) {
			die("ERROR: wrong format: prune_distances>composition");
		    } 
		    $comp{$comp_fact[0]} = $comp_fact[1]; 
		}
	    }
	    $prudis_args{composition} = \%comp;
	}
	if (exists $prudis_args{min_distance}) {
	    my @dis = ();
	    if ($prudis_args{min_distance} =~ m/(,|;)/) {
		my @dis_list = split(/$1/, $prudis_args{min_distance});
		foreach my $dis_element (@dis_list) {
		    $dis_element =~ s/\s+//g;
		    
		    if ($dis_element =~ m/(=|-)/) {
			my $sep = $1;
			my @dis_fact = split(/$sep/, $dis_element);
			if (scalar(@dis_fact) != 2) {
			    die("ERROR: wrong format: min_distance");
			} 
			push @dis, \@dis_fact; 
		    }
		}
		$prudis_args{min_distances} = \@dis;
	    }
	}
	if (exists $prudis_args{max_distance}) {
	    my @dis = ();
	    if ($prudis_args{max_distance} =~ m/(,|;)/) {
		my @dis_list = split(/$1/, $prudis_args{max_distance});
		foreach my $dis_element (@dis_list) {
		    $dis_element =~ s/\s+//g;
		    
		    if ($dis_element =~ m/(=|-)/) {
			my $sep = $1;
			my @dis_fact = split(/$sep/, $dis_element);
			if (scalar(@dis_fact) != 2) {
			    die("ERROR: wrong format: max_distance");
			} 
			push @dis, \@dis_fact; 
		    }
		}
		$prudis_args{max_distances} = \@dis;
	    }
	}
	    
	## Run the pruning
	
	my ($rem_o_cls, $rem_o_mbs) = $paphyg->prune_by_strains(\%prudis_args);
	my $rem_o_c = scalar(keys %{$rem_o_cls});
	my $rem_o_m = scalar(keys %{$rem_o_mbs});

	print STDERR "\n\t\t\t\tDONE. $rem_o_c clusters and $rem_o_m members ";
	print STDERR "have been removed.\n\n";

	## Rerun alignment and distances

	print STDERR "\t\t\t\tRERUNING alignments.\n\n";
	my @failed_pr_aln = $paphyg->run_alignments($align_args);
	
	print STDERR "\n\t\t\t\tRERUNING distances.\n\n";
	my @failed_pr_dist = $paphyg->run_distances($dist_args);

	if ($opt_O) {

	    ## Create a folder:
	    my $prdir = $pathoutdir . '/06_pruned_distances';
	    mkdir($prdir);
	
	    ## Print files
	    my $pr_alnbase = File::Spec->catfile( $prdir, 
						  'dist_pruned.alignments');
	    
	    my %pr_alnfiles = $paphyg->out_alignfile({ 
		'rootname'     => $pr_alnbase,
		'distribution' => 'single',
		'format'       => 'clustalw',
						     });
	
	    ## Print alignments stats
	    my $pralnstatbase = $prdir . '/dist_pruned.stats_alignments.tab';
	    my $prstats = $paphyg->out_alignstats({basename => $pralnstatbase});
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "a stats alignment file have been created";
	    print STDERR " with filename:\n\t\t\t$prstats\n\n";

	    my $pr_disbase = $prdir . '/dist_pruned.distances';
	    my %pr_disfiles = $paphyg->out_distancefile({ 
		'rootname'     => $pr_disbase,
		'distribution' => 'single',
		'format'       => 'phylip',
							});
	    my $pr_aln_n = scalar(keys %pr_alnfiles);
	    my $pr_dis_n = scalar(keys %pr_disfiles);
	    print STDERR "\t\t\tOPTION -O enabled: $pr_aln_n align. and ";
	    print STDERR "$pr_dis_n distance files (after prunning) have been ";
	    print STDERR "created with basename:\n\t\t\t$pr_alnbase\n\t\t\t";
	    print STDERR "$pr_disbase\n\n";
    
	    if (scalar(@failed_pr_aln) > 0) {
		open my $ffh_praln, '>', $prdir. '/dist_pruned.failed_alig.err';
		foreach my $f_praln (@failed_pr_aln) {
		    print $ffh_praln "$f_praln\n";
		}
		close($ffh_praln);
	    }
	    if (scalar(@failed_pr_dist) > 0) {
		open my $ffh_prdist, '>', $prdir.'/dist_pruned.failed_dist.err';
		foreach my $f_prdist (@failed_pr_dist) {
		    print $ffh_prdist "$f_prdist\n";
		}
		close($ffh_prdist);
	    }
	}

	$cls_n = scalar( keys %{$paphyg->get_clusters()});
	print STDERR "\t\t\t$cls_n clusters remain after prunning.\n\n";
    }
    else {
    	print STDERR "\t\t\tNO PRUNE DISTANCES ARGUMENT WAS USED.\n";
	print STDERR "\t\t\tSKIPPING PRUNE BY DISTANCES.\n\n";
    }
    
    ## Now it will continue, only if $cls_n > 0
    
    my %topotypes;

    if ($cls_n > 0) {
    
	## 2.7) Run tree methods ##########################################
	
	print STDERR "\t\t2.7) RUN TREES (" .  date() . "):\n\n";

	my @failed_trees = ();
	
	my %tree_args = %{$conf{path}->{$path_idx}->{tree}};
	my ($nj_tree, $ml_tree);

	if ($tree_args{method} =~ m/^(NJ|UPGMA)$/) {

	    my %nj_args = ( type => $tree_args{method} );
	    if (defined $tree_args{nj_parameters}) {
		foreach my $nj_arg (keys %{$tree_args{nj_parameters}}) {
		    $nj_args{$nj_arg} = $tree_args{nj_parameters}->{$nj_arg};
		}
	    }
	    if (defined $tree_args{outgroup}) {
		$nj_args{outgroup_strain} = $tree_args{outgroup};
	    }
	    if (defined $tree_args{quiet}) {
		$nj_args{quiet} = $tree_args{quiet};
	    }
	    
	    @failed_trees = $paphyg->run_njtrees(\%nj_args);
	    $nj_tree = \%nj_args;
	}
	else {
	    
	    my %ml_args = ( $tree_args{program} => {} );
	    if (defined $tree_args{ml_parameters}) {
		foreach my $ml_arg (keys %{$tree_args{ml_parameters}}) {
		    my $val = $tree_args{ml_parameters}->{$ml_arg};
		    $ml_args{$tree_args{program}}->{$ml_arg} = $val;
		}
	    }
	    
	    if (defined $tree_args{outgroup}) {
		$ml_args{outgroup_strain} = $tree_args{outgroup};
	    }
	    if (defined $tree_args{quiet}) {
		$ml_args{$tree_args{program}}->{quiet} = $tree_args{quiet};
	    }

	    @failed_trees = $paphyg->run_mltrees(\%ml_args);
	    $ml_tree = \%ml_args;
	}

	## Get number of trees

	my $trees_n = 0;
	my %t_cls = %{$paphyg->get_clusters()};
	foreach my $clid (keys %t_cls) {
	    if (defined $t_cls{$clid}->tree()) {
		$trees_n++;
	    }
	}
	
	print STDERR "\n\t\t\tDONE. $trees_n trees have been added.\n\n";

        if ($opt_O) {

	    ## Create a folder:
	    my $treedir = $pathoutdir . '/07_trees';
	    mkdir($treedir);
	
	    ## Print files
	    my $treebase = $treedir . '/tree';
	    my %treefiles = $paphyg->out_treefile({ 
	         'rootname'     => $treebase,
	         'distribution' => 'single',
	         'format'       => 'newick',
		   				    });

	    my $tree_n = scalar(keys %treefiles);
	    print STDERR "\t\t\tOPTION -O enabled: ";
	    print STDERR "$tree_n tree files have been created ";
	    print STDERR "with basename:\n\t\t\t$treebase\n\n";
        
	    open my $ffh_trees, '>', $treedir . '/failed_trees.err';
	    foreach my $f_tree (@failed_trees) {
		print $ffh_trees "$f_tree\n";
	    }
	    close($ffh_trees);
	}


	## 2.8) Reroot tree ###############################################

	print STDERR "\t\t2.8) REROOTING TREES (" .  date() . "):\n\n";
	if (defined $conf{path}->{$path_idx}->{reroot_tree}) {
		    
	    my $reroot_args = $conf{path}->{$path_idx}->{reroot_tree};
	    $paphyg->reroot_trees($reroot_args);

	    my $rtrees_n = 0;
	    my %rt_cls = %{$paphyg->get_clusters()};
	    foreach my $clid (keys %rt_cls) {
		if (defined $rt_cls{$clid}->tree()) {
		    $rtrees_n++;
		}
	    }
	
	    print STDERR "\t\t\tDONE. $rtrees_n trees have been rerooted.\n\n";

            if ($opt_O) {

	        ## Create a folder:
	        my $retreedir = $pathoutdir . '/08_trees_rerooted';
	        mkdir($retreedir);
	
	        ## Print files
	        my $retreebase = $retreedir . '/tree_rerooted';
	        my %retreefiles = $paphyg->out_treefile({ 
	            'rootname'     => $retreebase,
	            'distribution' => 'single',
	            'format'       => 'newick',
		      				      });

  	        my $retree_n = scalar(keys %retreefiles);
	        print STDERR "\t\t\tOPTION -O enabled: ";
	        print STDERR "$retree_n rerooted tree files have been created ";
	        print STDERR "with basename:\n\t\t\t$retreebase\n\n";
            }
	}
	else {
	    print STDERR "\t\t\tNo Arguments. SKIPPING REROOT TREES.\n\n";
	}


	## 2.9) Run bootstrapping ##########################################

	print STDERR "\t\t2.9) RUN BOOTSTRAPPING (" .  date() . "):\n\n";
	if (defined $conf{path}->{$path_idx}->{bootstrapping}) {
	    
	    my %bootarg = %{$conf{path}->{$path_idx}->{bootstrapping}};
	    
	    ## Process the arguments

	    my $boots_args = { 
		run_bootstrap => { datatype => 'Sequence' },
		run_distances => $dist_args,
		run_consensus => {}
	    };
	    if (defined $nj_tree) {
		$boots_args->{run_njtrees} = $nj_tree;
	    }
	    else {
		$boots_args->{run_mltrees} = $ml_tree;
	    }

	    if ($bootarg{replicates}) {
		$boots_args->{run_bootstrap}->{replicates} = 
		    $bootarg{replicates};
	    }
	    else {
		$boots_args->{run_bootstrap}->{replicates} = 100;
	    }
            
	    if (defined $bootarg{quiet}) {
		$boots_args->{run_bootstrap}->{quiet} = $bootarg{quiet};
		$boots_args->{run_consensus}->{quiet} = $bootarg{quiet};
	    }
	    
	    if (defined $bootarg{outgroup}) {
		$boots_args->{run_consensus}->{outgroup_strain} = 
		    $bootarg{outgroup};
	    }	    
	    if (defined $bootarg{midpoint}) {
		$boots_args->{run_consensus}->{root_by_midpoint} = 
		    $bootarg{midpoint};
	    }
	    if (defined $bootarg{normalized}) {
		$boots_args->{run_consensus}->{normalized} = 
		    $bootarg{normalized};
	    }
	    

	    $paphyg->run_bootstrapping($boots_args);
	    
	    my %boots = %{$paphyg->get_bootstrapping()};
	    my $bootsn = scalar(keys %boots);
	    
	    print STDERR "\n\t\t\tDONE. $bootsn bootstraps have been added\n\n";
	    
	    if (defined $bootarg{filter}) {
		
		print STDERR "\t\t\tPRUNING BY BOOTSTRAP VALUES.\n\n";
		my %rmcls_boots = $paphyg->prune_by_bootstrap($bootarg{filter});
	    
		my $rboots_n = scalar(keys %rmcls_boots);
	    	
		print STDERR "\n\t\t\tDONE. $rboots_n trees have been removed ";
		print STDERR "according the bootstrapping values.\n\n";
	    }
	    
            if ($opt_O) {

	        ## Create a folder:
	        my $bootsdir = $pathoutdir . '/09_bootstrap_consensus';
	        mkdir($bootsdir);
		
	        ## Print files
	        my $bootsbase = $bootsdir . '/consensus';
	        my %bootsfiles = $paphyg->out_bootstrapfile(
		    { 
			'rootname'     => $bootsbase,
			'distribution' => 'single',
			'format'       => 'newick',	
		    }
		    );

	        my $boots_n = scalar(keys %bootsfiles);
	        print STDERR "\t\t\tOPTION -O enabled: ";
	        print STDERR "$boots_n bootstrapping consensus files have been";
	        print STDERR " created with basename:\n\t\t\t$bootsbase\n\n";
             }	    
	}
	else {
	    print STDERR "\t\t\tNo Arguments. SKIPPING BOOTSTRAPPING.\n\n";
	}

	## 2.10) Run phygetopo ###########################################

	print STDERR "\t\t2.10) RUN TOPOANALYSIS (" .  date() . "):\n\n";
		    
	my $seqfams_href = $paphyg->get_clusters();
	my $strains_href = $paphyg->get_strains();

	my $phygetopo = PhyGeTopo->new({ seqfams => $seqfams_href,
					 strains => $strains_href,
				       });

	my $topo_args;
	if (defined $conf{path}->{$path_idx}->{topoanalysis}) {
	    my %topoan = %{$conf{path}->{$path_idx}->{bootstrapping}};
	    if (defined $topoan{branch_cutoffs}) {
		$topo_args->{branch_cutoffs} = {};

		my $sep = ',';
		if ($topoan{branch_cutoffs} =~ m/(,|;)/) {
		    $sep = $1;
		}

		my @bra = split(/$sep/, $topoan{branch_cutoffs});
		foreach my $bra (@bra) {
		    $bra =~ s/\s+//;
		 
		    my $sep2 = '-';
		    if ($bra =~ m/(=>|->|-|=)/) {
			$sep2 = $1;
		    }
		    my @intbra = split(/$sep2/, $bra);
		    $topo_args->{branch_cutoff}->{$intbra[0]} = $intbra[1];
		}		
	    }
	}

	%topotypes = $phygetopo->run_topoanalysis($topo_args);
	
	$topostats{$phname} = $phygetopo;

	my $topot_n = scalar(keys %topotypes);
	print STDERR "\t\t\tDONE. $topot_n topotypes have been created\n\n";

        if ($opt_O) {

	    ## Create a folder:
	    my $topodir = $pathoutdir . '/10_topologies';
            mkdir($topodir);
	
            ## Print files
            my $topobase1 = $topodir . '/topology_analysis';
            my $topofile1 = $phygetopo->out_topoanalysis( 
                                                      {basename => $topobase1 }
                                                       );

            my $topobase2 = $topodir . '/topology_composition';
            my $topofile2 = $phygetopo->out_topocomposition(
                                                     {basename => $topobase2}
                                                      );

            print STDERR "\t\t\tOPTION -O enabled: ";
            print STDERR "topoanalysis and topocomposition files have been";
            print STDERR " created with basename:\n\t\t\t$topobase1\n";
            print STDERR "\t\t\t$topobase2\n\n";

         }	 
    }

    ## 2.11) Run allele_identification  ####################################

    print STDERR "\t\t2.11) RUN ALLELE IDENTIFICATION (" .  date() . "):\n\n";

    my %seqfams = %{$paphyg->get_clusters()};
    my %strains = %{$paphyg->get_strains()};
    my %alleles_strains = ();

    if (defined $conf{path}->{$path_idx}->{allele_identification} ) {

	my %alleleident = %{$conf{path}->{$path_idx}->{allele_identification}};

	if (defined $alleleident{target} && defined $alleleident{parents}) {

	    my %parents = ();
	    if ($alleleident{parents} =~ m/(,|;)/) {
		my @parents = split(/$1/, $alleleident{parents});
		
		foreach my $parent (@parents) {
		    
		    if ($parent =~ m/(=>|->|=|-)/) {
			my @pair = split(/$1/, $parent);
			$parents{$pair[0]} = $pair[1];
		    }
		}
	    }
	    if (scalar(keys %parents) < 2) {
		die("ERROR: Less than 2 parents was used for identify_alleles");
	    }

	    if (scalar(keys %seqfams) > 0) { 

		## Create a folder and file:
		my $alleledir = $pathoutdir . '/11_allele_identification';
		mkdir($alleledir);

		my $allelefile = File::Spec->catfile($alleledir, 'alleles.tab');
		open my $alfh, '>', $allelefile;

		my $allele_n = 0;

		foreach my $clid (sort keys %seqfams) {
	    
		    my $famtree = $seqfams{$clid}->tree();
		    my $famalign = $seqfams{$clid}->alignment();

		    if ($opt_S) {
			print STDERR "\t\t\t\tProcessing cluster_id: $clid  \r";
		    }

		    
		    if (defined $famtree && defined $famalign) {
			
			my %alident = (
			    strains   => \%strains,
			    alignment => $famalign,
			    tree      => $famtree,
			    parents   => \%parents,
			    target    => $alleleident{target},
			    );

			if (defined $alleleident{filter_length}) {
			    $alident{filter}->{length} = 
				$alleleident{filter_length};
			}
			if (defined $alleleident{filter_identity}) {
			    $alident{filter}->{identity} = 
				$alleleident{filter_identity};
			}

			my %alleles = identify_alleles(\%alident);
			if (scalar( keys %alleles) > 0) {
			    foreach my $mb_id (sort keys %alleles) {
				
				$allele_n++;
				print $alfh "$clid\t$mb_id\t$alleles{$mb_id}\n";
				$alleles_strains{$mb_id} = $alleles{$mb_id};
			    }
			}
		    }
		}
		close($alfh);
		print STDERR "\n\n\t\t\t\t$allele_n alleles have been ";
		print STDERR "identified\n\n";
	    }	
	    else {
		print STDERR "\t\tAll the clusters have been removed in ";
		print STDERR "previous steps.\n";
		print STDERR "\t\tSKIPPING ALLELE IDENTIFICATION\n\n";
	    }
	}
	else {
	    print STDERR "\t\tNo target or/and parents were used for allele";
	    print STDERR "_identification.\n\t\tSKIPPING ALLELE IDENTIFICATION";
	    print STDERR "\n\n";
	}
    }
    else {
	print STDERR "\t\tNo allele_identification arguments were used.\n";
	print STDERR "\t\tSKIPPING ALLELE IDENTIFICATION\n\n";
    }

    ## 2.12) Run codeml for dN/dS  ####################################

    print STDERR "\t\t2.12) RUN CODEML ANALYSIS (" .  date() . "):\n\n";

    if (defined $conf{path}->{$path_idx}->{codeml_analysis} ) {

	my %codeml = %{$conf{path}->{$path_idx}->{codeml_analysis}};

	## It will replace the strains with the alleles when in the 
	## cases that they have been identified

	my %str_al = ();
	my %str = %{$paphyg->get_strains()};

	foreach my $member_id (sort keys %str) {
	    if (defined $alleles_strains{$member_id}) {
		$str_al{$member_id} = $alleles_strains{$member_id};
	    }
	    else {
		$str_al{$member_id} = $str{$member_id};
	    }
	}

	## Create the PhyGePaml object

	my $phygecodeml = PhyGePaml->new({ 
	    seqfams   => $paphyg->get_clusters(),
	    strains   => \%str_al,
	    topotypes => \%topotypes
					 });

	if (defined $opt_S) {
	    $phygecodeml->enable_reportstatus();
	}

	## Now it will feed the cds depending of the source

	if ($codeml{source_cds} eq 'file') {  ## It will get the CDS from file
	
	    if (defined $codeml{file_cds}) {
	    
		## Parse the file
		
		my %cds = ();
		my $cds_path = File::Spec->catfile($indir, $codeml{file_cds});
	    
		my $cdsio = Bio::SeqIO->new( -file => $cds_path, 
					     -format => 'fasta' );

		while( my $cdsseq = $cdsio->next_seq() ) {
		
		    $cds{$cdsseq->display_id()} = $cdsseq;
		}

		## Set the cds

		$phygecodeml->set_cds(\%cds);
	    }
	    else {
		die("ERROR: No cds_file was used for source_cds=file.\n\n");
	    }
	}
	else {  ## It will calculate the CDS using longest6frame
	    
	    $phygecodeml->predict_cds({ method => 'longest6frame' });
	}

	## It will set the align cds

	$phygecodeml->set_seqfam_cds();

	$phygecodeml->run_codeml($codeml{codeml_parameters});

	## Print the message of how many codeml has been run

	my %phygecml = %{$phygecodeml->get_codeml_results()};
	my $codeml_n = scalar(keys %phygecml);

	print STDERR "\n\n\t\t\t\t$codeml_n have been obtained\n\n";

	my $codemldir = $pathoutdir . '/12_codeml_analysis';
	mkdir($codemldir);
	my $codemlfile = File::Spec->catfile($codemldir, 'codeml_analysis.tab');

	$phygecodeml->out_codeml($codemlfile);
    }


    print STDERR "\t$path_idx ==> END. PATH=$phname (" .  date() . "):\n\n";
    print STDERR "\t\t\t\t<======+======>\n\n\n";
}

## 3) Analysis of phygetopos

print STDERR "\t3) RUNING STATISTICAL ANALISIS (" .  date() . "):\n\n";

## Check the number of phygetopo

my $phytopo_count = 0;
foreach my $method (keys %topostats) {
    my $phygtp = $topostats{$method};
    if (defined $phygtp) {
        $phytopo_count += scalar(keys (%{$phygtp->get_topotypes()}));
    }
}

if ($phytopo_count > 0) {

    ## Create the R base object

    my $rbase = R::YapRI::Base->new();
    
    ## UNCOMMENT TO TRACE ERRORS INTO R 
    ## $rbase->enable_debug();
    ## $rbase->enable_keepfiles();

    ## Create the PhyGeStat object.

    my $phygestats = PhyGeStats->new( { rbase     => $rbase, 
	  			        phygetopo => \%topostats,
				      } );

    ## Create the result dir

    my $statsdir = $outdir . '/12_Results';
    mkdir($statsdir);

    my $filetable = $statsdir . '/topology_composition.tab';
    my $barplot = $statsdir . '/topology_composition.bmp';
    my $treeplot = $statsdir . '/topology_trees.bmp';

    ## Create the matrix with the phygetopo data

    $phygestats->create_matrix(); 

    $phygestats->create_composition_table($filetable);
    $phygestats->create_composition_graph($barplot);
    $phygestats->create_tree_graph($treeplot);

    print STDERR "\t\tDONE. 3 files have been created:\n\n";
    print STDERR "\t\t\t$filetable\n\t\t\t$barplot\n\t\t\t$treeplot\n\n";

}
else {
    print STDERR "\t\tNO TOPOLOGIES WERE PRODUCED BY ANY PATH.\n";
    print STDERR "\t\tCHECK PRUNNING ARGUMENTS AND INITIAL DATASET.\n\n";
   
}

############################################
## 4) Print end message

print STDERR "\n\n";
print STDERR "*************************************************************";
print STDERR "\n** ENDING PhyGomic PIPELINE (" . date() . "):\n";
print STDERR "*************************************************************";
print STDERR "\n\nThanks for use PhyGomics.\n\n\n";



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
         2.1) Search homologous using blast                [optional]
         2.2) Run alignments, 'clustalw' by default        [mandatory]*
         2.3) Run distances, 'kimura' by default           [mandatory for NJ]*
         2.4) Prune members, it will rerun 2.2 and 2.3     [optional]
         2.5) Run trees, 'ML' & 'phyml' by default         [mandatory]*
         2.6) Tree rerooting                               [optional]
         2.7) Run bootstrapping                            [optional]*
         2.8) Run topoanalysis                             [mandatory]*

      3) Data integration and analysis [1 cycle]
         3.1) Create graphs and tables.

      (see -C, configuration file for more details about data)
      (using -O argument will create output files for each of these steps)


    Usage:
  
     phygomics.pl [-h] -c config.file -i input.dir -o output.dir [-C][-S][-O]
                              
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
      -O <create_steps_files>  Create the files for each step (alignments, 
                               distances, trees, consensus and topologies)

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
	source       => 'cluster_file',
	memberseq    => 'membersequence_file',
	memberstrain => 'memberstrain_file',
	);

    foreach my $req_k (keys %reqfiles) {
	unless (exists $confhref->{general}->{input_filenames}->{$req_k}) {
	    my $e = "ERROR: Required file:$reqfiles{$req_k} didnt was";
	    die("$e supplied to the configuration file.\n");
	}
    }

    my $indir = $confhref->{indir} || '';
    my %fileinput = %{$confhref->{general}->{input_filenames}};

    foreach my $confk (sort keys %fileinput) {

	if ($confk ne 'source_filetype') {

	    my $file = File::Spec->catfile($indir, $fileinput{$confk});
	    
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
		if ($1 eq 'max_cluster_members') {
		    $clval{$1} = $3;
		}
		else {
		    $clval{$1} = [$2, $3];
		}	
	    }
	}
	$clval_href = \%clval;
    }

    return $clval_href;
}


=head2 parse_homolog_args

  Usage: my $hom_val_href = parse_homolog_args(\%conf_for_path);

  Desc: parse homologous_search arguments and return a hash ref. usable by
        homologous_search function.

  Ret: $hom_val_href, a hash reference with homologous_search arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $hom_val_href = parse_homolog_args(\%conf_for_path);

=cut

sub parse_homolog_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_homolog_args()");

    my %homsearch = ();
    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied to parse_homolog_args() isnt a href.");
    }
    else {
	%homsearch = %{$pconf_href};
    }

    my %homsearch_form = ( blast => ['-p', 'blastn'] );

    my $blastdb = $homsearch{blast_database} ||
	die("ERROR: No blastdb was specified for homologous_search");
	
    my $path_blastdb = File::Spec->catfile($homsearch{indir}, $blastdb);

    push @{$homsearch_form{blast}}, ('-d', $path_blastdb);

    my $blastarg = $homsearch{blast_arguments} || '-e 1e-10'; 
    if ($blastarg =~ m/(,|;)/) {
	my @bargs = split(/$1/, $blastarg);
	foreach my $barg (@bargs) {

	    if ($barg !~ m/-\w/) {
		die("ERROR: ");
	    }

	    push @{$homsearch_form{blast}}, split(/\s+/, $barg);
	}
    }
    else {
	my @bargs = split(/\s+/, $blastarg);
	push @{$homsearch_form{blast}}, @bargs;
    }
    
    my $homstrain = $homsearch{homologous_strain} ||
	die("ERROR: No homologous_strain was specified homologous_search");
    
    $homsearch_form{strain} = $homstrain;
    
    if (defined $homsearch{blastmatch_filter}) {
	    
	my %filter = ();
	
	my $sep = ',';
	if ($homsearch{blastmatch_filter} =~ m/(;)/) {
	    $sep = $1;
	}
	my @filargs = split(/$sep/, $homsearch{blastmatch_filter});
	foreach my $fil (@filargs) {

	    $fil =~ s/^\s+//;
	    $fil =~ s/\s+$//;
	    
	    my @single_fil = split(/\s+/, $fil);
	    if (scalar(@single_fil) != 3) {
		die("ERROR: Wrong blastmatch_filter:$fil (format: x > 10)");
	    }
	    else {
		$filter{$single_fil[0]} = [$single_fil[1], $single_fil[2]];
	    }
	}
	$homsearch_form{filter} = \%filter;
    }

    return %homsearch_form;







    my %pargs = %{$pconf_href};
    
    my %hom_args = (
	blast  => [ -p => 'blastn', -d => $indir . '/' . $pargs{hs_dts} ],
	strain => $pargs{hs_str},
	);
	
    if (exists $pargs{hs_arg}) {
	my @bargs = split(/;/, $pargs{hs_arg});
	foreach my $barg (@bargs) {
	    if ($barg =~ m/\s*(-\w)\s+(.+)\s*/) {
		my $b_ar = $1;
		my $b_vl = $2;
		unless ($b_ar =~ m/-(o|i|p|d)/) {
		    push @{$hom_args{blast}}, ($b_ar => $b_vl); 
		}
	    }
	}
    }
    if (exists $pargs{hs_fil}) {
	my @bfils = split(/;/, $pargs{hs_fil});
	foreach my $bfil (@bfils) {
	    $bfil =~ s/\s+$//;
	    $bfil =~ s/^\s+//;
	    if ($bfil =~ m/\s*(.+)\s+(.+)\s+(.+)\s*/) {
		my ($arg, $ope, $val) = ($1, $2, $3);
		if (exists $hom_args{filter}) {
		    $hom_args{filter}->{$arg} = [$2, $3];
		}
		else {
		    $hom_args{filter} = { $arg => [$2, $3] };
		}
	    }
	}
    }
    return \%hom_args;
}

=head2 parse_align_args

  Usage: my $align_href = parse_align_args(\%conf_for_path);

  Desc: parse run_alignment arguments and return a hash ref. usable by
        run_alignment function.

  Ret: $align_href, a hash reference with run alignment arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $align_href = parse_align_args(\%conf_for_path);

=cut

sub parse_align_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_align_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied to parse_align_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    my %align_args = ( program => $pargs{al_prg} );
	
    if (exists $pargs{al_arg}) {
	
	$align_args{parameters} = {};

	my @aargs = split(/;/, $pargs{al_arg});
	foreach my $aarg (@aargs) {
	    $aarg =~ s/^\s+//;  ## Clean the spaces
	    $aarg =~ s/\s+$//;
	    if ($aarg =~ m/(.+)=(.+)/) {
		my $func = $1;
		my $val = $2;
		$func =~ s/\s+$//;  ## Clean the spaces
		$val =~ s/^\s+//;
		
		$align_args{parameters}->{$func} = $val; 
	    }
	    else {
		$align_args{parameters}->{$aarg} = undef;
	    }
	}
    }
    else {
	$align_args{parameters} = {};
    }
    
    return \%align_args;
}

=head2 parse_align_args

  Usage: my $prune_aln_href = parse_prune_aln_args(\%conf_for_path);

  Desc: parse prune_by_align arguments and return a hash ref. usable by
        prune_by_align function.

  Ret: $prune_aln_href, a hash reference with prune_by_align arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $prune_aln_href = parse_prune_aln_args(\%conf_for_path);

=cut

sub parse_prune_aln_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_prune_aln_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_prune_aln_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    my %pr_aln = ();

    if (exists $pargs{pr_aln}) {
	my @pr_alns = split(/;/, $pargs{pr_aln});
	foreach my $pr_aln (@pr_alns) {

	    $pr_aln =~ s/\s+$//;
	    $pr_aln =~ s/^\s+//;

	    if ($pr_aln =~ m/(.+)\s+(.+)\s+(.+)/) {
		my ($arg, $cond, $val) = ($1, $2, $3);
		$pr_aln{$arg} = [$cond, $val];
	    }
	}
    }

    return \%pr_aln;
}

=head2 parse_prune_str_args

  Usage: my $prune_str_href = parse_prune_str_args(\%conf_for_path);

  Desc: parse prune_by_strains arguments and return a hash ref. usable by
        prune_by_strains function.

  Ret: $prune_str_href, a hash reference with prune_by_strains arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $prune_str_href = parse_prune_str_args(\%conf_for_path);

=cut

sub parse_prune_str_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_prune_str_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_prune_str_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    my %pr_str = ();

    if (exists $pargs{pr_str}) {
	my @pr_strs = split(/;/, $pargs{pr_str});
	foreach my $pr_str (@pr_strs) {
	    if  ($pr_str =~ m/composition\s+=>\s+(.+)/) {
		$pr_str{composition} = {};

		my @comps = split(/,/, $1);
		foreach my $comp (@comps) {
		    $comp =~ s/^\s+//;
		    $comp =~ s/\s+$//;
		    if ($comp =~ m/(.+)\s*=\s*(\d+)/) {
			my ($arg, $val) = ($1, $2);
			$arg =~ s/\s+$//;
			$val =~ s/^\s+//;
			$pr_str{composition}->{$arg} = $val;
		    }
		}
	    }
	    elsif ($pr_str =~ m/(\w+_distance)\s+=>\s+(.+)/) {
		my $arg = $1;
		$pr_str{$arg} = [];

		my @dists = split(/,/, $2);
		foreach my $dist (@dists) {
		    $dist =~ s/^\s+//;
		    $dist =~ s/\s+$//;
		    if ($dist =~ m/(.+)\s*=\s*(.+)/) {
			my ($pair1, $pair2) = ($1, $2);
			$pair1 =~ s/\s+$//;
			$pair2 =~ s/^\s+//;
			push @{$pr_str{$arg}}, [$pair1, $pair2];
		    }
		}
	    }
	}
    }

    return \%pr_str;
}


=head2 parse_prune_ovl_args

  Usage: my $prune_ovl_href = parse_prune_ovl_args(\%conf_for_path);

  Desc: parse prune_by_overlaps arguments and return a hash ref. usable by
        prune_by_overlaps function.

  Ret: $prune_ovl_href, a hash reference with prune_by_overlaps arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $prune_str_href = parse_prune_ovl_args(\%conf_for_path);

=cut

sub parse_prune_ovl_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_prune_ovl_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_prune_ovl_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    my %pr_ovl = ();

    if (exists $pargs{pr_ovl}) {
	my @pr_strs = split(/;/, $pargs{pr_ovl});
	foreach my $pr_ovl (@pr_strs) {
	    
	    $pr_ovl =~ s/^\s+//;
	    $pr_ovl =~ s/\s+$//;

	    if  ($pr_ovl =~ m/composition\s+=>\s+(.+)/) {
		$pr_ovl{composition} = {};

		my @comps = split(/,/, $1);
		foreach my $comp (@comps) {
		    $comp =~ s/^\s+//;
		    $comp =~ s/\s+$//;
		    if ($comp =~ m/(.+)\s*=\s*(\d+)/) {
			my ($arg, $val) = ($1, $2);
			$arg =~ s/\s+$//;
			$val =~ s/^\s+//;
			$pr_ovl{composition}->{$arg} = $val;
		    }
		}
	    }
	    elsif ($pr_ovl =~ m/(.+)=>(.+)/i) {
		my $arg = $1;
		my $val = $2,
		
		$arg =~ s/\s+$//;
		$val =~ s/^\s+//;

		$pr_ovl{$arg} = $val;
	    }
	}
    }

    return \%pr_ovl;
}

=head2 parse_njtrees_args

  Usage: my $njtrees_href = parse_njtrees_args(\%conf_for_path);

  Desc: parse run_njtrees arguments and return a hash ref. usable by
        run_njtrees function.

  Ret: $njtrees_href, a hash reference with njtrees arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $njtrees_href = parse_njtrees_args(\%conf_for_path);

=cut

sub parse_njtrees_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_njtrees_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_njtrees_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    my %njtrees = ();

    if (exists $pargs{tr_arg}) {
	my @trargs = split(/;/, $pargs{tr_arg});
	foreach my $trarg (@trargs) {
	    $trarg =~ s/^\s+//;
	    $trarg =~ s/\s+$//;
	    if  ($trarg =~ m/(.+)\s*=>\s*(.+)/) {
		my $arg = $1;
		my $val = $2;
		$arg =~ s/\s+$//;
		$val =~ s/^\s+//;
		$njtrees{$arg} = $val;
	    }
	}
    }

    ## Overwrite type with $pargs{tr_met}

    $njtrees{type} = $pargs{tr_met};

    return \%njtrees;
}


=head2 parse_mltrees_args

  Usage: my $mltrees_href = parse_mltrees_args(\%conf_for_path);

  Desc: parse run_mltrees arguments and return a hash ref. usable by
        run_mltrees function.

  Ret: $mltrees_href, a hash reference with mltrees arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $mltrees_href = parse_mltrees_args(\%conf_for_path);

=cut

sub parse_mltrees_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_mltrees_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_mltrees_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    my %mltrees = ();

    my $prog;
    if (exists $pargs{tr_arg}) {
	my @trargs = split(/;/, $pargs{tr_arg});
	foreach my $trarg (@trargs) {
	    if ($trarg =~ m/phyml\s*=?>?\s*(.*)/) {
		$prog = 'phyml';
		$mltrees{phyml} = {};
		my $args1 = $1;
		if (defined $args1) {
		    my @args1 = split(/,/, $args1);
		    foreach my $arg1 (@args1) {
			$arg1 =~ s/^\s+//;
			$arg1 =~ s/\s+$//;
			if ($arg1 =~ m/(.+)\s*=\s*(.+)/) {
			    my ($arg, $val) = ($1, $2);
			    $arg =~ s/\s+$//;
			    $val =~ s/^\s+//;
			    $mltrees{phyml}->{$arg} = $val;
			}
		    }
		}
	    }
	    elsif ($trarg =~ m/dnaml\s*=>\s*(.+)/) {
		
		$prog = 'dnaml';
		$mltrees{dnaml} = {};
		my $args2 = $1;
		if (defined $args2) {
		    my @args2 = split(/,/, $args2);
		    foreach my $arg2 (@args2) {
			$arg2 =~ s/^\s+//;
			$arg2 =~ s/\s+$//;
			if ($arg2 =~ m/(.+)\s*=\s*(.+)/) {
			    my ($arg, $val) = ($1, $2);
			    $arg =~ s/\s+$//;
			    $val =~ s/^\s+//;
			    $mltrees{dnaml}->{$arg} = $val;
			}
		    }
		}
	    }
	    elsif ($trarg =~ m/outgroup_strain\s*=>?\s*(.+)/) {
		$mltrees{outgroup_strain} = $1;
	    }
	}
    }
    unless (defined $prog) {  ## If no program is defined add phyml by default
	$mltrees{phyml} = {};
    }


    return \%mltrees;
}


=head2 parse_reroot_args

  Usage: my $reroot_href = parse_reroot_args(\%conf_for_path);

  Desc: parse reroot_trees arguments and return a hash ref. usable by
        reroot_trees function.

  Ret: $reroot_href, a hash reference with reroot arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $reroot_href = parse_reroot_args(\%conf_for_path);

=cut

sub parse_reroot_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_reroot_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_reroot_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    my %reroot = ();

    if (exists $pargs{tr_rer}) {
	my @trargs = split(/;/, $pargs{tr_rer});
	foreach my $trarg (@trargs) {
	    $trarg =~ s/^\s+//;
	    $trarg =~ s/\s+$//;
	    if  ($trarg =~ m/(.+?)\s*=\s*(.+?)/) {
		my ($arg, $val) = ($1, $2);
		$arg =~ s/\s+$//;
		$val =~ s/^\s+//;
		$reroot{$arg} = $val;
	    }
	}
    }

    return \%reroot;
}


=head2 parse_boots_args

  Usage: my $boots_href = parse_boots_args(\%conf_for_path);

  Desc: parse run_bootstrapping arguments and return a hash ref. usable by
        run_bootstrapping function.

  Ret: $boots_href, a hash reference with bootstrapping arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $boots_href = parse_boots_args(\%conf_for_path);

=cut

sub parse_boots_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_boots_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_boots_args() isnt a href.");
    }

    my %pargs = %{$pconf_href};
    
    ## Define the bootstrapping args

    my %boots = ( 
	run_bootstrap => { datatype => 'Sequence', quiet => 1 },
	run_distances => { quiet => 1},
	## Tree will depend of the tree_method
	run_consensus => { normalized => 1, quiet => 1  }
	);

    ## Bootstrap values

    if (exists $pargs{bo_run}) {
	my @boargs = split(/;/, $pargs{bo_run});
	
	foreach my $bo (@boargs) {
	    $bo =~ s/^\s+//;
	    $bo =~ s/\s+$//;
	    if  ($bo =~ m/(.+)\s*=\s*(.+)/) {
		my ($arg, $val) = ($1, $2);
		$arg =~ s/\s+$//;
		$val =~ s/^s+//;
		$boots{run_bootstrap}->{$arg} = $val;
	    }
	}
    }

    ## Distances

    if (exists $pargs{di_fun}) {
	$boots{run_distances}->{method} = $pargs{di_fun};
    }
    else {
	$boots{run_distances}->{method} = 'Kimura';
    }

    ## Trees

    my $tree_mth = $pargs{tr_met} || 'ML';
    if ($tree_mth =~ m/^(NJ|UPGMA)$/) {

	my $njtree_args = parse_njtrees_args(\%pargs);
	$boots{run_njtrees} = $njtree_args;

	if (exists $njtree_args->{outgroup_strain}) {
	    $boots{run_consensus}->{outgroup_strain} = 
		$njtree_args->{outgroup_strain};
	}
    }
    else {
	
	my $mltree_args = parse_mltrees_args(\%pargs);
	$boots{run_mltrees} = $mltree_args;

	if (exists $mltree_args->{outgroup_strain}) {
	    $boots{run_consensus}->{outgroup_strain} = 
		$mltree_args->{outgroup_strain};
	}
    }

    ## Consensus

    if (defined $pargs{tr_rer}) {
		    
	my $reroot_args = parse_reroot_args(\%pargs);
	if (exists $reroot_args->{midpoint}) {
	    $boots{run_consensus}->{root_by_midpoint} = 1;
	}
    }

    return \%boots;
}


=head2 parse_topo_args

  Usage: my $topo_href = parse_topo_args(\%conf_for_path);

  Desc: parse run_topoanalysis arguments and return a hash ref. usable by
        run_topoanalysis function.

  Ret: $topo_href, a hash reference with topoanalysis arguments

  Args: \%conf_for_path, a hash ref. with the configuration data for a concrete
        path

  Side_Effects: None

  Example: my $boots_href = parse_topo_args(\%conf_for_path);

=cut

sub parse_topo_args {
    my $pconf_href = shift ||
	die("ERROR: No conf. hashref. was supplied to parse_topo_args()");

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied parse_topo_args() isnt a href.");
    }

    my %topo = ();

    my %pargs = %{$pconf_href};

    if (defined $pargs{tp_run}) {

	my @tpargs = split(/;/, $pargs{tp_run});

	foreach my $tparg (@tpargs) {

	    if  ($tparg =~ m/(.+?)\s*=>\s*(.+)/) {
		my ($arg, $val) = ($1, $2);

		if ($val =~ m/=/) {

		    $topo{$arg} = {};
		    my @equs = split(/,/, $val);
		    foreach my $eq (@equs) {

			if ($eq =~ m/(.+)\s*=\s*(.+)/) {
			    $topo{$arg}->{$1} = $2;
			}
		    }
		}
		else {

		    $topo{$arg} = $val;
		}
	    }
	}
    }

    return \%topo;

}

####
1; #
####
