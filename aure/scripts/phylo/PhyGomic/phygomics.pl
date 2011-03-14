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
use PhyGeTopo;
use PhyGeStats;

our ($opt_c, $opt_i, $opt_o, $opt_h, $opt_C, $opt_S, $opt_O);
getopts("c:i:o:hCSO");
if (!$opt_c && !$opt_i && !$opt_o && !$opt_C && !$opt_S && !$opt_O) {
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
    print STDERR ("\nMESSAGE: Outdir doesnt exist. Creating outdir=$outdir\n");
    mkdir($outdir);
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

#############################################################################
## 0) Parse configuration file, check files and PhyGeCluster object creation.

print STDERR "0) PARSING CONFIGURATION FILE (" .  date() . "):\n";

my %conf = parse_config($control);
my $cv = scalar(keys %conf);
my %cpaths = %{$conf{paths}};
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
print STDERR "\n\tPARSED. $cl_v clusters have been extracted.\n\n";

if ($opt_O) {

    ## Create a folder:
    my $clusterdir = $outdir . '/0_cluster_prepath';
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

print STDERR "\n\tPARSED. $seqn sequences have been extracted.\n\n";


###################################################
## 1.3) Parse the strain file and add to phygecluster

my $str_args = { strainfile => $indir . '/' . $conf{str_fn} };
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

my %paths = %{$conf{paths}};
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


    ## 2.0) Clone the PhyGeCluster

    print STDERR "\t\t2.0) CLONING CLUSTER DATA (" .  date() . "):\n\n";

    my $paphyg = $phyg->clone($opt_S);

    ## 2.1) Homologous search

    print STDERR "\t\t2.1) HOMOLOGOUS SEARCH (" .  date() . "):\n\n";

    if (defined $pargs{hs_dts} && defined $pargs{hs_str}) {
	
	my $hom_args_href = parse_homolog_args(\%pargs);

	## And run homologous search
	
	$paphyg->homologous_search($hom_args_href);
	
	## Count homologous
	my $hom_c = 0;

	my %cls = %{$paphyg->get_clusters()};
	my %str = %{$paphyg->get_strains()};
	foreach my $clid (keys %cls) {
	    my @membs = $cls{$clid}->get_members();
	    foreach my $memb (@membs) {
		my $m_id = $memb->id();
		if (defined $m_id) {
		    if ($str{$m_id} eq $pargs{hs_str}) {
			$hom_c++;
		    }
		}
	    }
	}

	print STDERR "\t\t\tDONE. $hom_c homologous have been assigned.\n\n";

	if ($opt_O) {

	    ## Create a folder:
	    my $homologdir = $pathoutdir . '/1_cluster_homologs';
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
	print STDERR "SKIP STEP.\n\n";
    }

    ## 2.2) Run alignments

    print STDERR "\t\t2.2) RUN ALIGNMENTS (" .  date() . "):\n\n";

    my $align_args;

    if (defined $pargs{al_prg}) {
	
	 $align_args = parse_align_args(\%pargs);		
    }
    else {
	 print STDERR "\t\t\t USING DEFAULT ALIGNMENT (clustalw)\n\n";
	 $align_args = { program    => 'clustalw', 
			 parameters => { quiet => undef, matrix => 'BLOSUM'},
	 };
    }
    
    if ($opt_S) {
	$align_args->{report_status} = 1;
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
	my $alndir = $pathoutdir . '/2_alignments';
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

    ## 2.3) Run distances

    print STDERR "\t\t2.3) RUN DISTANCES (" .  date() . "):\n\n";

    my $dist_args;
    if (defined $pargs{di_fun}) {
	
	$dist_args = { method => $pargs{di_fun}, quiet => 1};	
    }
    else {

	print STDERR "\t\t\t USING DEFAULT DISTANCE (Kimura)\n\n";
	$dist_args = { method => 'Kimura', quiet => 1};
    }

    if ($opt_S) {
	$dist_args->{report_status} = 1;
    }

    ## And now use the Run distances method

    my @failed_dist = $paphyg->run_distances($dist_args);

    my $dist_n = scalar( keys %{$paphyg->get_distances()});    
    print STDERR "\n\t\t\tDONE. $dist_n distance matrices have been added.\n\n";

    if ($opt_O) {

	## Create a folder:
	my $disdir = $pathoutdir . '/3_distances';
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
	

    ## 2.4) Run prune member methods

    print STDERR "\t\t2.4) PRUNING METHODS (" .  date() . "):\n\n";
		
    my $prune = 0;
    my @failed_praln = ();
    my @failed_prdst = ();

    if (defined $pargs{pr_aln}) {

	print STDERR "\t\t\t2.4.1) PRUNE BY ALIGNMENTS (" .  date() . "):\n\n";
	
	my $pr_aln_args = parse_prune_aln_args(\%pargs);
	my %rem_cls = $paphyg->prune_by_align($pr_aln_args);
	my $rem_c = scalar(keys %rem_cls);

	print STDERR "\t\t\t\tDONE. $rem_c clusters have been removed.\n\n";
	
	print STDERR "\t\t\t\tRERUNING alignments.\n\n";
	my @failed_pr_aln = $paphyg->run_alignments($align_args);
	push @failed_praln, @failed_pr_aln;
 
	print STDERR "\n\t\t\t\tRERUNING distances.\n\n";
	my @failed_pr_dist = $paphyg->run_distances($dist_args);
        push @failed_prdst, @failed_pr_dist;

	print STDERR "\n";

        $prune = 1;
    }
    else {
	print STDERR "\t\t\tNo Arguments. SKIPPING PRUNE BY ALIGNMENTS.\n\n";
    }

    if (defined $pargs{pr_str}) {
    
    	print STDERR "\t\t\t2.4.2) PRUNE BY STRAINS (" .  date() . "):\n\n";

	my $pr_str_args = parse_prune_str_args(\%pargs);

	my ($rem_p_cls, $rem_p_mbs)  = $paphyg->prune_by_strains($pr_str_args);
	my $rem_p_c = scalar(keys %{$rem_p_cls});
	my $rem_p_m = scalar(keys %{$rem_p_mbs});

	print STDERR "\t\t\t\tDONE. $rem_p_c clusters and $rem_p_m members ";
	print STDERR "have been removed.\n\n";
	
	print STDERR "\t\t\t\tRERUNING alignments.\n\n";
	my @failed_pr_aln = $paphyg->run_alignments($align_args);
	push @failed_praln, @failed_pr_aln;
 
	print STDERR "\n\t\t\t\tRERUNING distances.\n\n";
	my @failed_pr_dist = $paphyg->run_distances($dist_args);
        push @failed_prdst, @failed_pr_dist;
        print STDERR "\n";
        $prune = 1;
    }
    else {
	print STDERR "\t\t\tNo Arguments. SKIPPING PRUNE BY STRAINS.\n\n";
    }

    if (defined $pargs{pr_ovl}) {
 	
   	print STDERR "\t\t\t2.4.3) PRUNE BY OVERLAPINGS (" .  date() . "):\n\n";
    
	my $pr_ovl_args = parse_prune_ovl_args(\%pargs);
	my ($rem_o_cls, $rem_o_mbs)  = $paphyg->prune_by_overlaps($pr_ovl_args);
	my $rem_o_c = scalar(keys %{$rem_o_cls});
	my $rem_o_m = scalar(keys %{$rem_o_mbs});

	print STDERR "\t\t\t\tDONE. $rem_o_c clusters and $rem_o_m members ";
	print STDERR "have been removed.\n\n";
	
	unless (exists $pr_ovl_args->{trim}) {

	    print STDERR "\t\t\t\tRERUNING alignments.\n\n";
	    my @failed_pr_aln = $paphyg->run_alignments($align_args);
	    push @failed_praln, @failed_pr_aln;
	}
	else {
	     print STDERR "\t\t\t\tTRIM option ENABLED. ";
	     print STDERR "SKIPPING  rerun alignments.\n\n";
	}
 
	print STDERR "\n\t\t\t\tRERUNING distances.\n\n";
	my @failed_pr_dist = $paphyg->run_distances($dist_args);
        push @failed_prdst, @failed_pr_dist;
        print STDERR "\n";
        $prune = 1;
    }
    else {
	print STDERR "\t\t\tNo Arguments. SKIPPING PRUNE BY OVERLAPINGS.\n\n";
    }

    my $cls_n = scalar( keys %{$paphyg->get_clusters()});
    print STDERR "\t\t\t$cls_n clusters remain after prunning.\n\n";
    
    if ($opt_O && $prune == 1) {

	## Create a folder:
	my $prdir = $pathoutdir . '/4_afterprunning';
	mkdir($prdir);
	
	## Print files
	my $pr_alnbase = $prdir . '/ap.alignments';
	my %pr_alnfiles = $paphyg->out_alignfile({ 
	    'rootname'     => $pr_alnbase,
	    'distribution' => 'single',
	    'format'       => 'clustalw',
						 });
	
	## Print alignments stats
	my $pralnstatbase = $prdir . '/ap.stats_aln.tab';
	my $prstats = $paphyg->out_alignstats({ basename => $pralnstatbase});
	print STDERR "\t\t\tOPTION -O enabled: ";
	print STDERR "a stats alignment file have been created";
	print STDERR " with filename:\n\t\t\t$prstats\n\n";

        my $pr_disbase = $prdir . '/ap.distances';
	my %pr_disfiles = $paphyg->out_distancefile({ 
	    'rootname'     => $pr_disbase,
	    'distribution' => 'single',
	    'format'       => 'phylip',
						 });
	my $pr_aln_n = scalar(keys %pr_alnfiles);
        my $pr_dis_n = scalar(keys %pr_disfiles);
	print STDERR "\t\t\tOPTION -O enabled: $pr_aln_n align. and $pr_dis_n ";
	print STDERR "distance files (after prunning) have been created ";
	print STDERR "with basename:\n\t\t\t$pr_alnbase\n\t\t\t$pr_disbase\n\n";
    
	if (scalar(@failed_praln) > 0) {
	    open my $ffh_praln, '>', $prdir . '/ap.failed_distances.err';
	    foreach my $f_praln (@failed_praln) {
		print $ffh_praln "$f_praln\n";
	    }
	    close($ffh_praln);
	}
	if (scalar(@failed_prdst) > 0) {
	    open my $ffh_prdist, '>', $prdir . '/ap.failed_distances.err';
	    foreach my $f_prdist (@failed_prdst) {
		print $ffh_prdist "$f_prdist\n";
	    }
	    close($ffh_prdist);
	}
    }


    ## Now it will continue, only if $cls_n > 0

    if ($cls_n > 0) {
    
	## 2.5) Run tree methods
	
	print STDERR "\t\t2.5) RUN TREES (" .  date() . "):\n\n";

	my @failed_trees = ();
	my $tree_mth = $pargs{tr_met} || 'ML';
	if ($tree_mth =~ m/^(NJ|UPGMA)$/) {

	    my $njtree_args = parse_njtrees_args(\%pargs);

            if ($opt_S) {
	        $njtree_args->{report_status} = 1;
            }

	    @failed_trees = $paphyg->run_njtrees($njtree_args);
	}
	else {

	    my $mltree_args = parse_mltrees_args(\%pargs);

            if ($opt_S) {
	        $mltree_args->{report_status} = 1;
            }

	    @failed_trees = $paphyg->run_mltrees($mltree_args);
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
	    my $treedir = $pathoutdir . '/5_trees';
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


	## 2.6) Reroot tree

	print STDERR "\t\t2.6) REROOTING TREES (" .  date() . "):\n\n";
	if (defined $pargs{tr_rer}) {
		    
	    my $reroot_args = parse_reroot_args(\%pargs);
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
	        my $retreedir = $pathoutdir . '/6_trees_rerooted';
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


	## 2.7) Run bootstrapping

	print STDERR "\t\t2.7) RUN BOOTSTRAPPING (" .  date() . "):\n\n";
	if (defined $pargs{bo_run}) {
		    
	    my $boots_args = parse_boots_args(\%pargs);

            if ($opt_S) {
	        $boots_args->{report_status} = 1;
            }
            
	    $paphyg->run_bootstrapping($boots_args);

	    my %boots = %{$paphyg->get_bootstrapping()};
	    my $boots_n = scalar(keys %boots);

	    print STDERR "\t\t\tDONE. $boots_n bootstraps have been added\n\n";

	    if (defined $pargs{bo_fil}) {

		print STDERR "\t\t\tPRUNING BY BOOTSTRAP VALUES.\n\n";
		my %rmcls_boots = $paphyg->prune_by_bootstrap($pargs{bo_fil});
	    
		my $rboots_n = scalar(keys %rmcls_boots);
	    	
		print STDERR "\n\t\t\tDONE. $rboots_n trees have been removed ";
		print STDERR "according the bootstrapping values.\n\n";
	    }

            if ($opt_O) {

	        ## Create a folder:
	        my $bootsdir = $pathoutdir . '/7_bootstrap_consensus';
	        mkdir($bootsdir);
	
	        ## Print files
	        my $bootsbase = $bootsdir . '/consensus';
	        my %bootsfiles = $paphyg->out_bootstrapfile({ 
	             rootname     => $bootsbase,
	             distribution => 'single',
	             format       => 'newick',
		   				      });

	        my $boots_n = scalar(keys %bootsfiles);
	        print STDERR "\t\t\tOPTION -O enabled: ";
	        print STDERR "$boots_n bootstrapping consensus files have been";
	        print STDERR " created with basename:\n\t\t\t$bootsbase\n\n";
             }	    
	}
	else {
	    print STDERR "\t\t\tNo Arguments. SKIPPING BOOTSTRAPPING.\n\n";
	}

	## 2.8) Run phygetopo

	print STDERR "\t\t2.8) RUN TOPOANALYSIS (" .  date() . "):\n\n";
		    
	my $seqfams_href = $paphyg->get_clusters();
	my $strains_href = $paphyg->get_strains();

	my $phygetopo = PhyGeTopo->new({ seqfams => $seqfams_href,
					 strains => $strains_href,
				       });

	my $topo_args = parse_topo_args(\%pargs);
	my %topotypes = $phygetopo->run_topoanalysis($topo_args);
	
	$topostats{$phname} = $phygetopo;

	my $topot_n = scalar(keys %topotypes);
	print STDERR "\t\t\tDONE. $topot_n topotypes have been created\n\n";

        if ($opt_O) {

	    ## Create a folder:
	    my $topodir = $pathoutdir . '/8_topologies';
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

    my $statsdir = $outdir . '/9_Results';
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
## Boths can use max_cluster_members

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
##           [trim   => {1}]
## example:  [composition => Sly=1,Nta=2,Nto=1,Nsy=1]
##           [trim => 1 ]
## arguments: composition, random, trim, ovlscore


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
## format:  [{argument}{space}{=>}{space}{value}{semicolon}...]
## example: [quiet => 1; outgroup_strain => Sly] for NJ
##          [ phyml => -b=1000, -m=JC69]
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
## format:  [{argument}={integer}{semicolon}]
## example: [replicates = 1000]

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
## format:  [branch_cutoffs{=>}{integer}{=}{integer}{comma}{semicolon}...]
## example: [branch_cutoffs => 0.1=1]






## Configuration file end.
/;

    print $cfh $info;
    print STDERR "\nCONTROL FILE: $contrfile has been created (".date().")\n\n";
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

    unless (ref($pconf_href) eq 'HASH') {
	die("ERROR: $pconf_href supplied to parse_homolog_args() isnt a href.");
    }

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
	    elsif ($pr_ovl =~ m/(random|trim|ovlscore)/i) {
		$pr_ovl{lc($1)} = 1;
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
