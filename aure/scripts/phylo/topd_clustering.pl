#!/usr/bin/perl

=head1 NAME

 topd_clustering.pl
 Script to parse TOPD output and cluster the trees by types

=cut

=head1 SYPNOSIS

 topd_clustering.pl [-h] -i <topd_result_file> [-t <tree_file>] [-l <lnlikehood_file>] [-o <output>] [-m <topology_method>] [-c <cluster_variable>] [-v <cluster_cutoff_value>]

=head1 EXAMPLE 

 topd_clustering.pl -i topd_sol.txt -c disagree
                                  
=head2 I<Flags:>

=over


=item -i

B<topd_result_file>            TOPD result file (mandatory)

=item -t

B<tree_file>                   Tree file used as input for add + TREE to the output file (mandatory)

=item -l

B<lnlikehood_file>             Likehood file (-f1 ID, -f2 LNLIKEHOOD_VALUE) to add ln_likehood value to the output file (optional).

=item -o

B<output>                      Output basename (topd_clustering.output by default)

=item -m

B<topology_method>             Method used for topologies comparisson (method used by TOPD to compare trees: nodal, split, quartets, triplets, disagree)(disagree by default)

=item -v

B<method_variable>             Variable produced by the method used for clustering (distance, taxa, DC, EA...)(taxa by default)

=item -c

B<variable_cutoff>             Cutoff value to assign a relation as true for clustering (0 by default)

=item -h

B<help>                         Print the help

=back

=cut

=head1 DESCRIPTION

   TOPD/FMTS [Puigbo et al., 2007] is a program that 
  compares phylogenetic trees based in five different
  methods (nodal, split, quartets, triplets, disagree). 
  This script parse the result file and use these 
  comparisons to cluster the phylogenetic trees based in 
  their topologies.

   Depending of the method used for the comparisons, it is
  possible use different variables for clustering:

   -m nodal -v distance -c [0-1]
   -m split -v distance -c [0-1]
   -m quartets -v [s|d|r1|r2|u|DC|EA|SJA|SSJA|SD] -c [0-1]
   -m triplets -v [s|d|r1|r2|u|DC|EA|SJA|SSJA|SD] -c [0-1]
   
   -m disagree -v taxa -c 0 (by default)

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 topd_clusteing.pl


=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;


our ($opt_i, $opt_t, $opt_l, $opt_o, $opt_m, $opt_v, $opt_c, $opt_h);
getopts("i:t:l:o:m:v:c:h");
if (!$opt_i && !$opt_t && !$opt_l && !$opt_o && !$opt_m && !$opt_v && !$opt_c && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Mandatory arguments (die if they are not supplied)

my $topd = $opt_i 
    || die("\n\nMANDATORY ARGUMENT ERROR: -i <topd_result_file> argument was not supplied.\n\n");
my $tree = $opt_t
    || die("\n\nMANDATORY ARGUMENT ERROR: -t <tree_file> argument was not supplied.\n\n");

## Default arguments (use default values if they are not supplied)

my $output = $opt_o 
    || 'topd_clustering.output';
my $method = $opt_m
    || 'disagree';
my $variable = $opt_v
    || 'taxa';
my $cutoff = $opt_c
    || 0;

if ($method =~ m/^nodal|split$/) {
    unless ($variable =~ m/^distance$/) {
	die("\n\nARGUMENT INCOMPATIBILITY: -m $method only can be used with -v distance.\n\n");
    }
}
elsif ($method =~ m/^quartets|triplets$/) {
    unless ($variable =~ m/^s|d|r1|r2|u|DC|EA|SJA|SSJA|SD$/) {
	die("\n\nARGUMENT INCOMPATIBILITY: -m $method only can be used with -v [s|d|r1|r2|u|DC|EA|SJA|SSJA|SD].\n\n");
    }
}
elsif ($method =~ m/^disagree$/) {
    unless ($variable =~ m/^taxa$/) {
	die("\n\nARGUMENT INCOMPATIBILITY: -m $method only can be used with -v taxa.\n\n");
    }
}
else {
    die("\n\nARGUMENT VALUE ERROR: -m method only can be [nodal|split|quartets|triplets|disagree].\n\n");
}

## Optional arguments

my $lnlikehood = $opt_l;

##############################################
#### PARSING FILE  ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);

print STDERR "\n\nSCRIPT START POINT ($date)\n\n";

$date = `date`;
chomp($date);
print STDERR "\t1) Parsing tree file: $tree ($date)\n\n";

my %tree_id = ();
    
my $t = 0;
my $T = `wc -l $tree`;
chomp($T);
$T =~ s/$tree//;
    
open my $t_fh, '<', $tree;
    
while(<$t_fh>) {
    chomp($_);
    $t++;
	
    print STDERR "\t\tParsing line $t of $T lines in the file: $tree\r";
	
    my @t_data = split(/\t/, $_);
    $tree_id{$t_data[0]} = $t_data[1];	
}
print STDERR "\n\n";


my $max_tree = scalar(keys %tree_id);

$date = `date`;
chomp($date);
print STDERR "\t2) Parsing TOPD result file ($date): $topd\n\n";

open my $ifh, '<', $topd;
my $l = 0;
my $L = `wc -l $topd`;
chomp($L);
$L =~ s/$topd//;

my %cluster_tree = ();
my %cluster_types = ();

## Declare variables.
my ($id1, $id2, $overlap);
my $cl = 1;
my $tr = 0;

while (<$ifh>) {
    $l++;
    chomp($_);
    my $dif_cl = $cl - 1;    
    print STDERR "\t\tParsing line $l of $L lines in the file: $topd ($dif_cl different clusters, $tr trees)\r";

    if ($max_tree > 0 && $tr < $max_tree) {  ## This will skip the rest of the relations once there are not more trees to assign

	## First take relation and overwrite the relation variables
	my $relation = 'FALSE';

	if ($_ =~ m/^#+?\s+topd\s+(.+)\s+-\s+(.+)\s+#+/) {
	$id1 = $1;
	$id2 = $2;
	$overlap = 'TRUE';
	}
	if ($_ =~ m/The overlap in the data between .+? and .+? is too separate./) {
	    $overlap = 'FALSE';
	}
    
	if (defined $id1 && defined $id2 && $overlap eq 'TRUE') {
	    if ($method eq 'nodal') {
		if ($_ =~ m/^\*\sNodal\sDistance\s\(Pruned\/Unpruned\):\s(.+)\s\/\s(.+)$/) {
		    my $pru_dist = $1;
		    my $unpru_dist = $2;
		    
		    ## Now check if the topology is similar under the cutoff
		    
		    if ($pru_dist <= $cutoff) {
			$relation = 'TRUE';
		    }
		}
	    }
	    elsif ($method eq 'split') {
		if ($_ =~ m/^\*\sSplit\sDistance\s\[differents\/possibles\]:\s(.+)\s\[\s.+\s\/\s.+\s\]/) {
		    my $dist = $1;
		
		    if ($dist <= $cutoff) {
			$relation = 'TRUE';
		    }
		}
	    }
	    elsif ($method eq 'disagree') {
		if ($_ =~ m/^\*\sDisagreement/) {
		    my @dis_data = split(/,/, $_);
		    foreach my $disdata (@dis_data) {
			if ($disdata =~ m/^\*\sDisagreement\s\[\staxa\sdisagree\s\/\sall\staxa\s\]:\s\[\s?(.+?)\s/) {
			    my $taxa = $1;
			
			    if ($taxa <= $cutoff) {
				$relation = 'TRUE';
			    }
			}
		    }
		}
	    }
	    else {
		if ($_ =~ m/^\*\sQuartets|Triplets/) {
		    my @qt_data = split(/, /, $_);
		    foreach my $qt (@qt_data) {
			my @pair = split(/: /, $qt);

			if ($pair[0] eq $variable) {
			    if ($pair[1] =~ m/(\d+\.?\d*)/) {
				my $value = $1;
				if ($value <= $cutoff) {
				    $relation = 'TRUE';
				}
			    }		    
			}		
		    }			    
		}
	    }

	    ## If the relation is selected as true

	    if ($relation eq 'TRUE') {
		if (exists $cluster_tree{$id1}) {

		    ## Add new relation to cluster	    
		    unless (exists $cluster_tree{$id2}) {
			$cluster_tree{$id2} = $cluster_tree{$id1};

			my @test = @{$cluster_types{$cluster_tree{$id1}}};
			my $c_test = scalar(@test);
			push @{$cluster_types{$cluster_tree{$id1}}}, $id2;
			$tr++;
		    }
		}
		else {

		    ## Add a new cluster
		    unless (exists $cluster_tree{$id2}) {
			$cluster_tree{$id1} = $cl;
			$cluster_tree{$id2} = $cl;
			$cluster_types{$cl} = [$id1, $id2];
			$cl++;
			$tr += 2;
		    }
		    else {
			$cluster_tree{$id1} = $cluster_tree{$id2};
			push @{$cluster_types{$cluster_tree{$id2}}}, $id1;
			$tr++;
		    }
		}
		
		## Finally undef relation variables
		undef($id1);
		undef($id2);
	    }
	}
    }
}
print STDERR "\n\n";

$date = `date`;
chomp($date);
print STDERR "\t3) Parsing ln_likehood file:";

my %lnlikehood_id;
if (defined $lnlikehood) {
    print STDERR " $lnlikehood ($date)\n\n";

    my $h = 0;
    my $H = `wc -l $lnlikehood`;
    chomp($H);
    $H =~ s/$lnlikehood//;
    
    open my $h_fh, '<', $lnlikehood;
    
    while(<$h_fh>) {
	chomp($_);
	$h++;
	
	print STDERR "\t\tParsing line $h of $H lines in the file: $lnlikehood\r";
	
	my @h_data = split(/\t/, $_);
	$lnlikehood_id{$h_data[0]} = $h_data[1];	
    }
}
else {
    print STDERR "\n\n\t\tSKIPPING ($date). No tree file was supplied to this script.\n";
}
print STDERR "\n\n";


##############################################
#### DATA ANALYSIS ###########################
##############################################


my $outfilename = $output . '.cluster.tab';

open my $o_fh, '>', $outfilename;

$date = `date`;
chomp($date);
print STDERR "\t4) Printing report file: $outfilename ($date)\n\n";

## Define singlet
my $sgl = 0;

foreach my $id (sort keys %tree_id) {
    my $cluster = $cluster_tree{$id};
    my $treedata= $tree_id{$id};

    my $cluster_name;
    if (defined $cluster) {
	$cluster_name = 'CL-' . $cluster;
    }
    else {
	$sgl++;
	$cluster_name = 'SG-' . $sgl;
    }

    print $o_fh "$id\t$cluster_name\t$treedata";

    if (defined $lnlikehood) {
	my $lnlikehooddata = $lnlikehood_id{$id};
	
	print $o_fh "\t";
	if (defined $lnlikehooddata) {
	    print $o_fh "$lnlikehooddata";
	}
    }
    print $o_fh "\n";
}

my $total = 0;
print STDERR "\tREPORT:\n\n";
foreach my $cl_id (sort {$a <=> $b} keys %cluster_types) {
    my $cl_n = scalar(@{$cluster_types{$cl_id}});
    $total += $cl_n;
    print STDERR "\t\tCLUSTER_$cl_id:\t$cl_n\n";
}
print STDERR "\n\t\tTotal Trees:\t$total\n";
print STDERR "\n\n";


$date = `date`;
chomp($date);
print STDERR "\n\nEND OF THE SCRIPT ($date).\n\n";


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
 
       TOPD/FMTS [Puigbo et al., 2007] is a program that 
     compares phylogenetic trees based in five different
     methods (nodal, split, quartets, triplets, disagree). 
     This script parse the result file and use these 
     comparisons to cluster the phylogenetic trees based in 
     their topologies.

       Depending of the method used for the comparisons, it is
     possible use different variables for clustering:

      -m nodal -v distance -c [0-1]
      -m split -v distance -c [0-1]
      -m quartets -v [s|d|r1|r2|u|DC|EA|SJA|SSJA|SD] -c [0-1]
      -m triplets -v [s|d|r1|r2|u|DC|EA|SJA|SSJA|SD] -c [0-1]
   
      -m disagree -v taxa -c 0 (by default)
     
    Usage:
  
      topd_clustering.pl [-h] -i <topd_result_file> 
                             [-t <tree_file>] 
			     [-l <lnlikehood_file>]
                             [-o <output>] 
                             [-m <topology_method>] 
                             [-c <cluster_variable>] 
                             [-v <cluster_cutoff_value>]
                              
    Example:

      topd_clustering.pl -i topd_sol.txt -c disagree

    Flags:
     
     -i <topd_result_file>    TOPD result file (mandatory)

     -t <tree_file>           Tree file used as input for 
                              TOPD if a new file with ID,
                              CLUSTER,TREE want to be 
                              created (optional)

     -l <lnlikehood_file>     Likehood file (-f1 ID, -f2 
                              LNLIKEHOOD_VALUE) to add 
                              ln_likehood value to the 
                              output file (optional).

     -o <output>              Output basename (by default:
                              (topd_clustering.output)

     -m <topology_method>     Method used for topologies 
                              comparisson (method used by 
                              TOPD to compare trees: nodal, 
                              split, quartets, triplets, 
                              disagree)(disagree by default)

      -v <method_variable>    Variable produced by the 
                              method used for clustering 
                              (distance, taxa, DC, EA...)
                              (taxa by default)

      -c <variable_cutoff>    Cutoff value to assign a 
                              relation as true for 
                              clustering (0 by default)

      -h <help>               Print the help  


EOF
exit (1);
}




####
1; #
####
