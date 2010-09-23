#!/usr/bin/perl

=head1 NAME

 reclustering.pl
 Script to assign new clusters names to old ones

=cut

=head1 SYPNOSIS

 reclustering.pl [-h] -i <cluster_file> -c <cluster_equivalences> [-o <output>]

=head1 EXAMPLE 

 reclustering.pl -i cluster.tab -c new_clusters.tab
                                  
=head2 I<Flags:>

=over


=item -i

B<cluster_file>                File with three columns -f1 ID, -f2 TYPE_ID, -f3 TREE (mandatory)

=item -c

B<cluster_equivalences>        File with to columns -f1 OLD_TYPE_ID, -f2 NEW_TYPE_ID (mandatory)

=item -o

B<output>                      Output basename (reclustering.output by default)

=item -h

B<help>                         Print the help

=back

=cut

=head1 DESCRIPTION

   Reclustering is a simple script to reassign new cluster types to
   a file with multiple tree topology types. Also it will print one file
   per new type.

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 reclustering.pl

=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;


our ($opt_i, $opt_c, $opt_o, $opt_h);
getopts("i:c:o:h");
if (!$opt_i && !$opt_c && !$opt_o && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Mandatory arguments (die if they are not supplied)

my $old_cluster = $opt_i 
    || die("\n\nMANDATORY ARGUMENT ERROR: -i <cluster_file> argument was not supplied.\n\n");
my $cluster_eq = $opt_c
    || die("\n\nMANDATORY ARGUMENT ERROR: -c <cluster_equivalences> argument was not supplied.\n\n");

## Default arguments (use default values if they are not supplied)

my $output = $opt_o 
    || 'reclustering.output';


##############################################
#### PARSING FILE  ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);

print STDERR "\n\nSCRIPT START POINT ($date)\n\n";

$date = `date`;
chomp($date);
print STDERR "\t1) Parsing cluster file: $old_cluster ($date)\n\n";

my %old_cluster = ();
my %old_tree = ();
    
my $t = 0;
my $T = `wc -l $old_cluster`;
chomp($T);
$T =~ s/$old_cluster//;
    
open my $oc_fh, '<', $old_cluster;
    
while(<$oc_fh>) {
    chomp($_);
    $t++;
	
    print STDERR "\t\tParsing line $t of $T lines in the file: $old_cluster\r";
	
    my @t_data = split(/\t/, $_);
    $old_cluster{$t_data[0]} = $t_data[1];
    $old_tree{$t_data[0]} = $t_data[2]
}
print STDERR "\n\n";


$date = `date`;
chomp($date);
print STDERR "\t2) Parsing equivalence file ($date): $cluster_eq\n\n";

open my $ce_fh, '<', $cluster_eq;
my $l = 0;
my $L = `wc -l $cluster_eq`;
chomp($L);
$L =~ s/$cluster_eq//;

my %cluster_eq = ();
my %new_cluster = ();

while (<$ce_fh>) {
    $l++;
    chomp($_);
    print STDERR "\t\tParsing line $l of $L lines in the file: $cluster_eq\r";

    my @e_data = split(/\t/, $_);
    $cluster_eq{$e_data[0]} = $e_data[1];

    unless (exists $new_cluster{$e_data[1]}) {
	$new_cluster{$e_data[1]} = 1;
    }

}
print STDERR "\n\n";

$date = `date`;
chomp($date);
print STDERR "\t3) Combining files ($date):\n\n";

my $old_n = scalar(keys %cluster_eq);
my $new_n = scalar(keys %new_cluster);

print STDERR "\t\tThere are $old_n old clusters types and $new_n new cluster types.\n";

my $general_outfile = $output . '.general.txt';

print STDERR "\n\t\tPrinting general file...\n";

my %clustering = ();

open my $ou_fh, '>', $general_outfile;

foreach my $id (sort keys %old_cluster) {
    my $oldcl = $old_cluster{$id};
    my $oldtree = $old_tree{$id};
    my $new_cluster = $cluster_eq{$oldcl};

    unless (defined $new_cluster) {
	print STDERR "\n\t\tWARNING: id=$id with cluster=$oldcl, has not cluster equivalence. SKIPPING.\n";
    }
    else {
	print $ou_fh "$id\t$new_cluster\t$oldtree\n";
	if (exists $clustering{$new_cluster}) {
	    push @{$clustering{$new_cluster}}, $id;
	}
	else {
	    $clustering{$new_cluster} = [$id];
	}
    }
}

print STDERR "\n\t\tPrinting specific files by types...\n";

foreach my $type (sort keys %clustering) {

    my $type_filename = $output . '.' . "$type" . '.txt';
    $type_filename =~ s/\s+/_/g;
    open my $ty_fh, '>', $type_filename;

    my @treeids = @{$clustering{$type}};

    foreach my $tree_id (sort @treeids) {
	my $tree = $old_tree{$tree_id};
	print $ty_fh "$tree_id\t$type\t$tree\n";
    }

    close $ty_fh;
}

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
 
        Reclustering is a simple script to reassign new cluster types to
      a file with multiple tree topology types. Also it will print one file
      per new type.
    
    Usage:
  
     reclustering.pl [-h] -i <cluster_file> -c <cluster_equivalences> [-o <output>]
                              
    Example:

     reclustering.pl -i cluster.tab -c new_clusters.tab

    Flags:
     -i <cluster_file>            File with three columns -f1 ID, -f2 TYPE_ID, -f3 TREE (mandatory)
     -c <cluster_equivalences>    File with to columns -f1 OLD_TYPE_ID, -f2 NEW_TYPE_ID (mandatory)
     -o <output>                  Output basename (reclustering.output by default)
     -h <help>                    Print the help  


EOF
exit (1);
}




####
1; #
####
