#!/usr/bin/perl

=head1 NAME

 blast_clustering.pl
 Parse selfblast m8 format to create clusters (version.0.1).

=cut

=head1 SYPNOSIS

 blast_clustering.pl [-h] -i <input_file> -f <filter_conditions> > output.tab

=head2 I<Flags:>

=over


=item -i

B<input_file>             input file (mandatory)

=item -f 

B<filter_conditions>      filter conditions writen as field_name,condition,value (example: e_value,<,1e-100)
                          separated by semicolon

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script parse the results of a selfblast in m8 format and cluster the ids according the filtering conditions
 (if there are any).

 A self blast can made after the formating of the dataset (formatdb --help), using as -d and -i the same dataset

 The m8 column names are: 
  -f1: Query_id
  -f2: Subject_id
  -f3: Identity
  -f4: Align_length
  -f5: Mismatches
  -f6: Gaps_open
  -f7: Q_start
  -f8: Q_end
  -f9: S_start
  -f10: S_end
  -f11: E_value
  -f12: Hit_score
  
 Example:

  blast_clustering.pl -i test.selfblastn.m8 -f 'e_value,<,1e-100;identity,>,95'


=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 blast_clustering.pl

=cut

use strict;
use warnings;

use Getopt::Std;
use Math::BigFloat;


our ($opt_i, $opt_f, $opt_h);
getopts("i:f:h");
if (!$opt_i && !$opt_f && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}


## Take the input file data

my $infile = $opt_i || die("DATA ARGUMENT -i input_file WAS NOT SUPPLIED.\n");

## Parse the filtering and store in a hash key=field_name and value=condition+value
## The filter only can have specific names

my @filter_names = ('Query_id', 'Subject_id', 'Identity', 'Align_length', 'Mismatches', 'Gaps_open',
		    'Q_start', 'Q_end', 'S_start', 'S_end', 'E_value', 'Hit_score');

my %filters;

if (defined $opt_f) {

    print STDERR "\n\nFILTER OPTION: [Enabled].\n";
    my @filters = split(/;/, $opt_f);

    foreach my $fil (@filters) {
	
	my $match = 0;
	foreach my $perm_name (@filter_names) {
	    
	    my @filt_data = split(/,/, $fil);

	    if ($filt_data[0] =~ m/$perm_name/i) {
		$filters{$perm_name} = [$filt_data[1], $filt_data[2]];
		$match = 1;
	    }
	}
	if ($match == 0) {
	    print STDERR "\nWARNING: The filter condition:$fil does not match with the permited names. It will skip it.\n";
	}
    }
}
else {
    print STDERR "\n\nFILTER OPTION: [Disabled].\n";
}

## Now parse the blast
## It will define an hash with keys=id and value=array ref with the hits
## This will produce redundant clusters, one per id, so after that, it will
## remove the redundancy

open my $bfh, '<', $infile 
    || die("\nOPEN FILE ERROR: The file:$infile can not be openned (system error: $!).\n\n");

my %precluster;


print STDERR "\n\nPARSING SELFBLAST FILE.\n";

my ($l, $f) = (0, 0);

while (<$bfh>) {
    chomp($_);
    $l++;

    my @bl_dat = split(/\t/, $_);
   
    my $n = 0;
    my $filt = 0;
    
    print STDERR "Parsing line:$l        \r";

    if (defined $bl_dat[0] && defined $bl_dat[1]) {

	foreach my $bl (@bl_dat) {

	    my $field = $filter_names[$n];
	    $n++;
	    
	    ## It will filter according the condition

	    ## In a self blast it will contain the same id, it will added as first element
	    ## To remove the redundancy ( multiple hits with the same id, it will use an hash)
	    
	    if ($opt_f) {
		if (exists $filters{$field}) {
		    my $cond = $filters{$field}->[0];
		    my $val = $filters{$field}->[1];

		    if (compare_cond($bl, $cond, $val)) {
			$filt = 1;
			if (exists $precluster{$bl_dat[0]}) {
			    $precluster{$bl_dat[0]}->{$bl_dat[1]} = 1;
			}
			else {
			    $precluster{$bl_dat[0]} = { $bl_dat[1] => 1 };
			}
		    }
		}
	    }
	    else {
		
		if (exists $precluster{$bl_dat[0]}) {		    
		    unless (exists $precluster{$bl_dat[0]}->{$bl_dat[1]}) {
			$precluster{$bl_dat[0]}->{$bl_dat[1]} = 1;
		    }
		}
		else {
		    $precluster{$bl_dat[0]} = { $bl_dat[1] => 1 };		    
		}
	    }
	}
	if ($filt == 0) {
	    $f++;
	}
    }
}

print STDERR "\n\tDone.\n";
if ($opt_f) { 
    print STDERR "\t$f hits have been rejected according filter condition ($opt_f).\n";
}

## To filter the percluster it will remove the redundancy

my %cluster;
my %ord_cluster;

my ($cl, $s, $c, $max) = (0, 0, 0, 0);
my $max_name = '';

print STDERR "\n\nCLUSTERING.";

foreach my $id (sort keys %precluster) {
   
    if (defined $precluster{$id}) {
	$cl++;
	my %ids_c = %{$precluster{$id}};

	my $el_c = scalar(keys %ids_c);
		    
	my $name;
	if ($el_c == 1) {
	    $s++;
	    $name = 'singlet_' . $cl;
	}
	else {
	    $c++;
	    $name = 'cluster_' . $cl; 
	}
	if ($el_c > $max) {
	    $max = $el_c;
	    $max_name = $name;
	}

	if (exists $ord_cluster{$el_c}) {
	    push @{$ord_cluster{$el_c}}, $name;
	}
	else {
	    $ord_cluster{$el_c} = [$name];
	}

	foreach my $id_c (sort keys %ids_c) {
	  
	    if (exists $cluster{$name}) {
		push @{$cluster{$name}}, $id_c
	    }
	    else {
		$cluster{$name} = [$id_c];
	    }
	
	
	    ## And delete the id from precluster
	
	    delete $precluster{$id_c};
	}      
    }
}

## Finally print the results ordered by number of elements

print STDERR "\n\tDone.\n\t$c clusters and $s singlets has been created\n";
print STDERR "\tBiggest cluster: $max_name with $max elements.\n";

foreach my $cl (reverse sort {$a <=> $b} keys %ord_cluster) {
    my @cl_names = @{$ord_cluster{$cl}};
    foreach my $cl_name (@cl_names) {
	my @ids = @{$cluster{$cl_name}};
	foreach my $id (@ids) {
	    print STDOUT "$cl_name\t$id\n";
	}
    }
}
print STDERR "\n";

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
  
     This script parse the results of a selfblast in m8 format and cluster the ids according the filtering conditions
     (if there are any).

     A self blast can made after the formating of the dataset (formatdb --help), using as -d and -i the same dataset

     The m8 column names are: 
       -f1: Query_id
       -f2: Subject_id
       -f3: Identity
       -f4: Align_length
       -f5: Mismatches
       -f6: Gaps_open
       -f7: Q_start
       -f8: Q_end
       -f9: S_start
       -f10: S_end
       -f11: E_value
       -f12: Hit_score

       
    Usage:
        
      blast_clustering.pl [-h] -i <input_file> -f <filter_conditions> > output.tab
  
    Example:

      blast_clustering.pl -i test.selfblastn.m8 -f 'e_value,<,1e-100;identity,>,95' > output.tab
      
    Flags:

      -i <input_file>             input file (mandatory)
      -f <filter_conditions>      filter conditions writen as field_name,condition,value (example: e_value,<,1e-100)
                                  separated by semicolon
      -h <help>                   print the help
     
EOF
exit (1);
}

=head2 compare_cond

  Usage: compare_cond($val1, $cond, $val2)
  Desc: compare two values using a var as operator
  Ret: 1=true and 0=false
  Args: Three scalar, two values and a condition
        condition can be: <,>,=,<= or >=
  Side_Effects: Use equal by default
  Example: if (compare_cond($val1, $cond, $val2)) {
             ## Do something
           }

=cut

sub compare_cond {
    my $val1 = shift;
    my $cond = shift;
    my $val2 = shift;

    my $c = 0;

    ## It will use BigFloat to compare values as e-value
    ## with scientific notation

    my $v1 = Math::BigFloat->new($val1);
    my $v2 = Math::BigFloat->new($val2);
    my $comp = $v1->bcmp($v2);

    if ($cond eq '>' && $comp > 0) {
	$c = 1;
    }
    elsif ($cond eq '<' && $comp < 0) {	
	$c = 1;
    }
    elsif ($cond eq '=' && $comp == 0) {
	$c = 1;
    }
    elsif ($cond eq '>=' && $comp >= 0) {
	$c = 1;
    }
    elsif ($cond eq '<=' && $comp <= 0) {
	$c = 1;
    }
    else {
	if ($comp == 0) {
	    $c = 1;
	}
    }
    
    return $c;
}
