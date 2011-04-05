
package Strain::AlleleIdentification;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;

## To export some functions

use Exporter qw( import );

our @EXPORT_OK = qw( identify_alleles );


###############
### PERLDOC ###
###############

=head1 NAME

Strain::AlleleIdentification.pm
An object to identify alleles using tree structure.

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use Strain::AlleleIdentification qw/ identify_alleles /;

  my %allele_args = (
    strains      => \%strains,
    alignment    => $alignment,
    tree         => $tree,
    parents      => \%parents,
    target       => $target,
    filter       => \%filter,
  );

  my %alleles = identify_alleles(\%allele_args);


=head1 DESCRIPTION

 Function to identify the alleles in a polyploid/heterozigous species 
 knowing its parents.

 Steps:
              |1
       -------+-------
       |2            |3
    ---+---       ---+---
    |     |       |     |
    T1    A       T2    B

  1) Get the leaves in the tree of the target specie. (T1 and T2)
  
  2) Get the ancestor node for these leaves and get all the leaves for
     this ancestor (T1 => 2 => T1,A and T2 => 3 => T2,B),
 
  3) Compare the leaves under the same node-ancestor. If there are
     just only one parent and one target it will identify that target
     as the allele that comes from the parent (T1-A and T2-B)

  4) Calculate the identity and the possible SNPs between alleles and
     parents.

  5) Apply the filter according the following parameters:
     i)   alignment length
     ii)  alignment identity percentage

  If the allele is in another branch ir will not be identify as allele

                 |1
           ------+------                   only T3 will be identify as 
           |2          |                   allele of B
       ----+----       |
       |3      |       |4
    ---+---    |    ---+---
    |     |    |    |     |
    T1    T2   A    T3    B

  
  

=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 identify_alelles

  Usage: my %alleles = identify_alleles(\%allele_args);

  Desc: Analyze a phylogenetic tree with some members and its parents or
        close related ancestors (see description for more details)

  Ret: A hash with key=member_id, value=allele_name

  Args: A hash reference with the following key-value pairs: 
         strains      => a hashref. with key=member_id, value=strains
         alignment    => $alignment, a Bio::SimpleAlign object
         tree         => $tree, a Bio::Tree::Tree object
         parents      => a hashref with key=parent/ancestors strain
                                        value=allele name
         target       => $target, a strain-target name
         filter       => a hashref. with key=(length|identity), 
                                         value=cutoff value
        
  Side_Effects: Die if the no argument is used or it is not a hash ref.
                Die if no strains,alignment,tree,parents or target argument
                is used.
                Die if the format of the argument is not right.
                Return an empty hash if not alleles were identified.

  Example: my %alleles = identify_alleles(\%allele_args);

=cut

sub identify_alleles {
    my $args_href = shift ||
	croak("ERROR: No arguments were supplied to identify_alleles()");

    my %args = ();
    if (ref($args_href) ne 'HASH') {
	croak("ERROR: Argument supplied to identify_alleles() isnt a HASHREF.");
    }
    else {
	%args = %{$args_href};
    }

    ## Check that the arguments have the right format

    my %perms = (
	strains   => 'HASH',
	alignment => 'Bio::SimpleAlign',
	tree      => 'Bio::Tree::Tree',
	parents   => 'HASH',
	target    => undef,
	filter    => 'HASH',
	);

    foreach my $arg (keys %args) {
	unless (exists $perms{$arg}) {
	    croak("ERROR: $arg is a non-permited argum. for identify_alleles");
	}
	else {
	    if (defined $perms{$arg} && ref($args{$arg}) ne $perms{$arg}) {
		croak("ERROR: $arg doesnt have right value identify_alleles()");
	    } 
	}
    }

    ## Check the mandatory arguments (everything except filter)
    
    foreach my $marg (keys %perms) {
	
	if ($marg ne 'filter') {
	    unless (exists $args{$marg}) {
		croak("ERROR: Mandatory argument $marg wasnt used as argument");
	    }
	}
	else {
	    if (exists $args{filter}) {
		foreach my $fi (keys %{$args{filter}}) {
		    if ($fi !~ m/^(length|identity)$/) {
			croak("ERROR: Filter arg. has not the right arg ($fi)");
		    }
		}
	    }
	}
    }

    ## Declare possible alleles hash

    my %pos_alleles = ();
    my %pos_alleles_neigh = ();

    ## Everything has been checked, now it will start the allele identification

    my %strains = %{$args{strains}};
    my %parents = %{$args{parents}};
    
    if (scalar(keys %parents) < 2) {
	croak("ERROR: Less than one parent was specified as argument.");
    }
    
    my $tree = $args{tree};

    ## 1) Get the ancestors nodes for all the target leaves
    ##    Get the leaves for these nodes and check if they are parents
   
    my @target_members = ();
  
    foreach my $node ( $tree->get_nodes() ) {
	
	if ($node->is_Leaf()) {

	    my $node_id = $node->id();
	    if (exists $strains{$node_id}) {
		if ($strains{$node_id} eq $args{target}) {

		    push @target_members, $node_id;

		    my $ancestor = $node->ancestor();
		    my $parents_match = 0;
		    my $parent_allele;
		    
		    foreach my $descent ($ancestor->get_all_Descendents()) {
			
			if ($descent->is_Leaf()) {
			    
			    if (exists $strains{$descent->id()}) {

				my $desc_strain = $strains{$descent->id()};
				if (exists $parents{$desc_strain}) {
				    $parents_match++;
				    $parent_allele = $descent->id();
				}
			    }
			}
		    }

		    if ($parents_match == 1) { ## Add only of one of the parents
                                               ## was found

			my $parentallele = $parents{$strains{$parent_allele}};
			$pos_alleles{$node_id} = $parentallele;
			$pos_alleles_neigh{$node_id} =  $parent_allele;
		    }
		}
	    }
	}
    }

    my %alleles = ();

    ## 2) The alleles should be identified, so now it will get some 
    ##    values from the alignment

    my $align = $args{alignment};
    
    if (defined $args{filter}) {
    
	foreach my $target (keys %pos_alleles) {

	    my $newalign = Bio::SimpleAlign->new();
	    
	    my $seq = $align->get_seq_by_id($target);
	    unless (defined $seq) {
		croak("ERROR: No seq was found in $align for id=$target.")
	    }
	    else {
		$newalign->add_seq($seq);
	    }
	    
	    my $par =  $align->get_seq_by_id($pos_alleles_neigh{$target});
	    unless (defined $par) {
		croak("ERROR: No seq was found id=$pos_alleles_neigh{$target}.")
	    }
	    else {
		$newalign->add_seq($par);
	    }

	    my $cond = scalar(keys %{$args{filter}});
	    if (exists $args{filter}->{length}) {
		if ($newalign->length >= $args{filter}->{length}) {
		    $cond--;
		}
	    }
	    
	    if (exists $args{filter}->{identity}) {

		my $req_ident =  $args{filter}->{identity};
		if ($newalign->percentage_identity >= $req_ident){
		    $cond--;
		}
	    }
	    if ($cond == 0) {
		$alleles{$target} = $pos_alleles{$target};
	    }
	}
    }
    else {
	%alleles = %pos_alleles;
    }
     
    return %alleles;
}

####
1; #
####
