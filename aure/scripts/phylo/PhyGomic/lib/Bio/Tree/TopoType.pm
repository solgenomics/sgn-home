
package Bio::Tree::TopoType;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Try::Tiny;
use Math::BigFloat;

use Bio::Tree::Tree;

###############
### PERLDOC ###
###############

=head1 NAME

TopoType.pm
A class to handle a tree topology type.

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use Bio::Tree::TopoType;

  my $topotype = Bio::Tree::TopoType->new({ topology     => $tree,
                                            members      => $trees_aref,
                                            descriptions => $description_href, 
                                   });


  ## Accessors

  my $topology_tree = $topotype->get_topology();
  $topotype->set_topology();

  my $members_aref = $topotype->get_members();
  $topotype->set_members($members_aref);
  $topotype->add_members($new_members_aref);
  $topotype->delete_members($old_members_aref);

  my $topology_description_href = $topotype->get_description();
  my $tag_description = $topology->get_description($tag);
  $topotype->set_description(\%descriptions);
  $topology->add_description($tag, $description);
  $topology->delete_description($tag);

  ## Analysis

  my $member_count = $topotype->members_count();
  my @species = $topotype->species_extraction();
  my %homologous_count = $topotype->homologous_count();
  my %paralogous_count = $topotype->paralogous_count();
  my %gene_duplications = $topotype->gene_duplications();
  my %gene_loss = $topotype->gene_loss();


=head1 DESCRIPTION

 Object to handle the topologies for trees


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $topotype = Bio::Tree::TopoType->new($arguments_href);

  Desc: Create a topotype object with the specified parameters.

  Ret: a Bio::Tree::TopoType object

  Args: A hash reference with the following key-value pairs: 
          topology     => $tree, a Bio::Tree::Tree object to represent the 
                         topology
          members      => \@trees, an arrayref of Bio::Tree::Tree members with
                         the same topology
          descriptions => $descriptions, a hash reference with key=tag and
                          values=description
        
  Side_Effects: Die if the argument used is not a hash or its values arent 
                right.

  Example: my $topotype = Bio::Tree::TopoType->new({ 
                                                      topology => $treetype,
                                                      members  => \@trees, 
                                                  });

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %permargs = ( 
	topology     => 'Bio::Tree::Tree',
	members      => 'ARRAY',
	descriptions => 'HASH',
	);

    ## Check argument

    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ARGUMENT ERROR: $args_href used for new() is not HASH REF.");
	}
	else {
	    my %args = %{$args_href};
	    foreach my $argkey (keys %args) {
		my $exp = $permargs{$argkey};
		unless (defined $exp) {
		    croak("ARG. ERROR: $argkey isnt a permited arg. for new()");
		}
		else {
		    if (ref($args{$argkey}) ne ref($exp)) {
			croak("ARG. ERROR: $args{$argkey} isnt $exp for new()");
		    }
		}	
	    }
	}
    }

    my $topology = $args_href->{topology} || '';
    my $members_aref = $args_href->{members} || [];
    my $descriptions_href = $args_href->{descriptions} || {};

    ## Set vars in the object

    $self->set_topology($topology);
    $self->set_members($members_aref);
    $self->set_descriptions($descriptions_href);

    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get_topology

  Usage: my $topology = $topotype->get_topology(); 

  Desc: Get the topology from the object as Bio::Tree::Tree object

  Ret: A Bio::Tree::Tree object

  Args: None

  Side_Effects: None

  Example: my $topology = $topotype->get_topology(); 

=cut

sub get_topology {
    my $self = shift;
    return $self->{topology};
}

=head2 set_topology

  Usage: $phygetopo->set_topology($tree);

  Desc: Set the topology in a Bio::Tree::topoType object

  Ret: None

  Args: A Bio::Tree::Tree object

  Side_Effects: Die if no argument is used or if the object isnt Bio::Tree::Tree
                object

  Example: $phygetopo->set_topology($tree);

=cut

sub set_topology {
    my $self = shift;
    my $topology = shift;
   
    unless (defined $topology) {
	croak("ARG. ERROR: No arg. was used for set_topology() function");
    }
    else {
	if ($topology =~ m/\w+/) {
	    my $obj = 'Bio::Tree::Tree';
	    unless (ref($topology) eq $obj) {
		croak("ARG. ERROR: $topology at set_topology() isnt a $obj");
	    }
	}
	$self->{topology} = $topology;    
    }
}

=head2 get_members

  Usage: my $members_aref = $topotype->get_members(); 

  Desc: Get the members from the Bio::Tree::TopoType object

  Ret: An array reference of Bio::Tree::Tree objects.

  Args: None

  Side_Effects: None

  Example: my @members = @{$topotype->get_members()}; 

=cut

sub get_members {
    my $self = shift;
    return $self->{members};
}

=head2 set_members

  Usage: $phygetopo->set_members(\@members);

  Desc: Set the arrayref. of members in a Bio::Tree::TopoType object

  Ret: None

  Args: An array reference of Bio::Tree::Tree objects.

  Side_Effects: Die if no argument is used or if the object isnt an array of 
                Bio::Tree::Tree objects

  Example: $phygetopo->set_members( [ $tree1, $tree2, $tree3 ] );

=cut

sub set_members {
    my $self = shift;
    my $members_aref = shift;

    unless (defined $members_aref) {
	croak("ARG. ERROR: No arg. was used for set_members() function");
    }
    else {
	if ($members_aref =~ m/\w+/) {
	    my $obj = 'Bio::Tree::Tree';
	    unless (ref($members_aref) eq 'ARRAY') {
		croak("ARG. ERROR: $members_aref set_members() isnt ARRAYREF");
	    }
	    else {
		foreach my $elem (@{$members_aref}) {
		    unless (ref($elem) eq $obj) {
			croak("ARG. ERROR: $elem used set_members isnt $obj");
		    }
		}
	    }
	}
	$self->{members} = $members_aref;    
    }
}


=head2 add_members

  Usage: $phygetopo->add_members(\@new_members);

  Desc: Add members to the arrayref. of members in a Bio::Tree::TopoType 
        object

  Ret: None

  Args: An array reference of Bio::Tree::Tree objects.

  Side_Effects: Die if no argument is used or if the object isnt an array of 
                Bio::Tree::Tree objects

  Example: $phygetopo->add_members( [ $tree4, $tree5 ] );

=cut

sub add_members {
    my $self = shift;
    my $members_aref = shift;
    my $obj_members_aref = $self->get_members();   
   
    unless (defined $members_aref) {
	croak("ARG. ERROR: No arg. was used for add_members() function");
    }
    else {
	if ($members_aref =~ m/\w+/) {
	    my $obj = 'Bio::Tree::Tree';
	    unless (ref($members_aref) eq 'ARRAY') {
		croak("ARG. ERROR: $members_aref add_members() isnt ARRAYREF");
	    }
	    else {
		foreach my $elem (@{$members_aref}) {
		    unless (ref($elem) eq $obj) {
			croak("ARG. ERROR: $elem used add_members isnt $obj");
		    }
		    else {
			push @{$obj_members_aref}, $elem;
		    }
		}
	    }
	}
    }
}

=head2 delete_members

  Usage: my @deleted_members = $phygetopo->delete_members(\@members);

  Desc: Delete members of the arrayref. of members in a Bio::Tree::TopoType 
        object

  Ret: None

  Args: An array reference of Bio::Tree::Tree objects.

  Side_Effects: Die if no argument is used or if the object isnt an array of 
                Bio::Tree::Tree objects

  Example: my @deleted = $phygetopo->delete_members( [ $tree4 ] );

=cut

sub delete_members {
    my $self = shift;
    my $memb_aref = shift;
    my $obj_members_aref = $self->get_members();   
   
    my @removed = ();
    unless (defined $memb_aref) {
	croak("ARG. ERROR: No arg. was used for delete_members() function");
    }
    else {
	if ($memb_aref =~ m/\w+/) {
	    my $obj = 'Bio::Tree::Tree';
	    unless (ref($memb_aref) eq 'ARRAY') {
		croak("ARG. ERROR: $memb_aref delete_members() isnt ARRAYREF");
	    }
	    else {	
		foreach my $elem (@{$memb_aref}) {
		    unless (ref($elem) eq $obj) {
			croak("ARG. ERROR: $elem used set_members isnt $obj");
		    }
		    else {
			my @selected = ();
			while (my $old = shift(@{$obj_members_aref})) { 
			    if ($old ne $elem) {
				push @selected, $old;
			    }
			    else {
				push @removed, $old;
			    }
			}
			push @{$obj_members_aref}, @selected;
		    }
		}

	    }
	}
    }
    return @removed;
}

=head2 get_descriptions

  Usage: my $descriptions_href = $topotype->get_descriptions(); 
         my $tag_description = $topotype->get_descriptions($tag);         

  Desc: Get the descriptions from the Bio::Tree::TopoType object

  Ret: An hash reference with key=tag and value=description or a scalar
       if a tag is used as argument

  Args: $tag, a scalar [optional]

  Side_Effects: None

  Example: my %descriptions = %{$topotype->get_descriptions()}; 
           my $evol_description = $topotype->get_descriptions('evol');  

=cut

sub get_descriptions {
    my $self = shift;
    my $tag = shift;

    if (defined $tag) {
	return $self->{descriptions}->{$tag};
    }
    else {    
	return $self->{descriptions};
    }
}

=head2 set_descriptions

  Usage: $topotype->set_descriptions(\%descriptions);

  Desc: Set the hashref. of descriptions in a Bio::Tree::TopoType object

  Ret: None

  Args: An hash reference with key=tag and value=descrition

  Side_Effects: Die if no argument is used or if the object isnt an hashref

  Example: $topotype->set_descriptions({ expected_topology => 1 } );

=cut

sub set_descriptions {
    my $self = shift;
    my $desc_href = shift;

    unless (defined $desc_href) {
	croak("ARG. ERROR: No arg. was used for set_descriptions() function");
    }
    else {
	if ($desc_href =~ m/\w+/) {
	    unless (ref($desc_href) eq 'HASH') {
		croak("ARG. ERROR: $desc_href set_descriptions() isnt HASHREF");
	    }
	}
	$self->{descriptions} = $desc_href;    
    }
}


=head2 add_description

  Usage: $topotype->add_description($tag, $description);

  Desc: Add a new pair tag/description to a Bio::Tree::TopoType object

  Ret: None

  Args: $tag, a scalar as tag (for example: 'global', 'mytag', 'foo')
        $description, a scalar to describe something under an specific tag

  Side_Effects: Die if no argument is used.
                Edit a tag, if the tag used exists into the object.

  Example: $topotype->add_description('global', 'expected topology');

=cut

sub add_description {
    my $self = shift;
    my $tag = shift ||
	croak("ARG. ERROR: No tag arg. was used with add_description()");

    my $description = shift;   ## It doesnt use || because 0 is reconized 
                               ## as undef

    unless (defined $description) {
	croak("ARG. ERROR: No description arg was used with add_description()");
    }


    my $descr_href = $self->get_descriptions();
    $descr_href->{$tag} = $description;
}

=head2 delete_description

  Usage: my $deleted_description = $phygetopo->delete_description($tag);

  Desc: Delete tag/description pair from a Bio::Tree::TopoType object
        
  Ret: $deleted_description, a scalar

  Args: $tag, a scalar

  Side_Effects: Die if no argument is used.
                Return undef in the tag used does not exists into the object

  Example: my $deleted_evdescription = $phygetopo->delete_description('ev');

=cut

sub delete_description {
    my $self = shift;
    my $tag = shift ||
	croak("ARG. ERROR: No tag arg. was used with delete_description()");
   
    my $deleted;
    my $descr_href = $self->get_descriptions();
    
    if (exists $descr_href->{$tag}) {
	$deleted = delete($descr_href->{$tag});
    }
    return $deleted;
}


####
1; #
####
