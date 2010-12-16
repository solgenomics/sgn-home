
package Bio::Tree::TopoType;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Try::Tiny;
use Math::BigFloat;

use Bio::Tree::Tree;
use Bio::TreeIO::newick;

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
		    if (ref($args{$argkey}) ne $exp) {
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


###################
## OTHER METHODS ##
###################

=head2 is_same_tree

  Usage: my $true = Bio::Tree::TopoType::is_same_tree($tree1, $tree2, $diff);

  Desc: Compare two trees and return if is the same tree (1) or not (0)
        If return differences is used, it will return an array with the
        differences
        
  Ret: $true, a scalar 0 (false) or 1 (true)
       @differences, an array with differences

  Args: Two trees (Bio::Tree::Tree objects), 
        An array of strings describing differences (id diff is used)

  Side_Effects: Die if no argument is used or if they are no tree objects

  Example: if ( Bio::Tree::TopoType::is_same_tree($tree1, $tree2) ) {
              ## Do something
           }
           my @diff = Bio::Tree::TopoType::is_same_tree($tree1,$tree2,'getdiff')

=cut

sub is_same_tree {
    my $tree1 = shift ||
	croak("ARG. ERROR: First tree was not supplied to is_same_tree()");
    my $tree2 = shift ||
	croak("ARG. ERROR: Second tree was not supplied to is_same_tree()");
    my $get_diff = shift;


    if (ref($tree1) ne 'Bio::Tree::Tree') {
	croak("ARG. ERROR: 1st tree supplied to is_same_tree() isnt Tree Obj");
    }
    elsif (ref($tree2) ne 'Bio::Tree::Tree') {
	croak("ARG. ERROR: 2nd tree supplied to is_same_tree() isnt Tree Obj");
    }

    my $treename1 = $tree1->id || 'tree1';
    my $treename2 = $tree2->id || 'tree2';
    my %objstree = ( $treename1 => $tree1, $treename2 => $tree2 );
    my %hashtree = ();

    ## Decompose the trees in leaves, inode and ids

    foreach my $treeid (keys %objstree) {
	my $tree = $objstree{$treeid};
	
	my %parts = ( leaves => {}, inodes => {}, ids => {} );
	my $in_idx = 1;
	my $lf_idx = 1;
	my @nodes1 = $tree->get_nodes();
	foreach my $node (@nodes1) {
	    if ($node->is_Leaf()) {
		$parts{leaves}->{$lf_idx} = { id => $node->id(),
					      br => $node->branch_length() };
		$parts{ids}->{$node->id()} = $lf_idx;
		$lf_idx++;
	    }
	    else {
		my @desc_nodes = ();
		foreach my $desc ($node->each_Descendent()) {
		    my $dsc_id = $desc->id() || 'internalnode';
		    push @desc_nodes, $dsc_id;
		}
		my $desc_line = join(',', sort(@desc_nodes));

		$parts{inodes}->{$in_idx} = { br => $node->branch_length(), 
		                              ds => $desc_line };
		$in_idx++;
	    }
	}
	$hashtree{$treeid} = \%parts;
    }

    ## Now it will compare both trees

    ## 1) First if they have different leaves ids
    
    my @id_errors = ();
    my %tree1ids = %{$hashtree{$treename1}->{ids}};
    my %tree2ids = %{$hashtree{$treename2}->{ids}};

    foreach my $leaf1 (sort {$a cmp $b} keys %tree1ids) {
	unless (exists $tree2ids{$leaf1}) {
	    push @id_errors, "ABSENT ID:$leaf1 FOR:$treename2";
	}
	else {
	    if ($tree1ids{$leaf1} != $tree2ids{$leaf1}) {
		my $msg = "DIFFERENT POSITION ID:$leaf1 ";
		$msg .= "FOR:$treename1 POSITION:$tree1ids{$leaf1}";
		$msg .= "FOR:$treename2 POSITION:$tree2ids{$leaf1}";
		push @id_errors, $msg;
	    }
	}
    }
    foreach my $leaf2 (sort {$a cmp $b} keys %tree2ids) {
	unless (exists $tree1ids{$leaf2}) {
	    push @id_errors, "ABSENT ID:$leaf2 FOR:$treename1";
	}
    }

    ## 2) Second, if they have different descendent and branch length for 
    ##    internal nodes

    my %tree1in = %{$hashtree{$treename1}->{inodes}};
    my %tree2in = %{$hashtree{$treename2}->{inodes}};

    foreach my $idx1 (sort {$a <=> $b} keys %tree1in) {
	unless (exists $tree2in{$idx1}) {
	    push @id_errors, "WITHOUT INTERNAL NODE for POSITION: $idx1";
	}
	else {
	    my $desc1 = $tree1in{$idx1}->{ds};
	    my $desc2 = $tree2in{$idx1}->{ds};
	    if ($desc1 ne $desc2) {
		my $msg = "DIFFERENT DESCENDENTS FOR INT-NODE $idx1 BETWEEN ";
		$msg .= "TREES: $treename1 ($desc1) AND $treename2 ($desc2)";
		push @id_errors, $msg;
	    }

	    my $brl1 = $tree1in{$idx1}->{br};
	    my $brl2 = $tree2in{$idx1}->{br};
	    if (defined $brl1 && defined $brl2 && $brl1 ne $brl2) {
		my $msg = "DIFFERENT BRANCH LENGTH FOR INT-NODE $idx1 BETWEEN ";
		$msg .= "TREES: $treename1 ($brl1) AND $treename2 ($brl2)";
		push @id_errors, $msg;
	    }	    
	}
    }

    ## 3) Third, if they have different leaves and branch distances

    my %tree1lf = %{$hashtree{$treename1}->{leaves}};
    my %tree2lf = %{$hashtree{$treename2}->{leaves}};

    foreach my $idx1 (sort {$a <=> $b} keys %tree1lf) {
     	unless (exists $tree2lf{$idx1}) {
    	    push @id_errors, "ABSENT LEAF $tree1lf{$idx1} for TREE: $treename2";
    	}
	else {
     	    my $leaf1 = $tree1lf{$idx1}->{id};
     	    my $leaf2 = $tree2lf{$idx1}->{id};
     	    if ($leaf1 ne $leaf2) {
    		push @id_errors, "DIFFERENT LEAF NODES: $idx1 ($leaf1 $leaf2)";
    	    }
	    else {
		my $brl1 = $tree1lf{$idx1}->{br};
		my $brl2 = $tree2lf{$idx1}->{br};
		if ($brl1 ne $brl2) {
		    my $msg = "DIFFERENT BRANCH LENGTH FOR LEAF $leaf1 BETWEEN";
		    $msg .= "TREES: $treename1 ($brl1) AND $treename2 ($brl2)"
		}
	    }	    
     	}
    }

    ## Finally it will sum all the errors.

    my $issame = 1;
    if (defined $get_diff && $get_diff =~ m/diff|1/) {
	return @id_errors;
    }
    else {
	if (scalar(@id_errors) > 0) {
	    $issame = 0;
	}
	return $issame;
    }


}

=head2 _make_topotype

  Usage: my $topotype = Bio::Tree::TopoType::_make_topotype($args_href);

  Desc: Make a topology tree (node_id=strains and branch_length=0|1) from
        a tree argument. The branch length of these topologies can be 
        changed using branch_cutoffs argument. For example if branch_length
        => { 0.1 => 0, 0.5 => 1 } is used, all the branch length with values 
        less or equal to 0.1 will be replaced by 0, all the values less or 
        qual than 05 will be replaced by 1 and the rest by 2.
        
  Ret: $topotype, a Bio::Tree::Tree object

  Args: A hash reference with:
        tree           => $tree, a Bio::Tree::Tree object
        strains        => \%strains, a hash reference with key   = Leaves_ids 
                                                           value = strains
        branch_cutoffs => \%cutoffs, a hash reference with key   = upper_cutoff
                                                           value = new_value

  Side_Effects: Die if no argument is used or if the does not exist the 
                strain for a leaf id.

  Example: my $topotype = Bio::Tree::TopoType::_make_topotype($tree, \%strains);

=cut

sub _make_topotype {
    my $args_href = shift ||
	croak("ARG. ERROR: No args. hashref. was used with _make_topotype()");
    
    ## Check variables

    unless(ref($args_href) eq 'HASH') {
	croak("ARG. ERROR: $args_href added to _make_topotype() isnt HASHREF");
    }

    my %permargs = ( 
	tree           => 'Bio::Tree::Tree',
	strains        => 'HASH',
	branch_cutoffs => 'HASH',
	);

    foreach my $argkey (keys %{$args_href}) {
	unless (exists $permargs{$argkey}) {
	    croak("ARG. ERROR: $argkey isnt a permited key for _make_topotype");
	}
	else {
	    if (ref($args_href->{$argkey}) ne $permargs{$argkey}) {
		croak("ARG. ERROR: Value for $argkey isnt permited value");
	    }
	}
    }
 
    ## Check args and use default values

    my $tree = $args_href->{'tree'} ||
	croak("ARG. ERROR: No tree arg. was used with _make_topotype()");
    my $strains_href = $args_href->{'strains'} ||
	croak("ARG. ERROR: No strain arg. was used with _make_topotype()");
    my $branch_cfs_href = $args_href->{'branch_cutoffs'} ||
	{ 0.01 => 0 };

    my @new_nodes = ();
    my @nodes = $tree->get_nodes();

    ## Order nodes by leaf names or descendents if it is an internal node
    ## for example a tree like ((A:1,B:1):0.9,(D:1,C:1):0,7) will be order as
    ## A,B,C,D,inode:0,7 (redefined as ~0) and inode:0.9 (redefined as ~1).
    ## This means that the internal nodes in the new tree will be order by
    ## original distance also.

    my %nodes = ();
    my %onodes = ();  ## Key=$internalid and Value=$id or inode_id
    my $inode_idx = 0;
    my $glob_idx = 0;
    foreach my $nod (@nodes) {
	$glob_idx++;
	if ($nod->is_Leaf()) {
	    $nodes{$nod->id()} = $nod;
	    $onodes{$nod->internal_id()} = $nod->id() . '_' . $glob_idx;
	}
	else {
	    my @desc = ();
	    foreach my $desc ($nod->each_Descendent()) {
		if ($desc->is_Leaf()) {
		    push @desc, $desc->id();
		}
		else {
		    push @desc, '~' . $inode_idx; 
		    $inode_idx++;
		}
	    }
	    my $iname = join(',', sort @desc) . '_' . $glob_idx;
	    $nodes{$iname} = $nod;
	    $onodes{$nod->internal_id()} = $iname;
	}
    }

    
    ## It will reproduce the tree in two steps.
    ## 1) Get all the nodes and create a new ones. Store the relations
    ##    based in the internal id

    my %node_rel = ();
    my %new_nodes = ();

    foreach my $nod_name (sort {$a cmp $b} keys %nodes) {
	my $node = $nodes{$nod_name};
	my $node_id = $node->id();
	my $int_node_id = $node->internal_id();
	$node_rel{$int_node_id} = {};
	
	my $strain;

	if (defined $node_id) {
	    $strain = $strains_href->{$node_id};
	    unless (defined $strain) {
		croak("DATA ERROR: $node_id isnt have defined strain");
	    }
	}

	## Get the node value for the topology

	my $branch_val;
	my $nbr_length = $node->branch_length();
	foreach my $cutoff (sort {$a <=> $b} keys %{$branch_cfs_href}) {
	    if (defined $nbr_length && !$branch_val) {
		if ($nbr_length <= $cutoff) {
		    $branch_val = $branch_cfs_href->{$cutoff};
		}
	    }
	}
	## If it is not defined by cutoff it will take the max + 1

	my %brcf = %{$branch_cfs_href};
	my @vals = sort {$brcf{$b} <=> $brcf{$a} } keys %brcf;
	my $max = $brcf{$vals[0]};

	if (defined $nbr_length && !$branch_val) {	    
	    $branch_val = $max+1;	
	}

	## Create the new node and put into a hash with internal id as key

	my $newnode = Bio::Tree::Node->new( -branch_length => $branch_val,
					    -id            => $strain,
	    );
	$new_nodes{$int_node_id} = $newnode;

	## Get the descendents for the old one

	foreach my $desc ($node->each_Descendent()) {
	    my $iid = $desc->internal_id();
	    $node_rel{$int_node_id}->{$iid} = $onodes{$iid};
	}
    }

    ## 2) After la creation of all the new nodes it will add the relations
    ##    to each of them. Ordered by onode hash values

    foreach my $node_iid (keys %new_nodes) {

	my %relat = %{$node_rel{$node_iid}};
	foreach my $child_iid (sort {$relat{$a} cmp $relat{$b}} keys %relat) {
	    $new_nodes{$node_iid}->add_Descendent($new_nodes{$child_iid});
	}
    }

    ## 3) Create the root
    
    my @new_rootdesc = ();
    my $old_root = $tree->get_root_node();
    foreach my $rootdesc ($old_root->each_Descendent()) {
	push @new_rootdesc, $new_nodes{$rootdesc->internal_id()};
    } 

    my $new_root = Bio::Tree::Node->new( -descendents => \@new_rootdesc );
    
    ## 4) Finally create the tree

    my $topotype = Bio::Tree::Tree->new( -root => $new_root );

    return $topotype;
}


=head2 _tree2newick

  Usage: my $newick = Bio::Tree::TopoType::_tree2newick($tree);

  Desc: Create a tree string in newick format
        
  Ret: $newick, a string

  Args: $tree (a Bio::Tree::Tree object)

  Side_Effects: Die if no argument is used

  Example: my $newick = Bio::Tree::TopoType::_tree2newick($tree);

=cut

sub _tree2newick {
    my $tree = shift ||
	croak("ARG. ERROR: No tree was not supplied to _tree2newick()");

    if (ref($tree) ne 'Bio::Tree::Tree') {
	croak("ARG. ERROR: Tree supplied to _tree2newick() isnt Tree Obj");
    }
     
    my @data = Bio::TreeIO::newick::_write_tree_Helper($tree->get_root_node);
    my $string = join('', @data);

    return $string;
}







####
1; #
####
