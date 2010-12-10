#!/usr/bin/perl

=head1 NAME

  topotype.t
  A piece of code to test the Bio::Tree::TopoType module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl topotype.t
 prove topotype.t

=head1 DESCRIPTION

 Test Bio::Tree::TopoType module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 38;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1 to 4

BEGIN {
    use_ok('PhyGeCluster');
    use_ok('Bio::Tree::TopoType');
    use_ok('Bio::Tree::Node');
    use_ok('Bio::Tree::Tree');
}


## Create an empty object and test the possible die functions. TEST 5 to 8

my $topoty0 = Bio::Tree::TopoType->new();

is(ref($topoty0), 'Bio::Tree::TopoType', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { Bio::Tree::TopoType->new(['fake']) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { Bio::Tree::TopoType->new({ fake => {} }) } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a valid arg.';

throws_ok { Bio::Tree::TopoType->new({ topology => 'fk2'}) } qr/ERROR: fk2/, 
    'TESTING DIE ERROR for new() when run_distances is used without run_boots.';


###############
## ACCESSORS ##
###############

## 1) Create a topology tree like:
##    ((Sp1:1,Sp2:1):1(Sp1:1,Sp3:1):1)

my $leaf1 = Bio::Tree::Node->new( -id => 'Sp1', -branch_length => 1 );
my $leaf2 = Bio::Tree::Node->new( -id => 'Sp2', -branch_length => 1 );
my $dn1 = [$leaf1, $leaf2];

my $leaf3 = Bio::Tree::Node->new( -id => 'Sp1', -branch_length => 1 );
my $leaf4 = Bio::Tree::Node->new( -id => 'Sp3', -branch_length => 1 );
my $dn2 = [$leaf3, $leaf4];

my $node1 = Bio::Tree::Node->new( -branch_length => 1, -descendents => $dn1);
my $node2 = Bio::Tree::Node->new( -branch_length => 1, -descendents => $dn2);
my $node3 = Bio::Tree::Node->new( -descendents   => [$node1, $node2]);

my $topotree0 = Bio::Tree::Tree->new( -id => 'topology1', -root => $node3);

## 2) Create the tree members (2)

my %first_tree = ( 0.8 => { Sp1_01 => 0.5, Sp1_02 => 0.7 }, 
		   0.2 => { Sp2_01 => 0.3, Sp3_01 => 0.8 } );

my %secon_tree = ( 0.3 => { Sp1_03 => 0.3, Sp1_04 => 0.2 }, 
		   0.4 => { Sp2_02 => 0.8, Sp3_03 => 0.4 } );

my %third_tree = ( 0.1 => { Sp1_06 => 0.6, Sp1_05 => 0.1 }, 
		   0.2 => { Sp2_03 => 0.3, Sp3_02 => 0.3 } );

my %forth_tree = ( 0.3 => { Sp1_09 => 0.4, Sp1_08 => 0.3 }, 
		   0.1 => { Sp2_04 => 0.3, Sp3_04 => 0.3 } );

my %tree_data = ('tree1' => \%first_tree, 
		 'tree2' => \%secon_tree, 
		 'tree3' => \%third_tree,
		 'tree4' => \%forth_tree,
    );

my @trees_members = ();
foreach my $tree_id (keys %tree_data) {
    my %data = %{$tree_data{$tree_id}};
    my @root_desc = ();
    foreach my $node_branch (keys %data) {
	
	my %leaves = %{$data{$node_branch}};
	my @node_desc = ();
	foreach my $leaf_id (keys %leaves) {
	    my $leaf = Bio::Tree::Node->new( -id            => $leaf_id,
					     -branch_length => $leaves{$leaf_id}
		);
	    push @node_desc, $leaf;
	}

	my $node =  Bio::Tree::Node->new( -branch_length => $node_branch, 
					  -descendents   => \@node_desc );
	push @root_desc, $node;
    }
    my $root =  Bio::Tree::Node->new( -descendents   => \@root_desc );
    my $tree = Bio::Tree::Tree->new( -id => $tree_id, -root => $root );
    push @trees_members, $tree;
}
my $tree_4th = pop(@trees_members);

## 3) Create the descrition hash

my %descrp = ( 'global' => 'This is a test', 'preferred' => 1 );

## Check set/get_topology, TEST 9 to 11

$topoty0->set_topology($topotree0);
my $topotree0_1 = $topoty0->get_topology();

is($topotree0_1, $topotree0, 
    "Testing get/set_topology, checking tree object identity")
    or diag("Looks like this has failed");

throws_ok { $topoty0->set_topology() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_topology function';

throws_ok { $topoty0->set_topology('fake') } qr/ARG. ERROR: fake at set_/, 
    'TESTING DIE ERROR when arg. supplied to set_topology function isnt Tree';

## Check set/get_members, TEST 12 to 16

$topoty0->set_members(\@trees_members);
my $members_aref = $topoty0->get_members();

my $member_ident_count = 0;
foreach my $memb (@{$members_aref}) {
    unless (ref($memb) eq 'Bio::Tree::Tree') {
	$member_ident_count++;
    }
}

is($member_ident_count, 0, 
    "Testing get/set_members, checking tree object identity for the members")
    or diag("Looks like this has failed");

is(scalar(@{$members_aref}), 3,
    "Testing get/set_members, checking number of members")
    or diag("Looks like this has failed");

throws_ok { $topoty0->set_members() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_members function';

throws_ok { $topoty0->set_members('fake') } qr/ARG. ERROR: fake set_/, 
    'TESTING DIE ERROR when arg. supplied to set_members() isnt ArrayRef';

throws_ok { $topoty0->set_members(['fake2']) } qr/ARG. ERROR: fake2 used/, 
    'TESTING DIE ERROR when arg. element supplied to set_members() isnt Tree';

## Checking add_members(), TEST 17 to 20

$topoty0->add_members([$tree_4th]);
my $members_aref2 = $topoty0->get_members();

is(scalar(@{$members_aref2}), 4,
    "Testing add_members, checking number of members")
    or diag("Looks like this has failed");

throws_ok { $topoty0->add_members() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to add_members function';

throws_ok { $topoty0->add_members('fake') } qr/ARG. ERROR: fake add_/, 
    'TESTING DIE ERROR when arg. supplied to add_members() isnt ArrayRef';

throws_ok { $topoty0->add_members(['fake2']) } qr/ARG. ERROR: fake2 used/, 
    'TESTING DIE ERROR when arg. element supplied to add_members() isnt Tree';

## Checking delete_members(), TEST 21 to 25

my @removed_members = ($trees_members[0], $trees_members[2]);
my @removed_members_ids = ($trees_members[0]->id(), $trees_members[2]->id());
my @del_memb = $topoty0->delete_members(\@removed_members);
my $members_aref3 = $topoty0->get_members();

my @del_memb_ids = ();
foreach my $del_memb (@del_memb) {
    push @del_memb_ids, $del_memb->id();   
}
my $expected_del = join(',', sort @removed_members_ids);
my $obtained_del = join(',', sort @del_memb_ids);

is(scalar(@{$members_aref3}), 2,
    "Testing delete_members, checking number of members")
    or diag("Looks like this has failed");

is($expected_del, $obtained_del,
    "Testing delete_members, checking identity of the deleted members")
    or diag("Looks like this has failed");

throws_ok { $topoty0->delete_members() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to delete_members function';

throws_ok { $topoty0->delete_members('fake') } qr/ARG. ERROR: fake delete_/, 
    'TESTING DIE ERROR when arg. supplied to delete_members() isnt ArrayRef';

throws_ok { $topoty0->delete_members(['fake2']) } qr/ARG. ERROR: fake2 used/, 
    'TESTING DIE ERROR when arg.element supplied to delete_members() isnt Tree';


## Check set/get_description, TEST 26 to 30

$topoty0->set_descriptions(\%descrp);
my %descriptions0 = %{$topoty0->get_descriptions()};

is(scalar(keys %descriptions0), 2, 
   "Testing get/set_descriptions, checking number of descriptions")
    or diag("Looks like this has failed");

is(join(',', sort(keys %descriptions0)), join(',', sort(keys %descrp)), 
   "Testing get/set_descriptions, checking tags of the descriptions")
    or diag("Looks like this has failed");

my $preferred_descript = $topoty0->get_descriptions('preferred');

is($preferred_descript, 1, 
    "Testing get_descriptions using tags, checking value for a concrete tag")
    or diag("Looks like this has failed");

throws_ok { $topoty0->set_descriptions() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_descriptions function';

throws_ok { $topoty0->set_descriptions('fake') } qr/ARG. ERROR: fake set_desc/, 
    'TESTING DIE ERROR when arg. supplied to set_descriptions() isnt HashRef';


## Check add_description, TEST 31 to 35

$topoty0->add_description('gene duplication', 1);
my %descriptions1 = %{$topoty0->get_descriptions()};

is(scalar(keys %descriptions1), 3, 
   "Testing add_description, checking number of descriptions")
    or diag("Looks like this has failed");

my $gdupl_descript1 = $topoty0->get_descriptions('gene duplication');

is($gdupl_descript1, 1, 
   "Testing add_description, checking value for a new tag added")
    or diag("Looks like this has failed");

$topoty0->add_description('gene duplication', 0);

my $gdupl_descript2 = $topoty0->get_descriptions('gene duplication');

is($gdupl_descript2, 0, 
   "Testing add_description (editing), checking value for a new tag added")
    or diag("Looks like this has failed");

throws_ok { $topoty0->add_description() } qr/ARG. ERROR: No tag arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to add_description function';

throws_ok { $topoty0->add_description('fake') } qr/ARG. ERROR: No desc/, 
    'TESTING DIE ERROR when description arg. supplied to add_description()';


## Check delete_description, TEST 36 to 38

my $deleted_tag = $topoty0->delete_description('gene duplication');

my %descriptions2 = %{$topoty0->get_descriptions()};

is(scalar(keys %descriptions2), 2, 
   "Testing delete_description, checking number of descriptions")
    or diag("Looks like this has failed");

is($deleted_tag, 0, 
   "Testing delete_description, checking value of the tag deleted")
    or diag("Looks like this has failed");

throws_ok { $topoty0->delete_description() } qr/ARG. ERROR: No tag arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to delete_description()';



####
1; #
####
