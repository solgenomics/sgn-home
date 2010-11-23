#!/usr/bin/perl

=head1 NAME

  phygeboots.t
  A piece of code to test the PhyGeBoots module used for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl phygeboots.t
 prove phygeboots.t

 Note: This test use three external files that should be in the same 
       test folder.
       + testfiles/selfblast.test.m8
       + testfiles/seq.test.fasta
       + testfiles/strains.test.tab
       + testfiles/assembly_out.test.ace

=head1 DESCRIPTION

 Test PhyGeBoots.pm module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 66;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1 to 7

BEGIN {
    use_ok('PhyGeBoots');
    use_ok('PhyGeCluster');
    use_ok('Bio::Cluster::SequenceFamily');
    use_ok('Bio::Seq');
    use_ok('Bio::SimpleAlign');
    use_ok('Bio::Matrix::PhylipDist');
    use_ok('Bio::Tree::Tree');
}


## Create an empty object and test the possible die functions. TEST 8 to 13

my $phygeb0 = PhyGeBoots->new();

is(ref($phygeb0), 'PhyGeBoots', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { PhyGeBoots->new(['fake']) } qr/ARGUMENT ERROR:/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { PhyGeBoots->new({ run_bootstrap => {}}) } qr/on: run_bootstrap/, 
    'TESTING DIE ERROR for new() when run_bootstrap is used without seqfam';

my $fkhref0_1 = { seqfam => 'test', run_distances => {} };
throws_ok { PhyGeBoots->new($fkhref0_1) } qr/on: run_distances/, 
    'TESTING DIE ERROR for new() when run_distances is used without run_boots.';

my $fkhref0_2 = { seqfam => 'test', run_bootstrap => {}, run_njtrees => {} };
throws_ok { PhyGeBoots->new($fkhref0_2) } qr/on: run_njtrees/, 
    'TESTING DIE ERROR for new() when run_njtrees is used without run_distanc.';

my $fkhref0_3 = { seqfam        => 'test', 
		  run_bootstrap => {}, 
		  run_distances => {},
		  run_consensus => {},
};
throws_ok { PhyGeBoots->new($fkhref0_3) } qr/on: run_consensus/, 
    'TESTING DIE ERROR for new() when run_consensus is used without run_trees';


## To test get/set functions it will need:
## 1) Bio::Cluster::SequenceFamily object, TEST 14 to 17

my %seqtest = ( 
    'seqid1' => 'ACTCGCTGCGTA', 
    'seqid2' => 'ACTCCCTTCGTT',
    'seqid3' => 'ACACCCTGGGTT',
    );

my @seqs = ();
foreach my $testid (keys %seqtest) {
    my $seq = Bio::Seq::Meta->new( 
	-id    => $testid, 
	-seq   => $seqtest{$testid},
	-start => 1,
	-end   => 12,
	);
    push @seqs, $seq;
}

my $seqfam0 = Bio::Cluster::SequenceFamily->new( 
    -family_id => 'seqfam0',
    -members   => \@seqs,
    );

$phygeb0->set_seqfam($seqfam0);
my $seqfam0_1 = $phygeb0->get_seqfam();

is(ref($seqfam0_1), 'Bio::Cluster::SequenceFamily',
    "Testing set/get_seqfam, checking object identity")
    or diag("Looks like this has failed");

is(scalar($seqfam0_1->get_members()), 3,
    "Testing set/get_seqfam, checking number of members for interla object")
    or diag("Looks like this has failed");

throws_ok { $phygeb0->set_seqfam() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR for set_seqfam() when no arg. is supplied';

throws_ok { $phygeb0->set_seqfam('fake') } qr/ARG. ERROR: Arg=fake/, 
    'TESTING DIE ERROR for set_seqfam() when arg supplied is not the right obj';

## 2) Bio::SimpleAlign object, TEST 18 to 22

my $align0 = Bio::SimpleAlign->new( -seqs => \@seqs );

$phygeb0->set_aligns([$align0]);
my @align0_1 = $phygeb0->get_aligns();

is(scalar(@align0_1), 1,
    "Testing set/get_aligns, checking number of alignments")
    or diag("Looks like this has failed");

is(ref($align0_1[0]), 'Bio::SimpleAlign',
    "Testing set/get_aligns, checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phygeb0->set_aligns() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR for set_aligns() when no arg. is supplied';

throws_ok { $phygeb0->set_aligns('fake') } qr/ARG. ERROR: Arg=fake/, 
    'TESTING DIE ERROR for set_aligns() when arg. supplied isnt arrayref';

throws_ok { $phygeb0->set_aligns(['fake']) } qr/ARG. ERROR: Element=fake/, 
    'TESTING DIE ERROR for set_aligns() when arg. supplied doesnt have the obj';


## 3) Bio::Matrix::PhylipDist object, TEST 23 to 27

my $dist_mtx0 = Bio::Matrix::Generic->new();

$phygeb0->set_dists([$dist_mtx0]);
my @dist_mtx0_1 = $phygeb0->get_dists();

is(scalar(@dist_mtx0_1), 1,
    "Testing set/get_dists, checking number of matrix")
    or diag("Looks like this has failed");

is(ref($dist_mtx0_1[0]), 'Bio::Matrix::Generic',
    "Testing set/get_dists, checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phygeb0->set_dists() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR for set_dists() when no arg. is supplied';

throws_ok { $phygeb0->set_dists('fake') } qr/ARG. ERROR: Arg=fake/, 
    'TESTING DIE ERROR for set_dists() when arg. supplied isnt arrayref';

throws_ok { $phygeb0->set_dists(['fake']) } qr/ARG. ERROR: Element=fake/, 
    'TESTING DIE ERROR for set_dists() when arg. supplied doesnt have the obj';


## 4) Bio::Tree::Tree object, TEST 28 to 32

my $tree0 = Bio::Tree::Tree->new();

$phygeb0->set_trees([$tree0]);
my @trees0 = $phygeb0->get_trees();

is(scalar(@trees0), 1,
    "Testing set/get_trees, checking number of tree objects")
    or diag("Looks like this has failed");

is(ref($trees0[0]), 'Bio::Tree::Tree',
    "Testing set/get_trees, checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phygeb0->set_trees() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR for set_trees() when no arg. is supplied';

throws_ok { $phygeb0->set_trees('fake') } qr/ARG. ERROR: Arg=fake/, 
    'TESTING DIE ERROR for set_trees() when arg. supplied isnt arrayref';

throws_ok { $phygeb0->set_trees(['fake']) } qr/ARG. ERROR: Element=fake/, 
    'TESTING DIE ERROR for set_trees() when arg. supplied doesnt have the obj';


## 5) For consensus it also can use a Bio::Tree::Tree object, TEST 33 to 35

my $cons0 = Bio::Tree::Tree->new();

$phygeb0->set_consensus($cons0);
my $cons0_1 = $phygeb0->get_consensus();

is(ref($cons0_1), 'Bio::Tree::Tree',
    "Testing set/get_consensus, checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phygeb0->set_consensus() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR for set_consensus() when no arg. is supplied';

throws_ok { $phygeb0->set_consensus('fake') } qr/ARG. ERROR: Arg=fake/, 
    'TESTING DIE ERROR for set_consensus() when arg supplied isnt a right obj';

#####################
## RUNNING OPTIONS ##
#####################

## First it will get some real Bio::Cluster::SequenceFamily object. To do that
## it will create a new PhyGeCluster object and it will get one cluster of
## 6 members (faster).

my $blastfile = "$FindBin::Bin/testfiles/selfblast.test.m8";
my $seqfile = "$FindBin::Bin/testfiles/seq.test.fasta";
my $strainfile = "$FindBin::Bin/testfiles/strains.test.tab";
my $acefile = "$FindBin::Bin/testfiles/assembly_out.test.ace";
my $blastdbfile = "$FindBin::Bin/testfiles/blastref.test.fasta";

my $phygecluster = PhyGeCluster->new(
    { 
	acefile    => $acefile, 
	strainfile => $strainfile,
    }
    );

my %clusters = %{$phygecluster->get_clusters()};
my %cl_sizes = $phygecluster->cluster_sizes(3,6);
my @ord_clusters = sort { $cl_sizes{$b} <=> $cl_sizes{$a} } keys %cl_sizes;

my $seqfam1 = $clusters{$ord_clusters[0]};
my $member_n1 = scalar($seqfam1->get_members());

## 1) Create a new object with this seqfam and run_bootstrap, TEST 36 to 38

my $phygeb1 = PhyGeBoots->new({ seqfam => $seqfam1 });
my @aligns1 = $phygeb1->run_bootstrap(
    { 
	datatype   => 'Sequence', 
	replicates => 500, 
	quiet      => 1,
    }
    );

is(scalar(@aligns1), 500, 
    "Testing run_bootstrap, checking number of replicates (500)")
    or diag("Looks like this has failed");

my $wrong_align_objs1 = 0;
my $wrong_memb_count1 = 0;
foreach my $align1 (@aligns1) {
    unless (ref($align1) eq 'Bio::SimpleAlign') {
	$wrong_align_objs1++;
	$wrong_memb_count1++;
    }
    else {
	unless ($align1->num_sequences() == $member_n1) {
	    $wrong_memb_count1++;
	}
    }
}

is($wrong_align_objs1, 0, 
    "testing run_bootstrap, checking object identity")
    or diag("Looks like this has failed");

is($wrong_memb_count1, 0, 
    "testing run_bootstrap, checking number of members per alignment")
    or diag("Looks like this has failed");

## Test the right croak for run_bootstrap function, TEST 39 to 42

throws_ok { $phygeb1->run_bootstrap() } qr/ARG. ERROR: No args./, 
    'TESTING DIE ERROR for run_bootstrap() when no arg. is supplied';

throws_ok { $phygeb1->run_bootstrap(['fk']) } qr/ARG. ERROR: Arg. supplied/, 
    'TESTING DIE ERROR for run_bootstrap() arg. supplied isnt a HASHREF';

throws_ok { $phygeb1->run_bootstrap({ fk => 1}) } qr/ARG. ERROR: fk/, 
    'TESTING DIE ERROR for run_bootstrap() arg. supplied isnt permited';

throws_ok { $phygeb1->run_bootstrap({ quiet => 2}) } qr/ARG. ERROR: 2/, 
    'TESTING DIE ERROR for run_bootstrap() when arg. supplied has wrong value';

## 2) test run_distances, TEST 43 to 45

my @dists1 = $phygeb1->run_distances();

is(scalar(@dists1), 500, 
    "Testing run_distances, checking number of replicates (500)")
    or diag("Looks like this has failed");

my $wrong_dists_objs1 = 0;
my $wrong_memb_count2 = 0;
foreach my $dist1 (@dists1) {
    unless (ref($dist1) eq 'Bio::Matrix::PhylipDist') {
	$wrong_dists_objs1++;
	$wrong_memb_count2++;
    }
    else {	
	unless ($dist1->num_rows() == $member_n1) {
	    $wrong_memb_count2++;
	}
    }
}

is($wrong_dists_objs1, 0, 
    "testing run_distances, checking object identity")
    or diag("Looks like this has failed");

is($wrong_memb_count2, 0, 
    "testing run_distances, checking number of rows per matrix")
    or diag("Looks like this has failed");

## And check the croaks, TEST 46 to 48

throws_ok { $phygeb1->run_distances(['fk']) } qr/ARG. ERROR: Arg. supplied/, 
    'TESTING DIE ERROR for run_distances() arg. supplied isnt a HASHREF';

throws_ok { $phygeb1->run_distances({ fk => 1}) } qr/ARG. ERROR: No method/, 
    'TESTING DIE ERROR for run_distances() when no method arg is supplied';

throws_ok { $phygeb1->run_distances({method =>'fk'}) } qr/ARG. ERROR: fk/, 
    'TESTING DIE ERROR for run_distances() when method supplied isnt permited';


## 3) test run_njtrees, TEST 49 to 51

my @njtrees1 = $phygeb1->run_njtrees({ type => 'NJ', quiet => 1 });

is(scalar(@njtrees1), 500, 
    "Testing run_njtrees, checking number of replicates (500)")
    or diag("Looks like this has failed");

my $wrong_trees_objs1 = 0;
my $wrong_memb_count3 = 0;
foreach my $njtree1 (@njtrees1) {
    unless (ref($njtree1) eq 'Bio::Tree::Tree') {
	$wrong_trees_objs1++;
	$wrong_memb_count3++;
    }
    else {
	unless (scalar($njtree1->get_leaf_nodes()) == $member_n1) {
	    $wrong_memb_count3++;
	}
    }
}

is($wrong_trees_objs1, 0, 
    "testing run_njtrees, checking object identity")
    or diag("Looks like this has failed");

is($wrong_memb_count3, 0, 
    "testing run_njtrees, checking number of rows per matrix")
    or diag("Looks like this has failed");

## And check the croaks, TEST 52 to 54

throws_ok { $phygeb1->run_njtrees(['fk']) } qr/ARG. ERROR: Arg. supplied/, 
    'TESTING DIE ERROR for run_njtrees() arg. supplied isnt a HASHREF';

throws_ok { $phygeb1->run_njtrees({ fk => 1}) } qr/ARG. ERROR: fk/, 
    'TESTING DIE ERROR for run_njtrees() when arg supplied is not permited';

throws_ok { $phygeb1->run_njtrees({quiet =>'fk'}) } qr/ARG. ERROR: quiet/, 
    'TESTING DIE ERROR for run_njtrees() when value used is not permited';


## 4) test run_mltrees, TEST 55 to 57

my @mltrees1 = $phygeb1->run_mltrees();

is(scalar(@mltrees1), 500, 
    "Testing run_mltrees, checking number of replicates (500)")
    or diag("Looks like this has failed");

my $wrong_trees_objs2 = 0;
my $wrong_memb_count4 = 0;
foreach my $mltree1 (@mltrees1) {
    unless (ref($mltree1) eq 'Bio::Tree::Tree') {
	$wrong_trees_objs2++;
	$wrong_memb_count4++;
    }
    else {
	unless (scalar($mltree1->get_leaf_nodes()) == $member_n1) {
	    $wrong_memb_count4++;
	}
    }
}

is($wrong_trees_objs2, 0, 
    "testing run_mltrees, checking object identity")
    or diag("Looks like this has failed");

is($wrong_memb_count4, 0, 
    "testing run_mltrees, checking number of rows per matrix")
    or diag("Looks like this has failed");

## And check the croaks, TEST 58 to 60

throws_ok { $phygeb1->run_mltrees(['fk']) } qr/ARG. ERROR: Arg. supplied/, 
    'TESTING DIE ERROR for run_mltrees() arg. supplied isnt a HASHREF';

throws_ok { $phygeb1->run_mltrees({ fk => 1}) } qr/ARG. ERROR: fk/, 
    'TESTING DIE ERROR for run_mltrees() when arg supplied is not permited';

throws_ok { $phygeb1->run_mltrees({phyml_arg =>'fk'}) } qr/ARG. ERROR: phy/, 
    'TESTING DIE ERROR for run_mltrees() when value used is not permited';


## 5) Testing run_consensus, TEST 61 to 66

my $consensus1 = $phygeb1->run_consensus({ quiet => 1});

is(ref($consensus1), 'Bio::Tree::Tree', 
    "Testing run_consensus, checking object identity")
    or diag("Looks like this has failed");

is(scalar($consensus1->get_leaf_nodes()), $member_n1, 
    "Testing run_consensus, checking number of leaf nodes")
    or diag("Looks like this has failed");

my $consensus2 = $phygeb1->run_consensus({ quiet => 1, normalized => 1});

my $wrong_norm = 0;
my @nodes = $consensus2->get_nodes();
foreach my $node (@nodes) {
    if (defined $node->branch_length()) {
	if ($node->branch_length() > 100) {
	    $wrong_norm++;
	}
    }
}

is($wrong_norm, 0, 
    "Testing run_consensus (normalized), cehccking that branches < 100")
    or diag("Looks like this has failed");

throws_ok { $phygeb1->run_consensus(['fk']) } qr/ARG. ERROR: Arg. supplied/, 
    'TESTING DIE ERROR for run_consensus() arg. supplied isnt a HASHREF';

throws_ok { $phygeb1->run_consensus({ fk => 1}) } qr/ARG. ERROR: fk/, 
    'TESTING DIE ERROR for run_consensus() when arg supplied is not permited';

throws_ok { $phygeb1->run_consensus({ quiet =>'fk'}) } qr/ARG. ERROR: quiet/, 
    'TESTING DIE ERROR for run_consensus() when value used is not permited';


####
1; #
####
