#!/usr/bin/perl

=head1 NAME

  phygecluster.t
  A piece of code to test the Phygecluster module used for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl phygecluster.t
 prove phygecluster.t

 Note: This test use three external files that should be in the same 
       test folder.
       + testfiles/selfblast.test.m8
       + testfiles/seq.test.fasta
       + testfiles/strains.test.tab
       + testfiles/assembly_out.test.ace

=head1 DESCRIPTION

 Test Phygecluster.pm module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 175;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1,2,3

BEGIN {
    use_ok('PhyGeCluster');
    use_ok('Bio::Seq');
    use_ok('Bio::SeqIO');
}


## TEST 4, create an empty object

my $phygecluster0 = PhyGeCluster->new();

is(ref($phygecluster0), 'PhyGeCluster', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

## TEST 5 to 11 test add/get/remove
## Create an empty Bio::Seq object

my @seq_objs1 = ();
my @seq_ids1 = ('member_test1', 'member_test2', 'member_test3');
foreach my $seq_id1 (@seq_ids1) {
    my $seq = Bio::Seq->new( -id => $seq_id1);
    push @seq_objs1, $seq;
}

$phygecluster0->add_cluster('cluster_test1', \@seq_objs1);
my %clusters1 = %{$phygecluster0->get_clusters()};

is(scalar(keys %clusters1), 1, 
    "Test add/get_cluster function; Checking number of clusters")
    or diag("Looks like this has failed");

is(join(',', keys %clusters1), 'cluster_test1', 
    "Test add/get_cluster function; Checking cluster id")
    or diag("Looks like this has failed");
    
is(ref($clusters1{'cluster_test1'}), 'Bio::Cluster::SequenceFamily', 
    "Test add/get_cluster function; Cheking object ref. for hash value")
    or diag("Looks like this has failed");

my @members_ids1 = ();
my @members1 = $clusters1{'cluster_test1'}->get_members();
foreach my $member_obj (@members1) {
    push @members_ids1, $member_obj->display_id();
}

is(join(',', sort @members_ids1), join(',', sort @seq_ids1), 
    "Test add_cluster function; Cheking members for a cluster 1")
    or diag("Looks like this has failed");

my @seq_objs2 = ();
my @seq_ids2 = ('member_test4', 'member_test5');
foreach my $seq_id2 (@seq_ids2) {
    my $seq = Bio::Seq->new( -id => $seq_id2);
    push @seq_objs2, $seq;
}

$phygecluster0->add_cluster('cluster_test2', \@seq_objs2);
my %clusters2 = %{$phygecluster0->get_clusters()};

is(scalar(keys %clusters2), 2, 
    "Test add_cluster function; Checking number of clusters")
    or diag("Looks like this has failed");

is(join(',', sort(keys %clusters2)), 'cluster_test1,cluster_test2', 
    "Test add_cluster function; Checking cluster id list")
    or diag("Looks like this has failed");

my @members_ids2 = ();
my @members2 = $clusters2{'cluster_test2'}->get_members();
foreach my $member_obj (@members2) {
    push @members_ids2, $member_obj->display_id();
}

is(join(',', sort @members_ids2), join(',', sort @seq_ids2), 
    "Test add_cluster function; Cheking members for a cluster 2")
    or diag("Looks like this has failed");

## Testing die option for the wrong use of these functions, TEST 12 to 24

throws_ok { PhyGeCluster->new(['fake']) } qr/ARGUMENT ERROR:/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

my $fakehref1 = { sequencefile => 'test'};
throws_ok { PhyGeCluster->new($fakehref1) } qr/ARGUMENT INCOMPATIBILITY/, 
    'TESTING DIE ERROR when sequencefile arg. is used without blastfile';

my $fakehref2 = { fastblast_parser => 1};
throws_ok { PhyGeCluster->new($fakehref2) } qr/ARGUMENT INCOMPATIBILITY/, 
    'TESTING DIE ERROR when fastblast_parser arg. is used without blastfile';

my $fakehref3 = { strainfile => 'test'};
throws_ok { PhyGeCluster->new($fakehref3) } qr/ARGUMENT INCOMPATIBILITY/, 
    'TESTING DIE ERROR when strainfile arg. is used without blastfile';

my $fakehref4 = { blastfile => 'test', run_alignments => 'test'};
throws_ok { PhyGeCluster->new($fakehref4) } qr/ARGUMENT INCOMPATIBILITY/, 
    'TESTING DIE ERROR when run_alignments arg. is used without sequencefile';

my $fakehref5 = { blastfile     => 'test', 
		  sequencefile  => 'test', 
		  run_distances => 'test' };
throws_ok { PhyGeCluster->new($fakehref5) } qr/ARGUMENT INCOMPATIBILITY/, 
    'TESTING DIE ERROR when run_distances arg. is used without run_alignments';


throws_ok { $phygecluster0->set_clusters() } qr/ARGUMENT ERROR: No/, 
    'TESTING DIE ERROR when none arg. has been supplied set_cluster function';

throws_ok { $phygecluster0->set_clusters("WrongRef.") } qr/ARGUMENT ERROR: Wr/, 
    'TESTING DIE ERROR when arg. supplied to set_cluster function isnt hashref';

throws_ok { $phygecluster0->set_clusters({t=>'fake'})} qr/ARGUMENT ERROR: val/, 
    'TESTING DIE ERROR when hash supplied to set_cluster has not Bio::Cluster';


throws_ok { $phygecluster0->add_cluster() } qr/ARG. ERROR: No /, 
    'TESTING DIE ERROR when none argument was supplied to add_cluster function';

throws_ok { $phygecluster0->add_cluster("test1") } qr/ARG. ERROR: None/, 
    'TESTING DIE ERROR when only one argument was supplied to add_cluster';

throws_ok { $phygecluster0->add_cluster("test1", 1) } qr/ARG. ERROR: member a/, 
    'TESTING DIE ERROR when 2nd arg supplied to add_cluster is not a array ref';

throws_ok { $phygecluster0->add_cluster("t", [1]) } qr/ARG. ERROR: member=1/, 
    'TESTING DIE ERROR when elements in the 2nd arg. are not Bio::Seq objs';



## TEST 25 and 26 (find_cluster)

my $cluster3 = $phygecluster0->find_cluster('member_test2');

is($cluster3->family_id(), 'cluster_test1',
   "Test find_cluster function with defined value; checking cluster name")
    or diag("Looks like this has failed");

is($phygecluster0->find_cluster('fake_cluster_test1'), undef,
   "Test find_cluster function with undefined value; checking cluster name")
    or diag("Looks like this has failed");

## TEST 27 to 29 (remove_cluster)

$phygecluster0->remove_cluster('cluster_test1');
my %clusters3 = %{$phygecluster0->get_clusters()};

is(scalar(keys %clusters3), 1, 
    "Test remove_cluster function; Checking number of clusters")
    or diag("Looks like this has failed");

is(join(',', sort(keys %clusters3)), 'cluster_test2', 
    "Test remove_cluster function; Checking cluster id list")
    or diag("Looks like this has failed");

throws_ok { $phygecluster0->remove_cluster() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to remove_cluster function';

## TEST 30 to 33 (get/set strains)

my %strains1 = (
    'member_test4' => 'strain1', 
    'member_test5' => 'strain2' 
    ); 

$phygecluster0->set_strains(\%strains1);
my %get_strains1 = %{$phygecluster0->get_strains()};

foreach my $key_strain (keys %get_strains1) {
    is($get_strains1{$key_strain}, $strains1{$key_strain},
       "Test get/set_strains; Checking strain value for $key_strain")
	or diag("Looks like this has failed");
}


throws_ok { $phygecluster0->set_strains() } qr/ARGUMENT ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to set_strains function';

throws_ok { $phygecluster0->set_strains("Fake") } qr/ARGUMENT ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied to set_strains is not hashref';


##################
## PARSING TEST ##
##################

my $blastfile = "$FindBin::Bin/testfiles/selfblast.test.m8";
my $seqfile = "$FindBin::Bin/testfiles/seq.test.fasta";
my $strainfile = "$FindBin::Bin/testfiles/strains.test.tab";
my $acefile = "$FindBin::Bin/testfiles/assembly_out.test.ace";
my $blastdbfile = "$FindBin::Bin/testfiles/blastref.test.fasta";

## First get the data from the files

my @seqs = ();
my %input_seqs = ();
my $seq_io = Bio::SeqIO->new( -file => $seqfile, -format => 'fasta');

while( my $seq = $seq_io->next_seq()) {
    push @seqs, $seq->display_id();
    $input_seqs{$seq->display_id()} = $seq;
}
my $seq_n = scalar(@seqs);


## Test the function TEST 34 to 39

my %clusters4 = PhyGeCluster::parse_blastfile({ blastfile => $blastfile});

## So it will create them and compare with the obtained number

my $cl4_number = scalar(keys %clusters4);;
my $cl4_count = 1;
my @expected_cluster_names = ();
while ($cl4_count < $cl4_number+1) {
    push @expected_cluster_names, 'cluster_' . $cl4_count;
    $cl4_count++;
}

is(join(',', sort(keys %clusters4)), join(',', sort(@expected_cluster_names)), 
   "Test parse_blastfile; Checking cluster_id list ($cl4_number clusters)")
    or diag("Looks like this has failed");

my @obtained_members = ();
foreach my $cluster_id (keys %clusters4) {
    my @members4 =  $clusters4{$cluster_id}->get_members();
    foreach my $member_obj4 (@members4) {
	push @obtained_members, $member_obj4->display_id();
    }
}

is(join(',', sort(@obtained_members)), join(',', sort(@seqs)), 
   "Test parse_blastfile; Checking member_id list ($seq_n members)")
    or diag("Looks like this has failed");

## Compare with fastparse_blastfile

my %clusters4fast = PhyGeCluster::fastparse_blastfile(
    { 
	blastfile => $blastfile
    }
    );

my $diff1 = scalar(keys %clusters4fast);

foreach my $cl4f (sort keys %clusters4fast) {
    if (exists $clusters4{$cl4f}) {
	my @memb_ids4r = ();
	my @memb_obj4r = $clusters4{$cl4f}->get_members();
	foreach my $member4r (@memb_obj4r) {
	    push @memb_ids4r, $member4r->display_id();
	}
	my $reg_members4 = join(',', sort @memb_ids4r);

	my @memb_ids4f = ();
	my @memb_obj4f = $clusters4fast{$cl4f}->get_members();
	foreach my $member4f (@memb_obj4f) {
	    push @memb_ids4f, $member4f->display_id();
	}
	my $fast_members4 = join(',', sort @memb_ids4f);

	if($reg_members4 eq $fast_members4) {
	    $diff1--;
	}
	else {
	    print STDERR "Diff cluster: $cl4f\tR=$reg_members4 ";
	    print STDERR "vs F=$fast_members4\n";
	}
    }
    else {
	print STDERR "Cluster $cl4f do not exists in regular method\n";
    }
}
is($diff1, 0,
   "Testing fastparse_blastfile; Checking parse_blastfile comparison")
    or diag("Looks like this has failed");


## Now it will check other cluster conditions

my %clusters5 = PhyGeCluster::parse_blastfile(
    { 
	blastfile     => $blastfile,
	clustervalues => { percent_identity => ['>', 96], 
			   hsp_length       => ['>', 30] },
    }
    );

## So it will create them and compare with the obtained number

my $cl5_n = scalar(keys %clusters5);
my $cl5_count = 1;
my @expected_cluster_names5 = ();
while ($cl5_count < $cl5_n+1) {
    push @expected_cluster_names5, 'cluster_' . $cl5_count;
    $cl5_count++;
}

is(join(',', sort(keys %clusters5)), join(',', sort(@expected_cluster_names5)), 
   "Test parse_blastfile (id=95,l=30); Checking cl_id list($cl5_n clusters)")
    or diag("Looks like this has failed");

my @obtained_members5 = ();
foreach my $cluster_id (keys %clusters5) {
    my @memb_objs5 = $clusters5{$cluster_id}->get_members();
    foreach my $memb_obj5 (@memb_objs5) {
	push @obtained_members5, $memb_obj5->display_id();
    }
}

is(join(',', sort(@obtained_members5)), join(',', sort(@seqs)), 
   "Test parse_blastfile (id=95,l=30); Checking member_id list($seq_n members)")
    or diag("Looks like this has failed");

## Checking fastparse_blastfile other conditions
## Note: parse_blastfile use Bio::SeqIO parser. This parser calculate the
## percent_identity as frac_identical * 100 insetad to use the percent_identity 
## from the blast result file (as fastparse_blastfile does), so the cluster
## can be different

my %clusters5fast = PhyGeCluster::fastparse_blastfile(
    { 
	blastfile     => $blastfile,
	clustervalues => { percent_identity => ['>', 96], 
			   align_length     => ['>', 30] },
    }
    );

my $diff2 = scalar(keys %clusters5fast);

foreach my $cl5f (sort keys %clusters5fast) {
    if (exists $clusters5{$cl5f}) {

	my @memb_ids5r = ();
	my @memb_obj5r = $clusters5{$cl5f}->get_members();
	foreach my $member5r (@memb_obj5r) {
	    push @memb_ids5r, $member5r->display_id();
	}
	my $reg_members5 = join(',', sort @memb_ids5r);

	my @memb_ids5f = ();
	my @memb_obj5f = $clusters5fast{$cl5f}->get_members();
	foreach my $member5f (@memb_obj5f) {
	    push @memb_ids5f, $member5f->display_id();
	}
	my $fast_members5 = join(',', sort @memb_ids5f);

	if($reg_members5 eq $fast_members5) {
	    $diff2--;
	}
	else {
	    print STDERR "Diff cluster: $cl5f\tR=$reg_members5 ";
	    print STDERR "vs F=$fast_members5\n";
	}
    }
    else {
	print STDERR "Cluster $cl5f do not exists in regular method\n";
    }
}
is($diff2, 0,
   "Testing fastparse_blastfile (id=95,l=30); Checking parse_blastfile compar.")
    or diag("Looks like this has failed");


## Testing die functions for parse_blastfile, TEST 40 to 48

throws_ok { PhyGeCluster::parse_blastfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to parse_blastfile function';

throws_ok { PhyGeCluster::parse_blastfile('Fake') } qr/ARGUMENT ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { PhyGeCluster::parse_blastfile({t => 1}) } qr/ARGUMENT ERROR: 'bl/, 
    'TESTING DIE ERROR when href arg. does not contain blastfile key';

my $fake_href1 = { blastfile     => $blastfile,
		   clustervalues => { percent_identity => ['>', 95], 
				      align_length     => ['>', 30] },
};

throws_ok { PhyGeCluster::parse_blastfile($fake_href1) } qr/WRONG BLAST V/, 
    'TESTING DIE ERROR when href arg., clustervalues use wrong keys';


throws_ok { PhyGeCluster::fastparse_blastfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied fastparse_blastfile function';

throws_ok { PhyGeCluster::fastparse_blastfile('Fake') } qr/ARGUMENT ERROR: F/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { PhyGeCluster::fastparse_blastfile({t => 1}) } qr/ARGUMENT ERROR: /, 
    'TESTING DIE ERROR when href arg. does not contain blastfile key';

my $fake_href2 = { blastfile     => $blastfile,
		   clustervalues => { percent_identity => ['>', 95], 
				      hps_length       => ['>', 30] },
};

throws_ok { PhyGeCluster::fastparse_blastfile($fake_href2) } qr/WRONG BLAST V/, 
    'TESTING DIE ERROR when href arg., clustervalues use wrong keys';

my $fake_href3 = { blastfile => $blastfile, blastformat => 'megablast'};
throws_ok { PhyGeCluster::fastparse_blastfile($fake_href3) } qr/ARG. ERROR: f/, 
    'TESTING DIE ERROR when href arg., blastformat is not m8 or blasttable';


## Checking parsed sequence file, TEST 49

my %seqs1 = PhyGeCluster::parse_seqfile({ sequencefile => $seqfile });

is(join(',', sort(keys %seqs1)), join(',', sort(@seqs)),
   "Test parse_seqfile; Cheking sequence_ids")
    or diag("Looks like this failed");

## Testing die functions for parse_seqfile, TEST 50 to 52

throws_ok { PhyGeCluster::parse_seqfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to parse_seqfile function';

throws_ok { PhyGeCluster::parse_seqfile('Fake') } qr/ARG. ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { PhyGeCluster::parse_seqfile({t => 1}) } qr/ARG. ERROR: 'se/, 
    'TESTING DIE ERROR when href arg. does not contain seqfile key';

## Checking parsed strain file, TEST 53 and 54

my %expected_strains = ();
my $strain_seq_n = 0;
open my $str_io, '<', $strainfile;
while (<$str_io>) {
    chomp($_);
    my @data = split('\t', $_);
    unless ($expected_strains{$data[1]}) {
	$expected_strains{$data[1]} = 1;
    }
    else {
	$expected_strains{$data[1]}++;
    }
    $strain_seq_n++;
}
close $str_io;

my %strain = PhyGeCluster::parse_strainfile({ strainfile => $strainfile });
my %rev_strain = ();
foreach my $seqid (keys %strain) {
    unless (exists $rev_strain{$strain{$seqid}}) {
	$rev_strain{$strain{$seqid}} = 1
    }
    else {
	$rev_strain{$strain{$seqid}}++;
    }
}

is(scalar(keys(%strain)), $strain_seq_n, 
   "Test parse_strainfile; Checking number of sequences")
   or diag("Looks like this failed");

is(scalar(keys(%rev_strain)), scalar(keys(%expected_strains)),
   "Test parse_strainfile; Checking number of strains")
   or diag("Looks like this failed");

## Testing die for the parse_strainfile function, TEST 55 to 57

throws_ok { PhyGeCluster::parse_strainfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to parse_strainfile function';

throws_ok { PhyGeCluster::parse_strainfile('Fake') } qr/ARG. ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { PhyGeCluster::parse_strainfile({t => 1}) } qr/ARG. ERROR: 'str/, 
    'TESTING DIE ERROR when href arg. does not contain strainfile key';


## Test parsing_acefile()

my %clusters6 = PhyGeCluster::parse_acefile({ acefile => $acefile});

## So it will create them and compare with the obtained number

my $cl6_number = scalar(keys %clusters6);;
my $cl6_count = 1;
my @expected_cluster_names6 = ();
while ($cl6_count < $cl6_number+1) {
    push @expected_cluster_names6, 'Nic_selected_c' . $cl6_count;
    $cl6_count++;
}

is(join(',', sort(keys %clusters6)), join(',', sort(@expected_cluster_names6)), 
   "Test parse_acefile; Checking cluster_id list ($cl6_number clusters)")
    or diag("Looks like this has failed");

my @obtained_members6 = ();
foreach my $cluster_id (keys %clusters6) {
    my @members6 = $clusters6{$cluster_id}->get_members();
    foreach my $seqobj6 (@members6) {
	push @obtained_members6, $seqobj6->display_id();
    }
}

## Get the expected members parsing the file, TEST 58 to 62

my @exp_seq6 = ();
open my $acefh, '<', $acefile;
while(<$acefh>) {
    chomp($_);
    if ($_ =~ m/^AF\s+(.+?)\s+(\w)\s+(-?\d+)/) {
	push @exp_seq6, $1;
    }
}
close $acefh;

is(join(',', sort(@obtained_members6)), join(',', sort(@exp_seq6)), 
   "Test parse_acefile; Checking member_id list ($seq_n members)")
    or diag("Looks like this has failed");

throws_ok { PhyGeCluster::parse_acefile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to parse_acefile function';

throws_ok { PhyGeCluster::parse_acefile('Fake') } qr/ARG. ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { PhyGeCluster::parse_acefile({t => 1}) } qr/ARG. ERROR: 'ace/, 
    'TESTING DIE ERROR when href arg. does not contain acefile key';


#######################
## LOADING FUNCTIONS ##
#######################

## Testing load_seqfile, TEST 63 to 66

my $phygecluster3 = PhyGeCluster->new({ blastfile => $blastfile, 
					fastblast_parser => 1,  });
$phygecluster3->load_seqfile({ sequencefile => $seqfile });

my %clusters7 = %{$phygecluster3->get_clusters()};

my $wrong_sequences = 0;
foreach my $cluster_name7 (keys %clusters7) {
    my $cluster_obj7 = $clusters7{$cluster_name7};
    my @members7 = $cluster_obj7->get_members();
    
    foreach my $member7 (@members7) {
	if ($member7->seq() ne $input_seqs{$member7->display_id()}->seq()) {
	    $wrong_sequences++;
	}
    }
}
is($wrong_sequences, 0,
    "Test load_sequencefile; Checking sequences one by one")
    or diag("Looks like this has failed");

throws_ok { $phygecluster3->load_seqfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to load_seqfile function';

throws_ok { $phygecluster3->load_seqfile('Fake') } qr/ARG. ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { $phygecluster3->load_seqfile({t => 1}) } qr/ARG. ERROR: 'se/, 
    'TESTING DIE ERROR when href arg. does not contain seqfile key';


## Testing load_strainfile, TEST 67 to 71

$phygecluster3->load_strainfile({ strainfile => $strainfile });
my %get_strains2 = %{$phygecluster3->get_strains()};

my %uniq_strains2 = ();
foreach my $key_strain2 (keys %get_strains2) {
    my $strain2 = $get_strains2{$key_strain2};
    unless (exists $get_strains2{$key_strain2}) {
	$uniq_strains2{$strain2} = 1;
    }
    else {
	$uniq_strains2{$strain2}++;
    }
}

is(scalar(keys(%get_strains2)), $strain_seq_n, 
   "Test load_strainfile; Checking number of sequences")
   or diag("Looks like this failed");

is(scalar(keys(%uniq_strains2)), scalar(keys(%expected_strains)),
   "Test load_strainfile; Checking number of strains")
   or diag("Looks like this failed");

throws_ok { $phygecluster3->load_strainfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to load_strainfile function';

throws_ok { $phygecluster3->load_strainfile('Fake') } qr/ARG. ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { $phygecluster3->load_strainfile({t => 1}) } qr/ARG. ERROR: 'str/, 
    'TESTING DIE ERROR when href arg. does not contain strainfile key';


## Testing clone function ## TEST 72 to 75

my $phygecluster_cloned = $phygecluster3->clone();

is(ref($phygecluster_cloned), 'PhyGeCluster',
   "Testing clone(), checking object identity")
    or diag("Looks like this has failed");

is($phygecluster_cloned ne $phygecluster3, 1, 
   "Testing clone(), checking identity cloned object is dif. from source")
    or diag("Looks like this has failed");

is(scalar(keys %{$phygecluster_cloned->get_clusters()}), 
   scalar(keys %{$phygecluster3->get_clusters()}),
   "Testing clone(), same cluster number than original")
    or diag("Looks like this has failed");

is(scalar(keys %{$phygecluster_cloned->get_strains()}), 
   scalar(keys %{$phygecluster3->get_strains()}),
   "Testing clone(), same strains number than original")
    or diag("Looks like this has failed");

## Testing the creation of an object using acefile instead blastfile, TEST 76
## and 77

my $phygecluster_ace = PhyGeCluster->new(
    { 
	acefile    => $acefile,
	strainfile => $strainfile,
    }
    );

my %clusters_ace1 = %{$phygecluster_ace->get_clusters()};
is(scalar(keys %clusters_ace1), scalar(keys %clusters6),
    "Testing new for acefile argument, checking cluster number")
    or diag("Looks like this has failed");

throws_ok { PhyGeCluster->new({ acefile => 1, blastfile => 1}) } qr/ARGUMENT /, 
    'TESTING DIE ERROR when none arg. is supplied to load_strainfile function';

## testing adding a reference using blast (search_homolog function), 
## TEST 78 to 81

## It will clone the object before test the homologous_search function

my $phygecluster_ace2 = $phygecluster_ace->clone();
$phygecluster_ace2->homologous_search({ 
    blast  => [ -p => 'tblastx', -d => $blastdbfile, -e => '1e-10', -a => 2],
    strain => 'Sly',
				     });

my %clusters_ace2 = %{$phygecluster_ace2->get_clusters()};
my %strains_ace2 = %{$phygecluster_ace2->get_strains()};

## Now it will compare that it has at least one Sly strain for each cluster.

my $strain_hmlg_count = 0;
my $more_than_one_hmlg = 0;
foreach my $cl_aceid (keys %clusters_ace1) {
    my %acememb_1 = ();
    my @members1 = $clusters_ace1{$cl_aceid}->get_members();
    foreach my $memb1 (@members1) {
	my $id1 = $memb1->display_id();
	$acememb_1{$id1} = $memb1;
    }
    my $diff_strain_count = 0;
    my @members2 = $clusters_ace2{$cl_aceid}->get_members();
    foreach my $memb2 (@members2) {
	my $id2 = $memb2->display_id();

	unless (defined $acememb_1{$id2}) {
	    if ($strains_ace2{$id2} eq 'Sly') {
		$strain_hmlg_count++;
		$diff_strain_count++;
	    }
	}
    }
    if ($diff_strain_count > 1) {
	$more_than_one_hmlg++;
    }
}

is($strain_hmlg_count <=> 0, 1,
    "Testing homologous_search, checking strains homologous count")
    or diag("Looks like this has failed");

is($more_than_one_hmlg <=> 0, 0,
    "Testing homologous_search, checking more than one homologous per cluster")
    or diag("Looks like this has failed");



## Testing homologous_search with constrains using filter

my $phygecluster_ace3 = $phygecluster_ace->clone();
$phygecluster_ace3->homologous_search({ 
    blast  => [ -p => 'tblastx', -d => $blastdbfile, -e => '1e-10', -a => 2],
    strain => 'Sly',
    filter => { percent_identity    => ['>', 60],
		hsp_length          => ['>', 50],
    }
				      });

my %clusters_ace3 = %{$phygecluster_ace3->get_clusters()};
my %strains_ace3 = %{$phygecluster_ace3->get_strains()};

## Now it will compare that it has at least one Sly strain for each cluster.

my $strain_hmlg_count2 = 0;
my $more_than_one_hmlg2 = 0;
foreach my $cl_aceid2 (keys %clusters_ace1) {
    my %acememb_1 = ();
    my @members1 = $clusters_ace1{$cl_aceid2}->get_members();
    foreach my $memb1 (@members1) {
	my $id1 = $memb1->display_id();
	$acememb_1{$id1} = $memb1;
    }
    my $diff_strain_count2 = 0;
    my @members3 = $clusters_ace3{$cl_aceid2}->get_members();
    foreach my $memb3 (@members3) {
	my $id3 = $memb3->display_id();

	unless (defined $acememb_1{$id3}) {
	    if ($strains_ace3{$id3} eq 'Sly') {
		$strain_hmlg_count2++;
		$diff_strain_count2++;
	    }
	}
    }
    if ($diff_strain_count2 > 1) {
	$more_than_one_hmlg2++;
    }
}

is($strain_hmlg_count2 <=> 0, 1,
    "Testing homologous_search (using filter), checking strains homologous")
    or diag("Looks like this has failed");

is($more_than_one_hmlg2 <=> 0, 1,
    "Testing homologous_search (using filter), checking more than hmlg/cluster")
    or diag("Looks like this has failed");

## Check the die functions for homologous_search, TEST 82 to 90

throws_ok { $phygecluster_ace3->homologous_search() } qr/ARG. ERROR: None/, 
    'TESTING DIE ERROR when no arguments were supplied to homologous_search';

throws_ok { $phygecluster_ace3->homologous_search('fa') } qr/ARG. ERROR: Arg=/, 
    'TESTING DIE ERROR when args. supplied to homologous_search arent hashref';

throws_ok { $phygecluster_ace3->homologous_search({}) } qr/ARG. ERROR: No bl/, 
    'TESTING DIE ERROR when no blast args were supplied to homologous_search';

my $fhrf1 = { blast => 'fake'};
throws_ok { $phygecluster_ace3->homologous_search($fhrf1) } qr/ARG. ERROR: bl/, 
    'TESTING DIE ERROR when no blast args were supplied to homologous_search';

my $fhrf2 = { blast => [ -p => 'blastn']};
throws_ok { $phygecluster_ace3->homologous_search($fhrf2) } qr/ARG. ERROR: No/, 
    'TESTING DIE ERROR when no blast db were supplied to homologous_search';

my $fhrf3 = { blast  => [ -p => 'blastn', -d => $blastdbfile], 
	      filter => 'fake' };
throws_ok { $phygecluster_ace3->homologous_search($fhrf3) } qr/ARG. ERROR: fi/, 
    'TESTING DIE ERROR when filter arg. supplied homologous_search isnot hash';

my $fhrf4 = { blast  => [ -p => 'blastn', -d => $blastdbfile], 
	      filter => { gaps => 'fake'} };
throws_ok { $phygecluster_ace3->homologous_search($fhrf4) } qr/ARG. ERROR: Va/, 
    'TESTING DIE ERROR when filter arg. supplied homologous_search isnot hash';

my $fhrf5 = { blast  => [ -p => 'blastn', -d => $blastdbfile], 
	      filter => { gaps => ['>', 'fake']} };
throws_ok { $phygecluster_ace3->homologous_search($fhrf5) } qr/WRONG filterva/, 
    'TESTING DIE ERROR when filter arg. supplied homologous_search isnot int.';

my $fhrf6 = { blast  => [ -p => 'blastn', -d => $blastdbfile], 
	      filter => { gaps => ['??', 1]} };
throws_ok { $phygecluster_ace3->homologous_search($fhrf6) } qr/WRONG filterva/, 
    'TESTING DIE ERROR when filter arg. supplied homologous_search isnot perm.';



##########################
## ANALYTICAL FUNCTIONS ##
##########################

## Checking cluster_sizes function, TEST 91 to 96

my $phygecluster4 = PhyGeCluster->new(
     { 
 	blastfile        => $blastfile,
 	fastblast_parser => 1,
 	sequencefile     => $seqfile,
 	strainfile       => $strainfile,
     }
     );
my %cl_sizes = $phygecluster4->cluster_sizes();
my @cl_sizes_desc = sort({ $cl_sizes{$b} <=> $cl_sizes{$a}} keys %cl_sizes);
my @cl_sizes_asc = sort({ $cl_sizes{$a} <=> $cl_sizes{$b}} keys %cl_sizes);

## Checking some values that depends of the test files.
## Min. cluster size = 1 (singlets)
## Max. cluster size = 54
## Singlets number = 17
## Contigs number = 41
## Contigs number bigger than size 3 = 23

is($cl_sizes{$cl_sizes_asc[0]}, 1, 
   "Testing cluster_sizes; Checking min. value")
    or diag("Looks like this has failed");
is($cl_sizes{$cl_sizes_desc[0]}, 54, 
   "Testing cluster_sizes; Checking max. value")
    or diag("Looks like this has failed");

my %cl_singlets = $phygecluster4->cluster_sizes(1);
is(scalar(keys %cl_singlets), 17,
   "Testing cluster_sizes to get singlets, Checking singlet count")
    or diag("Looks like this has failed");

my %cl_contigs = $phygecluster4->cluster_sizes(2,54);
is(scalar(keys %cl_contigs), 41,
   "Testing cluster_sizes to get contigs, Checking contigs count")
    or diag("Looks like this has failed");

my %cl_contigs_bigger4 = $phygecluster4->cluster_sizes(4,54);
is(scalar(keys %cl_contigs_bigger4), 23,
   "Testing cluster_sizes to get contigs bigger than 4, Checking contigs count")
    or diag("Looks like this has failed");

## Check a cluster size that does not exist

is($phygecluster3->cluster_sizes(1000000), 0,
   "Testing cluster_sizes for non exists values, Checking 0")
    or diag("Looks like this has failed");


## Checking run_alignments
## The program used will be clustalw.
## Testing died functions, TEST 97 to 102

throws_ok { $phygecluster3->run_alignments() } qr/ARG. ERROR: None/, 
    'TESTING DIE ERROR when none args. were supplied to run_alignments';

throws_ok { $phygecluster3->run_alignments('fake') } qr/ARG. ERROR: Arg=/, 
    'TESTING DIE ERROR when arg. supplied to run_alignments is not hashref.';

throws_ok { $phygecluster3->run_alignments({ t => 1}) } qr/ARG. ERROR: Args/, 
    'TESTING DIE ERROR when arg. without program is supplied to run_alignments';

my $run_href1 = { program => 'fake'};
throws_ok { $phygecluster3->run_alignments($run_href1) } qr/ARG. ERROR: prog/, 
    'TESTING DIE ERROR when program supplied to run_alignments is not valid';

my $run_href2 = { program => 'clustalw' };
throws_ok { $phygecluster3->run_alignments($run_href2) } qr/ARG. ERROR: Args/, 
    'TESTING DIE ERROR when arg. without param. is supplied to run_alignments';

my $run_href3 = { program => 'clustalw', parameters => 'fake' };
throws_ok { $phygecluster3->run_alignments($run_href3) } qr/ARG. ERROR: 'par/, 
    'TESTING DIE ERROR when parameters supplied to run_alignments is not valid';


my @parameters1 = ('quiet' => 'yes', 'matrix' => 'BLOSUM');
$phygecluster3->run_alignments(
    {
	program    => 'clustalw',
	parameters => \@parameters1,
    }
    );

## It will check that all the singlets have not any Bio::SimpleAlign object
## and the clusters have one. TEST 103 to 105

my $singlets_aligns_expected = 0;
my $contigs_aligns_expected = 41;
my $diff_objects_expected = 0;
my $singlets_aligns_count = 0;
my $contigs_aligns_count = 0;
my $diff_objects_count = 0;

my %clusters8 = %{$phygecluster3->get_clusters()};
foreach my $cl_id (keys %clusters8) {
    my @membs = $clusters8{$cl_id}->get_members();
    my $memb_n = scalar(@membs);
    my $align = $clusters8{$cl_id}->alignment();

    if ($memb_n == 0) {
	if (defined $align) {
	    $singlets_aligns_count++;
	    unless (ref($align) eq 'Bio::SimpleAlign') {
		$diff_objects_count++;
	    }
	}
    }
    else {
	if (defined $align) {
	    $contigs_aligns_count++;
	    unless (ref($align) eq 'Bio::SimpleAlign') {
		$diff_objects_count++;
	    }
	}
    }
}

is($singlets_aligns_count, $singlets_aligns_expected,
   "Testing run_alignment, checking alignment absent for singlets")
    or diag("Looks like this has failed");

is($contigs_aligns_count, $contigs_aligns_expected,
   "Testing run_alignment, checking alignment objects for contigs")
    or diag("Looks like this has failed");

is($diff_objects_count, $diff_objects_expected,
   "Testing run_alignment, checking absent of other objects")
    or diag("Looks like this has failed");

## Check that the phygecluster_cloned has at least one align object, TEST 106

my $phygecluster_cloned2 = $phygecluster3->clone();
my $cloned_alignments = 0;
my %cloned_clusters = %{$phygecluster_cloned2->get_clusters()};
foreach my $cl_cluster_id (keys %cloned_clusters) {
    my $cloned_seqfam = $cloned_clusters{$cl_cluster_id};
    my $cloned_align = $cloned_seqfam->alignment();
    if (defined $cloned_align) {
	if(ref($cloned_align) eq 'Bio::SimpleAlign') {
	    $cloned_alignments++;
	}
    }
}

is($cloned_alignments <=> 0, 1,
   "Testing clone(), checking that exists more than one Bio::SimpleAlign obj.")
    or diag("Looks like this has failed");



## testing run_distances() function, and check that the identity of the objects 
## is Bio::Matrix::PhylipDist, TEST 107 to 109

$phygecluster3->run_distances();
my %dist_jc = %{$phygecluster3->get_distances()};

my $wrongobjs = 0;
foreach my $cl_id8 (keys %dist_jc) {
    my $biomat = $dist_jc{$cl_id8};
    
    unless(ref($biomat) eq 'Bio::Matrix::PhylipDist') {
	$wrongobjs++;
    }
}

throws_ok { $phygecluster3->run_distances('fake') } qr/ERROR METHOD: fake/, 
    'TESTING DIE ERROR when a non available method is used with run_distance';

is($wrongobjs, 0, 
    "Testing run_distance() with default method, checking wrong objects")
    or diag("Looks like this has failed");

$phygecluster3->run_distances('Kimura');
my %dist_kim = %{$phygecluster3->get_distances()};

## Kimura method should give a different number of distance matrix than 
## JukesCantor method (by default) because it has not the gap >= length rule
## It should be 41 for Kimura and 16 for JukesCantor

my $diff_dist_count = 0;
my $jc_count = scalar(keys %dist_jc);
my $kim_count = scalar(keys %dist_kim);
if ($jc_count != $kim_count) {
    $diff_dist_count = 1;
}

is($diff_dist_count, 1, 
    "Testing run_distance('Kimura'), checking different number of objects")
    or diag("Looks like this has failed");


##############################
## Checking prune functions ##
##############################

## Checking prune_by_align, TEST 110 to 112

my $cluster_count1 = scalar(keys( %{$phygecluster3->get_clusters}));

my $phygecluster_cl1 = $phygecluster3->clone();
my %rmclusters1 = $phygecluster_cl1->prune_by_align(
    {
	num_sequences => [ '<', 4 ],
	length        => [ '<', 100 ],
    }
    );

my $cluster_count2 = scalar(keys( %{$phygecluster_cl1->get_clusters}));
my $rm_count1 = scalar(keys %rmclusters1);


is($cluster_count2 > 0, 1,
    "Testing prune_by_align, counting removed clusters is not zero")
    or diag("Looks like this has failed");

is($cluster_count1 - $rm_count1, $cluster_count2,
    "Testing prune_by_align, counting removed clusters")
    or diag("Looks like this has failed");

my $wrong_rm_sequences = 0;
foreach my $rm_clusters (keys %rmclusters1) {
    my $num_sequences = $rmclusters1{$rm_clusters}->alignment()
	                                          ->num_sequences();
    if ($num_sequences >= 4) {
	$wrong_rm_sequences++;
    }
}

is($wrong_rm_sequences, 0, 
    "Testing prune_by_align, counting wrong removed clusters")
    or diag("Looks like this has failed");

## Checking croak functions for prune_by_align, TEST 113 to 120

throws_ok { $phygecluster3->prune_by_align() } qr/ARG. ERROR: None args./, 
    'TESTING DIE ERROR when none argument is used with prune_by_align';

throws_ok { $phygecluster3->prune_by_align('fake') } qr/ARG. ERROR: Arg/, 
    'TESTING DIE ERROR when argument used with prune_by_align is not hashref';

throws_ok { $phygecluster3->prune_by_align({ fake => 1}) } qr/ARG. ERROR: f/,
    'TESTING DIE ERROR when key argument used is not permited';

throws_ok { $phygecluster3->prune_by_align({score => 1}) } qr/ARG. ERROR: 1/,
    'TESTING DIE ERROR when value argument used is not array ref.';

throws_ok { $phygecluster3->prune_by_align({score => []}) } qr/ERROR: None/,
    'TESTING DIE ERROR when none condition was used for a function';

throws_ok { $phygecluster3->prune_by_align({score => [1]}) } qr/ERROR: None/,
    'TESTING DIE ERROR when none value was used for a function';

throws_ok { $phygecluster3->prune_by_align({score => ['>','t']})} qr/ERROR:/,
    'TESTING DIE ERROR when the value used in the condition is not an integer';

throws_ok { $phygecluster3->prune_by_align({score => ['*',1]}) } qr/ERROR:/,
    'TESTING DIE ERROR when the condition used is not an int. operator';


##########################################
## Check the function: prune_by_strains ##
##########################################

## Define first pruning arguments ##
## It will create a new object with the function clone to keep the old object

my %srcclusters = %{$phygecluster3->get_clusters()};

## 1) Get three members at random, TEST 121 to 126

my $phygecluster_cl2 = $phygecluster3->clone();
my $prune_args1 = { 
    composition => { 'Sly' => 1, 'Nsy' => 1, 'Nto' => 1},
};

my ($rmcluster2, $rmmemb2) = $phygecluster_cl2->prune_by_strains($prune_args1);
my %rmclusters2 = %{$rmcluster2};
my %rmmembers2 = %{$rmmemb2};

my %selclusters2 = %{$phygecluster_cl2->get_clusters()};

is(scalar(keys %selclusters2) != 0, 1,
    "Testing prune_by_strains, checking prune cluster > 0")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters2) != scalar(keys %srcclusters), 1,
    "Testing prune_by_strains, checking dif. between prune and unprune cluster")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters2) != scalar(keys %rmclusters2), 1,
    "Testing prune_by_strains, checking dif. between sel and rm clusters")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters2) + scalar(keys %rmclusters2) == 
   scalar(keys %srcclusters), 1,
    "Testing prune_by_strains, checking prune_rm+prune_select=unprune clusters")
    or diag("Looks like this has failed");

## Check if all the clusters have three members and if each of them is the
## specified strain

my $clt_diff3mb = 0;
my $clt_diffstr = 0;
my %strains2 = %{$phygecluster_cl2->get_strains()};


foreach my $fam_id (sort keys %selclusters2) {
    my $famseq2 = $selclusters2{$fam_id};
    my @memb2list = $famseq2->get_members();

    my %comp2 = ( 'Sly' => 1, 'Nsy' => 1, 'Nto' => 1 );
    
    if (scalar(@memb2list) != 3) {
	$clt_diff3mb++;
    }
    foreach my $member2 (@memb2list) {
	my $member_id2 = $member2->display_id();
	my $member2str = $strains2{$member_id2};
	$comp2{$member2str}--;
    }
    my $comp_errors = 0;
    foreach my $key (keys %comp2) {
	if ($comp2{$key} != 0) {
	    $clt_diffstr++;
	}
    }
}

is($clt_diff3mb, 0,
    "Testing prune_by_strains, checking selected clusters with 3 members")
    or diag("Looks like this has failed");

is($clt_diffstr, 0,
    "Testing prune_by_strains, checking clusters with different strain comp.")
    or diag("Looks like this has failed");


## 2) Get three members for two min_distance constraints, TEST 127 to 134

my $phygecluster_cl3 = $phygecluster3->clone();
my $prune_args3 = { 
    composition  => { 'Sly' => 1, 'Nsy' => 1, 'Nto' => 1 },
    min_distance => [ [ 'Sly', 'Nsy' ], [ 'Sly', 'Nto' ] ]
};

my ($rmcluster3, $rmmemb3) = $phygecluster_cl3->prune_by_strains($prune_args3);
my %rmclusters3 = %{$rmcluster3};
my %rmmembers3 = %{$rmmemb3};

my %selclusters3 = %{$phygecluster_cl3->get_clusters()};

is(scalar(keys %selclusters3) != 0, 1,
    "Testing prune_by_strains with min_distance, checking prune cluster > 0")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters3) != scalar(keys %srcclusters), 1,
    "Testing prune_by_strains with min_distance, checking dif. prune-unprune")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters3) != scalar(keys %rmclusters3), 1,
    "Testing prune_by_strains with min_distance, checking dif. sel-rm clusters")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters3) + scalar(keys %rmclusters3) == 
   scalar(keys %srcclusters), 1,
    "Testing prune_by_strains with min_distance, checking rm+sel=unprune")
    or diag("Looks like this has failed");

## To know if the distance is the min. distance it will take the distance
## matrix and compare

my $clt_diff3mb3 = 0;
my $clt_diffstr3 = 0;
my $clt_mindist3_ab = 0;
my $clt_mindist3_ac = 0;
my %strains3 = %{$phygecluster_cl3->get_strains()};
my %distance3 = %{$phygecluster_cl3->get_distances()};

foreach my $fam_id3 (sort keys %selclusters3) {
    my $famseq3 = $selclusters3{$fam_id3};
    my $unp_famseq3 = $srcclusters{$fam_id3};

    my $distmt3 = $distance3{$fam_id3};
    my @memb3list = $famseq3->get_members();
    my @unp_memb3list = $unp_famseq3->get_members();

    my $test1 = join(',', @memb3list);
    my $test2 = join(',', @unp_memb3list);

    my %comptypes = ( 'prune' => \@memb3list, 'unprune' => \@unp_memb3list );

    my @str_ab = ();    
    my @str_ac = ();
    my @unp_str_ab = ();
    my @unp_str_ac = ();
    foreach my $type (keys %comptypes) {

	my @members = @{$comptypes{$type}};
	foreach my $mb_aseq (@members) {

	    my $mb_a = $mb_aseq->display_id();
	    foreach my $mb_bseq (@members) {

		my $mb_b = $mb_bseq->display_id();
		if ($strains3{$mb_a} eq 'Sly') {
		    if ($strains3{$mb_b} eq 'Nsy') {

			my $dist_val = $distmt3->get_entry($mb_a, $mb_b);
			if ($type eq 'prune') {

			    push @str_ab, $dist_val;
			}
			else {
			     push @unp_str_ab, $dist_val;
			}
		    }
		    elsif ($strains3{$mb_b} eq 'Nto') {
			my $dist_val = $distmt3->get_entry($mb_a, $mb_b);
			if ($type eq 'prune') {
			    push @str_ac, $dist_val;
			}
			else {
			     push @unp_str_ac, $dist_val;
			}
		    }
		}
	    }
	}
    }

    ## Now it will get the min for each relation

    my @min_str_ab = sort { $a <=> $b } @str_ab;
    my @min_str_ac = sort { $a <=> $b } @str_ac;
    my @unpmin_str_ab = sort { $a <=> $b } @unp_str_ab;
    my @unpmin_str_ac = sort { $a <=> $b } @unp_str_ac;

    if ($min_str_ab[0] != $unpmin_str_ab[0]) {	
	$clt_mindist3_ab++;
    }
     if ($min_str_ac[0] != $unpmin_str_ac[0]) {
	$clt_mindist3_ac++;
    }
    

    my %comp3 = ( 'Sly' => 1, 'Nsy' => 1, 'Nto' => 1 );
    
    if (scalar(@memb3list) != 3) {
	$clt_diff3mb3++;
    }
    foreach my $member3 (@memb3list) {
	my $member_id3 = $member3->display_id();
	my $member3str = $strains3{$member_id3};
	$comp3{$member3str}--;
    }
    my $comp_errors = 0;
    foreach my $key (keys %comp3) {
	if ($comp3{$key} != 0) {
	    $clt_diffstr3++;
	}
    }
}

is($clt_diff3mb3, 0,
    "Testing prune_by_strains min_distance, checking selected with 3 members")
    or diag("Looks like this has failed");

is($clt_diffstr3, 0,
    "Testing prune_by_strains min_distance, checking with diff. strain comp.")
    or diag("Looks like this has failed");

is($clt_mindist3_ab, 0,
    "Testing prune_by_strains min_distance, checking all min_distance a-b")
    or diag("Looks like this has failed");

is($clt_mindist3_ac, 0,
    "Testing prune_by_strains min_distance, checking all min_distance a-c")
    or diag("Looks like this has failed");


## 3) Get five members for six min_distance constraints, TEST 135 to 144 

my $phygecluster_cl4 = $phygecluster3->clone();
my $prune_args4 = { 
    composition  => { 'Sly' => 1, 'Nsy' => 1, 'Nto' => 1, 'Nta' => 2 },
    min_distance => [ 
	[ 'Nta', 'Nta' ],
	[ 'Nta', 'Nsy' ], 
	[ 'Nta', 'Nto' ], 
	[ 'Nta', 'Sly' ],
	[ 'Nsy', 'Sly' ],
	[ 'Nto', 'Sly' ],		      
	],
};

my ($rmcluster4, $rmmemb4) = $phygecluster_cl4->prune_by_strains($prune_args4);
my %rmclusters4 = %{$rmcluster4};
my %rmmembers4 = %{$rmmemb4};

my %selclusters4 = %{$phygecluster_cl4->get_clusters()};

is(scalar(keys %selclusters4) != 0, 1,
    "Testing prune_by_strains with min_distance2, checking prune cluster > 0")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters4) != scalar(keys %srcclusters), 1,
    "Testing prune_by_strains with min_distance2, checking dif. prune-unprune")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters4) != scalar(keys %rmclusters4), 1,
    "Testing prune_by_strains with min_distance2, checking dif.sel-rm clusters")
    or diag("Looks like this has failed");

is(scalar(keys %selclusters4) + scalar(keys %rmclusters4) == 
   scalar(keys %srcclusters), 1,
    "Testing prune_by_strains with min_distance2, checking rm+sel=unprune")
    or diag("Looks like this has failed");

## To know if the distance is the min. distance it will take the distance
## matrix and compare

my $clt_diff4mb4 = 0;
my $clt_diffstr4 = 0;
my $clt_mindist4_aa = 0;
my $clt_mindist4_ab = 0;
my $clt_mindist4_ac = 0;
my $clt_mindist4_ad = 0;
my %strains4 = %{$phygecluster_cl4->get_strains()};
my %distance4 = %{$phygecluster_cl4->get_distances()};

foreach my $fam_id4 (sort keys %selclusters4) {
    my $famseq4 = $selclusters4{$fam_id4};
    my $unp_famseq4 = $srcclusters{$fam_id4};

    my $distmt4 = $distance4{$fam_id4};
    my @memb4list = $famseq4->get_members();
    my @unp_memb4list = $unp_famseq4->get_members();

    my %comptypes = ( 'prune' => \@memb4list, 'unprune' => \@unp_memb4list );

    my @str_aa = (); 
    my @str_ab = ();    
    my @str_ac = ();
    my @str_ad = (); 
    my @unp_str_aa = ();
    my @unp_str_ab = ();
    my @unp_str_ac = ();
    my @unp_str_ad = ();
    
    foreach my $type (keys %comptypes) {

	my @members = @{$comptypes{$type}};
	foreach my $mb_aseq (@members) {

	    my $mb_a = $mb_aseq->display_id();
	    foreach my $mb_bseq (@members) {

		my $mb_b = $mb_bseq->display_id();
		if ($strains4{$mb_a} eq 'Nta') {
		    if ($strains4{$mb_b} eq 'Nta') {

			my $dist_val = $distmt4->get_entry($mb_a, $mb_b);
			if ($type eq 'prune') {

			    push @str_aa, $dist_val;
			}
			else {
			     push @unp_str_aa, $dist_val;
			}
		    }
		    elsif ($strains4{$mb_b} eq 'Nsy') {
			my $dist_val = $distmt4->get_entry($mb_a, $mb_b);
			if ($type eq 'prune') {
			    push @str_ab, $dist_val;
			}
			else {
			     push @unp_str_ab, $dist_val;
			}
		    }
		    elsif ($strains4{$mb_b} eq 'Nto') {
			my $dist_val = $distmt4->get_entry($mb_a, $mb_b);
			if ($type eq 'prune') {
			    push @str_ac, $dist_val;
			}
			else {
			     push @unp_str_ac, $dist_val;
			}
		    }
		    elsif ($strains4{$mb_b} eq 'Sly') {
			my $dist_val = $distmt4->get_entry($mb_a, $mb_b);
			if ($type eq 'prune') {
			    push @str_ad, $dist_val;
			}
			else {
			     push @unp_str_ad, $dist_val;
			}
		    }
		}
	    }
	}
    }

    ## Now it will get the min for each relation
    
    my @min_str_aa = sort { $a <=> $b } @str_aa;
    my @min_str_ab = sort { $a <=> $b } @str_ab;
    my @min_str_ac = sort { $a <=> $b } @str_ac;
    my @min_str_ad = sort { $a <=> $b } @str_ad;
    my @unpmin_str_aa = sort { $a <=> $b } @unp_str_aa;
    my @unpmin_str_ab = sort { $a <=> $b } @unp_str_ab;
    my @unpmin_str_ac = sort { $a <=> $b } @unp_str_ac;
    my @unpmin_str_ad = sort { $a <=> $b } @unp_str_ad;

    if ($min_str_aa[0] != $unpmin_str_aa[0]) {	
	$clt_mindist4_aa++;
    }
    if ($min_str_ab[0] != $unpmin_str_ab[0]) {	
	$clt_mindist4_ab++;
    }
    if ($min_str_ac[0] != $unpmin_str_ac[0]) {
	$clt_mindist4_ac++;
    }
    if ($min_str_ad[0] != $unpmin_str_ad[0]) {	
	$clt_mindist4_ad++;
    }
    
    my %comp4 = ( 'Sly' => 1, 'Nsy' => 1, 'Nto' => 1, 'Nta' => 2 );
    
    if (scalar(@memb4list) != 5) {
	$clt_diff4mb4++;
    }
    foreach my $member4 (@memb4list) {
	my $member_id4 = $member4->display_id();
	my $member4str = $strains4{$member_id4};
	$comp4{$member4str}--;
    }
    my $comp_errors = 0;
    foreach my $key (keys %comp4) {
	if ($comp4{$key} != 0) {
	    $clt_diffstr4++;
	}
    }
}

is($clt_diff4mb4, 0,
    "Testing prune_by_strains min_distance2, checking selected with 5 members")
    or diag("Looks like this has failed");

is($clt_diffstr4, 0,
    "Testing prune_by_strains min_distance2, checking with diff. strain comp.")
    or diag("Looks like this has failed");

is($clt_mindist4_aa, 0,
    "Testing prune_by_strains min_distance2, checking all min_distance a-a")
    or diag("Looks like this has failed");

is($clt_mindist4_ab, 0,
    "Testing prune_by_strains min_distance2, checking all min_distance a-b")
    or diag("Looks like this has failed");

is($clt_mindist4_ac, 0,
    "Testing prune_by_strains min_distance2, checking all min_distance a-c")
    or diag("Looks like this has failed");

is($clt_mindist4_ad, 0,
    "Testing prune_by_strains min_distance2, checking all min_distance a-d")
    or diag("Looks like this has failed");

## Finally check the croak for prune_by_strains function, TEST 145 to 151

throws_ok { $phygecluster3->prune_by_strains() } qr/ARG. ERROR: None hash./, 
    'TESTING DIE ERROR when none argument is used with prune_by_strains';

throws_ok { $phygecluster3->prune_by_strains('fake') } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR when argument used with prune_by_align is not hashref';

throws_ok { $phygecluster3->prune_by_strains({}) } qr/ARG. ERROR: No comp/,
    'TESTING DIE ERROR when no composition argument was used';

my $fprune_href1 = { composition => {}, fake => 1};
throws_ok { $phygecluster3->prune_by_strains($fprune_href1) } qr/ERROR: Const/,
    'TESTING DIE ERROR when constraint used is not valid';

my $fprune_href2 = { composition => 'fake'};
throws_ok { $phygecluster3->prune_by_strains($fprune_href2) } qr/ERROR: Value/,
    'TESTING DIE ERROR when value used is not valid';

my $phygecluster3_f1 = $phygecluster3->clone();
$phygecluster3_f1->set_strains({});

my $prunehref3 = { composition => { 'Nsy' => 1, 'Nto' => 1 }};
throws_ok { $phygecluster3_f1->prune_by_strains($prunehref3) } qr/ERROR: No st/,
    'TESTING DIE ERROR when none strains were loaded into the object';

my $phygecluster3_f2 = $phygecluster3->clone();
$phygecluster3_f2->set_distances({});

throws_ok { $phygecluster3_f2->prune_by_strains($prunehref3) } qr/ERROR: No di/,
    'TESTING DIE ERROR when none distances were loaded into the object';


## Checking the outputs croak functions, TEST 152 to 159

my @outfunctions = ('out_clusterfile', 
		    'out_alignfile', 
		    'out_distancefile',
		    'out_bootstrapfile'
    );

foreach my $outfunc (@outfunctions) {
    throws_ok { $phygecluster3->$outfunc('fake') } qr/ARG.ERROR: fake/,
    "TESTING DIE ERROR when a non hash ref. arg is used with $outfunc";
    
    throws_ok { $phygecluster3->$outfunc({'fk' => 1 }) } qr/ARG.ERROR: fk/,
    "TESTING DIE ERROR when a non valid key is used in args. with $outfunc";
}


#############################
## BOOTSTRAPPING FUNCTIONS ##
#############################

## First test if get/set_bootstrapping died when should do it. TEST 160 to 163

throws_ok { $phygecluster3->set_bootstrapping() } qr/ARGUMENT ERROR: No arg/,
    "TESTING DIE ERROR when no argument is supplied with set_bootstrapping()";

throws_ok { $phygecluster3->set_bootstrapping('fake') } qr/ARGUMENT ERROR: fak/,
    "TESTING DIE ERROR when arg. supplied is not hash ref. set_bootstrapping()";

throws_ok { $phygecluster3->set_bootstrapping({ 1 => 't'}) } qr/VAL. ERROR: va/,
    "TESTING DIE ERROR when values for href. arent aref. set_bootstrapping()";

throws_ok { $phygecluster3->set_bootstrapping({1 =>['t']}) } qr/VAL. ERROR: t/,
    "TESTING DIE ERROR when values for aref. arent Bio::SimpleAlign. objects";


## Now it will test run_bootstrapping. TEST 164 to 166

$phygecluster3->run_bootstrapping(
    { 
	datatype   => 'Molecular sequences', 
	replicates => 100,
	quiet      => 'yes',
    }
    );

my %bootstr = %{$phygecluster3->get_bootstrapping()};

is(scalar(keys %bootstr) <=> 0, 1,
    "Testing run_bootstrapping, checking cluster count with boots > 0")
    or diag("Looks like this has failed");

my $wrong_boots_objs = 0;
my $wrong_boots_num = 0;

foreach my $bo_clid (keys %bootstr) {
    my @bo_aligns = @{$bootstr{$bo_clid}};
    
    if (scalar(@bo_aligns) != 100) {
	$wrong_boots_num++;
    }
    foreach my $bo_align (@bo_aligns) {
	unless (ref($bo_align) eq 'Bio::SimpleAlign') {
	    $wrong_boots_objs++;
	}
    }
}

is($wrong_boots_num, 0, 
    "Testing run_bootstrapping, checking number of align for cluster (100)")
    or diag("Looks like this has failed");

is($wrong_boots_objs, 0, 
    "Testing run_bootstrapping, checking align objects are Bio::SimpleAlign.")
    or diag("Looks like this has failed");

## Testing die options for run_bootstrapping

throws_ok { $phygecluster3->run_bootstrapping() } qr/ARG. ERROR: No args./,
    "TESTING DIE ERROR when no argument is supplied with run_bootstrapping()";

throws_ok { $phygecluster3->run_bootstrapping('fk') } qr/ARG. ERROR: Arg/,
    "TESTING DIE ERROR when arg. supplied to run_bootstrapping() isnt hashref";

throws_ok { $phygecluster3->run_bootstrapping({'fk'=>1}) } qr/ARG. ERROR: fk/,
    "TESTING DIE ERROR when arg. supplied to run_bootstrapping() isnt permited";

my $fk_hr1 = { datatype => 'fk'};

throws_ok { $phygecluster3->run_bootstrapping($fk_hr1) } qr/ARG. ERROR: fk/,
    "TESTING DIE ERROR when arg. supplied to run_bootst.() isnt permited value";

my $fk_hr2 = { datatype => 'Sequence', replicates => 'A lots'};

throws_ok { $phygecluster3->run_bootstrapping($fk_hr2) } qr/ARG. ERROR: A lot/,
    "TESTING DIE ERROR when arg. supplied to run_bootst.() isnt an integer";

my $fk_hr3 = { datatype => 'fk', quiet => 1};

throws_ok { $phygecluster3->run_bootstrapping($fk_hr3) } qr/ARG. ERROR: 1/,
    "TESTING DIE ERROR when arg. supplied to run_bootst.() isnt yes|no";


##################################
## TESTING OVERLAPPING REGIONS ###
##################################

my %overlaps = $phygecluster3->calculate_overlaps();

## Define some expected values (counted by other methods), TEST 173

my %ov_expval = ( 
    'cluster_1' => { 
	'Sly_01101' => {
	    'Sly_01101' => {'start' => 0,    'end' => 0,    'length'=> 0    },
	    'Nsy_31457' => {'start' => 1408, 'end' => 1490, 'length'=> 82   },
	    'Nsy_05034' => {'start' => 1195, 'end' => 2400, 'length'=> 1205 },
	    'Nto_10658' => {'start' => 1317, 'end' => 1778, 'length'=> 461  },
	},
	'Nsy_31457' => {
	    'Sly_01101' => {'start' => 1408, 'end' => 1490, 'length'=> 82   },
	    'Nsy_31457' => {'start' => 0,    'end' => 0,    'length'=> 0    },
	    'Nsy_05034' => {'start' => 1408, 'end' => 1490, 'length'=> 82   },
	    'Nto_10658' => {'start' => 1408, 'end' => 1490, 'length'=> 82   },
	},
	'Nsy_05034' => {
	    'Sly_01101' => {'start' => 1195, 'end' => 2400, 'length'=> 1205 },
	    'Nsy_31457' => {'start' => 1408, 'end' => 1490, 'length'=> 82   },
	    'Nsy_05034' => {'start' => 0,    'end' => 0,    'length'=> 0    },
	    'Nto_10658' => {'start' => 1317, 'end' => 1778, 'length'=> 461  },
	},
	'Nto_10658' => {
	    'Sly_01101' => {'start' => 1317, 'end' => 1778, 'length'=> 461  },
	    'Nsy_31457' => {'start' => 1408, 'end' => 1490, 'length'=> 82   },
	    'Nsy_05034' => {'start' => 1317, 'end' => 1778, 'length'=> 461  },
	    'Nto_10658' => {'start' => 0,    'end' => 0,    'length'=> 0    },
	},	
    },
    'cluster_2' => {
	'Sly_01219' => {
	    'Sly_01219' => {'start' => 0,    'end' => 0,    'length'=> 0    },
	    'Nta_08736' => {'start' => 136,  'end' => 1011, 'length'=> 875  },
	},
	'Nta_08736' => {
	    'Sly_01219' => {'start' => 136,  'end' => 1011, 'length'=> 875  },
	    'Nta_08736' => {'start' => 0,    'end' => 0,    'length'=> 0    },
	},
    }
    );

my $wrong_ovvalues = 0;
foreach my $clid_o (keys %overlaps) {
    my $mtx = $overlaps{$clid_o};
    my @rownames = $mtx->row_names();
    my @colnames = $mtx->column_names();
    foreach my $row (@rownames) {
	foreach my $col (@colnames) {
	    my $entry = $mtx->get_entry($row, $col);
	    my @vars = ('start', 'end', 'length');
	    foreach my $v (@vars) {
		if (defined $ov_expval{$clid_o}) {
		    my %ovlp = %{$ov_expval{$clid_o}};
		    if ($ovlp{$row}->{$col}->{$v} ne $entry->{$v}) {
			$wrong_ovvalues++;
		    }
		}
	    }
	}
    }
}

is($wrong_ovvalues, 0, 
    "Testing calculate_overlaps, checking know values for cluster_1 and _2")
    or diag("Looks like this has failed");

my %best_overlaps = $phygecluster3->best_overlaps();

my $best_ovcluster1 = join(',', sort @{$best_overlaps{'cluster_1'}});
my $best_ovcluster2 = join(',', sort @{$best_overlaps{'cluster_2'}});

is($best_ovcluster1, 'Nsy_05034,Sly_01101',
   "Testing best_overlaps, checking best overlap for cluster_1")
    or diag("Looks like this has failed");

is($best_ovcluster2, 'Nta_08736,Sly_01219',
   "Testing best_overlaps, checking best overlap for cluster_2")
    or diag("Looks like this has failed");

####
1; #
####
