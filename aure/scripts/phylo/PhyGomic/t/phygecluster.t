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
use Test::More tests => 75;
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

## Testing die option for the wrong use of these functions, TEST 12 to 18

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



## TEST 19 and 20 (find_cluster)

my $cluster3 = $phygecluster0->find_cluster('member_test2');

is($cluster3->family_id(), 'cluster_test1',
   "Test find_cluster function with defined value; checking cluster name")
    or diag("Looks like this has failed");

is($phygecluster0->find_cluster('fake_cluster_test1'), undef,
   "Test find_cluster function with undefined value; checking cluster name")
    or diag("Looks like this has failed");

## TEST 21 to 23 (remove_cluster)

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

## TEST 24 to 27 (get/set strains)

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

## First get the data from the files

my @seqs = ();
my %input_seqs = ();
my $seq_io = Bio::SeqIO->new( -file => $seqfile, -format => 'fasta');

while( my $seq = $seq_io->next_seq()) {
    push @seqs, $seq->display_id();
    $input_seqs{$seq->display_id()} = $seq;
}
my $seq_n = scalar(@seqs);


## Test the function TEST 28 to 33

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
    push @obtained_members, @{$clusters4{$cluster_id}};
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

foreach my $cl4fast (sort keys %clusters4fast) {
    if (exists $clusters4{$cl4fast}) {
	my $reg_members4 = join(',', sort @{$clusters4{$cl4fast}});
	my $fast_members4 = join(',', sort @{$clusters4fast{$cl4fast}});

	if($reg_members4 eq $fast_members4) {
	    $diff1--;
	}
	else {
	    print STDERR "Diff cluster: $cl4fast\tR=$reg_members4 ";
	    print STDERR "vs F=$fast_members4\n";
	}
    }
    else {
	print STDERR "Cluster $cl4fast do not exists in regular method\n";
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
    push @obtained_members5, @{$clusters5{$cluster_id}};
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

foreach my $cl5fast (sort keys %clusters5fast) {
    if (exists $clusters5{$cl5fast}) {
	my $reg_members5 = join(',', sort @{$clusters5{$cl5fast}});
	my $fast_members5 = join(',', sort @{$clusters5fast{$cl5fast}});

	if($reg_members5 eq $fast_members5) {
	    $diff2--;
	}
	else {
	    print STDERR "Diff cluster: $cl5fast\tR=$reg_members5 ";
	    print STDERR "vs F=$fast_members5\n";
	}
    }
    else {
	print STDERR "Cluster $cl5fast do not exists in regular method\n";
    }
}
is($diff2, 0,
   "Testing fastparse_blastfile (id=95,l=30); Checking parse_blastfile compar.")
    or diag("Looks like this has failed");


## Testing die functions for parse_blastfile, TEST 34 to 41

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


## Checking parsed sequence file, TEST 43

my %seqs1 = PhyGeCluster::parse_seqfile({ sequencefile => $seqfile });

is(join(',', sort(keys %seqs1)), join(',', sort(@seqs)),
   "Test parse_seqfile; Cheking sequence_ids")
    or diag("Looks like this failed");

## Testing die functions for parse_seqfile, TEST 44 to 46

throws_ok { PhyGeCluster::parse_seqfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to parse_seqfile function';

throws_ok { PhyGeCluster::parse_seqfile('Fake') } qr/ARG. ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { PhyGeCluster::parse_seqfile({t => 1}) } qr/ARG. ERROR: 'se/, 
    'TESTING DIE ERROR when href arg. does not contain seqfile key';

## Checking parsed strain file, TEST 47 and 48

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

## Testing die for the parse_strainfile function, TEST 49 to 51

throws_ok { PhyGeCluster::parse_strainfile() } qr/ARG. ERROR: No arg/, 
    'TESTING DIE ERROR when none arg. is supplied to parse_strainfile function';

throws_ok { PhyGeCluster::parse_strainfile('Fake') } qr/ARG. ERROR: Fake/, 
    'TESTING DIE ERROR when arg. supplied is not a hash reference';

throws_ok { PhyGeCluster::parse_strainfile({t => 1}) } qr/ARG. ERROR: 'str/, 
    'TESTING DIE ERROR when href arg. does not contain strainfile key';


#######################
## LOADING FUNCTIONS ##
#######################

## Testing load_seqfile, TEST 52 to 55

my $phygecluster3 = PhyGeCluster->new({ blastfile => $blastfile, 
					fastblast_parser => 1,  });
$phygecluster3->load_seqfile({ sequencefile => $seqfile });

my %clusters6 = %{$phygecluster3->get_clusters()};

my $wrong_sequences = 0;
foreach my $cluster_name6 (keys %clusters6) {
    my $cluster_obj6 = $clusters6{$cluster_name6};
    my @members6 = $cluster_obj6->get_members();
    
    foreach my $member6 (@members6) {
	if ($member6->seq() ne $input_seqs{$member6->display_id()}->seq()) {
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


## Testing load_strainfile, TEST 56 to 60

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


##########################
## ANALYTICAL FUNCTIONS ##
##########################

## Checking cluster_sizes function, TEST 53 to 57

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
## Testing died functions, TEST 67 to

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


my @parameters1 = ('quiet' => 1, 'matrix' => 'BLOSUM');
$phygecluster3->run_alignments(
    {
	program    => 'clustalw',
	parameters => \@parameters1,
    }
    );

## It will check that all the singlets have not any Bio::SimpleAlign object
## and the clusters have one.

my $singlets_aligns_expected = 0;
my $contigs_aligns_expected = 41;
my $diff_objects_expected = 0;
my $singlets_aligns_count = 0;
my $contigs_aligns_count = 0;
my $diff_objects_count = 0;


my %clusters7 = %{$phygecluster3->get_clusters()};
foreach my $cl_id (keys %clusters7) {
    my @membs = $clusters7{$cl_id}->get_members();
    my $memb_n = scalar(@membs);
    my $align = $clusters7{$cl_id}->alignment();
    
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

####
1; #
####
