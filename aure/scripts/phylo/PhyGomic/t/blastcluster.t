#!/usr/bin/perl

=head1 NAME

  blastcluster.t
  A piece of code to test the BlastCluster module used for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl blastcluster.t
 prove blastcluster.t

 Note: This test use three external files that should be in the same 
       test folder.
       + testfiles/selfblast.test.m8
       + testfiles/seq.test.fasta
       + testfiles/strains.test.tab



=head1 DESCRIPTION

 Test BlastCluster.pm module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;

use Data::Dumper;
use Test::More;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1,2

BEGIN {
    use_ok('BlastCluster');
    use_ok('Bio::Seq');
}


## TEST 3, create an empty object

my $blastcluster0 = BlastCluster->new();

is(ref($blastcluster0), 'BlastCluster', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

## TEST 4,5,6,7,8,9 and 10 test add/get/remove
## Create an empty Bio::Seq object

my @seq_objs1 = ();
my @seq_ids1 = ('member_test1', 'member_test2', 'member_test3');
foreach my $seq_id1 (@seq_ids1) {
    my $seq = Bio::Seq->new( -id => $seq_id1);
    push @seq_objs1, $seq;
}

$blastcluster0->add_cluster('cluster_test1', \@seq_objs1);
my %clusters1 = %{$blastcluster0->get_clusters()};

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

$blastcluster0->add_cluster('cluster_test2', \@seq_objs2);
my %clusters2 = %{$blastcluster0->get_clusters()};

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

## TEST 11 and 12 (find_cluster)

my $cluster3 = $blastcluster0->find_cluster('member_test2');

is($cluster3->family_id(), 'cluster_test1',
   "Test find_cluster function with defined value; checking cluster name")
    or diag("Looks like this has failed");

is($blastcluster0->find_cluster('fake_cluster_test1'), undef,
   "Test find_cluster function with undefined value; checking cluster name")
    or diag("Looks like this has failed");

## TEST 13 and 14 (remove_cluster)

$blastcluster0->remove_cluster('cluster_test1');
my %clusters3 = %{$blastcluster0->get_clusters()};

is(scalar(keys %clusters3), 1, 
    "Test remove_cluster function; Checking number of clusters")
    or diag("Looks like this has failed");

is(join(',', sort(keys %clusters3)), 'cluster_test2', 
    "Test remove_cluster function; Checking cluster id list")
    or diag("Looks like this has failed");

## TEST 15 and 16 (get/set strains)

my %strains1 = (
    'member_test4' => 'strain1', 
    'member_test5' => 'strain2' 
    ); 

$blastcluster0->set_strains(\%strains1);
my %get_strains1 = %{$blastcluster0->get_strains()};

foreach my $key_strain (keys %get_strains1) {
    is($get_strains1{$key_strain}, $strains1{$key_strain},
       "Test get/set_strains; Checking strain value for $key_strain")
	or diag("Looks like this has failed");
}


##################
## PARSING TEST ##
##################

my $blastfile = "$FindBin::Bin/testfiles/selfblast.test.m8";
my $seqfile = "$FindBin::Bin/testfiles/seq.test.fasta";
my $strainfile = "$FindBin::Bin/testfiles/strains.test.tab";

## First get the data from the files

my @seqs = ();
my $seq_io = Bio::SeqIO->new( -file => $seqfile, -format => 'fasta');

while( my $seq = $seq_io->next_seq()) {
    push @seqs, $seq->display_id();
}
my $seq_n = scalar(@seqs);

## Test the function

my %clusters4 = BlastCluster::parse_blastfile({ blastfile => $blastfile });

my @obtained_members = ();
foreach my $cluster_id (keys %clusters4) {
    push @obtained_members, @{$clusters4{$cluster_id}};
}
is(join(',', sort(@obtained_members)), join(',', sort(@seqs)), 
   "Test parse_blastfile; Checking member_id list ($seq_n members)")
    or diag("Looks like this has failed");


####
1; #
####
