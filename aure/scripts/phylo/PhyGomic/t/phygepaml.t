#!/usr/bin/perl

=head1 NAME

  phygepaml.t
  A piece of code to test the PhyGePaml module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl phygepaml.t
 prove phygepaml.t

=head1 DESCRIPTION

 Test PhyGePaml module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 29;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1 to 6

BEGIN {
    use_ok('PhyGeCluster');
    use_ok('PhyGeTopo');
    use_ok('Bio::Tree::TopoType');
    use_ok('Bio::Tree::Node');
    use_ok('Bio::Tree::Tree');
    use_ok('PhyGePaml');
}


## Create an empty object and test the possible die functions. TEST 7 to 10

my $phygepaml0 = PhyGePaml->new();

is(ref($phygepaml0), 'PhyGePaml', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { PhyGePaml->new(['fake']) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { PhyGePaml->new({ fake => {} }) } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a valid arg.';

throws_ok { PhyGePaml->new({ seqfams => 'fk2'}) } qr/ERROR: fk2/, 
    'TESTING DIE ERROR for new() when arg. doesnt have valid value';


###############
## ACCESSORS ##
###############


## Sequences

my %seq = (
    seq0 => 'ATGCGTGAGACTAGACAGTTGACCAGGTAGACGATGAGATCAGAGTCAGGGATTTTCAGGATTGA',
    seq1 => 'ATGCGTGAGTTTAGACAGTTGACCAGGTAGACACGGAGATCAGAGTCAGGGAATTTCACGATTGA',
    seq2 => 'ATGCGTGAGTTTAGACAGTAGACCAGGTAGACCCGGAGATAAGAGTCAGGGAATATCACGATTGA',
    seq3 => 'ATGCGTGAGTTTAGACAGTAGACCAGGAAGACCCGGAGATAAGAGCCAGGGAATATCAGGACTGA',
    seq4 => 'ATGCGGACGATGGCAGTTGGGTGGGTGCAGGACGAGAGAGAGCTGACGTGGACGATGGACGATGA',
    seq5 => 'ATGCGGACGATGGCAGTTCCGTGGGTGCAGGACGAGTGAGAGCTGACGTGGACGCTGGACGATGA',
    seq6 => 'ATGCGGACGTTGGCAGTTCCGTGCGTGCAGGACGAGTGCGAGCTGTCGTGGACGCCGGACGTTGA',
    seq7 => 'ATGCGGACGTTGGCAGTTCCGTGCGAGCAGGACGAGTGCGAGCTGACGTGGACGGCGGACGTTGA',
    seq8 => 'ATGCGACGTGACCGATGGACAAAAAAGCTAGGCACATCATTTACATTACGGGACAGGGATTGA',
    seq9 => 'ATGGGACGTGACCGATGGACAAATAAGCTAGGCACATCATTTACAATACGGCACAGGGAATGA',
    seq10 => 'ATGGGTCGTGACCGATGGACCAATAAGCTAGGCACATCATTTACACTACGGCAGAGGGAATGA'
    );

## 1) CHECK GET/SET_SEQFAMS, TEST 11 to 15
##    Create three empty seqfams objects

my %seqs0 = ( 
    seq0 => Bio::Seq->new(-id  => 'seq0', -seq => $seq{seq0}),
    seq1 => Bio::Seq->new(-id  => 'seq1', -seq => $seq{seq1}),
    seq2 => Bio::Seq->new(-id  => 'seq2', -seq => $seq{seq2}),
    seq3 => Bio::Seq->new(-id  => 'seq3', -seq => $seq{seq3}),
    seq4 => Bio::Seq->new(-id  => 'seq4', -seq => $seq{seq4}),
    seq5 => Bio::Seq->new(-id  => 'seq5', -seq => $seq{seq5}),
    seq6 => Bio::Seq->new(-id  => 'seq6', -seq => $seq{seq6}),
    seq7 => Bio::Seq->new(-id  => 'seq7', -seq => $seq{seq7}),
    seq8 => Bio::Seq->new(-id  => 'seq8', -seq => $seq{seq8}),
    seq9 => Bio::Seq->new(-id  => 'seq9', -seq => $seq{seq9}),
    seq10 => Bio::Seq->new(-id  => 'seq10', -seq => $seq{seq10}),
    );

my $seqmemb01 = [ $seqs0{seq0}, $seqs0{seq1}, $seqs0{seq2}, $seqs0{seq3} ];
my $seqmemb02 = [ $seqs0{seq4}, $seqs0{seq5}, $seqs0{seq6}, $seqs0{seq7} ];
my $seqmemb03 = [ $seqs0{seq8}, $seqs0{seq9}, $seqs0{seq10} ];

my %seqfams0 = (
    'cl_id1' => Bio::Cluster::SequenceFamily->new( -family_id => 'cl_id1',
						    -members   => $seqmemb01 ),
    'cl_id2' => Bio::Cluster::SequenceFamily->new( -family_id => 'cl_id2',
						    -members   => $seqmemb02 ),
    'cl_id3' => Bio::Cluster::SequenceFamily->new( -family_id => 'cl_id3',
						    -members   => $seqmemb03 ),
    );

## Test set/get_seqfams function

$phygepaml0->set_seqfams(\%seqfams0);
my %seqfams0g = %{$phygepaml0->get_seqfams()};

is(scalar(keys %seqfams0g), 3,
    "Testing get/set_seqfams, checking number of SeqFam objects")
    or diag("Looks like this has failed");

my $wrong_seqfam_objs = 0;
foreach my $fam_id0 (keys %seqfams0) {
    unless (ref($seqfams0{$fam_id0}) eq 'Bio::Cluster::SequenceFamily') {
	$wrong_seqfam_objs++;
    }
}

is($wrong_seqfam_objs, 0, 
    "Testing get/seq_seqfams, checking SeqFam objs. identity")
    or diag("Looks like this has failed");

throws_ok { $phygepaml0->set_seqfams() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_seqfams function';

throws_ok { $phygepaml0->set_seqfams('fake') } qr/ARG. ERROR: fake set_/, 
    'TESTING DIE ERROR when arg. supplied to set_seqfams() isnt HashRef';

throws_ok { $phygepaml0->set_seqfams({id => 'fake2'}) } qr/ARG. ERROR: fake2 u/, 
    'TESTING DIE ERROR when arg. element supplied to set_seqfams() isnt SeqFam';

## 2) CHECK GET/SET_STRAINS, TEST 16 to 18
##    Create strains

my %strains0 = ( 
    'seq0' => 'Str1',
    'seq1' => 'Str2',
    'seq2' => 'Str3',
    'seq3' => 'Str4',
    'seq4' => 'Str1',
    'seq5' => 'Str2',
    'seq6' => 'Str3',
    'seq7' => 'Str4',
    'seq8' => 'Str1',
    'seq9' => 'Str2',
    'seq10' => 'Str4',
    );

## Test set/get_strains function

$phygepaml0->set_strains(\%strains0);
my %strains0g = %{$phygepaml0->get_strains()};

is(scalar(keys %strains0g), 11,
    "Testing get/set_strains, checking number of strains objects")
    or diag("Looks like this has failed");

throws_ok { $phygepaml0->set_strains() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_strains function';

throws_ok { $phygepaml0->set_strains('fake') } qr/ARG. ERROR: When arg/, 
    'TESTING DIE ERROR when arg. supplied to set_strains() isnt HashRef';



## 3) CHECK GET/SET_CDS, TEST 19 to 24
##    Create cds hash
##    Prepare the cds as consensus of the seqfam

## Calculate the consensus sequence

my %cons = (
    cl_id1 => 
    'ATGCGTGAGTTTAGACAGTTGACCAGGTAGACCCGGAGATCAGAGTCAGGGAATATCACGATTGA',
    cl_id2 => 
    'ATGCGGACGATGGCAGTTCCGTGGGTGCAGGACGAGTGCGAGCTGACGTGGACGCCGGACGTTGA',
    cl_id3 => 
    'ATGGGACGTGACCGATGGACAAATAAGCTAGGCACATCATTTACATTACGGCACAGGGAATGA',
    );

my %consobj = (
    cl_id1 => Bio::Seq->new(-id  => 'cl_id1', -seq => $cons{cl_id1}),
    cl_id2 => Bio::Seq->new(-id  => 'cl_id2', -seq => $cons{cl_id2}),
    cl_id3 => Bio::Seq->new(-id  => 'cl_id3', -seq => $cons{cl_id3}),
    );


$phygepaml0->set_cds(\%consobj);
my %cons_g = %{$phygepaml0->get_cds()};

foreach my $cons_id (sort keys %cons_g) {
    is($cons_g{$cons_id}->seq(), $cons{$cons_id},
       "Testing get/set_cds, checking consensus seq for $cons_id")
	or diag("Looks like this has failed");
}

throws_ok { $phygepaml0->set_cds() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_cds function';

throws_ok { $phygepaml0->set_cds('fake1') } qr/ARG. ERROR: fake1/, 
    'TESTING DIE ERROR when arg. supplied to set_cds() isnt HashRef';

throws_ok { $phygepaml0->set_cds({ test1 => 'fake2'}) } qr/ARG. ERROR: fake2/, 
    'TESTING DIE ERROR when val. for arg. supplied to set_cds() arent Bio:Seq';

## 4) CHECK GET/SET_TOPOTYPES, TEST 25 to 29 

my $phygecl0 = PhyGeCluster->new();
$phygecl0->set_clusters(\%seqfams0);
$phygecl0->set_strains(\%strains0);

my %param1 = (quiet => undef, matrix => 'BLOSUM');
$phygecl0->run_alignments({ program    => 'clustalw', parameters => \%param1 });
$phygecl0->run_distances({ method => 'Kimura' });
$phygecl0->run_mltrees({ dnaml => {} });

my $phygetp0 = PhyGeTopo->new({ seqfams => $phygecl0->get_clusters(), 
				strains => \%strains0 });

my %topotypes0 = $phygetp0->run_topoanalysis();

$phygepaml0->set_topotypes(\%topotypes0);
my %topotypes0g = %{$phygepaml0->get_topotypes()};

foreach my $topo (keys %topotypes0g) {
    my $topology = $topotypes0{$topo}->get_topology_as_newick();
    my $topologyg = $topotypes0g{$topo}->get_topology_as_newick();
    is($topology, $topologyg,
       "Testing get/set_topotypes, checking topology for $topo")
	or diag("Looks like this has failed");
}

throws_ok { $phygepaml0->set_topotypes() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_topotypes function';

throws_ok { $phygepaml0->set_topotypes('fake1') } qr/ARG. ERROR: When arg/, 
    'TESTING DIE ERROR when arg. supplied to set_topotypes() isnt HashRef';

throws_ok { $phygepaml0->set_topotypes({ test1 => 'fake2'}) } qr/ARG. ERROR: V/,
    'TESTING DIE ERROR when values supplied to set_topologies() are wrong';








####
1; #
####
