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
use Test::More tests => 56;
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
    seq0 => 'CTGATGGAGACTAGACAGTTGACCAGGTATACGATGAGATCAGAGTCAGGGATTTTCTGACCG',
    seq1 => 'CTGATGGAGTTTAGACAGTTGACCAGGTATACACGGAGATCAGAGTCAGGGAATTTCTGACCG',
    seq2 => 'CTGATGGAGTTTAGACAGTTGACCAGGTATACCCGGAGATCAGAGTCAGGGAATATCTGACCG',
    seq3 => 'CTGATGGAGTTTAGACAGTTGACCAGGAATACCCGGAGATCAGAGCCAGGGAATATCTGACCG',
    seq4 => 'ATGCGGACGATGGCAGTTGGGTGGGTGCAGGACGAGAGAGAGCTGACGTGGACGATGGACTGA',
    seq5 => 'ATGCGGACGATGGCAGTTCCGTGGGTGCAGGACGAGTGAGAGCTGACGTGGACGCTGGACTGA',
    seq6 => 'ATGCGGACGTTGGCAGTTCCGTGCGTGCAGGACGAGTGCGAGCTGTCGTGGACGCCGGACTGA',
    seq7 => 'ATGCGGACGTTGGCAGTTCCGTGCGAGCAGGACGAGTGCGAGCTGACGTGGACGGCGGACTGA',
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
    'ATGGAGTTTAGACAGTTGACCAGGTATACCCGGAGATCAGAGTCAGGGAATATCTGA',
    cl_id2 => 
    'ATGCGGACGATGGCAGTTCCGTGGGTGCAGGACGAGTGCGAGCTGACGTGGACGCCGGACTGA',
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

## Test get/set_codeml_results, TEST 30 to 33

my %restest = (
    id1 => Bio::Matrix::Generic->new(),
    id2 => Bio::Matrix::Generic->new(),
    );

$phygepaml0->set_codeml_results(\%restest);

is(scalar(keys %{$phygepaml0->get_codeml_results()}), 2, 
    "Testing get/set_codeml_results, checking number of hash elements")
    or diag("Looks like this has failed");

throws_ok { $phygepaml0->set_codeml_results() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_codeml_results()';

throws_ok { $phygepaml0->set_codeml_results('fake1') } qr/ARG. ERROR: When ar/, 
    'TESTING DIE ERROR when arg. supplied to set_codeml_results() isnt HashRef';

throws_ok { $phygepaml0->set_codeml_results({ test1 => 'fake2'}) } qr/ERROR: V/,
    'TESTING DIE ERROR when values supplied to set_codeml_results() are wrong';


######################
## Analytical tools ##
######################

## Calculate the consensus sequence, TES 34 to 36

my %exp_prot = (
    cl_id1 => 'MEFRQLTRYTRRSESGNI*',
    cl_id2 => 'MRTMAVPWVQDECELTWTPD*',
    cl_id3 => 'MGRDRWTNKLGTSFTLRHRE*',
    );


my %proteins = $phygepaml0->translate_cds();
foreach my $prot_id (sort keys %proteins) {
    is ($proteins{$prot_id}->seq(), $exp_prot{$prot_id},
	"Testing translate_cds, checking protein seq for $prot_id")
	or diag("Looks like this has failed");
}


## Align member cds, TEST 37 and 38

my %seqfam_cds_objs = $phygepaml0->align_member_cds();

is(scalar(keys %seqfam_cds_objs), 3, 
    "Testing align_member_cds, checking number of new seqfams")
    or diag("Looks like this has failed");

throws_ok { $phygepaml0->align_member_cds('fk') } qr/ERROR: Parameters used/, 
    'TESTING DIE ERROR when parameters used for align_member_cds are not href';

## Test set_seqfam_cds, TEST 39 to 41

$phygepaml0->set_seqfam_cds();

my %newseqfams = %{$phygepaml0->get_seqfams()};
foreach my $sf_id (sort keys %newseqfams) {
    is(ref($newseqfams{$sf_id}), 'Bio::Cluster::SequenceFamily',
	"Testing set_seqfam_cds, checking ref. object for $sf_id")
	or diag("Looks like this has failed");
}

## Test _cds_by_longest6frame, TEST 42 to 46

my @cds0_data = PhyGePaml::_cds_by_longest6frame($seqs0{seq0}, 1);

is(length($cds0_data[0]->seq()), length($seqs0{seq0}->seq()) - 6, 
    "Testing _cds_by_longest6frame (forcing first met), cheking cds")
    or diag("Looks like this has failed");

throws_ok { PhyGePaml::_cds_by_longest6frame() } qr/ERROR: No seq./, 
    'TESTING DIE ERROR when no arg. was used for _cds_by_longest6frame';

throws_ok { PhyGePaml::_cds_by_longest6frame('fake') } qr/ERROR: Argument/, 
    'TESTING DIE ERROR when arg. supplied _cds_by_longest6frame isnt Bio::Seq';

throws_ok { PhyGePaml::_cds_by_longest6frame($seqs0{seq0}, 't') } qr/ERROR: W/, 
    'TESTING DIE ERROR when arg. supplied _cds_by_longest6frame isn boolean';

my @cds1_data = PhyGePaml::_cds_by_longest6frame($seqs0{seq0}, 0);

is(length($cds1_data[0]->seq()), length($seqs0{seq0}->seq()) - 3, 
    "Testing _cds_by_longest6frame (no first met), cheking cds")
    or diag("Looks like this has failed");


## Test predict_cds, TEST 47 to 51

$phygepaml0->predict_cds({ method    => 'longest6frame', 
			   arguments => { force_firstmet => 1}});

my %cds1 = %{$phygepaml0->get_cds()};

is($cds1{cl_id1}->seq(), 
   'ATGGAGWYTAGACAGTTGACCAGGWATACVMKGAGATCAGAGYCAGGG', 
   "Testing predict_cds, checking sequence of the consensus")
    or diag("looks like this has failed");

throws_ok { $phygepaml0->predict_cds() } qr/ERROR: No parameters/, 
    'TESTING DIE ERROR when no parameters were used for predict_cds';

throws_ok { $phygepaml0->predict_cds('fake') } qr/ERROR: Parameters spec/, 
    'TESTING DIE ERROR when parameters used for predict_cds are not hashref';

throws_ok { $phygepaml0->predict_cds({ fake => 1}) } qr/ERROR: fake/, 
    'TESTING DIE ERROR when parameters used for predict_cds is not permitted';

throws_ok { $phygepaml0->predict_cds({method => 1}) } qr/ERROR: method/, 
    'TESTING DIE ERROR when par. used for predict_cds has no-permitted value';

## Test run_codeml, TEST 52 to 55

$phygepaml0->run_codeml();

my %codeml = %{$phygepaml0->get_codeml_results()};

my $mlmtx1 = $codeml{cl_id1};


is(join(',', sort $mlmtx1->column_names()), 'seq0,seq1,seq2,seq3', 
    "Testing run_codeml, checking column names for mlmatrix (cl_id1)")
    or diag("Looks like this has failed");

is($mlmtx1->entry('seq3', 'seq2')->{'dS'}, 0.0006,
    "Testing run_codeml, checking dS value for mlmatrix (seq3,seq2)")
    or diag("Looks like this has failed");

is($mlmtx1->entry('seq3', 'seq3'), undef,
    "Testing run_codeml, checking undef value for selfentry mlmatrix")
    or diag("Looks like this has failed");

throws_ok { $phygepaml0->run_codeml('fake') } qr/ERROR: fake/, 
    'TESTING DIE ERROR when parameter supplied to run_codeml isnt a hashref.';

## TEST 56

throws_ok { $phygepaml0->out_codeml() } qr /ERROR: No filename/,
    'TESTING DIE ERROR when no filename is supplied to out_codeml';


####
1; #
####
