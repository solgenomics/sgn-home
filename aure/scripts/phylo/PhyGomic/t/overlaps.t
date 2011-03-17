#!/usr/bin/perl

=head1 NAME

  overlaps.t
  A piece of code to test the Bio::Align::Overlaps module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl overlaps.t
 prove overlaps.t

=head1 DESCRIPTION

 Test Bio::Align::Overlaps module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 57;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1 to 4

BEGIN {
    use_ok('Bio::Align::Overlaps');
    use_ok('Bio::LocatableSeq');
    use_ok('Bio::SimpleAlign');
}


## First create some sequences

############1---5----10---15---20---25---30---35---40---45---50---55---60---65
my $seq1 = 'ATGCCGCGTGCTGGCAGTTCGGAATCGGACGTAGCACCAGCTTGGTCACGTGGCATCCC------';
my $seq2 = '-----GCGTGCTGGCAGTTCGGATTCGGCCGTAGC------------------------------';
my $seq3 = '---------CCTGGCAGTTCGGATAGGGACGTCGCACCAGCTT---CACGTGGCATGACCCGATA';
my $seq4 = '----CGCGTGCTCCCAGTTCGGATTCGGACGTAGCACCAGCTTGGACACGTGGCATGACC-----';
my $seq5 = '---------------------GATTCGGACGTAGCACCAGCTTGGACACGTGGCATGACCCGATA';
my $seq6 = '----------CTGGCAGTTCGGATACGGACGTAGCACCAGCTTGGGCACCTGGCATGACCCGA--';
my $seq7 = '-------------------------------------TAGCTTGGACACCTCGCATGACCCGC--';
my $seq8 = '---CCGCGTGCTGCCAGTTCGGACTCGGACGTA--------------------------------';
my $seq9 = '------------------------------------TCAGCATGTTTTCGTGGAAAGACCCGA--';
############1---5----10---15---20---25---30---35---40---45---50---55---60---65

my %seqs = (
    seq1 => [ -id => 'seq1', -seq => $seq1, -start => 1,  -end => 59 ],
    seq2 => [ -id => 'seq2', -seq => $seq2, -start => 6,  -end => 35 ],
    seq3 => [ -id => 'seq3', -seq => $seq3, -start => 10, -end => 62 ],
    seq4 => [ -id => 'seq4', -seq => $seq4, -start => 5,  -end => 60 ],
    seq5 => [ -id => 'seq5', -seq => $seq5, -start => 22, -end => 65 ],
    seq6 => [ -id => 'seq6', -seq => $seq6, -start => 11, -end => 63 ],
    seq7 => [ -id => 'seq7', -seq => $seq7, -start => 38, -end => 63 ],
    seq8 => [ -id => 'seq8', -seq => $seq8, -start => 4,  -end => 33 ],
    seq9 => [ -id => 'seq9', -seq => $seq9, -start => 37, -end => 63 ],
    );

my %seqloc = ();
foreach my $seqid (sort keys %seqs) {
    my $seqloc = Bio::LocatableSeq->new(@{$seqs{$seqid}});
    $seqloc{$seqid} = $seqloc;
}

my @metaseqs = values %seqloc;
my $align = Bio::SimpleAlign->new(-seqs => \@metaseqs, -id => 'align1');


## test get_coordinates, TEST 4 to 23

my @seqs = $align->each_seq();
foreach my $seq (@seqs) {
  
    my $mb_id = $seq->id();
    my ($st1, $en1) = Bio::Align::Overlaps::get_coordinates($seq);

    if ($mb_id eq 'seq3') {
	$en1 -= 3;  ## It has an internal gap
    }
    
    is($st1, $seqs{$mb_id}[5], 
       "testing get_coordinates, checking start coord. for seq $mb_id")
	or diag("Looks like this has failed");
       
    is($en1, $seqs{$mb_id}[7], 
       "testing get_coordinates, checking end coord. for seq $mb_id")
	or diag("Looks like this has failed");

}

throws_ok { Bio::Align::Overlaps::get_coordinates() } qr/ERROR: No seq/,
    "TESTING DIE ERROR: when no argument is used with get_coordinates";

throws_ok { Bio::Align::Overlaps::get_coordinates('fk') } qr/ERROR: Object fk/,
    "TESTING DIE ERROR: when argument used for get_coordinates isnt seq. obj.";


## Test calculate_overlaps, TEST

my $mtx = Bio::Align::Overlaps::calculate_overlaps($align);

## Check some entries

## BBAB type, TEST 24 to 27

my %val_seq56 = %{$mtx->entry('seq5', 'seq6')};

is($val_seq56{start}, 22, 
   "testing calculate_overlaps, checking start coord. for ovl type BBAB")
    or diag("Looks like this has failed");

is($val_seq56{end}, 63, 
   "testing calculate_overlaps, checking end coord. for ovl type BBAB")
    or diag("Looks like this has failed");

is($val_seq56{length}, 42, 
   "testing calculate_overlaps, checking length coord. for ovl type BBAB")
    or diag("Looks like this has failed");

is($val_seq56{identity}, 92.8571428571429, 
   "testing calculate_overlaps, checking identity coord. for ovl type BBAB")
    or diag("Looks like this has failed");


## AAAB type, TEST 28 to 31

my %val_seq13 = %{$mtx->entry('seq1', 'seq3')};

is($val_seq13{start}, 10, 
   "testing calculate_overlaps, checking start coord. for ovl type AAAB")
    or diag("Looks like this has failed");

is($val_seq13{end}, 59, 
   "testing calculate_overlaps, checking end coord. for ovl type AAAB")
    or diag("Looks like this has failed");

is($val_seq13{length}, 50, 
   "testing calculate_overlaps, checking length coord. for ovl type AAAB")
    or diag("Looks like this has failed");

is($val_seq13{identity}, 85.1063829787234, 
   "testing calculate_overlaps, checking identity coord. for ovl type AAAB")
    or diag("Looks like this has failed");


## BAAB type, TEST 32 to 35

my %val_seq24 = %{$mtx->entry('seq2', 'seq4')};

is($val_seq24{start}, 6, 
   "testing calculate_overlaps, checking start coord. for ovl type BAAB")
    or diag("Looks like this has failed");

is($val_seq24{end}, 35, 
   "testing calculate_overlaps, checking end coord. for ovl type BAAB")
    or diag("Looks like this has failed");

is($val_seq24{length}, 30, 
   "testing calculate_overlaps, checking length coord. for ovl type BAAB")
    or diag("Looks like this has failed");

is($val_seq24{identity}, 90, 
   "testing calculate_overlaps, checking identity coord. for ovl type BAAB")
    or diag("Looks like this has failed");


## ABAB type, TEST 36 to 39

my %val_seq39 = %{$mtx->entry('seq3', 'seq9')};

is($val_seq39{start}, 37, 
   "testing calculate_overlaps, checking start coord. for ovl type ABAB")
    or diag("Looks like this has failed");

is($val_seq39{end}, 63, 
   "testing calculate_overlaps, checking end coord. for ovl type ABAB")
    or diag("Looks like this has failed");

is($val_seq39{length}, 27, 
   "testing calculate_overlaps, checking length coord. for ovl type ABAB")
    or diag("Looks like this has failed");

is($val_seq39{identity}, 75, 
   "testing calculate_overlaps, checking identity coord. for ovl type ABAB")
    or diag("Looks like this has failed");


## BBBB type, TEST 40 to 43

my %val_seq78 = %{$mtx->entry('seq7', 'seq8')};

foreach my $par (keys %val_seq78) {
    is($val_seq78{$par}, 0, 
       "testing calculate_overlaps, checking $par coord. for ovl type BBBB")
	or diag("Looks like this has failed");
}

## AAAA type, TEST 44 to 47

my %val_seq98 = %{$mtx->entry('seq9', 'seq8')};

foreach my $par (keys %val_seq98) {
    is($val_seq98{$par}, 0, 
       "testing calculate_overlaps, checking $par coord. for ovl type AAAA")
	or diag("Looks like this has failed");
}

## Throws test, TEST 48 and 49

throws_ok { Bio::Align::Overlaps::calculate_overlaps() } qr/ERROR: No align/,
    "TESTING DIE ERROR: when no argument is used with calculate_overlaps";

throws_ok { Bio::Align::Overlaps::calculate_overlaps('fk') } qr/ERROR: fk/,
    "TESTING DIE ERROR: when argument used for calculate_overlaps isnt alnobj.";


## test seed_list, TEST 50 to 57

my @seedlist1 = Bio::Align::Overlaps::seed_list($mtx);
my @seedlist2= Bio::Align::Overlaps::seed_list($mtx, 'length');
my @seedlist3 = Bio::Align::Overlaps::seed_list($mtx, 'identity');

is(join(':', @{$seedlist1[0]}), 'seq3:seq6', 
    "testing seed_list, checking first for ovlscore")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist1[-1]}), 'seq1:seq9', 
    "testing seed_list, checking last for ovlscore")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist2[0]}), 'seq1:seq4', 
    "testing seed_list, checking first for length")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist2[-1]}), 'seq5:seq8', 
    "testing seed_list, checking last for length")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist3[0]}), 'seq4:seq5', 
    "testing seed_list, checking first for identity")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist3[-1]}), 'seq7:seq9', 
    "testing seed_list, checking last for identity")
    or diag("Looks like this has failed");


throws_ok { Bio::Align::Overlaps::seed_list() } qr/ERROR: No argument/,
    "TESTING DIE ERROR: when no argument is used with seed_list";

throws_ok { Bio::Align::Overlaps::seed_list('fk') } qr/ERROR: fk/,
    "TESTING DIE ERROR: when argument used for seed_list isnt mtxobj.";



####
1; #
####
