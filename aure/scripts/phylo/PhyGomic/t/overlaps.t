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
use Test::More tests => 94;
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


## test seed_list, TEST 50 to 64

my @seedlist1 = Bio::Align::Overlaps::seed_list($mtx);
my @seedlist2= Bio::Align::Overlaps::seed_list($mtx, 'length');
my @seedlist3 = Bio::Align::Overlaps::seed_list($mtx, 'identity');

is(join(':', @{$seedlist1[0]}), 'seq3:seq6', 
    "testing seed_list, checking first for ovlscore")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist1[-1]}), 'seq1:seq9', 
    "testing seed_list, checking last for ovlscore")
    or diag("Looks like this has failed");

is(scalar(@seedlist1), 32, 
    "testing seed_list, checking seed number for ovlscore")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist2[0]}), 'seq1:seq4', 
    "testing seed_list, checking first for length")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist2[-1]}), 'seq5:seq8', 
    "testing seed_list, checking last for length")
    or diag("Looks like this has failed");

is(scalar(@seedlist2), 32, 
    "testing seed_list, checking seed number for length")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist3[0]}), 'seq4:seq5', 
    "testing seed_list, checking first for identity")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist3[-1]}), 'seq7:seq9', 
    "testing seed_list, checking last for identity")
    or diag("Looks like this has failed");

is(scalar(@seedlist3), 32, 
    "testing seed_list, checking seed number for identity")
    or diag("Looks like this has failed");

my @seedlist4 = Bio::Align::Overlaps::seed_list($mtx, undef, { length => 45 });

is(join(':', @{$seedlist4[0]}), 'seq3:seq6', 
    "testing seed_list, checking first for ovlscore with length filter")
    or diag("Looks like this has failed");

is(join(':', @{$seedlist4[-1]}), 'seq1:seq3', 
    "testing seed_list, checking last for ovlscore with length filter")
    or diag("Looks like this has failed");

is(scalar(@seedlist4), 6, 
    "testing seed_list, checking seed number for ovlscore with length filter")
    or diag("Looks like this has failed");


throws_ok { Bio::Align::Overlaps::seed_list() } qr/ERROR: No argument/,
    "TESTING DIE ERROR: when no argument is used with seed_list";

throws_ok { Bio::Align::Overlaps::seed_list('fk') } qr/ERROR: fk/,
    "TESTING DIE ERROR: when argument used for seed_list isnt mtxobj.";

throws_ok { Bio::Align::Overlaps::seed_list($mtx, undef, { fk1 => 1}) } qr/fk1/,
    "TESTING DIE ERROR: when filter param. used for seed_list isnt permited.";


## Checking calculate_overseeds, TEST 65 to 79

my $seed1vals = $mtx->entry($seedlist1[0]->[0], $seedlist1[0]->[1]);
my %overseeds = Bio::Align::Overlaps::calculate_overseeds( $align, 
							   $seedlist1[0], 
							   $seed1vals->{start}, 
							   $seed1vals->{end} );

is($overseeds{seq8}->{start}, 11,
    "testing calculate_overseeds, checking start for seq8")
    or diag("Looks like this has failed");

is($overseeds{seq1}->{end}, 59,
    "testing calculate_overseeds, checking end for seq1")
    or diag("Looks like this has failed");

is($overseeds{seq4}->{length}, 50,
    "testing calculate_overseeds, checking length for seq4")
    or diag("Looks like this has failed");

is($overseeds{seq5}->{identity}, 92.5,
    "testing calculate_overseeds, checking identity for seq5")
    or diag("Looks like this has failed");

is($overseeds{seq3}, undef, 
    "testing calculate overseeds, checking seed 1 absent")
    or diag("Looks like this has failed");

is($overseeds{seq6}, undef, 
    "testing calculate overseeds, checking seed 2 absent")
    or diag("Looks like this has failed");

throws_ok { Bio::Align::Overlaps::calculate_overseeds() } qr/ERROR: No align/,
    "TESTING DIE ERROR: when no argument is used with calculate_overseeds";

throws_ok { Bio::Align::Overlaps::calculate_overseeds('fk1') } qr/ERROR: fk1/,
    "TESTING DIE ERROR: when align used with calculate_overseeds isnt align";

throws_ok { Bio::Align::Overlaps::calculate_overseeds($align) } qr/ERROR: No s/,
    "TESTING DIE ERROR: when no seed pair is used with calculate_overseeds";

throws_ok { Bio::Align::Overlaps::calculate_overseeds($align, 'fk2') } qr/fk2/,
    "TESTING DIE ERROR: when seed used with calculate_overseeds isnt aref.";

my $fkaref = [1,2,3];

throws_ok { Bio::Align::Overlaps::calculate_overseeds($align, $fkaref) } qr/se/,
    "TESTING DIE ERROR: when seedaref used doesnt have 3 members";

my @no_start = ($align, $seedlist1[0]);

throws_ok { Bio::Align::Overlaps::calculate_overseeds(@no_start) } qr/No start/,
    "TESTING DIE ERROR: when no start is used with calculate_overseeds()";

my @fk_start = ($align, $seedlist1[0], 'fk3');

throws_ok { Bio::Align::Overlaps::calculate_overseeds(@fk_start) } qr/fk3/,
    "TESTING DIE ERROR: when wrong start is used with calculate_overseeds()";

my @no_end = ($align, $seedlist1[0], 5);

throws_ok { Bio::Align::Overlaps::calculate_overseeds(@no_end) } qr/No end/,
    "TESTING DIE ERROR: when no end is used with calculate_overseeds()";

my @fk_end = ($align, $seedlist1[0], 5, 'fk4');

throws_ok { Bio::Align::Overlaps::calculate_overseeds(@fk_end) } qr/fk4/,
    "TESTING DIE ERROR: when wrong start is used with calculate_overseeds()";


## Checking extension_list, TEST 80 to 94

my @extseed1_ovl = Bio::Align::Overlaps::extension_list(\%overseeds);
my @extseed1_len = Bio::Align::Overlaps::extension_list(\%overseeds,'length');
my @extseed1_ide = Bio::Align::Overlaps::extension_list(\%overseeds,'identity');

is($extseed1_ovl[0], 'seq4', 
    "testing extension_list, checking first member for overscore option")
    or diag("Looks like this has failed");

is($extseed1_ovl[-1], 'seq9', 
    "testing extension_list, checking last member for overscore option")
    or diag("Looks like this has failed");

is(scalar(@extseed1_ovl), 7, 
    "testing extension_list, checking number of member for extension (ovlsco.)")
    or diag("Looks like this has failed");

is($extseed1_len[0], 'seq4', 
    "testing extension_list, checking first member for length option")
    or diag("Looks like this has failed");

is($extseed1_len[-1], 'seq8', 
    "testing extension_list, checking last member for length option")
    or diag("Looks like this has failed");

is(scalar(@extseed1_len), 7, 
    "testing extension_list, checking number of member for extension (length)")
    or diag("Looks like this has failed");

is($extseed1_ide[0], 'seq5', 
    "testing extension_list, checking first member for identity option")
    or diag("Looks like this has failed");

is($extseed1_ide[-1], 'seq9', 
    "testing extension_list, checking last member for identity option")
    or diag("Looks like this has failed");

is(scalar(@extseed1_ide), 7, 
    "testing extension_list, checking number of member for extension (ident.)")
    or diag("Looks like this has failed");

my @extseed1_fil = Bio::Align::Overlaps::extension_list(\%overseeds, 
							undef,
							{ length => 40 },
    );

is($extseed1_fil[0], 'seq4', 
    "testing extension_list, checking first member for overscore and filter")
    or diag("Looks like this has failed");

is($extseed1_fil[-1], 'seq5', 
    "testing extension_list, checking last member for overscore and filter")
    or diag("Looks like this has failed");

is(scalar(@extseed1_fil), 3, 
    "testing extension_list, checking number of member for oversco. and filter")
    or diag("Looks like this has failed");


throws_ok { Bio::Align::Overlaps::extension_list() } qr/ERROR: No argument/,
    "TESTING DIE ERROR: when no argument is used with extension_list";

throws_ok { Bio::Align::Overlaps::extension_list('fk') } qr/ERROR: fk/,
    "TESTING DIE ERROR: when argument used for extension_list isnt hashref..";

throws_ok { Bio::Align::Overlaps::extension_list(\%overseeds, 
						 undef, 
						 { fk1 => 1}) } qr/fk1/,
    "TESTING DIE ERROR: when filter p. used for extension_list isnt permited.";


####
1; #
####
