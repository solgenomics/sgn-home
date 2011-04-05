#!/usr/bin/perl

=head1 NAME

  strain_alleleidentification.t
  A piece of code to test the Strain::AlleleIdentification module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl strain_alleleidentification.t
 prove strain_alleleidentification.t

=head1 DESCRIPTION

 Test Strain::AlleleIdentification module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 27;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1

BEGIN {
    use_ok('Strain::AlleleIdentification', qw( identify_alleles ));
    use_ok('Bio::LocatableSeq');
    use_ok('Bio::SimpleAlign');
    use_ok('Bio::Tree::Node');
    use_ok('Bio::Tree::Tree');
}

## Prepare the data

## Alignment

my %seqs = (
    sq1 => 'ATGCGCCATGACCACGATGGGCAGGATAAATCGGGATAGATAGTAGCCAGCATGACGTGGCAGTAG',
    sq2 => 'ATGCGCCACGACCACGATGGGCAGGATATATCGGGATAGACCGTAG-CAGCATGACGCGGCAGTAG',
    sq3 => 'ATGCCGGATGACCAGGATGGGCGGGATAAATCGAGATAGTTACTACCCAGCATGACCTCCCAGTAG',
    sq4 => 'ATGCCGGATGACCAGGATGCGCGGGATAAATCGAGATAGTTACAAC-CAGCATGACCACCCAGTAG',
    sq5 => 'ATGCCCGATGACCAGGATCCGCGGGATATATCGAGATAGATACAAC-CAGCATGACGAGGCAGTAG',
    sq6 => 'ATGCCGGACGACCAGGATGCGCGGGATAAAACGAGATAGTTACAACCCAGCAAGACCACCCAGTAG',
    );

my @locseqs = ();
foreach my $seqname (sort keys %seqs) {
    my $end = 66;
    if ($seqs{$seqname} =~ m/-/) {
	$end--;
    }

    push @locseqs, Bio::LocatableSeq->new( -id    => $seqname, 
					   -seq   => $seqs{$seqname},
					   -start => 1,
					   -end   => $end,
	);
}

my $align = Bio::SimpleAlign->new(-seqs => \@locseqs);

## Tree                              +----------- sq1 (A)
##                                   |
## +---------------------------------+ 
## |                                 |
## +                                 +---------------- sq2 (T)
## |   +----------- sq5 (C)
## |   |
## +---+                 +------------ sq3 (T) 
##     |                 |
##     +-----------------+        +--- sq4 (B)
##                       |        |
##                       +--------+
##                                |
##                                +--------------- sq6 (T)

my $node1 = Bio::Tree::Node->new(-id => 'sq1', -branch_length => 0.02413);
my $node2 = Bio::Tree::Node->new(-id => 'sq2', -branch_length => 0.05280);
my $node3 = Bio::Tree::Node->new(-id => 'sq3', -branch_length => 0.02552);
my $node4 = Bio::Tree::Node->new(-id => 'sq4', -branch_length => 0.00350);
my $node5 = Bio::Tree::Node->new(-id => 'sq5', -branch_length => 0.04289);
my $node6 = Bio::Tree::Node->new(-id => 'sq6', -branch_length => 0.04266);

my $node7 = Bio::Tree::Node->new( -descendents   => [$node1, $node2], 
				  -branch_length => 0.08957 );
my $node8 = Bio::Tree::Node->new( -descendents   => [$node4, $node6], 
				  -branch_length => 0.01993 );
my $node9 = Bio::Tree::Node->new( -descendents   => [$node3, $node8], 
				  -branch_length => 0.05839 );
my $node10 = Bio::Tree::Node->new( -descendents   => [$node5, $node9], 
	 			   -branch_length => 0.02139 );

my $root = Bio::Tree::Node->new( -descendents => [$node7, $node10]);

my $tree = Bio::Tree::Tree->new( -root => $root );


## Strains

my %strains = ( sq1 => 'A', 
		sq2 => 'T', 
		sq3 => 'T', 
		sq4 => 'B', 
		sq5 => 'C', 
		sq6 => 'T');


## Test the function

my %allele_args = ( strains   => \%strains,
		    alignment => $align,
		    tree      => $tree,
		    target    => 'T',
		    parents   => { 'A' => 'T-a', 'B' => 'T-b'},
    );

my %alleles = identify_alleles(\%allele_args);


## Test

is(scalar(keys %alleles), 3,
    "testing identify_alleles, checking number of alleles")
    or diag("Looks like this has failed");

is($alleles{sq1}, undef,
    "testing identify_alleles, checking allele for one of the parents(undef)")
    or diag("Looks like this has failed");

is($alleles{sq2}, 'T-a',
    "testing identify_alleles, checking allele for one of the targets (T-a)")
    or diag("Looks like this has failed");

is($alleles{sq3}, 'T-b',
    "testing identify_alleles, checking allele for one of the targets (T-b)")
    or diag("Looks like this has failed");

is($alleles{sq4}, undef,
    "testing identify_alleles, checking allele for the other parent (undef)")
    or diag("Looks like this has failed");

is($alleles{sq5}, undef,
    "testing identify_alleles, checking allele for one of the targets (undef)")
    or diag("Looks like this has failed");

is($alleles{sq6}, 'T-b',
    "testing identify_alleles, checking allele for one of the targets (T-b)")
    or diag("Looks like this has failed");

throws_ok { identify_alleles() } qr/ERROR: No arguments/,
    "TESTING DIE ERROR: when no argument are used for identify_alleles()";

throws_ok { identify_alleles([]) } qr/ERROR: Argument supplied/,
    "TESTING DIE ERROR: when arg. used for identify_alleles() isnt hashref";

throws_ok { identify_alleles({ fake => 1}) } qr/ERROR: fake is a non-permited/,
    "TESTING DIE ERROR: when arg. used for identify_alleles() isnt permited";

throws_ok { identify_alleles({ strains => 1})} qr/ERROR: strains doesnt have/,
    "TESTING DIE ERROR: when arg. used doesnt have right value";

throws_ok { identify_alleles({})} qr/ERROR: Mandatory argument/,
    "TESTING DIE ERROR: when mandatory argument is not used";

my %ff_allele_args = %allele_args;
$ff_allele_args{parents} = { 'A' => 'T-a'};

throws_ok { identify_alleles(\%ff_allele_args)} qr/ERROR: Less than one parent/,
    "TESTING DIE ERROR: when less tha one parent was used";

$node2->id('seq2');
$strains{seq2} = 'T';
$allele_args{filter} = {identity => 60};

throws_ok { identify_alleles(\%allele_args)} qr/ERROR: No seq was found/,
    "TESTING DIE ERROR: when target seqids in the tree are not in the align";

$node2->id('sq2');

$node1->id('seq1');
$strains{seq1} = 'A';

throws_ok { identify_alleles(\%allele_args)} qr/ERROR: No seq was found/,
    "TESTING DIE ERROR: when parents seqids in the tree are not in the align";

$node1->id('sq1');


## Test the function using filters

## Test the function

my %allele_args2 = ( strains   => \%strains,
		     alignment => $align,
		     tree      => $tree,
		     target    => 'T',
		     parents   => { 'A' => 'T-a', 'B' => 'T-b'},
		     filter    => { identity => 95},
    );

my %alleles2 = identify_alleles(\%allele_args2);


## Test

is(scalar(keys %alleles2), 2,
    "testing identify_alleles with filter, checking number of alleles")
    or diag("Looks like this has failed");

is($alleles2{sq1}, undef,
    "testing identify_alleles with filter, checking one of the parents(undef)")
    or diag("Looks like this has failed");

is($alleles2{sq2}, undef,
    "testing identify_alleles with filter, checking one of the targets (T-a)")
    or diag("Looks like this has failed");

is($alleles2{sq3}, 'T-b',
    "testing identify_alleles with filter, checking one of the targets (T-b)")
    or diag("Looks like this has failed");

is($alleles2{sq4}, undef,
    "testing identify_alleles with filter, checkingv the other parent (undef)")
    or diag("Looks like this has failed");

is($alleles2{sq5}, undef,
    "testing identify_alleles with filter, checking one of the targets (undef)")
    or diag("Looks like this has failed");

is($alleles2{sq6}, 'T-b',
    "testing identify_alleles with filter, checking one of the targets (T-b)")
    or diag("Looks like this has failed");



####
1; #
####
