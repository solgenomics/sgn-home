#!/usr/bin/perl

=head1 NAME

  strain_composition.t
  A piece of code to test the Strain::Composition module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl strain_composition.t
 prove strain_composition.t

=head1 DESCRIPTION

 Test Strain::Composition module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 47;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1

BEGIN {
    use_ok('Strain::Composition');
}

## Test a new empty object, TEST 2 to 5

my $comp = Strain::Composition->new();

is(ref($comp), 'Strain::Composition',
    "testing new(), checking object identity for an empty object")
    or diag("Looks like this has failed");

throws_ok { Strain::Composition->new([]) } qr/ERROR: ARRAY/,
    "TESTING DIE ERROR: when argument used for new() isnt hashref.";

throws_ok { Strain::Composition->new({ fake => 1 }) } qr/ERROR: fake/,
    "TESTING DIE ERROR: when key-argument used for new() isnt permitted.";

throws_ok { Strain::Composition->new({ strains => 'fk'}) } qr/ERROR: fk/,
    "TESTING DIE ERROR: when value-argument used for new() isnt permitted";

## Test accessors

my %strains = ( 'seq1' => 'str1', 
		'seq2' => 'str1',
		'seq3' => 'str1', 
		'seq4' => 'str2',
		'seq5' => 'str2', 
		'seq6' => 'str3',
		'seq7' => 'str3', 
		'seq8' => 'str4',
		'seq9' => 'str4', 
    );

## Test get/set_strains, TEST 6 to 9

$comp->set_strains(\%strains);
my %getstrains = %{$comp->get_strains()};

is(join(',', sort keys(%getstrains)), join(',', sort keys(%strains)), 
    "testing get/set_strains, checking strains keys")
    or diag("Looks like this has failed");

is(join(',', sort values(%getstrains)), join(',', sort values(%strains)), 
    "testing get/set_strains, checking strains values")
    or diag("Looks like this has failed");

throws_ok { $comp->set_strains() } qr/ERROR: No arg./,
    "TESTING DIE ERROR: when no argument was used for set_strains";

throws_ok { $comp->set_strains([]) } qr/ERROR: When/,
    "TESTING DIE ERROR: when argument used for set_strains isnt hashref.";


## Test get/set_composition, TEST 10 to 14

my %composition = ( 'str1' => 2, 'str2' => 1, 'str3' => 1, 'str4' => 1);

$comp->set_composition(\%composition);
my %getcomp = %{$comp->get_composition()};

is(join(',', sort keys(%getcomp)), join(',', sort keys(%composition)),
   "testing get/set_composition, checking composition keys")
    or diag("Looks like this has failed");

is(join(',', sort values(%getcomp)), join(',', sort values(%composition)),
   "testing get/set_composition, checking composition values")
    or diag("Looks like this has failed");

throws_ok { $comp->set_composition() } qr/ERROR: No arg./,
    "TESTING DIE ERROR: when no argument was used for set_composition";

throws_ok { $comp->set_composition([]) } qr/ERROR: When/,
    "TESTING DIE ERROR: when argument used for set_composition isnt hashref.";

throws_ok { $comp->set_composition({ 'str1' => 'fk'}) } qr/ERROR: value for/,
    "TESTING DIE ERROR: when value for hash isnt integer for set_composition";


## Test get/set_members, TEST 15 to 24

my %members = ( 'str1' => ['seq1', 'seq3'],
		'str2' => ['seq4'],
		'str3' => ['seq7'],
		'str4' => ['seq8'],
    );

$comp->set_members(\%members);
my %getmemb = %{$comp->get_members()};

is(join(',', sort keys(%getmemb)), join(',', sort keys(%members)),
   "testing get/set_members, checking members keys")
    or diag("Looks like this has failed");

my @getmemb = ();
my @expmemb = ();

foreach my $memb (keys %members) {
    push @expmemb, @{$members{$memb}};
    push @getmemb, @{$getmemb{$memb}};
}

is(join(',', sort @getmemb), join(',', sort @expmemb), 
    "testing get/set_members, checking member values aref.")
    or diag("Looks like this has failed");

throws_ok { $comp->set_members() } qr/ERROR: No arg./,
    "TESTING DIE ERROR: when no argument was used for set_members";

throws_ok { $comp->set_members([]) } qr/ERROR: When/,
    "TESTING DIE ERROR: when argument used for set_members isnt hashref.";

throws_ok { $comp->set_members({ 'fk1' => 1 }) } qr/ERROR: fk1 used/,
    "TESTING DIE ERROR: when strain used for set_members doesnt exist at comp.";

throws_ok { $comp->set_members({ 'str1' => 'fk2' }) } qr/ERROR: Value for/,
    "TESTING DIE ERROR: when value for strain used for set_members isnt AREF.";

throws_ok { $comp->set_members({ 'str2' => ['seq1', 'seq2'] }) } qr/ERROR: Mor/,
    "TESTING DIE ERROR: when more than comp. count members are used.";

throws_ok { $comp->set_members({ 'str2' => ['seq1', 'seq2'] }) } qr/ERROR: Mor/,
    "TESTING DIE ERROR: when more than comp. count members are used.";

throws_ok { $comp->set_members({ 'str2' => ['seqx'] }) } qr/ERROR: seqx/,
    "TESTING DIE ERROR: when member doesnt exists at strains";

throws_ok { $comp->set_members({ 'str2' => ['seq1'] }) } qr/ERROR: Strain/,
    "TESTING DIE ERROR: when strain used for member isnt same than strain acc.";


## Check is_complete, TEST 25 to 30

is($comp->is_complete(), 1, 
    "testing is_complete(), checking boolean true for all the strains")
    or diag("Looks like this has failed");

is($comp->is_complete('str1'), 1, 
    "testing is_complete(), checking boolean true for a specific strain")
    or diag("Looks like this has failed");

my %members2 = ( 'str1' => ['seq3'], 
		'str2' => ['seq4'],
		'str3' => ['seq7'],
		'str4' => ['seq8'],
    );

$comp->set_members(\%members2);

is($comp->is_complete(), 0, 
    "testing is_complete(), checking boolean false for all the strains")
    or diag("Looks like this has failed");

is($comp->is_complete('str1'), 0, 
    "testing is_complete(), checking boolean false for a specific strain")
    or diag("Looks like this has failed");

is($comp->is_complete('str2'), 1, 
    "testing is_complete(), checking boolean true for a specific strain 2")
    or diag("Looks like this has failed");

throws_ok { $comp->is_complete('strx') } qr/ERROR: Strain strx/,
    "TESTING DIE ERROR: when strain used for is_complete isnt in composition";


## Test add_member, TEST 31 to 36

is($comp->add_member('seq3'), undef, 
   "testing add_member(), checking that can not be added a member that exist")
    or diag("Looks like this has failed");

is($comp->add_member('seq2'), 'str1', 
   "testing add_member(), checking returning of the strain when it success")
    or diag("Looks like this has failed");

my %getmemb2 = %{$comp->get_members()};

is(join(',', sort @{$getmemb2{str1}}), 'seq2,seq3',
   "testing add_member(), checking member ids")
    or diag("Looks like this has failed");

is($comp->add_member('seq1'), undef, 
    "testing add_member(), checking undef when no more members can be added")
    or diag("Looks like this has failed");

my %getmemb3 = %{$comp->get_members()};

is(join(',', sort @{$getmemb3{str1}}), 'seq2,seq3',
   "testing add_member(), checking member ids has not changed")
    or diag("Looks like this has failed");

throws_ok { $comp->add_member('seqx') } qr/ERROR: seqx/,
    "TESTING DIE ERROR: when member used doesnt exist into the strain acc.";


## Check count_members, TEST 37 to 40

my %count_members = $comp->count_members();

foreach my $str1 (sort keys %count_members) {
    is($count_members{$str1}, $composition{$str1}, 
	"testing count_members(), checking members count for $str1")
	or diag("Looks like this has failed");
}

## Check total_members, TEST 41

is($comp->total_members(), 5, 
    "testing total_members(), checking total number of members")
    or diag("Looks like this has failed");


## Check delete_members, TEST 42 to 47

my @delmembers = $comp->delete_members('str2');

is($comp->total_members(), 4, 
    "testing delete_members(), checking total number of members after deletion")
    or diag("Looks like this has failed");

is(join(',', @delmembers), 'seq4', 
   "testing delete_members(), checking member deleted using strain")
    or diag("Looks like this has failed");

my %getmemb4 = %{$comp->get_members()};

is($getmemb4{str2}, undef, 
    "testing delete_members(), checking that the strain has been deleted")
    or diag("Looks like this has failed");

my %delmembers = $comp->delete_members();

is($comp->total_members(), 0, 
    "testing delete_members(), checking number of members for complete del.")
    or diag("Looks like this has failed");

is(join(',', sort keys %delmembers), 'str1,str3,str4', 
   "testing delete_members(), checking strain deleted using total deletion")
    or diag("Looks like this has failed");

throws_ok { $comp->delete_members('strx') } qr/ERROR: Strain strx/,
    "TESTING DIE ERROR: when strain used doesnt exist into the composition acc";



####
1; #
####
