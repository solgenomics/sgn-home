#!/usr/bin/perl

=head1 NAME

  phygestats.t
  A piece of code to test the PhyGeStats module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl phygestats.t
 prove phygestats.t

=head1 DESCRIPTION

 Test PhyGeStats module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 32;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1 to 4

BEGIN {
    use_ok('PhyGeStats');
    use_ok('PhyGeTopo');
    use_ok('PhyGeCluster');
    use_ok('YapRI::Base');
}


## Create an empty object and test the possible die functions. TEST 5 to 8

my $phstats0 = PhyGeStats->new({ rbase => '' });

is(ref($phstats0), 'PhyGeStats', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { PhyGeStats->new(['fake']) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { PhyGeStats->new({ fake => {} }) } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a valid arg.';

throws_ok { PhyGeStats->new({ rbase => 'fk2'}) } qr/ERROR: fk2/, 
    'TESTING DIE ERROR for new() when arg. doesnt have valid value';


#####################
## PREPARING FILES ##
#####################

## 0) Test files:

my $blastfile = "$FindBin::Bin/testfiles/selfblast.test.m8";
my $seqfile = "$FindBin::Bin/testfiles/seq.test.fasta";
my $strainfile = "$FindBin::Bin/testfiles/strains.test.tab";
my $acefile = "$FindBin::Bin/testfiles/assembly_out.test.ace";
my $blastdbfile = "$FindBin::Bin/testfiles/blastref.test.fasta";


## 1) PhyGeCluster:
## It will create two phygecluster for two methods to be able to compare them

my $phygecluster1 = PhyGeCluster->new({ blastfile => $blastfile, 
                                        fastblast_parser => 1,  });
$phygecluster1->load_seqfile({ sequencefile => $seqfile });
$phygecluster1->load_strainfile({ strainfile => $strainfile });

my @align0 = ('quiet' => 'yes', 'matrix' => 'BLOSUM');
$phygecluster1->run_alignments({program => 'clustalw', parameters => \@align0});

$phygecluster1->run_distances({ method => 'Kimura' });

my ($rm_clusters_href, $rm_members_href) = $phygecluster1->prune_by_strains({
    composition  => { 'Sly' => 1, 'Nsy' => 1, 'Nto' => 1, 'Nta' => 2 },
    min_distance => [ 
        [ 'Nta', 'Nta' ],
        [ 'Nta', 'Nsy' ], 
        [ 'Nta', 'Nto' ], 
        [ 'Nta', 'Sly' ],
        [ 'Nsy', 'Sly' ],
        [ 'Nto', 'Sly' ],                     
        ],
									    });

$phygecluster1->run_alignments({program => 'clustalw', parameters => \@align0});

$phygecluster1->run_distances({ method => 'Kimura' });

my $phygecluster2 = $phygecluster1->clone();

$phygecluster1->run_mltrees({ dnaml => {}, outgroup_strain => 'Sly' });

my %seqfams1 = %{$phygecluster1->get_clusters()};
my %strains1 = %{$phygecluster1->get_strains()};

$phygecluster2->run_njtrees({ quiet => 1, outgroup_strain => 'Sly' });

my %seqfams2 = %{$phygecluster2->get_clusters()};
my %strains2 = %{$phygecluster2->get_strains()};


## 2) PhyGeTopo:

my $phygetopo1 = PhyGeTopo->new({ seqfams => \%seqfams1, 
				  strains => \%strains1 });

my %topotypes1 = $phygetopo1->run_topoanalysis();

my $phygetopo2 = PhyGeTopo->new({ seqfams => \%seqfams2, 
				  strains => \%strains2 });

my %topotypes2 = $phygetopo2->run_topoanalysis();

## 3) Statistics::R

my $rdir = "$FindBin::Bin/..";
my $srh0 = YapRI::Base->new();



###############
## ACCESSORS ##
###############

## Get/Set_phygetopo, TEST 9 to 18

$phstats0->set_phygetopo({'ML' => $phygetopo1});
my %phygetopo1_r = %{$phstats0->get_phygetopo()};

is(ref($phygetopo1_r{'ML'}), 'PhyGeTopo',
   "Testing Get/Set_phygetopo, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phstats0->set_phygetopo() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_phygetopo function';

throws_ok { $phstats0->set_phygetopo('fake') } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to set_phygetopo isnt HASH obj.';

throws_ok { $phstats0->set_phygetopo({'name' => 'fake'}) } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to set_phygetopo isnt PhyGeTopo obj.';

$phstats0->add_phygetopo('NJ', $phygetopo2);
my %phygetopo_2r = %{$phstats0->get_phygetopo()};

is(scalar(keys %phygetopo_2r), 2, 
    "Testing Add_phygetopo, checking number of elements")
    or diag("Looks like this has failed");

throws_ok { $phstats0->add_phygetopo() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to add_phygetopo function';

throws_ok { $phstats0->add_phygetopo('fake') } qr/ERROR: No phygetopo/, 
    'TESTING DIE ERROR when no phygetopo arg. wassupplied to add_phygetopo';

throws_ok { $phstats0->add_phygetopo('name', 'fake') } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to add_phygetopo isnt PhyGeTopo obj.';

my %phygetopo_3r = %{$phstats0->delete_phygetopo()};

is(scalar(keys %phygetopo_3r), 2, 
    "Testing delete_phygetopo, checking number of elements deleted")
    or diag("Looks like this has failed");

my %phygetopo_4r = %{$phstats0->get_phygetopo()};

is(scalar(keys %phygetopo_4r), 0, 
    "Testing delete_phygetopo, checking empty hash")
    or diag("Looks like this has failed");


## Get/Set_rbase, TEST 19 to 21

$phstats0->set_rbase($srh0);
my $srh1 = $phstats0->get_rbase();

is(ref($srh1), 'YapRI::Base',
   "Testing Get/Set_rbase, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phstats0->set_rbase() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_rbase function';

throws_ok { $phstats0->set_rbase(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_rbase isnt YapRI::Base';


## Get/Set_matrix, TEST 22 to 25

my $mtx0 = YapRI::Data::Matrix->new({ name => 'mtx0' });
$phstats0->set_matrix($mtx0);
my $mtx1 = $phstats0->get_matrix();

is(ref($mtx1), 'YapRI::Data::Matrix',
   "Testing Get/Set_matrix, Checking object identity")
    or diag("Looks like this has failed");

is($mtx1->get_name(), 'mtx0',
   "Testing Get/Set_matrix, Checking matrix name")
    or diag("Looks like this has failed");


throws_ok { $phstats0->set_matrix() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_matrix function';

throws_ok { $phstats0->set_matrix(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_matrix isnt YapRI::Data::Matrix';



##################
## R connection ##
##################

## Test 26

my $phystats1 = PhyGeStats->new();
my $srh2 = $phystats1->get_rbase();

$srh2->create_block('PHYGEBLOCK1');
$srh2->add_command('print(1:30)', 'PHYGEBLOCK1');
$srh2->run_block('PHYGEBLOCK1');
my $result2 = $srh2->get_resultfiles('PHYGEBLOCK1');

my @results = ();
open my $rfh2, '<', $result2;
while(<$rfh2>) {
    chomp($_);
    $_ =~ s/^\s*\[\d+\]\s+//;
    push @results, split(/\s+/, $_);
}
close($rfh2);

my @exp_results = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
		   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		   21, 22, 23, 24, 25, 26, 27, 28, 29, 30 );

is(join(',', @results), join(',', @exp_results), 
   "Testing run_block for YapRI::Base object with an operation 1:30")
    or diag("Looks like this has failed");



#########################
## ANALYTICAL FUNCTION ##
#########################

$phystats1->set_phygetopo({ 'NJ' => $phygetopo2, 'ML' => $phygetopo1 });


## Test _compare_phygetopos, TEST 27

my %phycomp = $phystats1->_compare_phygetopos();

is(scalar(keys %phycomp), 3, 
    "testing _compare_phygetopos, checking number of new topologies")
    or diag("Looks like this has failed");


## Test _phygt2matrix, TEST 28 to 32

my $matrix = $phystats1->_phygt2matrix();

my @rownames = @{$matrix->get_rownames()};
my @colnames = @{$matrix->get_colnames()};
my @data = @{$matrix->get_data()};

is($matrix->get_coln(), 2, 
    "testing _phygt2matrix, checking column number (2)")
    or diag("Looks like this has failed");

is($matrix->get_rown(), 3, 
    "testing _phygt2matrix, checking row number (3)")
    or diag("Looks like this has failed");

is(join(',', @colnames), "ML,NJ", 
    "testing _phygt2matrix, checking colnames")
    or diag("Looks like this has failed");

is(join(',', @rownames), "topology_1,topology_2,topology_3",
    "testing _phygt2matrix, checking rownames")
    or diag("Looks like this has failed");

is(scalar(@data), 6, 
    "testing _phygt2matrix, checking data number (6)")
    or diag("Looks like this has failed");


################
## TEST GRAPH ##
################



## Clean the R dirs

#$srh0->cleanup();

####
1; #
####
