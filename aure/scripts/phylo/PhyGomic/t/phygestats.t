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
use Test::More tests => 23;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1 to 4

BEGIN {
    use_ok('PhyGeStats');
    use_ok('PhyGeTopo');
    use_ok('PhyGeCluster');
    use_ok('Statistics::R');
}


## Create an empty object and test the possible die functions. TEST 5 to 8

my $phstats0 = PhyGeStats->new({ r_connection => 0 });

is(ref($phstats0), 'PhyGeStats', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { PhyGeStats->new(['fake']) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { PhyGeStats->new({ fake => {} }) } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a valid arg.';

throws_ok { PhyGeStats->new({ r_connection => 'fk2'}) } qr/ERROR: fk2/, 
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
 

my $phygecluster1 = PhyGeCluster->new({ acefile => $acefile,
		 		        seqfile    => $seqfile,
				        strainfile => $strainfile,
				     }
    );

$phygecluster1->homologous_search(
    { blast  => [ -p => 'blastn', -d => $blastdbfile, -e => '1e-10', -a => 2],
      strain => 'Sly',
      filter => { hsp_length => ['>', 100], }
    }
    );

my @align0 = ('quiet' => 'yes', 'matrix' => 'BLOSUM');
$phygecluster1->run_alignments({program => 'clustalw', parameters => \@align0});

my $phygecluster2 = $phygecluster1->clone();

$phygecluster1->run_distances({ method => 'Kimura' });
$phygecluster1->run_mltrees({ dnaml => {}, outgroup_strain => 'Sly' });

my %seqfams1 = %{$phygecluster1->get_clusters()};
my %strains1 = %{$phygecluster1->get_strains()};

$phygecluster2->run_distances({ method => 'Kimura' });
$phygecluster2->run_njtrees({ quiet => 1 });

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

my $srh0 = Statistics::R->new();



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


## Get/Set_srh, TEST 19 to 21

$phstats0->set_r_connection($srh0);
my $srh1 = $phstats0->get_r_connection();

is(ref($srh1), 'Statistics::R',
   "Testing Get/Set_r_connection, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phstats0->set_r_connection() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_r_connection function';

throws_ok { $phstats0->set_r_connection(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_r_connection isnt Statistics::R';


##################
## R connection ##
##################

## Test 22 and 23

my $phystats1 = PhyGeStats->new();
my $srh2 = $phystats1->get_r_connection();

is($srh2->is_started(), 1, 
    "Testing startR over Statistics::R object when a new object is created")
    or diag("Looks like this has failed");

$srh2->send('print(1:30)');
my $result2 = $srh2->read();

## Parse the results removing different lines

$result2 =~ s/\s?\[\d+\]\s/ /g;
$result2 =~ s/\n/ /g;
$result2 =~ s/^\s+//g;
$result2 =~ s/\s+/ /g;

my @results = split(/ /, $result2);

my @exp_results = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
		   11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
		   21, 22, 23, 24, 25, 26, 27, 28, 29, 30 );

is(join(',', @results), join(',', @exp_results), 
   "Testing startR over Statistics::R object with an operation 1:30")
    or diag("Looks like this has failed");












## END THE R COMUNICATION ##

$srh2->stopR();

####
1; #
####
