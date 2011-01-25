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
use Test::More tests => 14;
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

my $phstats0 = PhyGeStats->new();

is(ref($phstats0), 'PhyGeStats', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { PhyGeStats->new(['fake']) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { PhyGeStats->new({ fake => {} }) } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a valid arg.';

throws_ok { PhyGeStats->new({ phygetopo => 'fk2'}) } qr/ERROR: fk2/, 
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

my $phygecluster = PhyGeCluster->new({ acefile => $acefile,
				       seqfile    => $seqfile,
				       strainfile => $strainfile,
				     }
    );

$phygecluster->homologous_search(
    { blast  => [ -p => 'blastn', -d => $blastdbfile, -e => '1e-10', -a => 2],
      strain => 'Sly',
      filter => { hsp_length => ['>', 100], }
    }
    );

my @align0 = ('quiet' => 'yes', 'matrix' => 'BLOSUM');
$phygecluster->run_alignments({program => 'clustalw', parameters => \@align0});
$phygecluster->run_distances({ method => 'Kimura' });
$phygecluster->run_mltrees({ dnaml => {}, outgroup_strain => 'Sly' });

my %seqfams = %{$phygecluster->get_clusters()};
my %strains = %{$phygecluster->get_strains()};


## 2) PhyGeTopo:

my $phygetopo0 = PhyGeTopo->new({ seqfams => \%seqfams, 
				  strains => \%strains });

my %topotypes = $phygetopo0->run_topoanalysis();


## 3) Statistics::R

my $srh0 = Statistics::R->new();



###############
## ACCESSORS ##
###############

## Get/Set_phygetopo, TEST 9 to 11

$phstats0->set_phygetopo($phygetopo0);
my $phygetopo1 = $phstats0->get_phygetopo();

is(ref($phygetopo1), 'PhyGeTopo',
   "Testing Get/Set_phygetopo, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phstats0->set_phygetopo() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_phygetopo function';

throws_ok { $phstats0->set_phygetopo(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied to set_phygetopo isnt PhyGeTopo obj.';


## Get/Set_srh, TEST 12 to 14

$phstats0->set_srh($srh0);
my $srh1 = $phstats0->get_srh();

is(ref($srh1), 'Statistics::R',
   "Testing Get/Set_srh, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phstats0->set_srh() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_srh function';

throws_ok { $phstats0->set_srh(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied to set_srh isnt Statistics::R obj.';

####
1; #
####
