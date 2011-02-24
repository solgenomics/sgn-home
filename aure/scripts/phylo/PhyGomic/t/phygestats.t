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
use Test::More tests => 58;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

use Image::Size;

## TEST 1 to 4

BEGIN {
    use_ok('PhyGeStats');
    use_ok('PhyGeTopo');
    use_ok('PhyGeCluster');
    use_ok('R::YapRI::Base');
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

my $srh0 = R::YapRI::Base->new();



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

is(ref($srh1), 'R::YapRI::Base',
   "Testing Get/Set_rbase, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phstats0->set_rbase() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_rbase function';

throws_ok { $phstats0->set_rbase(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_rbase isnt R::YapRI::Base';


## Get/Set_matrix, TEST 22 to 25

my $mtx0 = R::YapRI::Data::Matrix->new({ name => 'mtx0' });
$phstats0->set_matrix($mtx0);
my $mtx1 = $phstats0->get_matrix();

is(ref($mtx1), 'R::YapRI::Data::Matrix',
   "Testing Get/Set_matrix, Checking object identity")
    or diag("Looks like this has failed");

is($mtx1->get_name(), 'mtx0',
   "Testing Get/Set_matrix, Checking matrix name")
    or diag("Looks like this has failed");


throws_ok { $phstats0->set_matrix() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_matrix function';

throws_ok { $phstats0->set_matrix(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_matrix isnt R::YapRI::Data::Matrix';



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
   "Testing run_block for R::YapRI::Base object with an operation 1:30")
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


## create_matrix, TEST 33 and 34

## phstats0 should be empty

throws_ok { $phstats0->create_matrix() } qr/ERROR: There isnt/, 
    'TESTING DIE ERROR when there isnt any phygt. set before run create_matrix';

$phystats1->create_matrix('topomtx');

is($phystats1->get_matrix()->get_name(), 'topomtx', 
    "testing create_matrix, checking matrix name set into the accessor")
    or diag("Looks like this has failed");




################
## TEST GRAPH ##
################

## create_composition_graph, TEST 35 and 41

my $tempdir = $phystats1->get_rbase()->get_cmddir();
my $graphfile = $tempdir . '/TopoComp.bmp';


$phystats1->create_composition_graph($graphfile);

my ($img_x, $img_y) = Image::Size::imgsize($graphfile);

is($img_x, 800, 
    "testing create_composition_graph, checking image width")
    or diag("Looks like this has failed");

is($img_y, 600, 
    "testing create_composition_graph, checking image height")
    or diag("Looks like this has failed");

throws_ok { $phstats0->create_composition_graph() } qr/ERROR: No filename/, 
    'TESTING DIE ERROR when no filename was supplied create_composition_graph';

throws_ok { $phstats0->create_composition_graph('test', []) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied create_composition_graph isnt HREF';

$phstats0->set_rbase('');
throws_ok { $phstats0->create_composition_graph('test') } qr/ERROR: No rbase/, 
    'TESTING DIE ERROR when rbase is empty before run create_composition_graph';
$phstats0->set_rbase($srh0);

$phstats0->set_matrix('');
throws_ok { $phstats0->create_composition_graph('test') } qr/ERROR: No matrix/, 
    'TESTING DIE ERROR when rbase is empty before run create_composition_graph';
$phstats0->set_matrix($matrix);

throws_ok { $phstats0->create_composition_graph('test', { fk => {} }) } qr/fk/, 
    'TESTING DIE ERROR when no-valid graph arg. used run create_composition_gr';



################
## TEST TABLE ##
################

## create_composition_table, TEST 42 to 46

my $tabfile = $tempdir . '/TopoTable.tab';
$phystats1->create_composition_table($tabfile);

open my $tabfh, '<', $tabfile;

my $match = 0;

my %exptabfile = (
    1 => "\tML\tNJ",
    2 => "cluster_03\ttopology_1\ttopology_2",
    3 => "cluster_09\ttopology_1\ttopology_3",
    );

my %obttabfile = ();
my $l = 0;
while(<$tabfh>) {
    chomp($_);
    $l++;
    $obttabfile{$l} = $_;
}
close($tabfh);

foreach my $idx (sort {$a <=> $b} keys %exptabfile) {
    is($obttabfile{$idx}, $exptabfile{$idx},
	"testing create_composition_table, checking tab file line $idx")
	or diag("Looks like this has failed");
}

throws_ok { $phstats0->create_composition_table() } qr/ERROR: No filename/, 
    'TESTING DIE ERROR when no filename was supplied create_composition_table';

throws_ok { $phstats0->create_composition_table('fake') } qr/ERROR: There isn/, 
    'TESTING DIE ERROR when there isnt any set before create_composition_table';


#######################
## TEST TREE METHODS ##
#######################

## _tree_list, TEST 47 to 49

my %trees = $phystats1->_tree_list();
my %expected_trees = (
    topology_1	=> "(((Nsy:1,Nto:1):1,Nta:1):1,Nta:1,Sly:1)",
    topology_2	=> "((Nsy:1,(Nta:1,Nto:1):1):1,Nta:1,Sly:1)",
    topology_3	=> "((Nsy:1,Nto:1):1,(Nta:1,Nta:1):1,Sly:1)",
    );


foreach my $tpid (sort keys %expected_trees) {
    is($trees{$tpid}, $expected_trees{$tpid},
	"testing _tree_list function, checking trees and keys for $tpid")
	or diag("Looks like this has failed");
}

## create_tree_graph, TEST 50 to 58

my %treefile = (
    vertical   => { filename => $tempdir . '/TopoTreesVertical.bmp',
		    width    => 400,
		    height   => 900,
    },
    horizontal => { filename => $tempdir . '/TopoTreesHorizontal.bmp',
		    width    => 900,
		    height   => 400,
    },
    matrix     => { filename => $tempdir . '/TopoTreesMatrix.bmp',
		    width    => 900,
		    height   => 300,
    },
    );


foreach my $stack (sort keys %treefile) {

    $phystats1->create_tree_graph( $treefile{$stack}->{filename}, 
				   { stack => $stack } 
	);
    my ($w_img, $h_img) = Image::Size::imgsize($treefile{$stack}->{filename});
    
    is($w_img, $treefile{$stack}->{width}, 
       "testing create_tree_graph for $stack stack, checking width")
	or diag("Looks likke this has failed");

    is($h_img, $treefile{$stack}->{height}, 
       "testing create_tree_graph for $stack stack, checking height")
	or diag("Looks likke this has failed");
}


throws_ok { $phstats0->create_tree_graph() } qr/ERROR: No filename/, 
    'TESTING DIE ERROR when no filename was supplied create_tree_graph';

throws_ok { $phstats0->create_tree_graph('test', []) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied create_tree_graph isnt HREF';

$phstats0->set_rbase('');
throws_ok { $phstats0->create_tree_graph('test') } qr/ERROR: No rbase/, 
    'TESTING DIE ERROR when rbase is empty before run create_tree_graph';
$phstats0->set_rbase($srh0);



## Clean the R dirs

$srh0->cleanup();

####
1; #
####
