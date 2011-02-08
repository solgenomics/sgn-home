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
use Test::More tests => 43;
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


## Get/Set_srh, TEST 19 to 21

$phstats0->set_rbase($srh0);
my $srh1 = $phstats0->get_rbase();

is(ref($srh1), 'YapRI::Base',
   "Testing Get/Set_rbase, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phstats0->set_rbase() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_rbase function';

throws_ok { $phstats0->set_rbase(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_rbase isnt Statistics::R';


##################
## R connection ##
##################

## Test 22

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

my $file0 = $phystats1->_r_infile_tm();

## It should not create any file because the phygetopo is empty, TEST 23

is(-s $file0, undef, 
    "Testing _r_infile_topomemb internal function when no file is created")
    or diag("Looks like this has failed");

## Testing the file creation (_r_infile_tm), TEST 24 to 28

$phystats1->set_phygetopo({ 'NJ' => $phygetopo2, 'ML' => $phygetopo1 });
my $file1 = $phystats1->_r_infile_tm();

my $wrong_coln = 0;
my $wrong_headers = 0;
my $wrong_rownames = 0;
my $wrong_topoformat = 0;
my $wrong_data_format = 0;

open my $fh, '<', $file1;
my $l = 0;
while (<$fh>) {
    my @data = split(/\t/, $_);
    if (scalar(@data) != 4) {
	$wrong_coln++;
    }

    ## Check headers for the first line
    if ($l == 0) {
	unless ($_ =~ m/\t"\w+"\t"\w+"\t"Topology"/) {
	    $wrong_headers++;
	}
    }
    else {
	my $row_names = shift(@data);
	my $topologies = pop(@data);
	
	unless ($row_names =~ m/"\w+"/) {
	    $wrong_rownames++;
	}
	unless ($topologies =~ m/"\(.+\)"/) {
	    $wrong_topoformat++;
	}
	foreach my $data (@data) {
	    unless ($data =~ m/^\d+$/) {
		$wrong_data_format++;
	    }
	}
    }
    $l++;
}

is($wrong_coln, 0, 
    "Testing _r_infile_topomemb internal function, checking col. number")
    or diag("Looks like this has failed");

is($wrong_headers, 0, 
    "Testing _r_infile_topomemb internal function, checking headers")
    or diag("Looks like this has failed");

is($wrong_rownames, 0, 
    "Testing _r_infile_topomemb internal function, checking row names")
    or diag("Looks like this has failed");

is($wrong_topoformat, 0, 
    "Testing _r_infile_topomemb internal function, checking topology format")
    or diag("Looks like this has failed");

is($wrong_data_format, 0, 
    "Testing _r_infile_topomemb internal function, checking data format")
    or diag("Looks like this has failed");

## Check the die functions for _r_infile_tm() function, TEST 29 to 32

throws_ok { $phstats0->_r_infile_tm([]); } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied _r_infile_tm isnt HASH REF.';

throws_ok { $phstats0->_r_infile_tm({ 'fake' => 1}); } qr/ERROR: fake/, 
    'TESTING DIE ERROR when key arg supplied _r_infile_tm isnt valid';

throws_ok { $phstats0->_r_infile_tm({ tempfile => 2}); } qr/ERROR: 2/, 
    'TESTING DIE ERROR when value arg supplied _r_infile_tm isnt valid';

throws_ok { $phstats0->_r_infile_tm({ phygetopo => 'HASH'}); } qr/ERROR: phyg/, 
    'TESTING DIE ERROR when phygetopo arg supplied _r_infile_tm isnt HASHREF';


## Checking the _r_loadfile function, TEST 33 to 35

my $r_obj_name = $phystats1->_r_loadfile($file1, 'TopoMembTable');


## Check that the object contains the data

$srh2->create_block('CHECK_' . $r_obj_name, $r_obj_name);
$srh2->add_command("print($r_obj_name)", 'CHECK_' . $r_obj_name);
$srh2->run_block('CHECK_' . $r_obj_name);
my $rfile3 = $srh2->get_resultfiles('CHECK_' . $r_obj_name);

my $r_col_n = 0;
my $r_row_n = 0;
my $r_data_n = 0;


open my $rfh3, '<', $rfile3;
while(<$rfh3>) {
    chomp($_);
    $_ =~ s/^\s+//;
    $_ =~ s/\s+$//;
    $_ =~ s/\s+/ /g;
    
    $r_row_n++;
    my @data = split(/\s+/, $_);
    $r_col_n = scalar(@data);
    $r_data_n += $r_col_n;
}


is($r_row_n, 4, 
    "Testing _r_loadfile internal function, checking row number (4)")
    or diag("Looks like this has failed");

is($r_col_n, 4, 
    "Testing _r_loadfile internal function, checking col number (4)")
    or diag("Looks like this has failed");

is($r_data_n, 15, 
    "Testing _r_loadfile internal function, checking data number (4x4)")
    or diag("Looks like this has failed");

## Check if it dies properly, TEST 36 and 37

throws_ok { $phstats0->_r_loadfile(); } qr/ERROR: No file/, 
    'TESTING DIE ERROR when no file arg was supplied _r_loadfile()';

throws_ok { $phstats0->_r_loadfile('test'); } qr/ERROR: No R basename/, 
    'TESTING DIE ERROR when no basename arg was supplied _r_loadfile()';


################
## TEST GRAPH ##
################

## Test _initR_grDevices, TEST 38 to 43

my ($grfile, $grblock) = $phystats1->_initR_grDevices();

is($grfile =~ m/rGraph/, 1, 
    "testing _initR_grDevices, checking filename")
    or diag("Looks like this has failed");

$srh2->create_block('TEST_' . $grblock, $grblock);
$srh2->add_command('dev.list()', 'TEST_' . $grblock);
$srh2->run_block('TEST_' . $grblock);

my $gr_checkfile1 = $srh2->get_resultfiles('TEST_' . $grblock);

my $init_bmp = 0;
open my $grfh1, '<', $gr_checkfile1;
while (<$grfh1>) {
    chomp($_);
    if ($_ =~ m/bmp/) {
	$init_bmp = 1;
    }
}
close($grfh1);

is($init_bmp, 1, 
    "testing _initR_grDevices, checking dev.list R command over a new block")
    or diag("Looks like this has failed");

throws_ok { $phstats0->_initR_grDevices('fake'); } qr/ERROR:fake isnt/, 
    'TESTING DIE ERROR when no valid device is used for _initR_grDevices()';

throws_ok { $phstats0->_initR_grDevices(undef, { fk => 1}); } qr/ERROR: fk/, 
    'TESTING DIE ERROR when no valid grarg key is used for _initR_grDevices()';

throws_ok { $phstats0->_initR_grDevices(undef,{width => 'big'}); } qr/OR: w/, 
    'TESTING DIE ERROR when no valid grarg val is used for _initR_grDevices()';

$phstats0->set_rbase('');

throws_ok { $phstats0->_initR_grDevices(); } qr/ERROR: rbase/, 
    'TESTING DIE ERROR when rbase is not set for _initR_grDevices()';


## Clean the R dirs

#$srh0->cleanup();

####
1; #
####
