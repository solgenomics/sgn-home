#!/usr/bin/perl

=head1 NAME

  phygeannot.t
  A piece of code to test the PhyGeAnnot module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl phygeannot.t
 prove phygeannot.t

=head1 DESCRIPTION

 Test PhyGeAnnot module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 45;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

use Image::Size;

## TEST 1 to 4

BEGIN {
    use_ok('PhyGeAnnot');
    use_ok('PhyGeTopo');
    use_ok('PhyGeCluster');
    use_ok('YapRI::Base');
}


## Create an empty object and test the possible die functions. TEST 5 to 8

my $phannot0 = PhyGeAnnot->new({ rbase => '' });

is(ref($phannot0), 'PhyGeAnnot', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { PhyGeAnnot->new(['fake']) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { PhyGeAnnot->new({ fake => {} }) } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a valid arg.';

throws_ok { PhyGeAnnot->new({ rbase => 'fk2'}) } qr/ERROR: fk2/, 
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
my $gofile = "$FindBin::Bin/testfiles/go.test.tab";


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

my $srh0 = YapRI::Base->new();



###############
## ACCESSORS ##
###############

## Get/Set_phygetopo, TEST 9 to 18

$phannot0->set_phygetopo({'ML' => $phygetopo1});
my %phygetopo1_r = %{$phannot0->get_phygetopo()};

is(ref($phygetopo1_r{'ML'}), 'PhyGeTopo',
   "Testing Get/Set_phygetopo, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phannot0->set_phygetopo() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_phygetopo function';

throws_ok { $phannot0->set_phygetopo('fake') } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to set_phygetopo isnt HASH obj.';

throws_ok { $phannot0->set_phygetopo({'name' => 'fake'}) } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to set_phygetopo isnt PhyGeTopo obj.';

$phannot0->add_phygetopo('NJ', $phygetopo2);
my %phygetopo_2r = %{$phannot0->get_phygetopo()};

is(scalar(keys %phygetopo_2r), 2, 
    "Testing Add_phygetopo, checking number of elements")
    or diag("Looks like this has failed");

throws_ok { $phannot0->add_phygetopo() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to add_phygetopo function';

throws_ok { $phannot0->add_phygetopo('fake') } qr/ERROR: No phygetopo/, 
    'TESTING DIE ERROR when no phygetopo arg. wassupplied to add_phygetopo';

throws_ok { $phannot0->add_phygetopo('name', 'fake') } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to add_phygetopo isnt PhyGeTopo obj.';

my %phygetopo_3r = %{$phannot0->delete_phygetopo()};

is(scalar(keys %phygetopo_3r), 2, 
    "Testing delete_phygetopo, checking number of elements deleted")
    or diag("Looks like this has failed");

my %phygetopo_4r = %{$phannot0->get_phygetopo()};

is(scalar(keys %phygetopo_4r), 0, 
    "Testing delete_phygetopo, checking empty hash")
    or diag("Looks like this has failed");


## Get/Set_rbase, TEST 19 to 21

$phannot0->set_rbase($srh0);
my $rbase1 = $phannot0->get_rbase();

is(ref($rbase1), 'YapRI::Base',
   "Testing Get/Set_rbase, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phannot0->set_rbase() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_rbase function';

throws_ok { $phannot0->set_rbase(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_rbase isnt YapRI::Base';

## Get/Set_gene_annot, TEST 22 to 27

my %gene_annot0 = (
    Sly_01101 => { 'nr' => { subject_id  => 'NP_187488.1', 
			     description => 'ATSIK protein kinase',
			     evalue      => '2e-156',
		   },
    },
    Sly_01219 => { 'nr' => { subject_id  => 'XP_002517269.1',
			     description => 'Auxin-induced protein 5NG4',
			     evalue      => '1e-153',
		   },
    },
    Nta_00407 => { 'nr' => { subject_id  => 'AAV44205.1',
			     description => 'unknow protein',
			     evalue      => '4e-37',
		   },	
    },
    Nto_40803 => { 'nr' => { subject_id  => 'AAV44205.1',
			     description => 'unknow protein',
			     evalue      => '2e-49',
		   },
    },	
    Nsy_35378 => { 'nr' => { subject_id  => 'YP_358636.1',
			     description => 'hypothetical protein PhapfoPp090',
			     evalue      => '4e-47',
		   },
    }
    );

$phannot0->set_gene_annot(\%gene_annot0);
my %gene_annot1 = %{$phannot0->get_gene_annot()};

is(scalar(keys %gene_annot1), scalar(keys %gene_annot0),
    "Testing get/set_gene_annot, checking number of annotations")
    or diag("Looks like this has failed");

is(join(',', sort keys %gene_annot1), join(',', sort keys %gene_annot0),
    "Testing get/set_gene_annot, checking list of members")
    or diag("Looks like this has failed");

throws_ok { $phannot0->set_gene_annot() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_gene_annot function';

throws_ok { $phannot0->set_gene_annot('fake') } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to set_gene_annot isnt HASH obj.';

throws_ok { $phannot0->set_gene_annot({'name' => 'fake'})} qr/ERROR: href mem/, 
    'TESTING DIE ERROR when arg. supplied to set_gene_annot isnt HREF of HREF.';

throws_ok { $phannot0->set_gene_annot({ AB => {C => 'f'}})} qr/ERROR: href bl/, 
    'TESTING DIE ERROR when arg.supplied to set_gene_annot isnt HREF-HREF-HREF';


## Test add_gene_annot, TEST 28 to 32

$phannot0->add_gene_annot('Nto_07827', 'nr', 
			  { subject_id  => 'CBI21007.3',
			    description => 'unnamed protein produc',
			    evalue      => '8e-72',
			  }
    );

my %gene_annot2 = %{$phannot0->get_gene_annot()};

is(scalar(keys %gene_annot2), 6,
    "Testing add_gene_annot, checking number of annotations")
    or diag("Looks like this has failed");

throws_ok { $phannot0->add_gene_annot() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to add_gene_annot function';

throws_ok { $phannot0->add_gene_annot('A') } qr/ERROR: No blastdb/, 
    'TESTING DIE ERROR when no blastdb arg. was supplied to add_gene_annot';

throws_ok { $phannot0->add_gene_annot('A', 'B')} qr/ERROR: No blast/, 
    'TESTING DIE ERROR when no blast arg. was supplied to add_gene_annot.';

throws_ok { $phannot0->add_gene_annot('A','B','C')} qr/ERROR: C for add/, 
    'TESTING DIE ERROR when third arg. supplied to add_gene_annot isnt href.';


## Test get/set_go_annot, TEST 33 to 37

my %go_annot0 = ( Sly_01101 => { 'GO:0004672' => 'protein kinase activity'},
		  Sly_01219 => { 'GO:0009734' => 'response to auxin stimulus'},
		  Nta_00407 => { },
		  Nto_40803 => { },	
		  Nsy_35378 => { }
    );

$phannot0->set_go_annot(\%go_annot0);
my %go_annot1 = %{$phannot0->get_go_annot()};

is(scalar(keys %go_annot1), scalar(keys %go_annot0),
    "Testing get/set_go_annot, checking number of annotations")
    or diag("Looks like this has failed");

is(join(',', sort keys %go_annot1), join(',', sort keys %go_annot0),
    "Testing get/set_go_annot, checking list of members")
    or diag("Looks like this has failed");

throws_ok { $phannot0->set_go_annot() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_go_annot function';

throws_ok { $phannot0->set_go_annot('fake') } qr/ERROR: fake/, 
    'TESTING DIE ERROR when arg. supplied to set_go_annot isnt HASH obj.';

throws_ok { $phannot0->set_go_annot({'name' => 'fake'})} qr/ERROR: fake for/, 
    'TESTING DIE ERROR when arg. supplied to set_go_annot isnt HREF of HREF.';

## Add_go_annot, TEST 38 to 41

$phannot0->add_go_annot('Nto_07827', {});

my %go_annot2 = %{$phannot0->get_go_annot()};

is(scalar(keys %go_annot2), 6,
    "Testing add_go_annot, checking number of annotations")
    or diag("Looks like this has failed");

throws_ok { $phannot0->add_go_annot() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to add_go_annot function';

throws_ok { $phannot0->add_go_annot('A') } qr/ERROR: No go href/, 
    'TESTING DIE ERROR when no blastdb arg. was supplied to add_go_annot';

throws_ok { $phannot0->add_go_annot('A', 'B')} qr/ERROR: B for add_/, 
    'TESTING DIE ERROR when 2nd arg. supplied to add_go_annot. isnt HREF';


## parse_go_file, TEST 42 to 45

my $go_href = PhyGeAnnot::parse_go_file($gofile);

my $gocount = 0;
foreach my $memb (keys %{$go_href}) {
    $gocount += scalar(keys %{$go_href->{$memb}});
}

is(scalar(keys %{$go_href}), 39, 
   "Testing parse_go_file, checking number of members")
   or diag("Looks like this has failed");

is($gocount, 40, 
   "Testing parse_go_file, checking number of GO terms")
   or diag("Looks like this has failed");

throws_ok { PhyGeAnnot::parse_go_file() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to parse_go_file function';

throws_ok { PhyGeAnnot::parse_go_file('A', 'B') } qr/ERROR: B for parse/, 
    'TESTING DIE ERROR when no arg. href supplied to parse_go_file isnt HREF.';




## Clean the R dirs

$srh0->cleanup();

####
1; #
####
