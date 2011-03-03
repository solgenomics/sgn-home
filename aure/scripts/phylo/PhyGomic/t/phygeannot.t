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
use Test::More tests => 64;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

use Image::Size;

## TEST 1 to 4

BEGIN {
    use_ok('PhyGeAnnot');
    use_ok('PhyGeTopo');
    use_ok('PhyGeCluster');
    use_ok('R::YapRI::Base');
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
my $annotfile = "$FindBin::Bin/testfiles/seq.annot.blastx.ath.m8";
my $deflinefile = "$FindBin::Bin/testfiles/defline.test.tab";

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

is(ref($rbase1), 'R::YapRI::Base',
   "Testing Get/Set_rbase, Checking object identity")
    or diag("Looks like this has failed");

throws_ok { $phannot0->set_rbase() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_rbase function';

throws_ok { $phannot0->set_rbase(['fake']) } qr/ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg supplied set_rbase isnt R::YapRI::Base';

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

is(scalar(keys %{$go_href}), 63, 
   "Testing parse_go_file, checking number of members")
   or diag("Looks like this has failed");

is($gocount, 134, 
   "Testing parse_go_file, checking number of GO terms")
   or diag("Looks like this has failed");

throws_ok { PhyGeAnnot::parse_go_file() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to parse_go_file function';

throws_ok { PhyGeAnnot::parse_go_file('A', 'B') } qr/ERROR: B for parse/, 
    'TESTING DIE ERROR when no arg. href supplied to parse_go_file isnt HREF.';


## parse_go_file, TEST 46 to 48

## First set empty the go_terms

$phannot0->set_go_annot({});

$phannot0->load_go_file($gofile);

my %go_annot3 = %{$phannot0->get_go_annot()};

is(scalar(keys %go_annot3), 63,
    "Testing load_go_file, checking number of annotations")
    or diag("Looks like this has failed");

throws_ok { $phannot0->load_go_file() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to load_go_file function';

throws_ok { $phannot0->load_go_file('A', 'B') } qr/ERROR: B for load/, 
    'TESTING DIE ERROR when no arg. href supplied to load_go_file isnt HREF.';


## parse_blast_file, TEST 49 to 57


my %parse_bl_args = ( 
    blastdb       => 'ath',
    defline       => $deflinefile,
    );

my $blast_href = PhyGeAnnot::parse_blast_file( $annotfile, \%parse_bl_args );

my %blannot = %{$blast_href};

is(scalar(keys %blannot), 107, 
    "Testing parse_blast_file, checking number of annotations")
    or diag("Looks like this has failed");

my $right_bldb = 0;
my $right_params = 0;

foreach my $id (keys %blannot) {
    my %dbs = %{$blannot{$id}};
    foreach my $bldb (keys %dbs) {
	if ($bldb eq 'ath') {
	    $right_bldb++;
	}

	my %bldata = %{$dbs{$bldb}};
	my $test = join(',', keys %bldata);
	if (scalar(keys %bldata) == 13) {
	    $right_params++;
	}
    }
}

is($right_bldb, 107, 
    "Testing parse_blast_file, checking number of right blastdb names")
    or diag("Looks like this has failed");

is($right_params, 107, 
    "Testing parse_blast_file, checking number of right blast parameters")
    or diag("Looks like this has failed");

throws_ok { PhyGeAnnot::parse_blast_file() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to parse_blast_file function';

throws_ok { PhyGeAnnot::parse_blast_file('A') } qr/ARG. ERROR: Blast result/, 
    'TESTING DIE ERROR when blast result was supplied to parse_blast_file is 0';

my $bl = $annotfile;

throws_ok { PhyGeAnnot::parse_blast_file($bl, 'B') } qr/ERROR: B for parse/, 
    'TESTING DIE ERROR when arg. href supplied to parse_blast_file isnt HREF.';

throws_ok { PhyGeAnnot::parse_blast_file($bl, {}) } qr/ERROR: No blastdb/, 
    'TESTING DIE ERROR when no blastdb supplied to parse_blast_file';

throws_ok { PhyGeAnnot::parse_blast_file($bl,{blastdb =>'A'})} qr/ERROR: No de/,
    'TESTING DIE ERROR when no defline supplied to parse_blast_file';

my $wrong_args = {blastdb => 'A', defline => 'fake' };
throws_ok { PhyGeAnnot::parse_blast_file($bl, $wrong_args) } qr/ERROR: Defline/,
    'TESTING DIE ERROR when defline file supplied to parse_blast_file is 0';


## load_blast_file, TEST 58 TO 60

$phannot0->set_gene_annot({});

$phannot0->load_blast_file($annotfile, \%parse_bl_args);

my %gene_annot3 = %{$phannot0->get_gene_annot()};

is(scalar(keys %gene_annot3), 107,
    "Testing load_gene_file, checking number of annotations")
    or diag("Looks like this has failed");

throws_ok { $phannot0->load_blast_file() } qr/ARG. ERROR: No arg./, 
    'TESTING DIE ERROR when no arg. was supplied to load_gene_file function';

throws_ok { $phannot0->load_blast_file('A', 'B') } qr/ERROR: B for load/, 
    'TESTING DIE ERROR when no arg. href supplied to load_go_file isnt HREF.';


## go_conservative_annotation, TEST 61 to 64

$phannot0->set_phygetopo({'ML' => $phygetopo1});

my $phyg_annot = 0;
my %phygetopo2 = %{$phannot0->get_phygetopo()};
foreach my $method (keys %phygetopo2) {
    my $phygtp = $phygetopo2{$method};
    my %phygtp_annot = %{$phygtp->get_annotations()};
    $phyg_annot += scalar(%phygtp_annot);
}

is($phyg_annot, 0, 
    "Testing go_conservative_annotation, check annotations previous method.")
    or diag("Looks like this has failed");


$phannot0->go_conservative_annotation();

my %phygetopo3 = %{$phannot0->get_phygetopo()};
foreach my $method (keys %phygetopo3) {
    my $phygtp = $phygetopo3{$method};
    my %phygtp_annot = %{$phygtp->get_annotations()};
    $phyg_annot += scalar(keys %phygtp_annot);

    foreach my $id (keys %phygtp_annot) {
	my %go = %{$phygtp_annot{$id}->{'GO'}};
    }
}

is($phyg_annot, 2, 
    "Testing go_conservative_annotation, check annotations after method.")
    or diag("Looks like this has failed");

$phannot0->set_go_annot({});

throws_ok { $phannot0->go_conservative_annotation() } qr/No GO /, 
    'TESTING DIE ERROR when no GO was set for go_conservative_annotat.';

$phannot0->delete_phygetopo();

throws_ok { $phannot0->go_conservative_annotation() } qr/No phygetopo /, 
    'TESTING DIE ERROR when no phygetopo was set for go_conservative_annotat.';



####
1; #
####
