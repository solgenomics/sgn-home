#!/usr/bin/perl

=head1 NAME

  phygetopo.t
  A piece of code to test the PhyGeTopo module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl phygetopo.t
 prove phygetopo.t

=head1 DESCRIPTION

 Test PhyGeTopo module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 41;
use Test::Exception;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1 to 5

BEGIN {
    use_ok('PhyGeCluster');
    use_ok('PhyGeTopo');
    use_ok('Bio::Tree::TopoType');
    use_ok('Bio::Tree::Node');
    use_ok('Bio::Tree::Tree');
}


## Create an empty object and test the possible die functions. TEST 6 to 9

my $phygetop0 = PhyGeTopo->new();

is(ref($phygetop0), 'PhyGeTopo', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

throws_ok { PhyGeTopo->new(['fake']) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { PhyGeTopo->new({ fake => {} }) } qr/ARG. ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a valid arg.';

throws_ok { PhyGeTopo->new({ seqfams => 'fk2'}) } qr/ERROR: fk2/, 
    'TESTING DIE ERROR for new() when arg. doesnt have valid value';


###############
## ACCESSORS ##
###############

## 1) CHECK GET/SET_SEQFAMS, TEST 10 to 14
##    Create three empty seqfams objects

my %seqs0 = ( 
    'sid1' => Bio::Seq->new(-id => 'sid1', -seq => 'AACCTGTGCGATTAGACGGCAGTTT'),
    'sid2' => Bio::Seq->new(-id => 'sid2', -seq => 'TGGCTGTGGGATTAGACGGTTTTTT'),
    'sid3' => Bio::Seq->new(-id => 'sid3', -seq => 'ACTCTGTGCGATTAGACGGCAAATT'),
    'sid4' => Bio::Seq->new(-id => 'sid4', -seq => 'TGGCTGTGCGATTAGACCCTTTTTT'),
    'sid5' => Bio::Seq->new(-id => 'sid5', -seq => 'AACCTGTGCGAAAAGACGGCAGTTT'),
    'sid6' => Bio::Seq->new(-id => 'sid6', -seq => 'AACCTGTGCGAAAAGACGGCAGTTT'),
    'sid7' => Bio::Seq->new(-id => 'sid7', -seq => 'TGGCTGTGCGATTAGACGGTTTTTT'),
    'sid8' => Bio::Seq->new(-id => 'sid7', -seq => 'ACTCTGTGCGATTAGACCCCAAATT'),
    'sid9' => Bio::Seq->new(-id => 'sid7', -seq => 'ACTCTGTGCGATTTATCGGCAAATT'),
    );

my $seqmemb01 = [ $seqs0{sid1}, $seqs0{sid5}, $seqs0{sid6} ];
my $seqmemb02 = [ $seqs0{sid4}, $seqs0{sid2}, $seqs0{sid7} ];
my $seqmemb03 = [ $seqs0{sid8}, $seqs0{sid9}, $seqs0{sid3} ];

my %seqfams0 = (
    'cl_id1' => Bio::Cluster::SequenceFamily->new( -family_id => 'cl_id1',
						    -members   => $seqmemb01 ),
    'cl_id2' => Bio::Cluster::SequenceFamily->new( -family_id => 'cl_id2',
						    -members   => $seqmemb02 ),
    'cl_id3' => Bio::Cluster::SequenceFamily->new( -family_id => 'cl_id3',
						    -members   => $seqmemb03 ),
    );

## Test set/get_seqfams function

$phygetop0->set_seqfams(\%seqfams0);
my %seqfams0g = %{$phygetop0->get_seqfams()};

is(scalar(keys %seqfams0g), 3,
    "Testing get/set_seqfams, checking number of SeqFam objects")
    or diag("Looks like this has failed");

my $wrong_seqfam_objs = 0;
foreach my $fam_id0 (keys %seqfams0) {
    unless (ref($seqfams0{$fam_id0}) eq 'Bio::Cluster::SequenceFamily') {
	$wrong_seqfam_objs++;
    }
}

is($wrong_seqfam_objs, 0, 
    "Testing get/seq_seqfams, checking SeqFam objs. identity")
    or diag("Looks like this has failed");

throws_ok { $phygetop0->set_seqfams() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_seqfams function';

throws_ok { $phygetop0->set_seqfams('fake') } qr/ARG. ERROR: fake set_/, 
    'TESTING DIE ERROR when arg. supplied to set_seqfams() isnt HashRef';

throws_ok { $phygetop0->set_seqfams({id => 'fake2'}) } qr/ARG. ERROR: fake2 u/, 
    'TESTING DIE ERROR when arg. element supplied to set_seqfams() isnt SeqFam';

## 2) CHECK GET/SET_STRAINS, TEST 15 to 17
##    Create strains

my %strains0 = ( 
    'sid1' => 'Str1',
    'sid2' => 'Str2',
    'sid3' => 'Str3',
    'sid4' => 'Str1',
    'sid5' => 'Str2',
    'sid6' => 'Str3',
    'sid7' => 'Str1',
    'sid8' => 'Str1',
    'sid9' => 'Str2',
    );

## Test set/get_strains function

$phygetop0->set_strains(\%strains0);
my %strains0g = %{$phygetop0->get_strains()};

is(scalar(keys %strains0g), 9,
    "Testing get/set_strains, checking number of strains objects")
    or diag("Looks like this has failed");

throws_ok { $phygetop0->set_strains() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_strains function';

throws_ok { $phygetop0->set_strains('fake') } qr/ARG. ERROR: When arg/, 
    'TESTING DIE ERROR when arg. supplied to set_strains() isnt HashRef';


## 3) CHECK GET/SET_ANNOTATIONS, TEST 18 to 24
##    Create annotation hash

my %annot0 = (
    'cl_id1' => { db1 => 'description 1-1', db2 => 'description 2-1' },
    'cl_id2' => { db1 => 'description 1-2', db2 => 'description 2-2' },
    'cl_id3' => { db1 => 'description 1-3', db2 => 'description 2-3' },
    );

$phygetop0->set_annotations(\%annot0);
my %annot0g = %{$phygetop0->get_annotations()};

is(scalar(keys %annot0g), 3,
    "Testing get/set_annotations, checking number of annotations")
    or diag("Looks like this has failed");

is(scalar(keys %{$annot0g{'cl_id1'}}), 2,
    "Testing get/set_annotations, checking number of tags for first element")
    or diag("Looks like this has failed");

my %db1_annot0g = %{$phygetop0->get_annotations('db1')};

is(scalar(keys %db1_annot0g), 3,
    "Testing get/set_annotations for a single tag, checking number of annot.")
    or diag("Looks like this has failed");

is($db1_annot0g{'cl_id2'}, 'description 1-2',
    "Testing get/set_annotations for a single tag, checking that it has 1 tag")
    or diag("Looks like this has failed");

throws_ok { $phygetop0->set_annotations() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_annotations function';

throws_ok { $phygetop0->set_annotations('fake') } qr/ARG. ERROR: When arg./, 
    'TESTING DIE ERROR when arg. supplied to set_annotations() isnt HashRef';

throws_ok { $phygetop0->set_annotations({id => 'fk2'}) } qr/ ERROR: When val./, 
    'TESTING DIE ERROR when arg.value supplied to set_annotations() isnt HASH';

## 3) CHECK GET/SET_TOPOTYPES, TEST 25 to 30
##    Create topotype hash with the trees.

my %trees_data0 = (
    topo01 => { x => [{ Str1 => 1, Str2 => 1 }, 1], Str3 => 1 },
    cl_id1 => { x => [{ sid1 => 0.5, sid5 => 0.2 }, 0.1], sid6 => 0.3 },
    topo02 => { x => [{ Str1 => 1, Str3 => 1 }, 1], Str2 => 1 },
    cl_id2 => { x => [{ sid4 => 0.3, sid7 => 0.4 }, 0.5], sid2 => 0.1 },
    cl_id3 => { x => [{ sid8 => 0.6, sid3 => 0.3 }, 0.2], sid9 => 0.2 },
    );

my %trees_obj0 = ();
foreach my $tid (keys %trees_data0) {
    my %data = %{$trees_data0{$tid}};
    
    my @root_desc = ();
    foreach my $nid (keys %data) {
	if ($nid ne 'x') {
	    my $node = Bio::Tree::Node->new( -id            => $nid, 
					     -branch_length => $data{$nid} );
	    push @root_desc, $node;
	}
	else {
	    my %descdata = %{$data{$nid}->[0]};
	    my $brl = $data{$nid}->[1];
	    my @node_desc = ();
	    foreach my $nnid (keys %descdata) {
		my $nbrl = $descdata{$nnid};
		my $nnode = Bio::Tree::Node->new( -id            => $nnid, 
						  -branch_length => $nbrl );
		push @node_desc, $nnode;
	    }
	    my $node = Bio::Tree::Node->new( -descendents   => \@node_desc, 
					     -branch_length => $brl         );
	    
	}
    }
    my $root = Bio::Tree::Node->new( -descendents => \@root_desc );
    my $tree = Bio::Tree::Tree->new( -id => $tid, -root => $root);
    $trees_obj0{$tid} = $tree;
}

my $topoty01 = Bio::Tree::TopoType->new( { topology => $trees_obj0{topo01},
			 	           members  => [ $trees_obj0{cl_id1} ],
			      }
			      );
my $topoty02 = Bio::Tree::TopoType->new( { topology => $trees_obj0{topo02},
			 	           members  => [ $trees_obj0{cl_id2},
					                 $trees_obj0{cl_id3},
				           ],
			                 }
			               );

my %topotypes0 = ( cl_id1 => $topoty01, 
		   cl_id2 => $topoty02, 
		   cl_id3 => $topoty02 );

## Test set/get_topotypes function

$phygetop0->set_topotypes(\%topotypes0);
my %topotypes0g = %{$phygetop0->get_topotypes()};

is(scalar(keys %topotypes0g), 3, 
    "Testing get/set_topotypes, checking number of keys")
    or diag("Looks like this has failed");

my $wrong_topoobj = 0;
my %topodiff = ();
foreach my $clid0 (keys %topotypes0g) {
    my $topo0 = $topotypes0g{$clid0};
    if (ref($topo0) ne 'Bio::Tree::TopoType') {
	$wrong_topoobj++;
    }
    else {
	$topodiff{$topo0} = 1;
    }
}

is($wrong_topoobj, 0, 
    "Testing get/set_topotypes, checking identity of the objects")
    or diag("Looks like this has failed");

is(scalar(keys %topodiff), 2,
    "Testing get/set_topotypes, checking identity of the objects")
    or diag("Looks like this has failed");

throws_ok { $phygetop0->set_topotypes() } qr/ARG. ERROR: No arg. was/, 
    'TESTING DIE ERROR when no arg. was supplied to set_topotypes function';

throws_ok { $phygetop0->set_topotypes('fake') } qr/ARG. ERROR: When arg./, 
    'TESTING DIE ERROR when arg. supplied to set_topotypes() isnt HashRef';

throws_ok { $phygetop0->set_topotypes({id => 'fk2'}) } qr/ARG. ERROR: Values/, 
    'TESTING DIE ERROR when arg.value supplied to set_topotypes() isnt HASH';



##########################
## ANALYTICAL FUNCTIONS ##
##########################

## First it will prepare a new phygecluster object to get all the data
## to phygetopo object, TEST 31 to 36

my $blastfile = "$FindBin::Bin/testfiles/selfblast.test.m8";
my $seqfile = "$FindBin::Bin/testfiles/seq.test.fasta";
my $strainfile = "$FindBin::Bin/testfiles/strains.test.tab";
my $acefile = "$FindBin::Bin/testfiles/assembly_out.test.ace";
my $blastdbfile = "$FindBin::Bin/testfiles/blastref.test.fasta";

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

my @parameters1 = ('quiet' => 'yes', 'matrix' => 'BLOSUM');
$phygecluster->run_alignments({ program    => 'clustalw', 
				parameters => \@parameters1 });

$phygecluster->run_distances({ method => 'Kimura' });

$phygecluster->run_mltrees({ dnaml => {}, outgroup_strain => 'Sly' });

my %seqfams = %{$phygecluster->get_clusters()};
my %strains = %{$phygecluster->get_strains()};

my %revstrains = ();
foreach my $seqid (keys %strains) {
    unless (exists $revstrains{$strains{$seqid}}) {
	$revstrains{$strains{$seqid}} = 1;
    }
}


## Create the phygetopo object and run the topoanalysis

my $phygetopo1 = PhyGeTopo->new({ seqfams => \%seqfams, 
				  strains => \%strains });

my %topotypes = $phygetopo1->run_topoanalysis();


is(scalar(keys %topotypes) <=> 0, 1,
    "Testing run_topoanalysis, checking topotypes > 0")
    or diag("Looks like this has failed");

my $wrong_branch = 0;
my $wrong_ids = 0;

foreach my $type (sort keys %topotypes) {
    my $topology = $topotypes{$type}->get_topology();
    my $nwcks = Bio::Tree::TopoType::_tree2newick($topology);
   
    my @nodes = $topology->get_nodes();
    foreach my $node (@nodes) {
	my $branch = $node->branch_length();
	my $id = $node->id();
	if (defined $branch) {
	    if ($branch != 1) {
		$wrong_branch++;
	    }
	}
	if (defined $id) {
	    unless (exists $revstrains{$id}) {
		$wrong_ids++;
	    }
	}
    }
}

is($wrong_branch, 0, 
   "testing run_topoanalysis, checking wrong branch length") 
    or diag("Looks like this has failed");

is($wrong_ids, 0, 
   "testing run_topoanalysis, checking wrong node ids (all should be strains)") 
    or diag("Looks like this has failed");


throws_ok { $phygetop0->run_topoanalysis('fake') } qr/ERROR: fake used/, 
    'TESTING DIE ERROR when arg. supplied to run_topoanalysis() isnt HashRef';

throws_ok { $phygetop0->run_topoanalysis({ fake2 => 1}) } qr/ERROR: fake2 isn/, 
    'TESTING DIE ERROR when arg. supplied to run_topoanalysis() isnt permited';

throws_ok { $phygetop0->run_topoanalysis({ base_toponame => ''}) } qr/ERROR: /, 
    'TESTING DIE ERROR when value supplied to run_topoanalysis() isnt permited';


## Testing out_topoanalysis, TEST 37 to 41

my $filename = $phygetopo1->out_topoanalysis(
    {
	tempfile => 1, 
	r_in     => 1,
	headers  => ['member_n'],
    });

my $wrong_outid = 0;
my $right_header = 0;

open my $ffh, '<', $filename;
while(<$ffh>) {
    chomp($_);
    unless ($_ =~ m/member_n/) {
	my @data = split(/\t/, $_);
	
	## It is formated to be used with R so it added a couple of "
	if ($data[0] =~ m/^"(.+)"$/) {
	    unless (exists $topotypes{$1}) {
		$wrong_outid++;
	    }
	}
    }
    else {
	$right_header++;
    }
}

is($wrong_outid, 0, 
    "Testing out_topoanalysis, checking row names")
    or diag("Looks like this has failed");

is($right_header, 1, 
    "Testing out_topoanalysis, checking header names")
    or diag("Looks like this has failed");

throws_ok { $phygetop0->out_topoanalysis('fake') } qr/ERROR: fake used/, 
    'TESTING DIE ERROR when arg. supplied to out_topoanalysis() isnt HashRef';

throws_ok { $phygetop0->out_topoanalysis({ fake2 => 1}) } qr/ERROR: fake2 isn/, 
    'TESTING DIE ERROR when arg. supplied to out_topoanalysis() isnt permited';

throws_ok { $phygetop0->out_topoanalysis({ headers => 'fk'}) } qr/ERROR: /, 
    'TESTING DIE ERROR when value supplied to out_topoanalysis() isnt permited';


####
1; #
####
