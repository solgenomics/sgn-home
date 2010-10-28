#!/usr/bin/perl
=head1 NAME

load_protein_data.pl - Load protein sequences into the SGN database

=head1 USAGE

load_protein_data.pl [-e | -l ] [database] [protein_fasta] [cds_fasta]

=head1 DESCRIPTION

This script load protein data for SGN unigenes into the sgn.cds table of the SGN database. Proteins can either be predicted with ESTScan or with the longest 6-frame translation. 

The ids in the cds and protein files need to correspond to SGN unigenes identifiers.

When loading ESTScan data, use the -e option, and ESTScan specific data will be parsed from the description line.

When loading longest six frame translations, use the -l option to parse the direction information that is appended to the identifier.

The corresponding sequences in the protein and cds files need to be in the file in the same order.

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut

use strict;

use Getopt::Std;
use CXGN::DB::InsertDBH;
use CXGN::DB::Connection;
use CXGN::Transcript::Unigene;
use CXGN::Transcript::CDS;
use Bio::SeqIO;

our ($opt_e, $opt_l);
getopts('el');


if (!$opt_e && !$opt_l) { 
    print STDERR "You need to give either the -e option for ESTScan data, or the -l option for the longest protein sequence data\n";
    exit(-1);
}

if ($opt_e && $opt_l) { 
    print STDERR "-e and -l are mutually exclusive. Please choose -e for ESTScan and -l for longest frame translations.\n";
    exit(-1);
}

my $db = shift;
my $protein_fasta = shift;
my $cds_fasta = shift;

if (!$cds_fasta) { 
    print "USAGE: $0 [-e | -l] [ test | production ] protein_fasta_file cds_fasta_file\n";
    print "(use -e for estscan data [not yet implemented])\n\n";

    exit(0);
}

my $dbhost = "db-devel";
my $dbname = "sandbox";


if ($db eq "production") {
    $dbhost="db";
    $dbname="cxgn";
}

print STDERR "USING DBHOST $dbhost DBNAME $dbname. Type RETURN to continue.\n";
my $continue = <>;
if ($continue =~/n/i) { print STDERR "\naborted\n"; exit(); }

my $dbh = CXGN::DB::InsertDBH::connect({dbhost=>$dbhost,dbname=>$dbname});

my $io = Bio::SeqIO->new( -format=>"fasta", -file=>$protein_fasta);
my $cds_io = Bio::SeqIO->new( -format=>"fasta", -file=>$cds_fasta);

my $old_unigene;


my $store_count=0;
my $skipped_count = 0;
while (my $seq = $io->next_seq())  {
    my $id = $seq->id();
    my $prot = $seq->seq();

    my $cds_obj = $cds_io ->next_seq();
    my $cds_id = $cds_obj->id();
    my $cds_seq = $cds_obj ->seq();
    my $cds_desc = $cds_obj ->desc();

    if ($cds_id ne $id ) { 
	print STDERR "Warning! CDS and protein ids don't match ($cds_id vs $id)\n";
    }

    if ($id !~ /SGN-U\d+/i) { die "Need a unigene id for each sequence\n"; }
#    my $cds = CXGN::Transcript::CDS->new_with_unigene_id($dbh, $unigene_id);
    my $cds = CXGN::Transcript::CDS->new($dbh);


    my $unigene_id;
    my $frame;
    my $score;
    my $seq_edits;

    my $forward_reverse = "F"; 
    

    if ($opt_e) { 

	if ($id=~/SGN-U(\d+)/i) { 
	    $unigene_id=$1;
	}
	else { 
	    print STDERR "Skipping $id -- illegal id format [SGN-U\d+].\n";
	    next();
	}

	my $unigene = CXGN::Transcript::Unigene->new($dbh, $unigene_id);
	if (!$unigene) { print STDERR "Skipping $unigene_id... Not found in database.\n"; 
			 $skipped_count++;
			 next();
		     }
	my $seq_text = $cds_seq;
	my $seq_edits = $cds_seq; 
	$seq_edits =~ s/[a-z]//g; # remove lower case characters- these represented deleted nucleotides.
	
	if ($cds_desc=~/minus strand/i) { 
	    $forward_reverse = "R";
	}

	my ($score, $start, $end, @rest) = split /\s+/, $cds_desc;
#	print STDERR "SCORE: $score, START: $start, END: $end\n";
	$cds ->set_unigene_id($unigene_id);
	$cds ->set_seq_text($seq_text);
	$cds ->set_seq_edits($seq_edits);
	$cds ->set_cds_seq($seq_edits);
	$cds ->set_direction($forward_reverse);
	$cds ->set_begin($start);
	$cds ->set_end($end);
	$cds ->set_score($score);
	$cds ->set_protein_seq($prot);
	$cds->set_method('estscan');
    }

    elsif ($opt_l) { 
	if ($id =~ /(\d+)([-+]\d)/) { 
	    $unigene_id = $1;
	    $frame = $2;
	}
	if (!$unigene_id) { die "Need a unigene id"; }

	my $direction = "F";
	$cds->set_protein_seq($prot);
	$cds->set_unigene_id($unigene_id);
	$cds->set_method("longest6frame");
	$cds->set_cds_seq($cds_seq);
	if ($frame <0) { 
	    $direction = "R";
	}
	$cds->set_direction($direction);
	$cds->set_frame($frame);
    }

    

    if ($db !~ /test/) { 
	print STDERR "Storing cds for $id...";
	my $new_id = $cds->store();
	#print STDERR "Done\n";
	print "[$new_id]              \n";
	if ($new_id) { 
	    $store_count++;
	}
	else { 
	    die "An error occurred during store() for seq $id.";
	}
    }
    else { print STDERR "test mode, not storing...\n"; }
    $cds = undef;
    $old_unigene = $unigene_id;
}

print STDERR "Store $store_count proteins. Skipped $skipped_count proteins.\n";
print STDERR "Commit?\n";
my $yes = <>;
if ($yes !~ /n/i) { 
    print STDERR "Committing load...\n";
    $dbh->commit();
}
else { 
    print STDERR "Rolling back...\n";
    $dbh->rollback();
}

