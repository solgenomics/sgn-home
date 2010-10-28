#!/usr/bin/perl
use strict;

use CXGN::DB::InsertDBH;
use CXGN::Transcript::Unigene;

if (!@ARGV) { 
    print "USAGE: unigene_dump.pl <dbname> <build_id|organism name> [min seq len] [dump type]\n";
	print "Dump type can be 'sequence' (default) or 'scores'\n";
	exit 0;
}

my $dbname = shift;
my $build_id = shift;
my $min_seq_len = shift || 50;
my $dump_type = shift || "sequence";
validate_dump_type($dump_type);

print STDERR "NOTE: Dumping sequences longer than $min_seq_len.\n\n";

my $dbhost = "";
if ($dbname eq "sandbox") { 
    $dbhost = "scopolamine";
}
elsif ($dbname eq "cxgn") { 
    $dbhost = "db.sgn.cornell.edu";
}

my $dbh = CXGN::DB::InsertDBH -> new( { 
    dbname=>$dbname,
    dbhost=>$dbhost,
}
				       );

unless($build_id =~ /^\d+$/){
	#Resolve current build_id FROM case-insensitive internal match to groups.comment
	print STDERR "Resolving build_id FROM provided organism name (fragment) '$build_id'";
	my $sth = $dbh->prepare("
		SELECT unigene_build_id, groups.comment AS name
		FROM sgn.unigene_build
		LEFT JOIN sgn.groups ON (organism_group_id=group_id)
		WHERE groups.comment ILIKE ?
		AND status='C'
		");
	$sth->execute('%' . $build_id . '%');
	my $name = "";
	($build_id, $name) = $sth->fetchrow_array();
	print STDERR "\nOrganism: $name, Build ID: $build_id\n";
}

unless($build_id =~ /^\d+$/){
	print STDERR "Build id not provided, or could not be resolved from name\n";
	exit 1;
}


my @unigene_ids = CXGN::Transcript::Unigene::get_unigene_ids_by_build_id($dbh, $build_id);
foreach my $uid (@unigene_ids) { 
    print STDERR "Processing $uid...              \r";
    my $unigene = CXGN::Transcript::Unigene->new($dbh, $uid);
    
    my $sequence = $unigene->get_sequence();
	my $scores = $unigene->get_scores();
	my $sgn_id = $unigene->get_sgn_id();
	my $annotation = $unigene->get_annotation_string(1e-10, 300); #limit to 300 chars
    $annotation =~ s/\|//g; 
    if (length($sequence) >= $min_seq_len) { 
		my $content = $sequence;
		$content = $scores if $dump_type eq "scores";
		print ">$sgn_id $annotation\n$content\n"; 
	}
    else { 
	print STDERR "$sgn_id has a length of ".length($sequence)." and was therefore omitted from the dump (min sequence length is $min_seq_len).\n\n";
    }
}

sub validate_dump_type {
	my $dt = shift;
	my @dts = qw/ sequence scores /;
	my %dts = ();
	$dts{$_} = 1 foreach(@dts);
	die ("Unrecognized dump type: $dt\nValid dump types: " . join(", ", @dts) . "\n")
		unless $dts{$dt};
}	
