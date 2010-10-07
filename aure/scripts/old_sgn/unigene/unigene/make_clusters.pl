#!/usr/bin/perl -w
use strict;
use CXGN::BioTools::FastaParser;

if (!defined($ARGV[0]) || $ARGV[0] eq "help") {
    print STDERR <<EOF;

    Quick program to build FASTA cluster files (sequence and quality) given
    three files: FASTA soruce files sequence and quality, and a cluster file
    indicating which sequence ids belong in which cluster.

    order of arguments: <source seqfile> <source qualfile>

    Expects clustering information on STDIN, in format like this:
    >Cluster <n> COMMENT
    white space separated sequence ids on single line

    Dumps many cluster-<n>.seq and cluster-<n>.qual files in current directory

EOF
    exit 0;
}


my $seqfile = shift;
my $qualfile = shift;
print STDERR "Building clusters from $seqfile\n";
my $fasta_obj = CXGN::BioTools::FastaParser->load_file($seqfile,$qualfile);
my ($cluster_no, $inputline, %clusters, @singletons);
print STDERR "Loading cluster blocks\n";
while(<>) {
    chomp;
    next if m/^$/;

    if (m/^>Cluster ([0-9]+)/) {
	$cluster_no = $1;

	$inputline = <>;
	chomp $inputline;
	@{$clusters{$cluster_no}} = split/\s+/,$inputline;
	next;
    }
    if (m/^>Singletons/) {
	$inputline = <>;
	chomp $inputline;
	@singletons = split/\s+/,$inputline;
	next;
    }

    print "Unrecognized input: $_\n";
}

print STDERR "Sorting clusters by size\n";
my @cluster_ids = sort { @{$clusters{$b}} <=> @{$clusters{$a}} } keys %clusters;
my @clusters = map {$clusters{$_}} @cluster_ids;

my ($i, $seq, $qual, @qual, $header, $seqid);
for($i=0;$i<@clusters;$i++) {

    print STDERR "Writing ",scalar(@{$clusters[$i]})," sequences to cluster #$i\n";
    open OUTFILE, ">cluster-$i.seq"
	or die "Can't open cluster sequence file ($!)";
    open OUTFILE_QUAL, ">cluster-$i.seq.qual"
	or die "Can't open cluster quality file ($!)";

    foreach $seqid ( @{$clusters[$i]} ) {
	$seq = $fasta_obj->get_sequence($seqid, \@qual);
		if (! $seq) {
	    print "Error on query $seqid\n";
	    print $fasta_obj->last_error(), "\n";
	}
	if (! @qual) {
	    print "Error on query $seqid\n";
	    print $fasta_obj->last_error(), "\n";
	}
	$header = $fasta_obj->get_header($seqid);
	if (! $header) {		
	    print "Error retrieving header for $seqid\n";
	    print $fasta_obj->last_error(), "\n";
	}
	$seq = CXGN::BioTools::FastaParser::format_sequence($header, $seq);
	$qual = CXGN::BioTools::FastaParser::format_qualdata($header, \@qual);
	print OUTFILE "$seq";
	print OUTFILE_QUAL "$qual";
    }
    close OUTFILE;
    close OUTFILE_QUAL;
}

print STDERR "Outputting singletons to singletons.seq\n";
open OUTFILE, ">singletons.seq"
    or die "Can't open singletons sequence file ($!)";

open OUTFILE_QUAL, ">singletons.seq.qual"
    or die "Can't open singletons quality file ($!)";

foreach $seqid ( @singletons ) {
    $seq = $fasta_obj->get_sequence($seqid, \@qual);
    if (! $seq) {
	print "Error on query $seqid\n";
	print $fasta_obj->last_error(), "\n";
    }
    if (! @qual) {
	print "Error on query $seqid\n";
	print $fasta_obj->last_error(), "\n";
    }
    $header = $fasta_obj->get_header($seqid);
    if (! $header) {		
	print "Error retrieving header for $seqid\n";
	print $fasta_obj->last_error(), "\n";
    }
    print OUTFILE ">$seqid\n$seq\n";
    print OUTFILE_QUAL ">$seqid\n",join(" ",@qual),"\n";
}


close OUTFILE;
close OUTFILE_QUAL;
