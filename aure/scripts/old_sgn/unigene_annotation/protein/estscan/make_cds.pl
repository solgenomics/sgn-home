#!/usr/bin/perl
use strict;

usage() unless (@ARGV);

my $input_fasta = shift;
fatal("Input file not found") unless -f $input_fasta;

my $estpath = shift || "/data/local/website/cxgn/core/sgn-tools/estscan-3.0.2";
fatal("Estscan program directory not found") unless -d $estpath;

my $matrixfile = shift || $estpath . "/Arabidopsis_thaliana.smat";
fatal("Arabidopsis matrix file not found") unless -f $matrixfile;

system "$estpath/estscan -M $matrixfile $input_fasta";


sub usage {
	print "Usage:\n$0 <unigene_fasta> [non-svn estscan program dir] [Arabidopsis Matrix File]\n";
	print "Prints predicted cds to STDOUT\n";
	exit 0;
}

sub fatal {
	my $mesg = shift;
	chomp $mesg;
	print STDERR "Error: $mesg\n";
	exit 1;
}
