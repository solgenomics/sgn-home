#!/usr/bin/perl
use strict;

=head1 NAME

cleanup_signalp_output.pl

=head1 SYNOPSIS

Cleans up output by removing all lines beginning w/ #, except the first two in the file.  Also reformats file to tab-delimited form, w/ 2 tabs before the HMM data.  This is how secretary_loader.pl parses the SignalP output and also how the output is provided on the CBS website.

=head1 USAGE

IMPORTANT: input file should be formatted as the STDOUT from signalp running in "short" format mode:  ./signalp -t euk -f short [protein.fasta]
 ./cleanup_signalp_output.pl [file]

=cut

my $file = shift;
my $TEMP = $file . "_clean";
if(!(-f $file)) {die "Must send a file to clean"}
open(WF, ">$TEMP");
open(RF, $file);
my %seen;
my $iter = 0;
while(<RF>){
	if(/^\s*#/ && $iter>4) { next }	
	elsif(/^\s*#/) { print WF $_; next }
	$iter++;
	my ($unigene_id) = /SGN-U(\d+)/i;
	unless($unigene_id) { next }
	s/^\s+//;
	s/\s+/\t/g;
	s/\tSGN/\t\tSGN/;
	unless(/\n$/) { $_ .= "\n"; }
	print WF $_ unless($seen{$unigene_id}++)
}
close(WF);
system("cp $TEMP $file");
system("rm $TEMP");
