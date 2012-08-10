#!/usr/local/bin/perl

use warnings;
use strict;
no strict 'refs';

## Consensus sequence file must be in fasta format, with entries beginning with >
open (TOMATO, "tomato_Blaster_Rec_Map_consensus.fa");
my $consensusfile = "";
while (<TOMATO>) {
    $consensusfile=$consensusfile.$_;
}
close (TOMATO);

$consensusfile =~ s/\n//g;

## Separate consensus sequences into files, with label as filename and sequence as file contents
my @consensusseqs = split(/>/,$consensusfile);
shift(@consensusseqs);

my @labels;
my @seq;

foreach my $string (@consensusseqs) {
	$string =~ s/([0-9]{1})([A-Za-z]{1})/$1\t$2/;
	my @array = split(/\t/,$string);
	push(@labels, $array[0]);
	push(@seq, $array[1]);
};

$numberofelements = @labels;

for my $i (0 .. ($numberofelements-1)) {
	my $filehandler = "INPUT";
	$filehandler=$filehandler."$i";
	open $filehandler, ">>$labels[$i]";
	print $filehandler "$seq[$i]";
	close($filehandler);
};


## Blast all consensus sequences with database (S_lycopersicu.fa)
foreach my $label (@labels) {
	system("blastall -p blastn -d S_lycopersicum.fa -i $label -o $label.output -m 8");
};

## Filter output files by criteria (%identity >= 98%, aligned length >= 200bp), save as .hp files
foreach my $label (@labels) {
	my $filehandler = "output";
	$filehandler = $label.".".$filehandler;
	open $filehandler, "<$filehandler";
	while (<$filehandler>) {
		my $newfilehandler = $filehandler."hp";
		my @blastresults = split(/\t/,$_);
		if ($blastresults[2] >= 98 && ($blastresults[3]>=200)) {
			open $newfilehandler, ">>$filehandler.hp";
			print $newfilehandler "$_";
			close($newfilehandler);
		}
	};
	close($filehandler);
};

## Count number of hits that made it to the .hp files
my $filelist = `ls *.hp`;
my @files = split(/\n/, $filelist);
my %hash = ();
foreach my $file (@files) {
	my $filehandler2 = "file";
	my $count = 0;
	$filehandler2 = $filehandler2.$file;
	open $filehandler2, "<", $file;
	while (<$filehandler2>) {
		if ($_ =~ m/1_0/) {
			$count=$count+1;
		};
	};
	close($filehandler2);
	$hash{$file}=$count;
};

## Print to results file.
open RESULTS, ">resultsfile";
foreach my $key (keys %hash) {
	print RESULTS $key."\t".$hash{$key}."\n";
};
close(RESULTS);

