#! /usr/bin/perl

use strict;
use warnings;

# Author: Paul Vaneck

# This script parses the ortholog files outputted from Tom's ortholog finding script: ortholog_support.pl.
# It creates a file containing each gene family with only the orthologs of each gene with an NJ bootstrap value
# of a given value $min out of 100 or above. 

# This is not very modular and has lots of hardcoding.

# 60 was chosen as a default minimum value.
my $min = 60;
my $gene;
my $count = 0;

my @files = `ls ~/ReferenceGenomeProject/data/tax-fa/famfiles/alignedfiles/orthologs/*.orthologs`;
open(OUTPUT, ">ortholist60.out");

foreach my $temp (@files){
	open(FILE, $temp);
	$temp =~ m/orthologs\/(.*\.orthologs)/;
	print OUTPUT ">$1\n";

	while( my $line = <FILE>)  {

		if ($line =~ m/^\s*$/){
			next;
		}

		if ($line =~ m/bootstrap/){
			$line =~ m/^([a-zA-Z0-9.-]+)\s+/;
			$gene = $1;
			$count = 0;
			
			#print OUTPUT "+".$1."\n";
		}
		elsif ($line =~ m/\s+([a-zA-Z0-9.-]+)\s+[0-9]\s+[0-9]\s+([0-9]+)/){
			if ($2 >= $min) {
				if ($count == 0){
					print OUTPUT "+".$gene."\n";
					print OUTPUT "\t".$1."\t".$2."\n";
					$count++;
				}
				else {
					print OUTPUT "\t".$1."\t".$2."\n";
				}
				
			}
		}

	}
	
	close(FILE);
}

print "Done.\n";

close(OUTPUT);
