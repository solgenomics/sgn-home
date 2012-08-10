#! /usr/bin/perl

use strict;
use warnings;
use Getopt::Std;

# Author: Paul Vaneck

# This script generates gff3 track files for gene families to be used on GBrowse.

use vars qw($opt_s $opt_f $opt_o);
# -s <species> (tomato, arabidopsis, rice, potato, grape)
# -f <species gff3 file path>
# -o <gene family file> Should be orthoMCL output file.

# get options
getopts("s:f:o:");


die "Usage: perl famtrackgen.pl -s <species> -f <species gff3 file path> -o <orthoMCL gene family file> > Output.gff3\n" unless($opt_s && $opt_f && $opt_o);

open(SPECIES, $opt_f) || die("Could not open species gff3 file."); 
open(FAMILY, $opt_o) || die("Could not open gene family file.");

my @species_array = <SPECIES>;
my $species_gff3 = join('', @species_array);

my $key;
my $type = "gene";
my $max_family_amount = 20;

if (($opt_s ne "tomato") && ($opt_s ne "arabidopsis") && ($opt_s ne "rice") && ($opt_s ne "potato") && ($opt_s ne "grape")){
	die "Species must be 'tomato', 'arabidopsis', 'rice', 'potato', or 'grape'.\n";

}
else {
	if ($opt_s eq "tomato") {
		$key = "Soly";
	}
	elsif ($opt_s eq "arabidopsis")	{
		$key = "AT";
	}
	elsif ($opt_s eq "rice") {
		$key = "Os"
	}
	elsif ($opt_s eq "potato") {
		$key = "PGSC";
	}
	elsif ($opt_s eq "grape") {
		$key = "GSVIV";
	}

}

my $gene_count = 1;
my $grepline;
open(OUTPUT, ">$opt_s-gen2.gff3");

print "Processing...\n";
print OUTPUT "##gff-version 3\n";

# Go through gene family file looking for the existence of a given species in each family.
while( my $line = <FAMILY>)  {
	chomp($line);

	# If the given species exists in the given family.
	if ($line =~ m/$key/){
		$line =~ m/([0-9]+) genes/;
		$gene_count = $1;
		$line =~ s/ORTHO.*:\s*//;
		my @genearray = split(/ /,$line);

		# for each gene in the given gene family
		foreach my $temp (@genearray){

			# Parse out (../data/tomato/) type stuff, leaving just the name.
			$temp =~ s/\(.*\)//;
		
			# Grape can have multiple representations, namely GSVIVT and GSVIVG
			if (($temp =~ m/GSVIV12X_T/) || ($temp =~ m/GSVIVT/)){
				$type = "mRNA";
			}

			# If it is the correct species.
			if ($temp =~ m/$key/){

				# These grep lines get the corresponding gff3 line for the species name.
				my @grepline = grep /$temp/,@species_array;
				my @greplinegene = grep /$type\t/,@grepline;

				if ($greplinegene[0]){
				
					my $gene_line = $greplinegene[0];
					$gene_line =~ m/([a-zA-Z0-9-._]+)\t[a-zA-Z0-9-._]+\t$type\t([0-9]+)\t([0-9]+)\t[^\t]+\t[^\t]+\t[^\t]+\t(.+)\n/;
					#print $1."\t".$2."\t".$3."\t".$4;
					my $chr_id = $1;
					my $start = $2;
					my $end = $3;
					my $attributes = $4;
					if ($gene_count <= $max_family_amount) {

						# Add the given species coordinate to all other members of the family.
						foreach my $member (@genearray){

							if ($member eq $temp) {
								next;					
							}
							else {
								my $spec;

								# Determines the species of the given subject
								if ($member =~ m/AT/){
									$spec = "Arabidopsis thaliana";
								}
								elsif ($member =~ m/Soly/){
									$spec = "Solanum lycopersicum";
								}
								elsif ($member =~ m/Os/){
									$spec = "Oryza sativa";
								}
								elsif ($member =~ m/PGSC/){
									$spec = "Solanum tuberosum";
								}
								elsif ($member =~ m/GSVIV/){
									$spec = "Vitis vinifera";
								}
		
								$member =~ s/\(.*\)//;
								print OUTPUT $chr_id."\tOrthoMCL\tgene_family_member\t".$start."\t".$end."\t.\t-\t.\tName=$member;Note=$spec;\n";
							}
				
						}
					}
					else {

						print OUTPUT $chr_id."\tOrthoMCL\tgene_family_member\t".$start."\t".$end."\t.\t-\t.\tName=Large Family (Over $max_family_amount Members);\n";
						last;
					}
				}

			}

		}
		

	}

}

close(OUTPUT);

close(SPECIES);
close(FAMILY);

 







