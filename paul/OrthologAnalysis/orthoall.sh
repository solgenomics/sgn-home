#!/bin/bash
# Simple script to run Tom's ortholog finder on all gene-family aligned fasta files

# Change the path to the location of the aligned files.
FILES=/home/paul/ReferenceGenomeProject/data/tax-fasta/famfiles/alignedfiles/*

for f in $FILES
do

perl ./ortholog_support.pl -i $f -N 100 -k -m .015 -T ML

done


