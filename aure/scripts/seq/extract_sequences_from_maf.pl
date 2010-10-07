#!/usr/bin/perl

=head1 NAME

 extract_sequences_from_maf.pl
 Script to parse and get sequences from maf assembly file from MIRA (version.0.1).

=cut

=head1 SYPNOSIS

 extract_sequences_from_maf.pl [-h] -i <maf_file> [-o <output_basename>] [-U]

=head2 I<Flags:>

=over


=item -i

B<maf_file>                     assembly file in maf format format (mandatory)

=item -o

B<output_basename>              output basename (input_file.output by default)

=item -U

B<unpadded_sequences>           extract the unpadded sequences

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

 This script parse the assembly file in maf format that mira can produce.
 and extract the sequences (contigs and reads) producing the following
 outputs:

  + fasta
  + qual

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 extract_sequences_from_maf.pl

=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use Bio::Seq;
use Bio::Seq::Quality;
use Bio::SeqIO;
use Math::BigFloat;

our ($opt_i, $opt_o, $opt_U, $opt_h);
getopts("i:o:Uh");
if (!$opt_i && !$opt_o && !$opt_U && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check the mandatory ones

my $maf = $opt_i
    || die("MANDATORY ARGUMENT ERROR: -i <input_maffile> was not supplied.\n");
my $baseout = $opt_o 
    || $opt_i . '_output';

## First, parse the file and get the data

print STDERR "\n\nPARSING MAF FILE AND WRITTING OUTPUT FILES.\n\n";

open my $maf_fh, '<', $maf 
    || die("OPEN FILE ERROR: file=$maf could not be openned (system error=$!).\n");

my $l = 0;
my $n_co = 0;
my $n_rd = 0;

## To reduce the amount of memory used it will store only the data that it needs
## It will not use Bioperl to parse the ace because it parse and 
## everything (and for big files have memory problems)

## Define the catch variables

my ($contig_id, $contig_seq, $contig_qual);
my ($read_id, $read_seq, $read_qual);

## Define the files and open the IOs

my $contig_seq_file = $baseout . '.contigs.fasta';
my $contig_qual_file = $baseout . '.contigs.qual';
my $read_seq_file = $baseout . '.reads.fasta';
my $read_qual_file = $baseout . '.reads.qual';

my $contig_seqio = Bio::SeqIO->new( -file => ">$contig_seq_file", -format => 'fasta' );
my $contig_qualio = Bio::SeqIO->new( -file => ">$contig_qual_file", -format => 'qual' );
my $read_seqio = Bio::SeqIO->new( -file => ">$read_seq_file", -format => 'fasta' );
my $read_qualio = Bio::SeqIO->new( -file => ">$read_qual_file", -format => 'qual' );

while (<$maf_fh>) {
    $l++;

    print STDERR "\tParsing line=$l for file=$maf. ($n_co contigs and $n_rd reads have been parsed).      \r";

    chomp($_);

    ## Maf file works with tags and tabs. To catch a data type only needs to catch the tag
    ## and the data after the tag

    ## First pass the data to the objects

    if ($_ =~ m/^\\\\$/) { 
	if (defined $contig_id) {

	    ## Write the data into the seq object

	    my $co_seqobj = Bio::Seq::Quality->new( -id => $contig_id, -seq => $contig_seq, -qual => $contig_qual );

	    ## And write the object into the files
	    
	    $contig_seqio->write_seq($co_seqobj);
	    $contig_qualio->write_seq($co_seqobj);

	    $n_co++;

	    ## Finally clean the data undef them
	    
	    ($contig_id, $contig_seq, $contig_qual) = (undef, undef, undef);
	}
    }
    if (defined $read_id && $_ =~ m/^TN\t$read_id$/) {

	## Write the data into the seq object

	my $rd_seqobj = Bio::Seq::Quality->new( -id => $read_id, -seq => $read_seq, -qual => $read_qual );

	## And write the object into the files
	    
	$read_seqio->write_seq($rd_seqobj);
	$read_qualio->write_seq($rd_seqobj);

	$n_rd++;

	## Finally clean the data undef them
	    
	($read_id, $read_seq, $read_qual) = (undef, undef, undef);	
    }


    ## Catch the contig_id

    if ($_ =~ m/^CO\t(.+?)$/) { 
	$contig_id = $1;
    }

    if ($_ =~ m/^CS\t(.+?)$/) { 
	$contig_seq = $1;

	if ($opt_U) {
	    $contig_seq =~ s/\*//g;
	}
    }

    ## It needs from fastq to qual

    if ($_ =~ m/^CQ\t(.+?)$/) {

	my @qual;
    
	if (defined $1) {

	    my @ascii = split(//, $1);	    
	    foreach my $ascii (@ascii) {
		my $q = ord($ascii) - 33;
		if ($opt_U) {
		    if ($q != 1) {
			push @qual, $q;
		    }
		}
		else {
		    push @qual, $q;
		}
	    }
	}
	
	$contig_qual = join(' ', @qual);	
    }

    ## Catch the read_id

    if ($_ =~ m/^RD\t(.+?)$/) { 
	$read_id = $1;
    }

    if ($_ =~ m/^RS\t(.+?)$/) { 
	$read_seq = $1;
	if ($opt_U) {
	    $read_seq =~ s/\*//g;
	}
    }

    ## It needs from fastq to qual

    if ($_ =~ m/^RQ\t(.+?)$/) {

	my @qual;
    
	if (defined $1) {

	    my @ascii = split(//, $1);	    
	    foreach my $ascii (@ascii) {
		my $q = ord($ascii) - 33;
		if ($opt_U) {
		    if($q != 1) {
			push @qual, $q;
		    }
		}
		else {
		    push @qual, $q;
		}
	    }
	}
	
	$read_qual = join(' ', @qual);	
    }
}

print STDERR "\n\n";

## Now print the results

print STDERR "\n=======================================================================\n";
if ($opt_U) {
    print STDERR "\nGet Unpadded Sequence Option ENABLED.\n\n";
}
print STDERR "$n_co contigs and $n_rd reads\n\thave been extracted from the $maf file.\n";
print STDERR "\nThe outfiles are:\n\t+ $contig_seq_file\n\t+ $contig_qual_file\n\t+ $read_seq_file\n\t+ $read_qual_file\n\n";
print STDERR "\n=======================================================================\n\n\n";

=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script parse the assembly file in maf format that mira can produce.
     and extract the sequences (contigs and reads) producing the following
     outputs

       + fasta
       + qual
         
    Usage:

     extract_sequences_from_maf.pl [-h] -i <maf_file> [-o <output_basename>] [-U]
       
    Flags:
      
     -i <maf_file>              assembly file in ace format format (mandatory)
     -o <output_basename>       output basename (input_file.output by default)
     -U <unpadded_sequences>    extract the unpadded sequences
     -h <help>                  print the help
      

EOF
exit (1);
}
