#!/usr/bin/perl

=head1 NAME

 aln2fasta.pl
 Tool to extract sequences from an alignment file

=cut

=head1 SYPNOSIS

 aln2fasta.pl [-h] -i <input_alignment_file> -f <input_format> -o <output>
                   -G <remove_gaps>


=head2 I<Flags:>

=over

=item -i 

B<input_aln_file>         input alignment file (mandatory)

=item -f

B<file_format>            file format (clustalw by default)

=item -o

B<output_basename>        output (mandatory)

=item -G

B<remove gaps>            remove gaps from output

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script extract fasta sequences from an alignmnet file using bioperl.

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 aln2fasta.pl


=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::SeqIO;
use Bio::AlignIO;

our ($opt_i, $opt_f, $opt_o, $opt_G, $opt_h);
getopts("i:f:o:Gh");
if (!$opt_i && !$opt_f && !$opt_o && !$opt_G && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check them

my $input = $opt_i || 
    die("INPUT ARGS ERROR: Input arg. was not supplied (-i <aln_file>).\n");

my $format = $opt_f || 'clustalw';

my $output = $opt_o ||
    die("OUTPUT ARGS ERROR: Output arg. was not supplied (-o <fasta_file>).\n");

## Open alignment

my $in  = Bio::AlignIO->new(-file   => $input,
			    -format => $format);

## Get the sequences and put into a Sequence Object

my $out = Bio::SeqIO->new(-file => ">$output", -format => 'fasta' );

## Remove the gaps if $opt_G

my ($aln_count, $seq_count) = (0, 0);

print STDERR "\n\nSTART PARSING\n\n";

while ( my $aln = $in->next_aln() ) {

    $aln_count++;
    foreach my $seqobj ($aln->each_seq() ) {
	
	$seq_count++;
	if ($opt_G) {
	    my $seq = $seqobj->seq();
	    $seq =~ s/(-|\.)//g;
	    $seqobj->seq($seq);
	}
	$out->write_seq($seqobj);
    }
}


print STDERR "\n$aln_count alignments have been parsed.\n";
print STDERR "\n$seq_count sequences have been extracted.\n";
print STDERR "\nDONE.\n";




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

      This script extract fasta sequences from an alignmnet file using bioperl.

    Usage:
     
      aln2fasta.pl [-h] -i <input_alignment_file> -f <input_format> -o <output>
                        -G <remove_gaps>

    Examples:

      aln2fasta.pl [-h] -i tftomato.aln -o tomato.fasta -G

    Flags:

      -i <input_aln_file>         input alignment file (mandatory)
      -f <file_format>            file format (clustalw by default)
      -o <output_basename>        output (mandatory)
      -G <remove_gaps>            remove gaps from output
      -h <help>                   print the help
     

EOF
exit (1);
}

