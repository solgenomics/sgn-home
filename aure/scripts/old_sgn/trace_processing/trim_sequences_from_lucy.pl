#!/usr/bin/perl

=head1 NAME

 trim_sequences_from_lucy.pl.
 A script to trim the sequences in fasta format using the output of LUCY program.(version.1.0.).

=cut

=head1 SYPNOSIS

 trim_sequences_from_lucy.pl -f <fastasequence>  
    
=head2 I<Flags:>

=over

=item -f 

B<fastasequence>        sequence lucy output in fasta format (this means that have clean region coordenates after id (mandatory)
  
=back

=cut

=head1 DESCRIPTION

This script trim the sequences if add the coordeantes after them
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

trim_sequences_from_lucy.pl


=cut

use strict;
use Bio::SeqIO;
use File::Basename;
use Getopt::Std;

our ($opt_f, $opt_h);
getopts("f:h");
if (!$opt_f && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}
check_input();

print "\nProcessing ... $opt_f fasta file.\n\n";

my $fastaseq=$opt_f;

my $input = Bio::SeqIO->new(	-file 	=> 	$fastaseq,
  	              		-format => 	'fasta'
			);
my $output = 'trimmed'.$fastaseq;
open my $out, '>', $output || die "I can not open the output file $output\n";

while ( my $seq_obj = $input->next_seq() ) {
	my $id=$seq_obj->display_id;
        my $desc=$seq_obj->desc;
        my @coord=split (/ /, $desc);
        my $begin=$coord[3];
        my $length=$coord[4];
        my $seq=$seq_obj->seq;
        my $trimseq=substr($seq, $begin, $length);
        print $out ">$id\n$trimseq\n";
	
}


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
     A script to trim the sequences in fasta format that comes from LUCY program
    
    Usage: 
      trim_sequences_from_lucy.pl [-h] -f <fasta_file_from_lucy_output>

    Examples:
      trim_sequences_from_lucy.pl [-h] -f tobacco_seq_trimmed.fasta

    Flags:
      -f fasta file	      file with the sequences in fasta format.
      -h this help

EOF
exit (1);
}

=head2 check_input

  Usage: check_input()
  Desc: check if exists some mandatory inputs. Also check the compatibility of them (for example -r and -R).
  Ret: none
  Args: none
  Side_Effects: die if there are something wrong
  Example: check_input();

=cut

sub check_input {
  if (!$opt_f) {	
	die "Sorry, you do not specify a fasta file in the -f option.\n";
  }
  
}
