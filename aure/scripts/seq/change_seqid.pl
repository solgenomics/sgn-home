#!/usr/bin/perl

=head1 NAME

 change_seqid.pl
 
 A script to change the sequence id in a file (version.0.1).

=cut

=head1 SYPNOSIS

 change_seqif.pl -f <fasta_input_file> -r <rootname> -s <sequence_number>

=head2 I<Flags:>

=over


=item -f

B<fasta_file>           fasta file (mandatory)

=item -s

B<sequence_number>      sequence number (mandatory)

=item -r

B<rootname>             rootname for the new sequences

=item -h

B<help>                 print the help

=back

=cut

=head1 DESCRIPTION

 This script replace sequence id for rootname+sequence_number

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 change_seqid.pl

=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::SeqIO;
use Math::BigFloat;


our ($opt_f, $opt_s, $opt_r, $opt_h);
getopts("f:s:r:h");
if (!$opt_f && !$opt_s && !$opt_r && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}


## Check arguments

my $fastafile = $opt_f || 
    die("DATA ARGUMENT ERROR: -f <fasta_file> WAS NOT SUPPLIED.\n");

my $sequence_number = $opt_s ||
    die("DATA ARGUMENT ERROR: -s <sequence_number> WAS NOT SUPPLIED");

my $rootname = $opt_r ||
    die("DATA ARGUMENT ERROR: -r <rootname> WAS NOT SUPPLIED");


## Create the output files

my $seqout_filename = $fastafile . '.seqidreplaced.fasta';
my $equiv_filename = $fastafile . '.seqidreplaced.equiv.tab';

my $l = length($sequence_number);

my $seqin = Bio::SeqIO->new( -format => 'fasta', -file => "$fastafile");
my $seqout = Bio::SeqIO->new( -format => 'fasta', -file => ">$seqout_filename");

open my $equiv_io, '>', $equiv_filename;

my $n = 0;

print STDERR "Starting replacing...\n\n";

while( my $seq = $seqin->next_seq()) {
    $n++;
    my $seq_id = $seq->display_id();

    print STDERR "\tProcessing sequence $n ($seq_id)\r";

    my $format = '%0' . $l . "s";
    my $new_seqid = $rootname . '_' . sprintf($format, $n);
    print $equiv_io "$seq_id\t$new_seqid\n";
    my $new_seq = Bio::Seq->new( -id => $new_seqid, -seq => $seq->seq());
    $seqout->write_seq($new_seq);
} 
print STDERR "\n\nDone.\n\n";



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

      This script replace sequence id for rootname+sequence_number
    
    Usage:
        
      change_seqif.pl -f <fastafile> -r <rootname> -s <sequence_number>

    Example:
	
      change_seqif.pl -f tomato.fasta -r contig -s 20000
      
    Flags:

      -f <fasta_file>           fasta file (mandatory)
      -s <sequence_number>      sequence number (mandatory)
      -r <rootname>             rootname for the new sequences
      -h <help>                 print the help

     
EOF
exit (1);
}




###
1;#
###
