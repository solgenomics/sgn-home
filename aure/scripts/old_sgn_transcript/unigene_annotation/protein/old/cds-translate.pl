#!/usr/bin/perl

=head1 NAME

cds-translate.pl

=head1 USAGE

perl cds-translate.pl < cds_fasta > protein_fasta

=head1 DESCRIPTION

A script that translates cds into protein, being considerate of special characters that ESTScan inserts into the cds sequence.

=head1 AUTHOR(S)

Marty Kreuter, Chris Carpita.

=cut


use Bio::SeqIO;
#use CXGN::Phylo::Alignment::Member;
use  CXGN::Tools::Sequence;

my $s = Bio::SeqIO->new(-fh=>\*STDIN,-format=>'fasta');

while (my $seq = $s->next_seq) {
  my $seq2 =  CXGN::Tools::Sequence->new(
                  {type => "cds",
                  id => $seq->id,
                  seq => $seq->seq});

  eval { $seq2->translate_cds(); };
#  $seq2->translate_cds();
  if ($@) {
    print STDERR ">>".$@;
  } else {
    #Dies if there are problems
    #Does not return translated sequence, sets a hash value:
    my $protein_seq = $seq2->{translated_protein_seq};
    printf ">%s\n%s\n", $seq->id, $protein_seq;
  }
}
