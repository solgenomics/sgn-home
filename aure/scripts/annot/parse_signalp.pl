#!/usr/bin/perl -w
use strict;

=head1 Name

parse_signalp.pl

=cut

=head1 Description

Parse the result of SignalP 3.0.

NN:  if for parameter D, it says "YES" in the result, the parser accepts the query sequence as a signal peptide and takes the predicted cleavage site.

HMM:  if it says "signal peptide" in the result, the parser accepts the query sequence as a signal peptide and takes the predicted cleavage site.

The result is a tab-delimit file.  id, NN prediction, HMM prediction.

=cut

=head1 Usage

The input file is standard input.

=cut

=head1 Author

Chenwei Lin (cl295@cornell.edu)

=cut

=head1 References

Bendtsen JD etal, Improved Prediction od SIgnal Peptides: SIgnalP 3.0, J. Mol. Biol., 340:783-795 (2004)

=cut


while (<>){
  if (/^Using/) {
    my ($id, $nn, $hmm, $d_param, $prediction) = ('','','','','');
    while (<>) {
      if (/^>/) {
	($id) = $_ =~ />(\S+)/;
      }
      if (/^\s+D/) {
	my @items = split /\s+/;
	$d_param = $items[5];
      }
      if (/Most likely .+ (\d+) and \d+/){
	$nn = $1;
      }
      if (/Prediction/){
	($prediction) = $_ =~ /(\S+\s\S+)$/;
      }
      if (/Max cleavage .+\s+(\S+)\s+and\s+\d+$/) {
	$hmm = $1;  
	last;
      }
    }
    ($d_param ne 'YES') and $nn = '';
    ($prediction ne 'Signal peptide') and $hmm = '';
    print "$id\t$nn\t$hmm\n";
  }
}
