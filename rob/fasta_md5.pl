#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

use Bio::SeqIO;

use Digest::MD5 'md5_hex';

#use Data::Dumper;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script  file file file ...

  Given one or more fasta files on standard input or as arguments,
  print a tab-delimited list of the sequence name, and the md5
  checksum of the sequence.

  Options:

    none yet

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('',\%opt) or usage();


@ARGV = ('-') unless @ARGV;

foreach my $file (@ARGV) {
  my $i = Bio::SeqIO->new( -format => 'fasta', -file => $file );
  while( my $seq = $i->next_seq ) {
    print join "\t", $seq->id, md5_hex( $seq->seq );
    print "\n";
  }
}
