#!/usr/bin/perl

=head1 NAME

primers2fasta.pl

=head1 SYNOPSIS

primers2fasta.pl -i tab_delimited_primer_name_and_sequences.txt > multi-fasta_file.fa

=head1 DESCRIPTION

This script returns a multi-fasta sequence from a tab-delimited file containing a marker name, forward primer, and reverse primer on each line

=cut

use strict;
use Getopt::Std;
#use Text::Format;
our ($opt_i, $opt_n, $opt_r, $opt_h);
getopts('i:h');
if ($opt_h){
  help();
  exit;
}
if (!$opt_i) {
    print STDERR "\nFilename required\n\n\n";
    help();
}
my $filename = $opt_i;

open FILE, "<", $filename or die "No such file $filename";
while (<FILE>) {
  chomp $_;
  my @row =  split('\t', $_);
  if (scalar(@row) == 3) {
      print '>'.$row[0]."\n".$row[1]."nnnnnnnnnnnnnnnnnnnn".$row[2]."\n";
  }
  elsif (scalar(@row) == 2) {
      print '>'.$row[0]."\n".$row[1]."\n";
  }  
  else {
      die "Wrong number of elements in line:\n$_\n";
  }
}

sub help {
  print STDERR <<EOF;
  $0:

    Description:
     This script returns a multi-fasta sequence from a tab-delimited file containing 
     a marker name, forward primer, and reverse primer (optional) on each line. 
     When forward and reverse primers are supplied, they are concatenated with 20 Ns

    Usage:
     
      primers2fasta.pl -i <input_tab_delimited_primers_file> > output.fa

    Flags:

      -i <input_file>         input tab delimited primers file  (mandatory)
      -h <help>                     help

EOF
exit (1);
}

=pod

  =back

  =head1 LICENSE

  Same as Perl.

  =head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

  =cut


