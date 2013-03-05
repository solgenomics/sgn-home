#!/usr/bin/perl

=head1 NAME

seq_from_fasta.pl

=head1 SYNOPSIS

seq_from_fasta.pl -i multi-fasta_file.fa -n header_of_fasta_to_get -r 1-100

=head1 DESCRIPTION

This script returns a fasta sequence from a multi-fasta file, and optionally returns a subsequence.

=cut

use strict;
use Getopt::Std;
#use Text::Format;
our ($opt_i, $opt_n, $opt_r, $opt_h);
getopts('i:n:r:h');
if ($opt_h){
  help();
  exit;
}
if (!$opt_i || !$opt_n) {
    print STDERR "\nFilename and fasta sequence name required\n\n\n";
    help();
}
my $filename = $opt_i;
my $search_id = $opt_n;
my $range_string = $opt_r;
my $range_start;
my $range_end;
my $sequence_found = 0;

open FILE, "<", $filename or die "No such file $filename";
my $collect_sequence = 0;
my $sequence_string;
my $header_found_count = 0;
while (<FILE>) {
  chomp $_;
  if (substr($_,0,1) eq '>') {
    $collect_sequence = 0;
    my $header_line = substr($_,1);
    my @delim_header =  split(/\|| |\t/, $header_line);
    my $fasta_id = $delim_header[0];
    if ($fasta_id eq $search_id) {
      $sequence_found = 1;
      $header_found_count++;
      if ($header_found_count > 1){
	die "More than one header matching $search_id\n";
      }
      $collect_sequence = 1;
    }
  }
  elsif ($collect_sequence) {
    $sequence_string = $sequence_string.$_;
  }
}
if ($sequence_found==0) {
  die "No sequence found matching: $search_id\n";
}
if ($range_string) {
  my @range = split(/-|,/,$range_string);
  if (scalar(@range) != 2) {
    die "Invalid range: $range_string\nRange must be two numeric values delimited by a dash or comma.\n";
  }
  else {
    if ($range[0]<$range[1]){
      $range_start=$range[0];
      $range_end=$range[1];
      print ">".$search_id."_from_position_".$range_start."_to_".$range_end."\n";
      $sequence_string = substr($sequence_string,$range_start-1,$range_end-($range_start-1));
    }
    else {
      $range_start=$range[1];
      $range_end=$range[0];
      print ">".$search_id."_from_position_".$range_start."_to_".$range_end."_reverse_complement\n";
      $sequence_string = substr($sequence_string,$range_start-1,$range_end-($range_start-1));
      my $rev_sequence_string = reverse($sequence_string);
      $rev_sequence_string =~ tr/ACGTacgt/TGCAtgca/;
      $sequence_string=$rev_sequence_string;
    }
  }
}
else {
  print ">$search_id\n";
}
#limit output to 60 bases per line
my @line_matches = ($sequence_string =~ /.{1,60}/g);
foreach my $output_line (@line_matches) {
  print "$output_line\n";
}
#print "$sequence_string\n";
sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This script returns a fasta sequence from a multi-fasta file, and optionally returns a subsequence.

    Usage:
     
      seq_from_fasta.pl -i <input_fasta_file> -n <header_name_of_sequence> -r <range-begin>,<range-end>

    Flags:

      -i <input_fasta_file>         input fasta file (mandatory)
      -n <header_name_of_sequence>  name of sequence to return
      -r <rage-begin>,<range-end>   start and end (inclusive) of sequence range to return (optional)
                                    sequence will be reverse complement when end coordinate is smaller than start
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


