#!/usr/bin/perl

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

use strict;
use warnings;

use English;
use Carp;
use FindBin;

use Getopt::Std;
use Pod::Usage;

#use Data::Dumper;

our %opt;
getopts('v',\%opt) or pod2usage(1);
my $invert = $opt{v} ? 1 : 0;

my $pat = shift @ARGV;
$pat = qr/$pat/;


$/ = "\n>";

my $first = <>;
chomp $first;
print "$first\n" if $first =~ $pat xor $invert;

while(my $seq = <>) {
  chomp $seq;
  print ">$seq\n" if $seq =~ $pat xor $invert;
}


__END__

=head1 NAME

fasta_grep.pl - like 'grep', but for use on fasta-format sequence
files.

=head1 DESCRIPTION

Matches the perl regexp given on the command line to the
ENTIRE fasta record (ident, defline, sequence) and prints it if it
matches.

Idiosyncracy: Puts an extra newline at the end of the file if the last
seq in your input matches.

=head1 SYNOPSIS

  fasta_grep.pl  perl_regular_expression  list_of_files

  Options:

    -v invert match, show seqs that do NOT match

=head1 MAINTAINER

Robert Buels

=head1 AUTHOR(S)

Robert Buels, E<lt>rmb32@cornell.eduE<gt>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut
