#!/usr/bin/env perl
use strict;
use warnings;

use Pod::Usage;

my $index = shift;
$index && $index =~ /^\d+$/
    or pod2usage();

my %cnt;
my $printing = 1;
while(<>) {
    if( />\s*(\S+)/ ) {
        my $cnt = ++$cnt{$1};
        if( $cnt == $index ) {
            $printing = 1;
        } else {
            $printing = 0;
        }
    }

    print if $printing;
}


=head1 NAME

fasta_dedup.pl - extract first, second, etc occurrences of a given
ident in a fasta file

=head1 USAGE

fasta_dedup.pl index_num  file file file (or stdin)

=cut

