#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';

use Bio::SeqIO;

use File::Temp;

my ( $file ) = @ARGV;

my $dir = $file.'_split';
mkdir $dir;

my $in = Bio::SeqIO->new( -file => $file, -format => 'fasta' );
my $ct = 0;
while( my $s = $in->next_seq ) {
    Bio::SeqIO->new(
	-file   => sprintf('>%s/%05d.%s.fasta',$dir,$ct++,$s->id),
	-format => 'fasta',
	-width  => 80,
	)->write_seq( $s );
}
