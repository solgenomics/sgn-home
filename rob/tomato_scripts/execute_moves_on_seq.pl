#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use autodie ':all';

use Bio::SeqIO;

my ( $seqfile, $movesfile ) = @ARGV;

my %moves = read_moves( $movesfile );

my $seq_in  = Bio::SeqIO->new( -format => 'fasta', -file => $seqfile );
my $seq_out = Bio::SeqIO->new( -format => 'fasta', -fh => \*STDOUT, -width => 80 );

while( my $seq = $seq_in->next_seq ) {

    my $id = $seq->id;
    my $seqstr = $seq->seq;

    # execute the moves on it
    my $moves = $moves{$id} or die "no moves for ".$seq->id;
    foreach my $move ( @$moves ) {
        do_move( $move, \$seqstr );
    }

    $id =~ s/2\.3/2.4/;

    $seq->id( $id );
    $seq->seq( $seqstr );

    $seq_out->write_seq( $seq );
}

exit;

############################3

sub do_move {
    my ( $move, $seq ) = @_;

    my ( $location, $offset ) = @$move;

    if( $offset > 0 ) {
        substr( $$seq, $location-1, 0, 'N'x$offset );
    }
    elsif( $offset < 0 ) {
        substr( $$seq, $location-1+$offset, -$offset, '' );
    }
}

sub read_moves {
    my ( $f ) = @_;
    open( my $fh, $f );
    my %m;
    while( <$fh> ) {
        chomp;
        my @f = split;
        my $move = shift @f;
        $move eq 'move' or die;
        my $id = shift @f;
        push @{$m{$id}}, \@f;
    }

    #use Smart::Comments;
    ### moves: %m

    return %m;
}
