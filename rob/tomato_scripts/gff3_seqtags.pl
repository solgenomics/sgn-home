#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use autodie ':all';

use File::Basename;
use File::stat;

use Bio::Index::Fasta;

my $tag_halfsize = 40;

my ( $gff3_file, $ref_seqs_file ) = @ARGV;
die 'must give gff3 file and ref seqs file' unless @ARGV == 2;

my $seq_index = open_index( $ref_seqs_file );

open my $gff_fh, '<', $gff3_file;

$| = 1;

my $line_number = 0;
my $current_seq_id;
my $current_seq;
my $current_seq_length;
while( my $line = <$gff_fh> ) {
    $line_number++;
    next if $line =~ /^#/;
    my @f = split /\t/, $line, 9;
    chomp $f[-1];

    for my $coord ( @f[3,4] ) {
        print  ">$f[1] $f[2] $f[8]\n";
        my $seq;
        if( $current_seq_id && $current_seq_id eq $f[0] ) {
            $seq = $current_seq;
        } else {
            my $s = $seq_index->fetch( $f[0] )
                or die "ref seq not found for $f[0]";

            $current_seq_id = $s->id;
            $current_seq    = $s->seq;
            $current_seq_length = length $current_seq;
        }

        my $tag_start  = clamp( $current_seq_length, $coord-1-$tag_halfsize );
        my $tag_length = clamp( $current_seq_length, $coord-1+$tag_halfsize ) - $tag_start;
        #print "== $coord => $tag_start, $tag_length\n";
        print substr( $current_seq, $tag_start, $tag_length ), "\n";
    }
}

# 0-based string coordinate clamp
sub clamp {
    my ( $seqlen, $coord ) = @_;

    return 0 if $coord < 0;

    return $seqlen - 1 if $coord >= $seqlen;

    return $coord;
}


sub open_index {
    my ( $file ) = @_;

    my $index_file = basename( $file ).'.index';

    my $inx = Bio::Index::Fasta->new( -filename   => $index_file,
                                      -write_flag => 1 );
    $inx->make_index( $file );

    return $inx;
}
