#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use POSIX 'ceil';

use lib '/home/rob/dev/bioperl/Bio-FeatureIO/lib';
use Bio::SeqIO;
use Bio::FeatureIO;

my $cdna_seqs = Bio::SeqIO->new( -file => shift, -format => 'fasta' );
my $features_in  = Bio::FeatureIO->new( -file => shift, -format => 'gff', -version => 3 );
my $features_out = Bio::FeatureIO->new( -fh => \*STDOUT, -format => 'gff', -version => 3 );

my $trims = find_trim_regions( $cdna_seqs );

while( my @f = $features_in->next_feature_group ) {
    for my $feature ( @f ) {
        my ( $name ) = eval { $feature->get_tag_values('Name') };
        if( my $trim = $name && $trims->{$name} ) {
            #use Data::Dump 'dump';
            #print "# TRIMMING $name: ".dump( $trim )."\n";
            $feature = undef unless trim_feature( $feature, $trim );
        }
        $features_out->write_feature( $feature ) if $feature;
    }
}

############## subs ###############

sub trim_feature {
    my ( $feature, $trim ) = @_;

    $feature->primary_tag eq 'gene'
        or die "don't know how to handle ".$feature->primary_tag." feature ".$feature->get_tag_values('Name');

    my ( $start, $end ) = ( $feature->start, $feature->end );

    my ( $trim_start, $trim_end ) = ( $trim->{"5'"}, $trim->{"3'"} );
    # round the trim start up to the nearest codon to avoid frameshifts
    $trim_start = ceil($trim_start/3)*3 if $trim_start;
    ( $trim_start, $trim_end ) = ( $trim_end, $trim_start ) if $feature->strand == -1;

    $start += $trim_start if $trim_start;
    $end   -= $trim_end   if $trim_end;

    #print "# TRIMMING to bounds $start, $end\n";

    return recursive_trim_to_bounds( $feature, $start, $end );
}

sub recursive_trim_to_bounds {
    my ( $feature, $start, $end ) = @_;

    # trim the subfeatures, deleting any that disappear
    $feature->add_SeqFeature( $_ )
        for grep recursive_trim_to_bounds( $_, $start, $end ),
            $feature->remove_SeqFeatures;

    $feature->start( $start ) if $start && $feature->start < $start;
    $feature->end( $end )     if $end   && $feature->end   > $end;

    return $feature->start < $feature->end;
}

sub find_trim_regions {
    my ( $cdna_seqs ) = @_;
    my %trims;
    while ( my $cdna = $cdna_seqs->next_seq ) {
        my $seq = $cdna->seq;
        ( my $gene_id = $cdna->id ) =~ s/\.\d+$//;
        # leading Ns
        if ( $seq =~ /^(N+)/i ) {
            $trims{$gene_id}{"5'"} = length $1;
        }

        # trailing Ns
        if ( $seq =~ /(N+)$/i ) {
            $trims{$gene_id}{"3'"} = length $1;
        }
    }

    return \%trims;
}
