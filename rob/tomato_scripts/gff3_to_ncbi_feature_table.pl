#!/usr/bin/env perl
use strict;
use warnings;

use lib '/home/rob/dev/bioperl/Bio-FeatureIO/lib';

use IO::File;

use Bio::FeatureIO;

my $in = Bio::FeatureIO->new(
    -format  => 'gff',
    -version => 3,
    -file    => $ARGV[0],
    );

while( my @f = $in->next_feature_group ) {
    for my $gene ( @f ) {
        my $fh = get_output_fh( $gene );

        $gene->primary_tag eq 'gene'
            or die "top-level feature not a gene?";

        $fh->print_tab(
            [ $gene->start, $gene->end, 'gene' ],
            format_gene_attributes( $gene ),
        );

        for my $mrna ( $gene->get_SeqFeatures ) {
            $mrna->primary_tag eq 'mRNA'
                or die "second-level feature not an mRNA?";

            my @exon_tabs =
                map [ $_->start, $_->end, undef ],
                grep $_->primary_tag eq 'exon',
                $mrna->get_SeqFeatures;
            $exon_tabs[0][2] = 'mRNA';

            push @exon_tabs,
                format_mrna_attributes( $mrna );

            $fh->print_tab( @exon_tabs );

            my @cds_tabs =
                map [ $_->start, $_->end, undef ],
                grep $_->primary_tag eq 'CDS',
                $mrna->get_SeqFeatures;
            $cds_tabs[0][2] = 'CDS';

            push @cds_tabs,
                format_cds_attributes( $mrna );

            $fh->print_tab( @cds_tabs );
        }
    }
}

sub format_mrna_attributes {
    my ( $mrna ) = @_;
    my ( $id ) = $mrna->get_tag_values( 'Name' );

    my @attributes = (
        [ protein_id    => $id ],
        [ transcript_id => $id ],
    );

    if( my ( $desc ) = eval { $mrna->get_tag_values( 'Note' ) }) {
        if( my ( $product ) = $desc =~ /^(.+)(?=\(A[A-Z]+)/ ) {
            push @attributes, [ product => $product ];
        }
        if( my @interpro = $desc =~ /\b(IPR\d+)\b/g ) {
            push @attributes, map [ dbxref => "InterPro:$_" ], @interpro;
        }

        push @attributes, [ note => $desc ];
    }

    if( my @terms = eval { $mrna->get_tag_values( 'Ontology_term' ) }) {
        push @attributes, map [ dbxref => $_ ], @terms;
    }

    return format_attributes( \@attributes );
}

sub format_cds_attributes {
    my ( $mrna ) = @_;
    return format_mrna_attributes( $mrna );
}
sub format_gene_attributes {
    my ( $gene ) = @_;
    return format_attributes( [ [ locus_id => $gene->get_tag_values('Alias') ]  ]);
}

sub format_attributes {
    my ( $attributes ) = @_;

    return
        map [ undef, undef, undef, @$_ ],
        sort { $a->[0] cmp $b->[0] }
        @$attributes;
}

my %filehandle_cache;
sub get_output_fh {
    my ( $feat ) = @_;
    my $seq_id = $feat->seq_id or die "no seq_id on feature!";
    return $filehandle_cache{$seq_id} ||= do {
        my $fh = mytabfile->new( "$seq_id.tbl", '>' )
            or die "$! writing $seq_id.tbl";

        $fh->print( ">Feature $seq_id\n" );
        $fh
    };
}


package mytabfile;
use base 'IO::File';

sub print_tab {
    my $self = shift;
    no warnings 'uninitialized';
    $self->print( join( "\t", @$_ ), "\n" ) for @_;
}

