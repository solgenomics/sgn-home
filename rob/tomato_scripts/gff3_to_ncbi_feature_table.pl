#!/usr/bin/env perl
use strict;
use warnings;
use 5.10.0;

use lib '/home/rob/dev/bioperl/Bio-FeatureIO/lib';

use IO::File;

use Bio::FeatureIO;

my $in = Bio::FeatureIO->new(
    -format  => 'gff',
    -version => 3,
    -file    => $ARGV[0],
    );

my %handlers = (
    gene => \&handle_gene,
);
while( my @f = $in->next_feature_group ) {
    for my $feature ( @f ) {
        my $fh = get_output_fh( $feature )
            or next;

        no strict 'refs';

        my $type = $feature->primary_tag;
        my $handler = $handlers{$type}
            or die "no handler defined for top-level $type features";
        $handler->( $fh, $feature );
    }
}

exit;

################ subroutines ################

sub handle_gene {
    my ( $fh, $gene ) = @_;

    $gene->primary_tag eq 'gene'
        or die "not a gene?";

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


sub format_mrna_attributes {
    my ( $mrna ) = @_;
    my ( $locus_tag ) = $mrna->get_tag_values( 'Name' );
    $locus_tag =~ s/(\.\d+)/.x/;

    my $prefix = "gnl|WGS:AEKE|";
    my @attributes = (
        [ protein_id    => $prefix.$locus_tag ],
        [ transcript_id => $prefix."mRNA:$locus_tag" ],
    );

    if( my ( $desc ) = eval { $mrna->get_tag_values( 'Note' ) }) {
        if( my ( $product ) = $desc =~ /^(.+)(?=\(A[A-Z]+)/ ) {
            $product =~ s/^\s+|\s+$//g;
            push @attributes, [ product => $product ];
        }

        if( my @interpro = $desc =~ /\b(IPR\d+)\b/g ) {
            push @attributes, map [ dbxref => "InterPro:$_" ], @interpro;
            push @attributes, map [ inference => "protein motif:InterPro:$_" ], @interpro;
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
    return format_attributes([
        map {
            ( my $tag = $_ ) =~ s/Solyc/Solyc_/;
            [ locus_tag => $tag ]
        } $gene->get_tag_values('Alias')
    ]);
}

sub format_attributes {
    my ( $attributes ) = @_;

    return
        map [ undef, undef, undef, @$_ ],
        sort { $a->[0] cmp $b->[0] }
        @$attributes;
}


use Data::Dumper;

my %filehandle_cache;
sub get_output_fh {
    my ( $feat ) = @_;
    my $seq_id = $feat->seq_id or die "no seq_id on feature!";

    return $filehandle_cache{$seq_id} //= do {

        my $accession = {qw(
            SL2.40ch01	CM001064.1
            SL2.40ch02	CM001065.1
            SL2.40ch03	CM001066.1
            SL2.40ch04	CM001067.1
            SL2.40ch05	CM001068.1
            SL2.40ch06	CM001069.1
            SL2.40ch07	CM001070.1
            SL2.40ch08	CM001071.1
            SL2.40ch09	CM001072.1
            SL2.40ch10	CM001073.1
            SL2.40ch11	CM001074.1
            SL2.40ch12	CM001075.1
        )}->{$seq_id};

        unless( $accession ) {
            warn "I don't know the genbank accession for $seq_id, skipping all features on that sequence";
            0
        } else {
            my $fh = mytabfile->new( "$accession.tbl", '>' )
                or die "$! writing $accession.tbl";

            $fh->print( ">Feature $accession\n" );
            $fh
        }
    };
}


package mytabfile;
use base 'IO::File';

sub print_tab {
    my $self = shift;
    no warnings 'uninitialized';
    $self->print( join( "\t", @$_ ), "\n" ) for @_;
}

