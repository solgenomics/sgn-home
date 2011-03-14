#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use CXGN::BioTools::AGP qw/ agp_parse agp_format_part /;
use Bio::FeatureIO;

my ( $chr_from_scaffolds, $chr_from_contigs ) = @ARGV;

# read in all the features and sort by start coordinate
my @features = make_features_from_files( $chr_from_scaffolds, $chr_from_contigs );

# aggregate the features to make the contigs subfeatures of the scaffolds
my @aggregated_features;
for my $f ( @features ) {
    if( $f->primary_tag eq 'supercontig' ) {
        push @aggregated_features, $f;
    }
    elsif( $f->primary_tag eq 'contig' ) {
        my $scaff = $aggregated_features[-1];
        $scaff->add_SeqFeature( $f );

        $f->seq_id eq $scaff->seq_id
           && $f->start >= $scaff->start
           && $f->end <= $scaff->end
           or die "aggregation failed, contig ".$f->display_name." cannot be added as a subfeature of ".$scaff->display_name;
    }
    else {
        die "unknown feature";
    }
}

my $gff3_out = Bio::FeatureIO->new( -format => 'gff', -version => 3, -fh => \*STDOUT );
$gff3_out->write_feature( $_ ) for @aggregated_features;

##########  subroutines ########

sub make_features_from_files {
    my ( $scaffold_features, $contig_features ) =
        map agp_to_features( agp_parse( $_ ) ),
        @_;

    $_->primary_tag('supercontig') for @$scaffold_features;
    $_->primary_tag('contig')      for @$contig_features;

    return sort {
        $a->seq_id cmp $b->seq_id
          ||
        $a->start <=> $b->start
    } ( @$scaffold_features, @$contig_features );
}

sub agp_to_features {
    my $agp_lines = shift;
    my @features;
    for my $line ( @$agp_lines ) {
        push @features, agp_to_feature( $line );
    }

    return \@features;
}
sub agp_to_feature {
    my $line = shift;

    return if $line->{comment};

    if ( $line->{is_gap} ) {
        return;
    } else {
        my $overall_strand = ( $line->{orient} eq '-' ? -1 : 1 );
        return Bio::SeqFeature::Generic->new(
            -display_name => $line->{ident},
            -seq_id => $line->{objname},
            -start  => $line->{ostart},
            -end    => $line->{oend},
            -strand => $overall_strand,
            -primary_tag => 'contig',
            -source => 'SL2.40_assembly',
            -tag => {
                ID     => $line->{ident},
                Name   => $line->{ident},
                Target => "$line->{ident} $line->{cstart} $line->{cend} +",
                reliably_oriented => ( $line->{orient} eq '0' ? 0 : 1 ),
            },
          );
    }
}
