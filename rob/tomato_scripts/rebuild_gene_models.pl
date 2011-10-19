#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';
use feature 'say';

use lib '/home/rob/dev/bioperl/Bio-FeatureIO/lib';

use List::MoreUtils qw/ minmax any /;

use Bio::Range;
use Bio::FeatureIO;


my $rebuild_list = read_rebuild_list( $ARGV[1] );
#my $rebuild_list = { ('Solyc07g054880.2') x 2 };

process_gff3( $ARGV[0], \*STDOUT, $rebuild_list );

sub process_gff3 {
    my ( $gff3_file, $out_fh, $rebuild_list ) = @_;

    my $in = Bio::FeatureIO->new(
        -format  => 'gff',
        -version => 3,
        -file    => $gff3_file,
        );
    my $out = Bio::FeatureIO->new( -format => 'gff', -version => 3, -fh => $out_fh );

    while( my @f = $in->next_feature_group ) {
        for my $f ( @f ) {
            my $gene_name =  ($f->get_tag_values('Name'))[0];
            if( $gene_name && $rebuild_list->{$gene_name} ) {
                $f = rebuild_gene( $f );
            }
            $out->write_feature( $f );
        }
    }
}

sub rebuild_gene {
    my ( $gene ) = @_;
    my ( $mrna ) = my @mrnas = $gene->get_SeqFeatures;
    die 'multi mrnas??' if @mrnas > 1;
    die 'no mrna??' unless $mrna;

    my @features = $mrna->get_SeqFeatures;

    # recalculate the UTRs
    recalculate_utrs( $mrna, $gene );

    # just remove the intron subfeatures rather than fix them
    remove_introns( $mrna );

    # check mutual exclusivity of ( UTRs, CDS ), (exons, introns)
    # check containment of exon->UTR, exon->CDS
    return $gene;
}

sub range(@) {
    my ( $start, $end ) = minmax map { $_->start, $_->end } @_;
    return Bio::Range->new( -start => $start, -end => $end );
}

sub remove_introns {
    $_[0]->add_SeqFeature( $_ ) for grep $_->primary_tag ne 'intron', $_[0]->remove_SeqFeatures;
}

sub recalculate_utrs {
    my ( $mrna, $gene ) = @_;

    my @features = grep $_->primary_tag !~ /_UTR$/i,# && $_->primary_tag ne 'intron',
                   $mrna->remove_SeqFeatures;

    my @exons =
        sort { $a->start <=> $b->start }
        grep { $_->primary_tag eq 'exon' }
        @features;

    my @cds =
        grep { $_->primary_tag eq 'CDS' }
        @features;

    # merge any exons that place a CDS in an intron
    my @new_exons;
    for( my $i = 0; $i<@exons; $i++ ) {
        next unless $exons[$i];
        if( any { $_->contains( $exons[$i]->end + 1 ) || $exons[$i+1] && $_->contains( $exons[$i+1]->start - 1 ) } @cds ) {
            #warn "contains, merging";
            $exons[$i]->end( $exons[$i+1]->end );
            $exons[$i+1] = undef;
        }
        push @new_exons, $exons[$i];
    }

    # find transcription range
    my $cds_range =
       range
       grep $_->primary_tag eq 'CDS',
       @features;

    my @utrs =
       # split UTRs on introns if necessary
       map {
           my $u = $_;

           my @utrs =
             # upgrade the ranges to features
             map {
                 Bio::SeqFeature::Generic->new(
                     -start  => $_->start,
                     -end    => $_->end,
                     -strand => $u->strand,
                     -seq_id => $mrna->seq_id,
                     -primary_tag => $u->primary_tag,
                     -source_tag  => $mrna->source_tag,
                     )
             }
             # intersect with each exon (makes one or more Ranges)
             map  {
                 #warn "intersecting ".$u->gff_string." with ".$_->gff_string."\n";
                 my ( $start, $stop, $strand ) = $u->intersection($_);
                 if( defined $start ) {
                     #warn "makes $start, $stop, $strand\n";
                     Bio::Range->new(
                         -start  => $start,
                         -end    => $stop,
                         -strand => $strand,
                         );
                 } else {
                     ()
                 }
             }
             @new_exons
       }
       # make UTR features
       map Bio::SeqFeature::Generic->new( %$_ ),
       {
           -start  => $mrna->start,
           -end    => $cds_range->start - 1,
           -strand => $mrna->strand,
           -primary_tag => ( $mrna->strand > 0 ? 'five_prime_UTR' : 'three_prime_UTR' ),
       },
       {
           -start  => $cds_range->end + 1,
           -end    => $mrna->end,
           -strand => $mrna->strand,
           -primary_tag => ( $mrna->strand > 0 ? 'three_prime_UTR' : 'five_prime_UTR' ),
       };

    for my $feature( @utrs, @new_exons, @cds, $mrna, $gene ) {
        for my $tagname ('ID','Name') {
            if( $feature->has_tag( $tagname ) ) {
                my @ids = $feature->remove_tag( $tagname );
                s/\.(\d+)/'.'.($1+1)/e for @ids;
                $feature->add_tag_value( $tagname, $_ ) for @ids;
            }
        }
    }

    $mrna->add_SeqFeature($_)
       for sort {$a->start <=> $b->start || $a->primary_tag cmp $b->primary_tag }
           @utrs, @new_exons, @cds;

    return $mrna;
}

sub read_rebuild_list {
    my $file = shift;
    open my $f, '<', $file or die "$! reading $file";
    my %list;
    while(<$f>) {
        chomp;
        s/\.\d+$//;
        $list{$_} = 1;
    }
    return \%list;
}
