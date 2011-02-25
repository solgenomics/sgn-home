#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use feature 'state';

use Getopt::Std;
use Pod::Usage;

use Bio::FeatureIO;
use Bio::Index::Fasta;
use Bio::SeqIO;
use Bio::PrimarySeq;

use Smart::Comments;

@ARGV == 3 or pod2usage();

my $gff3_in = Bio::FeatureIO->new(
    -format  => 'gff',
    -version => 3,
    -file    => $ARGV[0],
  );

my $seq_out = Bio::SeqIO->new(
    -format => 'fasta',
    -fh     => \*STDOUT,
  );

my $tempfile = File::Temp->new;
my $seq_index = Bio::Index::Fasta->new(
    -filename   => "$tempfile",
    -write_flag => 1 );
$seq_index->make_index( $ARGV[1] );

my $output_type = $ARGV[2];

while( my @toplevel_features = $gff3_in->next_feature_group ) {
    for my $feature (@toplevel_features) {
        if( $feature->primary_tag eq 'gene' ) {
            # get the mRNA subfeatures
            my @mrnas =
                grep $_->primary_tag eq 'mRNA',
                $feature->get_SeqFeatures;

            for my $mrna ( @mrnas ) {
                my $ref_seq = _fetch( $mrna->seq_id );
                $seq_out->write_seq(
                    $output_type eq 'protein'
                        ? protein_seq( $mrna, $ref_seq )
                        : splice_feature_seqs( $mrna, $ref_seq, $output_type )
                  );
            }
        }
    }
}


sub _fetch {
    my $id = shift;
    state %seq_cache;
    return $seq_cache{ $id } if $seq_cache{$id};

    %seq_cache = ();

    my $seq = $seq_index->fetch( $id )
        or die "reference sequence $id not found in index\n";
    return $seq_cache{$id} = $seq;
}

sub splice_feature_seqs {
    my ( $feature, $ref_seq, $subfeature_type ) = @_;

    my @sub = grep $_->primary_tag eq $subfeature_type, $feature->get_SeqFeatures
        or die "no $subfeature_type subfeatures for ".($feature->get_tag_values('Name'))[0];

    my $seq = Bio::PrimarySeq->new(
        -id => ($feature->get_tag_values('Name'))[0],
        -seq => join( '', map $ref_seq->subseq( $_->start, $_->end ), @sub ),
      );
    $seq = $seq->revcom if $feature->strand < 0;
    return $seq;
}

sub protein_seq {
    my ( $feature, $ref_seq ) = @_;

    my @cds = grep $_->primary_tag eq 'CDS', $feature->get_SeqFeatures
        or die "no CDS subfeatures for ".($feature->get_tag_values('Name'))[0];

    my @nuc_seqs = map {
        my $f = $_;
        my $s = $ref_seq->trunc( $f->start, $f->end );
        if( my $p = $f->phase ) {
            if( $f->strand < 0 ) {
                $s->subseq( 1, $s->length - $p );
            } else {
                $s->subseq( $p+1, $s->length );
            }
        } else {
            $s->seq
        }
    } @cds;

    my $seq = Bio::PrimarySeq->new(
        -id => ($feature->get_tag_values('Name'))[0],
        -seq => join( '', @nuc_seqs ),
      );

    $seq = $seq->revcom if $feature->strand < 0;

    return $seq->translate;
}

__END__

=head1 USAGE

genemodels_to_seq.pl genes.gff3 reference_seqs.fasta CDS|exon|protein

=cut
