#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use feature 'state';

use Getopt::Std;
use Pod::Usage;

use Bio::FeatureIO;
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
my $seq_index = my_seq_index->new(
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
                my $feature_type = $output_type eq 'protein' ? 'CDS' : $output_type;
                my $seq = splice_feature_seqs( $mrna, $ref_seq, $feature_type );
                $seq = $seq->translate if $output_type eq 'protein';
                $seq_out->write_seq( $seq );
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

    my $id = $feature->has_tag('Name') ? ($feature->get_tag_values('Name'))[0] : $feature->primary_id;
    $id =~ s/^[^:]+://;

    my $seq = Bio::PrimarySeq->new(
        -id  => $id,
        -seq => join( '', map $ref_seq->subseq( $_->start, $_->end ), @sub ),
      );
    $seq = $seq->revcom if $feature->strand < 0;
    return $seq;
}

BEGIN {
    package my_seq_index;
    use base 'Bio::Index::Fasta';

    sub _get_SeqIO_object {
        my( $self, $i ) = @_;

        unless ($self->{'_seqio_cache'}[$i]) {
            my $fh = $self->_file_handle($i);
            # make a new SeqIO object
            my $seqio = Bio::SeqIO->new(
                -format => $self->_file_format,
                -fh     => $fh,
                -seqfactory => Bio::Seq::SeqFactory->new( -type => 'Bio::Seq::LargePrimarySeq' ),
              );
            $self->{'_seqio_cache'}[$i] = $seqio;
        }
        return $self->{'_seqio_cache'}[$i];
    }
}

__END__

=head1 USAGE

genemodels_to_seq.pl genes.gff3 reference_seqs.fasta CDS|exon|protein

=cut
