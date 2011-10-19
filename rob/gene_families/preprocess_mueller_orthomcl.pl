#!/usr/bin/env perl
use strict;
use warnings;

use Memoize;
memoize('cvterm');

use Bio::Chado::Schema;

our $schema = Bio::Chado::Schema->connect(shift);

while(<>) {
    s/ (PGSC\d{4}DMP\d+) (\([^\)]+\) \s+)/pgsc_protein_to_gene($1,$2)/egx;
    print;
    last;
}


sub pgsc_protein_to_gene {
    my ( $protein_name, $trailing_text ) = @_;

    my ( $gene ) = my @genes = $schema->resultset('Sequence::Feature')
           ->search({
               'me.name' => $protein_name,
               'me.type_id' => cvterm(sequence => 'polypeptide')->cvterm_id,
              })
           ->search_related('feature_relationship_subjects', {
               'feature_relationship_subjects.type_id' => cvterm( relationship => 'derives_from')->cvterm_id,
             })
           ->search_related('object')
           ->search_related('feature_relationship_subjects', {
               'feature_relationship_subjects_2.type_id' => cvterm( relationship => 'part_of' )->cvterm_id,
             })
           ->search_related('object', {
               'object_2.type_id' => cvterm( sequence => 'gene' )->cvterm_id,
             });

    if( @genes > 1 ) {
        die "multiple genes found for potato polypeptide $protein_name: ".join(", ", map $_->name.' ('.$_->feature_id.')', @genes)."\n";
    }
    elsif( !@genes ) {
        warn "WARNING, no gene found for $protein_name, removing this member\n";
        return '';
        #die "no genes found for potato polypeptide $protein_name\n";
    }

    return $gene->name.$trailing_text;
}


sub cvterm {
    my ( $cv, $term ) = @_;
    my $t = $schema->resultset('Cv::Cv')
                   ->search( { 'me.name' => $cv })
                   ->search_related( 'cvterms', { 'cvterms.name' => $term } )
                   ->single;
    $t or die "Cannot find required cvterm in $cv ontology called $term.  Is this term in the $cv ontology, and is that ontology loaded?";
    return $t;
}

