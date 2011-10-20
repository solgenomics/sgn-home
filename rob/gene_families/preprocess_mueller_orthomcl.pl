#!/usr/bin/env perl
use strict;
use warnings;

use Memoize;
memoize('cvterm');

use List::MoreUtils 'uniq';
use Bio::Chado::Schema;

our $schema = Bio::Chado::Schema->connect(shift);

while(<>) {

    # convert peptide names to gene names, with varying degrees of
    # difficulty depending on the organism.  also remove any genes
    # that are not already found in the db
    s/ (PGSC\d{4}DMP\d+) (\([^\)]+\) \s+) /gene_for_polypeptide($1,$2)/egx; #potato
    s/ GSVIVT(\d+) (\([^\)]+\) \s+)       /gene_for_polypeptide("GSVIV12X_T$1",$2)/egx; #grape
    s/ (Os\d+)t(\d+)\-\d+                 /confirm_gene($1.'g'.$2)/exg; #rice
    s/ (Solyc\d+g\d+\.\d+)(\.\d+)+        /confirm_gene($1)/exig; # tomato
    s/ (AT\d+G\d+)(\.\d+)+                /confirm_gene($1)/exig; # arabidopsis

    # correct the gene and taxa counts in case we removed some.
    my ( $cluster, @genes ) = split;
    my @taxa = uniq( /\(([^\)]+)\)/g );
    s/\(\d+ genes,\d+ taxa\)/'('.@genes.' genes,'.(@taxa-1).' taxa)'/e or die;

    print;
}

sub confirm_gene {
    my $gene_name = shift;

    return $gene_name
      if $schema->resultset('Sequence::Feature')
                ->find({
                    name => $gene_name,
                    type => cvterm( sequence => 'gene' ),
                  });

    die "gene $gene_name not found\n";
}

sub gene_for_polypeptide {
    my ( $protein_name, $trailing_text ) = @_;
    $trailing_text ||= '';

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
        die "multiple genes found for polypeptide $protein_name: ".join(", ", map $_->name.' ('.$_->feature_id.')', @genes)."\n";
    }
    elsif( !@genes ) {
        warn "WARNING, no gene found for $protein_name, removing this member\n";
        return '';
        #die "no genes found for polypeptide $protein_name\n";
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

