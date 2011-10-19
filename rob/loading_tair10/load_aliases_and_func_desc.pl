#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dumper;
use Memoize;
use List::MoreUtils 'uniq';

use Bio::Chado::Schema;

my $dsn = shift;
my $schema = Bio::Chado::Schema->connect($dsn);
my $ara = $schema->resultset('Organism::Organism')
                         ->find( { species => 'Arabidopsis thaliana'})
     or die "could not find A.thaliana organism in database";

$schema->txn_do( sub {
    while(<>) {
        next unless /^>/;

        my $data = parse_line( $_ );
        load( $data );
    }

});


sub load {
    my ( $data ) = @_;
    my $feature = $ara->search_related('features', {
                      name => $data->{id},
                      type_id => get_cvterm( $schema, 'sequence', $data->{type} )->cvterm_id,
                  })->single
        or die "could not find feature with type $data->{type} and name $data->{id}";
    $feature->add_to_featureprops({
        type  => get_cvterm( $schema, undef, 'Note' ),
        value => $data->{description},
    });
    for my $alias ( uniq @{$data->{aliases}} ) {
	#warn "inserting alias $alias for $data->{id}\n";
        $schema->resultset('Sequence::Synonym')
               ->find_or_create({
                   type => get_cvterm( $schema, 'synonym_type', 'exact' ),
                   name => $alias,
		   synonym_sgml => $alias,
                 })
               ->add_to_feature_synonyms({ feature => $feature, pub_id => 1 });
    }
}

memoize('get_cvterm');
sub get_cvterm {
    my ( $schema, $cvname, $type_name ) = @_;
    my $rs = $schema->resultset('Cv::Cv');
    if( defined $cvname ) {
        $rs = $rs->search({ 'me.name' => $cvname })
    }
    return   $rs->search_related('cvterms', {
                    'cvterms.name'        => $type_name,
                    'cvterms.is_obsolete' => 0,
                   })
                ->single or die "cannot find cvterm $cvname:$type_name";
}

sub parse_line {
    local $_ = shift;

    chomp;
    # parse ID and SO type
    my ( $id, $type ) = /^>(\S+)\s*type:(\S+)/g or die;

    # parse the aliases
    my @aliases = map { split /[\s,]+/, $_ } ( /Symbols: ([^\|]+)/ );

    # parse the description
    s/\|[^\|]+$//;
    my ( $description ) = / Symbols: [^\|]+ \| \s* (.+)/x or die;
    $description =~ s/\s*$//;

    return {
        id => $id,
        type => $type,
        aliases => \@aliases,
        description => $description,
    };
}
