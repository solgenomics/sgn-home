use strict;
use warnings;

use Bio::Chado::Schema;

my $schema = Bio::Chado::Schema->connect(shift);

my $features = $schema->resultset('Sequence::Feature');
my $synonyms = $schema->resultset('Sequence::Synonym');

my $polypeptide_type =
    $schema->resultset('Cv::Cv')
    ->search({ 'me.name' => 'sequence' })
    ->search_related('cvterms', { 'cvterms.name' => 'polypeptide' } )
    ->single or die;

open( my $names_fh, '<', shift ) or die;

$schema->txn_do(sub {

    while( <$names_fh> ) {
	chomp;
	my ( $mrna, $poly ) = split;
	#print "$mrna => $poly\n";
	my @poly_feats = 
	    $features
	    ->search({ name => $poly,
		       type_id => $polypeptide_type->cvterm_id });
	unless( @poly_feats ) {
	    print "$mrna polypeptide not found\n";
	    next;
	}
	if( @poly_feats > 1 ) {
	    print "multiple $mrna polys found: ",( join ", ", map $_->feature_id, @poly_feats),"\n";
	}
	for my $poly_feat (@poly_feats) {
	    my $synonym = $synonyms->find_or_create({ name => $mrna, type_id => 5, synonym_sgml => $mrna });
	    
	    $poly_feat->create_related('feature_synonyms', { pub_id => 1, synonym => $synonym});
	}
    }
    #die 'dry run';
});	     
