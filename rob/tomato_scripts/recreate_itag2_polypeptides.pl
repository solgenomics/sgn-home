#!/usr/bin/env perl
use strict;
use warnings;

use DBI;
use Bio::Chado::Schema;

my $schema = Bio::Chado::Schema->connect(@ARGV);
#$dbh->{RaiseError} = 1;
#$dbh->begin_work;

sub echo(@) {
    print @_;
    return @_;
}

# make polypeptide data for all the mRNAs
my $polys = $schema->storage->dbh->prepare( echo <<'');
select f.organism_id, f.name, replace(f.uniquename,'mRNA','polypeptide') as uniquename, 22027 as type_id, f.is_analysis, f.is_obsolete, fl.srcfeature_id,fl.fmin,fl.fmax,fl.strand,fl.phase,fl.residue_info,fl.locgroup,fl.rank, f.feature_id as mrna_id from feature f join featureloc fl using(feature_id) join analysisfeature using(feature_id) join analysis using(analysis_id) where analysis.name='ITAG_renaming' and f.type_id =  (select cvterm_id from cvterm where name = 'mRNA')

$polys->execute;

$schema->txn_do(sub {

    while( my $p = $polys->fetchrow_hashref ) {
        my $f = $schema->resultset('Sequence::Feature')->create({
            organism_id => $p->{organism_id},
            name => $p->{name},
            uniquename => $p->{uniquename},
            type_id => $p->{type_id},
            is_analysis => $p->{is_analysis},
            is_obsolete => $p->{is_obsolete},
        });

        #     # insert new feature
        #     my ($new_poly) = $dbh->selectrow_array( echo(<<''), undef, @{$poly_data}[0..5] );
        # insert into feature ( organism_id, name, uniquename, type_id, is_analysis, is_obsolete ) values ( ?, ?, ?, ?, ?, ? )

        $f->create_related('featureloc_features', {
            ( map { $_ => $p->{$_} }
          qw(
                srcfeature_id
                fmin
                fmax
                strand
                phase
                residue_info
                locgroup
                rank
            )
          )
            });

        #     # insert featureloc
        #     $dbh->do( echo(<<''), undef, $new_poly, @{$poly_data}[6..13] );
        # insert into featureloc (feature_id, srcfeature_id, fmin, fmax, strand, phase, residue_info, locgroup, rank) values ( ?, ?, ?, ?, ?, ?, ?, ? )

        $f->create_related('feature_relationship_subjects', {
            object_id => $p->{mrna_id},
            type_id => 21910,
            rank => 0,
        });

        #     # insert feature_relationship
        #     $dbh->do( echo(<<''), undef, $poly_data->[14], $new_poly, 21910, 0 );
        # insert into feature_relationship ( subject_id, object_id, type_id, rank )
        #                           values ( ?, ?, ?, ? )


    }
});

#$dbh->rollback;
#$dbh->disconnect;
