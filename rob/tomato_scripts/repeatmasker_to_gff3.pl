#!/usr/bin/env perl
use strict;
use warnings;

use Bio::FeatureIO;
use Bio::Tools::RepeatMasker;

my $g = Bio::FeatureIO->new(
    -format => 'gff',
    -version => 3,
    -fh => \*STDOUT,
  );

for my $file (@ARGV) {
    my $parser = Bio::Tools::RepeatMasker->new(-file => $file);
    while ( my $result = $parser->next_result ) {
        my $f = $result->feature1;

        # normalize the primary tags;
        my $type = $f->primary_tag;
        $f->primary_tag( 'repeat_region' );
        $f->add_tag_value('repeat_class',$type);

        #fix the Target attributes
        my @tgt = $f->remove_tag('Target');
        if ( $tgt[1] > $tgt[2] ) {
            @tgt[1,2] = @tgt[2,1];
        }
        $f->add_tag_value( 'Target', join ' ', @tgt, '+' );

        $g->write_feature( $f );
    }
}
