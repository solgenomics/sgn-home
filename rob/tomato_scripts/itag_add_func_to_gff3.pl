#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use autodie ':all';

use FindBin;
use lib "$FindBin::RealBin";

use load_functional 'load_functional_files';

use Bio::GFF3::LowLevel 'gff3_parse_attributes','gff3_format_attributes';

my ( $specific_csv, $general_csv, $descriptions_file, $gff3_file ) = @ARGV;

my $func = load_functional_files( $specific_csv, $general_csv, $descriptions_file );
# hashref of protein_id => {
#   specific => [ go, go, ...],
#   general  => [ go, go, ...],
#   combined => [ go, go, ...],
#   description => 'functional description',
# }

open my $gff3_fh, '<', $gff3_file;
while( <$gff3_fh> ) {
    next if /^#/;
    next unless m{ \t ITAG_eugene \t mRNA \t }x;
    s/([^\t]+)$/add_func_attrs($func,$1)/e;
} continue {
    print;
}

sub add_func_attrs {
    my ( $func, $attr_string ) = @_;

    my $attr = gff3_parse_attributes( $attr_string )
        or die "failed to parse attrs '$attr_string'";

    ( my $transcript_no_version = $attr->{Name}[0] ) =~ s/\.\d+//;

    if( my $terms = $func->{ $transcript_no_version }{combined} ) {
        $attr->{Ontology_term} = $terms;
    }
    if( my $terms = $func->{ $transcript_no_version }{general} ) {
        $attr->{interpro2go_term} = $terms;
    }
    if( my $terms = $func->{ $transcript_no_version }{specific} ) {
        $attr->{sifter_term} = $terms;
    }


    $attr->{Note} = [ $func->{ $transcript_no_version }{description} ||die ];

    return gff3_format_attributes( $attr )."\n";
}
