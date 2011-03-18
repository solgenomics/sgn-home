#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';

use FindBin;
use lib "$FindBin::RealBin";

use load_functional 'load_functional_files';

my ( $specific_csv, $general_csv, $descriptions_file, $fasta_file ) = @ARGV;

my $func = load_functional_files( $specific_csv, $general_csv, $descriptions_file );
# hashref of protein_id => {
#   specific => [ go, go, ...],
#   general  => [ go, go, ...],
#   combined => [ go, go, ...],
#   description => 'functional description',
# }


open my $fasta_in, '<', $fasta_file;
while( <$fasta_in> ) {
    if( my ($id) = /^>\s*(\S+)/ ) {
        $id =~ s/\.\d//;
        if( my $terms = $func->{ $id }{combined} ) {
            chomp;
            s/\s+go_terms:\S+//;
            $_ .= " go_terms:".join(',', @$terms )."\n";
        }
        if( my $desc = $func->{ $id }{description} ) {
            chomp;
            s/\s+(functional_)?description:"[^"]+"//;
            $desc =~ s/"/'\"'/ge;
            $_ .= qq| functional_description:"$desc"\n|;
        }
    }
    print;
}

exit;
