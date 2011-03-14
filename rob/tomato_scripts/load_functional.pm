package load_functional;
use strict;
use warnings;
use autodie ':all';

use base 'Exporter';

our @EXPORT = qw/
    load_functional_files
/;

sub load_functional_files {
    my ( $specific_csv, $general_csv, $descriptions_file ) = @_;

    my %func;
    # load the go terms
    load_kklee_go_file( $specific_csv, \%func, 'specific' );
    load_kklee_go_file( $general_csv,  \%func, 'general'  );
    combine_go( \%func );

    # load the human-readable description
    load_descriptions( $descriptions_file, \%func );

    return \%func;
}

sub load_descriptions {
    my ( $desc_file, $func ) = @_;
    open my $desc, '<', $desc_file;
    while( my $line = <$desc> ) {
        chomp $line;
        my @cols = split /";"/, $line;
        my $prot_id = $cols[1];
        $prot_id =~ s/\.\d//;
        my $desc = $cols[4];
        $func->{$prot_id}{description} = $desc;
    }

}

sub combine_go {
    my ( $go ) = @_;

    for my $prot ( keys %$go ) {
        my $rec = $go->{$prot};
        $rec->{combined} = $rec->{specific} || $rec->{general};
    }
}

sub load_kklee_go_file {
    my ( $file, $go, $type ) = @_;

    open my $fh, '<', $file;
    while ( my $line = <$fh> ) {
        chomp $line;
        next if $line =~ /No GOs/;
        my ( $protein_id, @terms) = split /\s+/, $line;
        $protein_id =~ s/\.\d//;
        s/,+$// for @terms;
        push @{ $go->{$protein_id}{$type} ||= [] }, @terms
    }
}

1;

