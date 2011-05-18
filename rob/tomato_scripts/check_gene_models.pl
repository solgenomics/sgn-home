#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';

use Bio::Range;
use Bio::GFF3::LowLevel 'gff3_parse_feature';

# check for CDS's that are not entirely contained within exons
my %genes;
while(<>) {
    next if /^#/ || ! /\t(exon|CDS)\t/;

    my ( $gene ) = /(Solyc\d\dg\d{6})/ or die "no gene ident found in $_";
    my $f = gff3_parse_feature($_);
    my $r = Bio::Range->new( -start => $f->{start}, -end => $f->{end}, -strand => $f->{strand} eq '-' ? -1 : 1 );

    push @{ $genes{$gene}{ $f->{type}||die } }, $r;
}

#use Data::Dump 'dump';
#warn dump( \%genes );

for my $gene_name ( sort keys %genes ) {
    my $gene = $genes{$gene_name};
    for my $cds_range ( @{ $gene->{CDS} } ) {
        grep $_->contains( $cds_range ), @{ $gene->{exon} }
           or warn "gene $gene_name, CDS ".$cds_range->toString." not contained by any exon\n"
    }
}
