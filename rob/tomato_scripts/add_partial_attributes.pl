#!/usr/bin/env perl

# adds partial_5prime and partial_3prime boolean attributes in-place
# to one or more files

use strict;
use warnings;
use autodie ':all';

use File::Copy;
use File::Temp;

use Bio::GFF3::LowLevel qw( gff3_parse_feature gff3_format_feature );

my %partial_5 = read_list( shift );
my %partial_3 = read_list( shift );

for my $file (@ARGV) {
    open my $fh, '<', $file;
    my $out = File::Temp->new;

    while( my $line = <$fh> ) {
        last if $line =~ /^##FASTA/;
        next if $line =~ /^#/;

        my $f = gff3_parse_feature( $line ) or die;
        next unless $f->{type} =~ /^(mRNA|gene)$/;

        my ($model_name) = $f->{attributes}{ID}[0] =~ /(Solyc\d+g\d+)/ or die;
        my $modified;

        if( $partial_5{$model_name} ) {
            push @{$f->{attributes}{partial}}, "5'";
            $modified++;
        }

        if( $partial_3{$model_name} ) {
            push @{$f->{attributes}{partial}}, "3'";
            $modified++;
        }

        $line = gff3_format_feature( $f ) if $modified;
    } continue {
        $out->print( $line );
    }

    $out->close;

    move( "$out", $file ) or die "$! moving $out -> $file";
}

######## subroutines ########

sub read_list {
    my $file = shift;
    open my $f, '<', $file;
    my @l;
    while( my $x = <$f> ) {
        chomp $x;
        push @l, ( $x => 1 );
    }
    return @l;
}
