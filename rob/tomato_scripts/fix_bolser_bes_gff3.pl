#!/usr/bin/perl
use strict;
use warnings;

use List::MoreUtils 'uniq';
use Bio::GFF3::LowLevel qw/ gff3_parse_feature gff3_format_feature /;

while(<>) {
    next if /^#/;
    unless( /^#/ ) {
        my $f= gff3_parse_feature( $_ );
        $f->{source} =~ s/^dund/DBolser_Dundee_BES_SSAHA/;
        $f->{source} =~ s/ /_/g;

        $f->{type} = {
            clone       => 'BAC_clone',
            'clone end' => 'BAC_end',
        }->{$f->{type}} || $f->{type};

        @{$f}{'start','end'} = sort { $a <=> $b } @{$f}{'start','end'};

        if( $f->{type} eq 'BAC_end' ) {
            $f->{attributes}{Name} = delete $f->{attributes}{ID};
            $f->{attributes}{Name}[0] =~ s/^_//;
        }
        elsif( $f->{type} eq 'BAC_clone' ) {
            # make the name of the BAC_clone features be the
            # SGN-standard names, with the other ones as aliases
            my @names = (
                @{ delete $f->{attributes}{Name}  || [] },
                @{ delete $f->{attributes}{Alias} || [] },
            );
            @names = uniq map {
                s/_[RF]$//;
                s/^_//;
                if( my ($lib) = /^([a-z_]+)/i ) {
                    $lib = {
                        SLm => 'SL_MboI',
                        SLe => 'SL_EcoRI',
                        HBa => 'LE_HBa',
                        FOS => 'SL_FOS',
                        SLs => 'SL_s',
                        SLt => 'SL_MT',
                        LPc => 'LpenCOS',
                        LPb => 'LpenBAC',
                        SLf => 'SL_FOS',
                        RH  => 'RH',
                    }->{$lib} || $lib;
                    s/^([a-z_]+)/$lib/;
                }
                $_
             } @names;

            $f->{attributes}{Name} = [ shift @names ];
            $f->{attributes}{Alias} = \@names if @names;
        }

        if( my $t = $f->{attributes}{Target} ) {
            if( $f->{attributes}{Alias} ) {
                $t->[0] =~ s/^x\s+/$f->{attributes}{Alias}[0] /;
            } elsif( $f->{attributes}{Name} ) {
                $t->[0] =~ s/^x\s+/$f->{attributes}{Name}[0] /;
            }

            my ( $id, $s, $e, $strand ) = split /\s+/, $t->[0];
            $strand ||= '+';
            ($s,$e) = sort { $a <=> $b } ($s,$e);
            $t->[0] = join ' ', $id, $s, $e, $strand;
        }

        $_ = gff3_format_feature( $f );
    }
    print;
}
