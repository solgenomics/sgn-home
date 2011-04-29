#!/usr/bin/perl
use strict;
use warnings;

use Bio::GFF3::LowLevel qw/ gff3_parse_feature gff3_format_feature /;

while(<>) {
    next if /^##genome-build/;
    unless( /^#/ ) {
        my $f= gff3_parse_feature( $_ );
        $f->{source} = 'DBolser_Dundee_BES_SSAHA';

        $f->{type} = {
            clone       => 'BAC_clone',
            'clone end' => 'BAC_end',
        }->{$f->{type}} || $f->{type};

        @{$f}{'start','end'} = sort { $a <=> $b } @{$f}{'start','end'};

        if( $f->{type} eq 'BAC_end' ) {
            $f->{attributes}{Name} = delete $f->{attributes}{ID};
        }
        elsif( $f->{type} eq 'BAC_clone' ) {
            # make the name of the BAC_clone features be the
            # SGN-standard names, with the other ones as aliases
            $f->{attributes}{Alias} = [ @{ $f->{attributes}{Name} } ];
            $f->{attributes}{Name} = [ map {
                my $lib = substr($_,0,3);
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
                }->{$lib} or die "unrecognized in $_";
                substr($_,0,3) = $lib;
                $_
               } @{$f->{attributes}{Name} || []}
             ];
        }

        if( my $t = $f->{attributes}{Target} ) {
            if( $f->{attributes}{Alias} ) {
                $t->[0] =~ s/^x\s+/$f->{attributes}{Alias}[0] /;
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
