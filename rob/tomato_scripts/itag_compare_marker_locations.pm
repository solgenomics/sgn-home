#!/usr/bin/env perl
package itag_compare_marker_locations;
use Moose;
#use warnings FATAL => 'all';
use feature 'say';

use Data::Dump 'dump';

use Path::Class;

use CXGN::Marker;
use CXGN::Map;

use Bio::GFF3::LowLevel 'gff3_parse_feature';

use List::MoreUtils qw/ uniq any /;

with 'MooseX::Getopt';
with 'MooseX::Runnable';
with 'MooseX::Role::DBIx::Connector' => {
    accessor_options => {
        db_conn  => [ traits => ['NoGetopt']],
        db_attrs => [ traits => ['NoGetopt']],
    },
};


has 'map_name' => (
    is  => 'ro',
    isa => 'Str',
    required => 0,
    documentation => 'restrict map names (regex)',
);

sub run {
    my ( $self, $markers_gff3 ) = @_;

    $self->db_conn->dbh->do('set search_path = sgn,public');

    my %mapping_data;
    $self->read_itag_locations( \%mapping_data, $markers_gff3 );
    $self->fetch_marker_locations( \%mapping_data );

    # find markers with chromosome disagreements
    for my $id ( keys %mapping_data ) {
        my $locs = $mapping_data{$id};
        my @chrs = uniq map $_->[0], @$locs;
        delete $mapping_data{$id} unless @chrs > 1;
    }

    $_ = [ sort { $a->[0] <=> $b->[0] } @$_ ] for values %mapping_data;
    print dump( \%mapping_data );
}

sub match_map_name {
    my ( $self, $name ) = @_;
    my $n = $self->map_name or return 1;
    return $name =~ /$n/i;
}

sub fetch_marker_locations {
    my ( $self, $mapping_data ) = @_;

    for my $id ( keys %$mapping_data ) {
        my $m = CXGN::Marker->new( $self->db_conn->dbh, $id ) or die "$id not found ";

        my @locs =
            map {
                my $loc = $_;
                die dump $loc unless $loc->{lg_name};
                my $map = CXGN::Map->new( $self->db_conn->dbh, {map_version_id=>$loc->{map_version_id}});
                my $mapname = $map->get_short_name;
                if( $self->match_map_name( $mapname ) ) {
                    [ $loc->{lg_name}, $map->get_short_name, $loc->{position}, $loc->{confidence} ]
                } else {
                    ()
                }
            }
            grep $_,
            map $_->{location},
            @{ $m->current_mapping_experiments };

        push @{ $mapping_data->{$id} ||= [] }, @locs;
    }
}


sub read_itag_locations {
    my ( $self, $mapping_data, $markers_gff3 ) = @_;
    $markers_gff3 = file($markers_gff3);
    my $fh = $markers_gff3->openr;

    while( <$fh> ) {
        next unless /\tITAG_sgn_markers\t/i && /^SL\d\.\d\dch(\d\d)\s/ && /\tmatch\t/;

        my $f = gff3_parse_feature( $_ );
        #die dump $f;
        my ( $id ) = $f->{attributes}{Alias}[0] =~ /^SGN-M(\d+)/
          or die dump $f;
        my ( $chr ) = $f->{seq_id} =~ /^SL\d\.\d\dch(\d\d)$/;
        $chr += 0;
        next if $chr == 0;
        my ( $s, $e ) = sort { $a <=> $b } @{$f}{'start','end'};
        push @{ $mapping_data->{$id+0} ||= [] }, [$chr, 'ITAG2.3', "$s..$e", $f->{score}];
    }
}


1;
