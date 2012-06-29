
=head1 NAME

update_solcap_marker_position.pl

=head1 SYNOPSIS

perl update_solcap_marker_position.pl [options] -H hostname -D dbname -i in_gff3_file


=head1 DESCRIPTION

Update the position of the solcap markers.
This was not done after running patch #21 , which altered the marker_location.position field from numeric(8,5) to (9,6).
As a result, the positions are now wrong.


=head1 AUTHOR

 Naama Menda<nm249@cornell.edu>

=head1 COPYRIGHT & LICENSE

Copyright 2012 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut


#!/usr/bin/perl

use strict;
use Getopt::Std;
use File::Slurp;
use CXGN::Marker;

use CXGN::DB::InsertDBH;
use Try::Tiny;

use Carp qw /croak/ ;


our ($opt_H, $opt_D, $opt_i);

getopts('H:i:D:');

my $dbhost = $opt_H;
my $dbname = $opt_D;
my $infile = $opt_i;

my $dbh = CXGN::DB::InsertDBH->new( { dbhost=>$dbhost,
				      dbname=>$dbname,
				      dbargs => {AutoCommit => 0,
						 RaiseError => 1}
				    }
    );


my @gff = read_file($infile);

try {
    foreach my $line (@gff) {
        chomp $line;
        #SL2.40ch07	ITAG_sgn_markers	match	1890767	1890867	.	.	.	Alias=SGN-M18753;ID=solcap_snp_sl_11171;Name=solcap_snp_sl_11171

        my ($chr_string, undef, undef, $pos1, $pos2, undef, undef, undef, $description) = split (/\t/, $line);
        my $marker_name ;
        if ( $description =~/^Alias=(SGN-M.*);ID=(solcap.*);Name=(.*)/ ) {
            $marker_name = $2;
        } else {
            croak("**** something is wrong with marker name . ($description) \n\n");
        }
        $chr_string =~ m/SL2.40ch0?(\d+)/;
        my $file_chr = $1;
        my $check = $pos2-$pos1;
        if ($check> 200) {
            croak("***something looks wrong with position for marker $marker_name, span = $check\n");
        }elsif ($check < 100) {
            warn("!!!marker span for $marker_name is less than 100 bp ($check) \n");
        }
        my $new_position = ($pos2+$pos1)/2;
        $new_position = sprintf "%.0f", $new_position; #round up if midpoint is at .5
        $new_position = $new_position/1000000; #divide by 10^6 to fit into the database numeric(9,6);
        my $marker = CXGN::Marker->new_with_name($dbh, $marker_name);
        if (!$marker) { warn "**No marker exists for $marker_name!!!\n\n"; next; }
        my $exps=$marker->experiments();
        my ($chr, $position);
        if($exps) {
            for my $exp(@{$exps}) {
                if($exp->{location}) {
                    if ( $exp->{location}->map_version_id == 105 ) {
                        $chr = $exp->{location}->lg_name;
                        my $location_id = $exp->{location}->location_id;
                        if ($chr != $file_chr) { die "***marker chromosome $chr does not match file $file_chr (marker name = $marker_name) \n "; }
                        $position = $exp->{location}->position;
                        print "Updating location for marker $marker_name. chr = $chr, position = $position. New position = $new_position  \n";
                        my $q = "UPDATE sgn.marker_location SET position = $new_position WHERE location_id = ?";
                        my $sth = $dbh->prepare($q);
                        $sth->execute($location_id) ;
                    }
                }
            }
        }
    }
    $dbh->commit;
} catch {
    $dbh->rollback;
    die "An error occured! Rolling back!! $_ \n" ;
};


