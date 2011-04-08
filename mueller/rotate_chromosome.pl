
use Modern::Perl;

use Getopt::Std;
use CXGN::DB::InsertDBH;

our ($opt_D, $opt_H);

# this script rotates a given chromosome on a map. Specifically written to turn around chromosome 12 on the potato QTL map.

my $map_version_id = shift;
my $chr_name = shift;

my $dbh = CXGN::DB::InsertDBH->new( { dbname=> $opt_D,
				      dbhost=> $opt_H,
				    });


my $sth = $dbh->prepare("SELECT lg_id FROM sgn.linkage_group WHERE map_version_id=? AND lg_name = ?");

$sth->execute($map_version_id, $chr_name);

my ($lg_id) = $sth->fetchrow_array();

print STDERR "map_version_id $map_version_id chr_name $chr_name lg_id = $lg_id\n";

print STDERR "Determining size of linkage group... ";
$sth = $dbh->prepare("SELECT max(position) FROM sgn.marker_location where map_version_id =? AND lg_id = ?");

$sth->execute($map_version_id, $lg_id);

my ($max_position) = $sth->fetchrow_array();

print STDERR "$max_position. \n";

print STDERR "Get marker data...\n";
$sth = $dbh->prepare("SELECT location_id, position, position_north, position_south FROM sgn.marker_location WHERE map_version_id =? AND lg_id=? order by position");

$sth->execute($map_version_id, $lg_id);

my $update_h = $dbh->prepare("UPDATE marker_location set position=?, position_north=?, position_south=? WHERE location_id=?");

print STDERR "Updating positions...\n";
my $count = 0;
while (my ($location_id, $position, $position_north, $position_south) = $sth->fetchrow_array()) { 

    my $new_position = $max_position - $position;
    my $new_position_north = undef;
    if ($position_north) { $new_position_north = $max_position - $position_north; }
    my $new_position_south = undef;
    if ($position_south) { $new_position_south = $max_position - $position_south; }
    my $result = $update_h -> execute($new_position, $new_position_north, $new_position_south, $location_id);

    print STDERR "OLD: $position. NEW: $new_position.\n";

    if ($result) { $count++; }
}

print STDERR "COMMITTING...\n";
$dbh->commit;


print STDERR "DONE.\n";
