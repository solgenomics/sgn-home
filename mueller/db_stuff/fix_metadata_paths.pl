
=head1 NAME

fix_metadata_paths.pl - fix dirnames in metadata table

=head1 DESCRIPTION

perl fix_metadata_paths.pl -h localhost -u postgres -d cxgn

Will replace /data in the paths with /export, as required for the installation in the BTI server room.

Note: it will prompt for a password.

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut

use strict;

use Getopt::Std;
use Term::ReadKey;
use CXGN::Metadata::Schema;

our($opt_h, $opt_d, $opt_u);

getopts('h:d:u:');

my $user = $opt_u || die "need a username with -u option";
my $dbname = $opt_d || die "need a dbname with -d option";
my $dbhost = $opt_h || die "need host with -h option";

my $dsn = "dbi:Pg:dbname=$dbname;host=$dbhost";



ReadMode("noecho");

print "password for $user: ";
my $password = <STDIN>;
chomp($password);
ReadMode("original");

my $schema = CXGN::Metadata::Schema->connect($dsn, $user, $password, { on_connect_do => 'SET search_path TO metadata, public' });

my @rows = $schema->resultset("MdFiles")->search();

my $count = 0;
foreach my $row (@rows) { 
    if ($row->dirname() =~ m|^/data|) { 
	my $dirname = $row->dirname();
	$dirname =~ s|^(/data)|/export|;
	
	print STDERR "NEW DIRNAME $dirname\n";
	$row->update( { dirname=>$dirname });

	
    }
}
print STDERR "Done.\n";
