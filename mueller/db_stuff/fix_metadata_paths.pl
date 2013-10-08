
use CXGN::Metadata::Schema;
use Term::ReadKey;

my $user = shift || die "need a username as first argument";
my $dbname = shift || die "need a dbname as second argument";

my $dsn = "dbi:Pg:dbname=$dbname";



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
