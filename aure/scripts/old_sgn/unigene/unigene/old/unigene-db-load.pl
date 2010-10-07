#!/usr/bin/perl -w
use strict;
use DBI;
use CXGN::BioTools::FastaParser;


if (!$ARGV[0] || $ARGV[0] eq "help") {
  print <<EOF;

  Script to load output of a new unigene build into new SGN database.

  Actions performed:

  - New unigene build entry created in unigene_build table. Organism group
    must already be created.
  - Source data group created from input sequence file. FASTA ids are
    expected to be the est_id from EST table
  - unigene output FASTA file is read and unigene table loaded.
  - Output of acefile_membership.pl is read and loaded into unigene_member
  - ACE/ output directory is traversed again, .ace files gzipd and loaded
    into unigene_ace table.

  Notes: 
     - Organism group must already be created.
     - unigene.seq FASTA output file must have corresponding .qual file, and
       FASTA header names should be cluster_id-contig_id for contigs and
       cluster-S-est_id for singlets, and X-S-est_id for preclustering singlets
     - acefile_membership.pl must have been run and output stored in a file
     - "working directory" should have unigene/ and ACE/ and basecall/ 
       directories.

  Usage: <host> <user> <pass> <organism group id> <unigene working directory> 
         <filename of acefile_membership.pl output>

EOF
  exit -1;
}

my ($host, $user, $pass, $organism_group_id, $unigene_dir, 
    $membership_filename) = @ARGV;

my $dbh = DBI->connect("dbi:mysql:host=$host;database=sgn",$user,$pass,
                          { RaiseError => 1, PrintError => 0 })
  or die "Failed to connect to database ($DBI::errstr)";

if ( ! -d $unigene_dir ) {
  print STDERR "Can't find unigene output directory \"$unigene_dir\" ($!)\n";
  exit -1;
}

if ( ! -f $membership_filename ) {
  print STDERR "Can't find output file of acefile_membership.pl " .
    "\"$membership_filename\" ($!)\n";
  exit -1;
}

# Check that unigene input file can be found
my $ifh = CXGN::BioTools::FastaParser->load_file("$unigene_dir/basecall/qc-passed.seq")
  or die "Failed to open \"$unigene_dir/basecall/qc-passed.seq\" ($!)";

# Check that FASTA output file can be found
my $ufh = CXGN::BioTools::FastaParser->load_file("$unigene_dir/unigene/unigene.seq","$unigene_dir/unigene/unigene.qual");


# Check that the organism group already exists...
my $sth;
$sth = $dbh->prepare("SELECT o.organism_name from organism as o, group_linkage as l where l.member_id=o.organism_id and l.group_id=?");

$sth->execute($organism_group_id);

if ($sth->rows == 0) {
  print "Organism group \"$organism_group_id\" not found in database\n";
  exit -1;
}

print "Organisms in group:\n";
while(($_) = $sth->fetchrow_array()) {
  print "\t$_\n";
}
print "\n";

$sth->finish();

# Find the build number in this series, create a new build.
my ($build) = $dbh->selectrow_array("SELECT MAX(build_nr) from unigene_build where organism_group_id=$organism_group_id");

if (!defined($build)) {
  $build = 1;
} else {
  $build++;
}

$dbh->do("INSERT into unigene_build set organism_group_id=$organism_group_id,build_nr=$build,method_id=2");
my $unigene_build_id = $dbh->{'mysql_insertid'};

print "Unigene build id is $unigene_build_id, build_nr = $build\n";

# Create an input data group and set it for the unigene build
$dbh->do("INSERT into groups set type=?,comment=?",undef,2,"Input sequence data group for unigene build (unigene id $unigene_build_id), organism group $organism_group_id");

my $input_group_id = $dbh->{'mysql_insertid'};

$dbh->do("UPDATE unigene_build set source_data_group_id=? where unigene_build_id=?",undef,$input_group_id,$unigene_build_id);


# Load input data group
my $idg_iq = $dbh->prepare("INSERT into group_linkage set group_id=$input_group_id,member_id=?,member_type=7");

foreach my $est_id ( $ifh->get_seqids() ) {
  $idg_iq->execute($est_id);
}

$idg_iq->finish();


open AP, "<$membership_filename"
  or die "Failed to open file \"$membership_filename\" ($!)";

open TF, ">member-load-data.tdv"
  or die "Failed to open temporary file for member load data ($!)";

my $iu_q = $dbh->prepare("INSERT into unigene set unigene_build_id=?,consensi_id=?,cluster_no=?,contig_no=?,nr_members=?");

my $ic_q = $dbh->prepare("INSERT into unigene_consensi set seq=?,qscores=?");

print STDERR "Loading unigenes and consensi...\n";
while(<AP>) {
  if (m/^Contig:/) {
    my (undef, $cluster, $contig, $n_members) = split;

    # Find/Load the consensus sequence and quality values
    my ($consensus_id, $unigene_id);
    my ($seq, @qual);
    my $seq_id = "$cluster-$contig";
    $seq = $ufh->get_sequence($seq_id,\@qual);

    if (!$seq) {
      print STDERR "Can't find consensus sequence for $cluster-$contig\n";
      $consensus_id = undef;
    } else {
      $ic_q->execute($seq,join(" ",@qual));
      $consensus_id = $dbh->{'mysql_insertid'};
      $iu_q->execute($unigene_build_id,$consensus_id,$cluster,
			       $contig,$n_members);
      $unigene_id = $dbh->{'mysql_insertid'};
    }

    # Read in the members ests
    for(my $i=0;$i<$n_members;$i++) {
      $_ = <AP>;
      my ($est_id,$dir,$start,$end,$start_trim,$end_trim) = split;

      if ($est_id =~ m/(SGN-E|E|sgn\|E)([0-9]+)/) {
	$est_id = $2;
      }
      print TF join("\t",$unigene_id,
		    $est_id, $start, $end, $start + $start_trim,
		    $end - $end_trim, $dir),"\n";
    }
  } elsif (m/^Singlet:/) {
    my (undef,$cluster) = split;
    $_ = <AP>;
    my ($est_id, undef, $start, $stop) = split;
    if ($est_id =~ m/(SGN-E|E|sgn\|E)([0-9]+)/) {
      $est_id = $2;
    }

    if ($cluster == -1) {
      $cluster = "NULL";
    }

    $iu_q->execute($unigene_build_id,undef,$cluster,undef,1);
    my $unigene_id = $dbh->{'mysql_insertid'};
    print TF join("\t",$unigene_id,$est_id,
			    $start,$stop,$start,$stop,"+"),"\n";
  } else {
    print STDERR "Following output from membership file not understood:\n";
    print STDERR $_;
  }
}

close AP;
close TF;

$ic_q->finish();

print STDERR "Loading membership data...\n";
$dbh->do("LOAD DATA LOCAL INFILE \"member-load-data.tdv\" into table unigene_member (unigene_id,est_id,start,stop,qstart,qend,dir)");

print STDERR "Loading ACE file data...\n";

my $ai_q = $dbh->prepare("INSERT into unigene_ace (unigene_build_id, cluster_no, acedata) VALUES (?,?,?)");

open FP, "find $unigene_dir/ACE -name \"*.ace\" |"
  or die "Failed to open pipe find ($!)";

my $cluster;
while(<FP>) {
  chomp;
  if (m/cluster-([0-9]+)\.seq/) {
    $cluster = $1;
  } else {
    print "Can't parse cluster number from ace filename \"$_\" -- ".
      "skipping file\n";
    next;
  }

  open AF, "gzip -c $_ |"
    or die "Failed to gzip pipe for ace file \"$_\" ($!)";

  my @data = <AF>;

  $ai_q->execute($unigene_build_id,$cluster,join("",@data));

  close AF;
}

close FP;
$ai_q->finish();
