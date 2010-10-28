#!/usr/bin/perl -w
use strict;
use DBI;

if (!$ARGV[0] || $ARGV[0] eq "help") {
  print <<EOF;

  Program to load the raw results of basecalling to the new sgn database. This
  program reads the ".phd.1" files produced by phred. Note that phred will
  produce one file for each chromatogram processed, so this may be an enormous
  number of files.

  This program inserts the information into database "sgn" as follows.

  For each identifier read from STDIN:
     - open the .phd.1 file from the phred output directory
     - Parse out the bases (and ucase them), quality scores, and basecalling
       positions.
     - Get phred version, trace name, etc, from .phd file
     - Select for read_id from table seqread, where trace_name eq the name of
       this trace
     - Insert information into table "est" linking to seqread.

  It is not the job of this script to create the clone table and seqread 
  table entries. The user should create them manually for each new plate 
  sent to a sequencing facility - or - migrate the relevant entries from the
  sgn_submit database if this is a data submission from outside.

  Usage: <host> <user> <pass> <phd output directory>

  NOTE: This script was modified from an previous form which required 
        information for creating seqread table entries. This information is
        now assumed to be already loaded. Likewise, this script will 
        query the database for each identifier on STDIN and make sure an
        entry exists in seqread before proceeding with any modifications.

EOF
  exit -1;
}

my ($host, $user, $pass, $phd_dir) = @ARGV;

my $dbh = DBI->connect("dbi:mysql:host=$host;database=sgn",$user,$pass, { RaiseError => 1 });

my $trace_nameq = $dbh->prepare("select read_id from seqread where trace_name=?");

my @db_notfound;
my @fs_notfound;
my %trace_namemap;

my @trace_names;
while(<STDIN>) {
  chomp;
  s/^.*\///;
  s/\.gz$//;
  push @trace_names, $_;
}

foreach ( @trace_names ) {
  # This extra qname business is to get rid of chromatogram extension names
  # names like .esd, .scf, and .ab1, because they are typically not entered
  # into the database by scripts that preceed this. Since this is not handled
  # consistently either by us (SGN) or by submitters and facilities, some
  # work arounds are needed pretty much anywhere that handles the files.
  my $qname = $_;
  $qname =~ s/\.(esd|ab.|scf)$//i;
  $trace_nameq->execute($qname);

  if ($trace_nameq->rows==0) {
    push @db_notfound, $_;
  } else {
    ($trace_namemap{$_}) = $trace_nameq->fetchrow_array();
  }
  if (! -f "$phd_dir/$_.phd.1") {
    push @fs_notfound, $_;
  }
}

if ( @fs_notfound>0 || @db_notfound>0) {
  print STDERR "Following trace names not found in database:\n";
  print STDERR join("\n",@db_notfound),"\n";
  print STDERR "Following trace names do not have .phd.1 files:\n";
  print STDERR join("\n",@fs_notfound),"\n";
  print STDERR scalar(@db_notfound)," not found in database, ",scalar(@fs_notfound)," not found on filesystem.\n";
  print STDERR "ERRORS DETECTED, NO INSERTS EXECUTED.\n";
  exit -1;
} else {
  print STDERR "All seqread entries and .phd.1 files found. Running inserts...\n";
}


# status and flags columns: Set status to indicate that vector and quality have
# not been assessed, contaminants have not been assessed, chimeric sequence 
# has not been assessed.
#
# Set flags to default 0 setting.
my $insert_estq = $dbh->prepare("insert into est (read_id,version,basecaller,seq,qscore,call_positions,date,status,flags) values (?,?,?,?,?,?,?,0x70,0)");

my $est_versionq = $dbh->prepare("SELECT MAX(version)+1 from est where read_id=?");

my $dateq = $dbh->prepare("select NOW()");
$dateq->execute;
my ($date) = $dateq->fetchrow_array();


foreach my $id ( @trace_names ) {
  my $phd_file = "$id.phd.1";

  if (open PHDFILE, "<$phd_dir/$phd_file") {
    my @seq = ();
    my @qual = ();
    my @pos = ();
    my ($phred_version, $version, $clone_id, $read_id);

    while(<PHDFILE>) {
      chomp;
      if (m/^BEGIN_SEQUENCE\ (.+)/) {
	my $sequence_name = $1;
	if ($sequence_name ne $id) {
	  print STDERR "Warning: Sequence name in $phd_file \"$sequence_name\" does ".
	    "not match file name for $id.\n";
	}
      }
      if (m/^PHRED_VERSION: (.+)/) {
	$phred_version = $1;
      }
      if (m/^BEGIN_DNA/) {
	while(<PHDFILE>) {
	  chomp;
	  last if m/^END_DNA/;
	  my ($base,$qual,$pos) = split;
	  push @seq, uc $base;
	  push @qual, $qual;
	  push @pos, $pos;
	}
      }
      last if m/^END_DNA/;
    }

    close PHDFILE;

    $read_id = $trace_namemap{$id};
    $est_versionq->execute($read_id);
    if ($est_versionq->rows > 1) {
      ($version) = $est_versionq->fetchrow_array();
    } else {
      $version = 1;
    }
    $insert_estq->execute($read_id,$version,$phred_version,join("",@seq),
			  join(" ",@qual),join(" ",@pos),$date);
	
  } else {
    print STDERR "Can't open PHD file \"$phd_file\" ($!)\n";
  }

}
