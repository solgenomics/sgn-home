#!/usr/bin/perl -w

=head1 NAME

query-unigene-input.pl

=head1 DESCRIPTION

Usage: query-unigene-input.pl [-h] -H dbhost -D dbname -o organism_ids  BASENAME

=head2 Notes

=over 5

=item 

Neither BASENAME.seq nor the BASENAME.qual may exist
already. Typically, "qc-passed.seq" should be used for usage
with unigene_assembly.pl.

=item 

several organism_ids can be specified as a comma separated list (no spaces).

=back

Program to create input FASTA files for running a unigene assembly. At time
of writing, the unigene assembly was performed by running 
"unified_assembly.pl" with FASTA input files "qc-passed.seq" and 
"qc-passed.qual" located in a directory called "basecall/" that is contained
in the directory the unigene_assembly.pl script is being run. 

This program functions to create these input files by querying the database.
Note that subsequent database loading scripts expect the FASTA headers to be
the est_id identifier, as these FASTA headers are passed along through the
assembly process to identify the assembled ESTs.

For input, the program expects a list of libraries fed as STDIN. This should
be formatted as a tab-delimited file with one line per library, with the db
library_id in the first column. This should be created using the mysql
monitor.



=head1 AUTHORS

unknown (lucky guys)

=cut

use strict;
use CXGN::DB::InsertDBH;
use Getopt::Std;

our ($opt_h, $opt_H, $opt_D, $opt_o);

getopts ("hH:D:o:");

if ($opt_h) {
  help ();
}

my ($basename) = @ARGV;
unless ($basename) {
  help ();
}
my $output_fname = $basename . ".seq";
if ( -f $output_fname . ".seq" ) {
  print STDERR "Output file \"$output_fname\" already exists\n";
  exit (1);
}
my $qual_fname = $basename . ".qual";
if ( -f $qual_fname ) {
  print STDERR "Output quality file $qual_fname already exists\n";
  exit (1);
}

my $dbh = CXGN::DB::InsertDBH->new( { dbname => $opt_D,
				       dbhost => $opt_H,
				       dbuser => "postgres" });

my $organism_q = "SELECT organism_name FROM sgn.organism WHERE organism_id=?";
my $organism_h = $dbh->prepare($organism_q);
$organism_h->execute($opt_o);
my ($organism_name) = $organism_h->fetchrow_array();

if (!$organism_name) { die "Can't find organism $opt_o!!!\n"; }

print STDERR "$opt_o IS ORGANISM $organism_name.\n"; sleep(1);

my $libraryq = "SELECT library_id FROM sgn.library WHERE type NOT IN (SELECT type_id FROM sgn.types WHERE comment = 'demethylated genomic sequence library') AND organism_id IN ($opt_o)";

my @libraries = @{$dbh->selectcol_arrayref($libraryq)};
print "".(join ", ", @libraries)."\n";
if (@libraries==0) {
  print STDERR "No libraries selected.\n";
  exit -1;
}

my $lq = $dbh->prepare("SELECT library_shortname from library where library_id = ?");
my $cursorname = "cursor";
my $sq = $dbh->prepare ("DECLARE $cursorname CURSOR FOR SELECT clone.clone_name, direction, est.est_id, seq, qscore, hqi_start, hqi_length from clone inner join seqread using (clone_id) inner join est using (read_id) inner join qc_report using (est_id) where clone.library_id=? and est.status=0 and est.flags=0");
my $sth = $dbh->prepare ("FETCH 100 FROM $cursorname");

open SF, ">$output_fname"
  or die "Can't open sequence output file \"$output_fname\" ($!)";

open QF, ">$qual_fname"
  or die "Can't open quality output file \"$qual_fname\" ($!)";

my $total = 0;
foreach my $library_id ( @libraries ) {
  $lq->execute($library_id);
  if ($lq->rows == 0) {
    print STDERR "No library found for library_id \"$library_id\"\n";
    next;
  }
  my ($shortname) = $lq->fetchrow_array();
  print STDERR "Querying library $shortname... ";

  my $n = 0;
  $sq->execute($library_id);
  while ($sth->execute) {
    if ($sth->rows == 0) {
      last;
    }
    while(my ($clone_name, $dir, $est_id, $seq, $qscore, $hqi_start, $hqi_length) =
	  $sth->fetchrow_array()) {

      my $header = ">SGN-E$est_id [$clone_name]";
      $seq = substr($seq,$hqi_start,$hqi_length);
      my @qual = split/\s+/,$qscore;
      @qual = splice(@qual,$hqi_start,$hqi_length);

# don't do this because the sequence needs to be the same as in the database.
   #    if ($dir eq "3") {
# 	$seq =~ tr/ACGT/TGCA/;
# 	$seq = join("",reverse split//,$seq);
# 	@qual = reverse @qual;
# 	$header .= "(-)";
#       }
      print SF "$header\n";
      print SF "$seq\n";

      print QF "$header\n";
      print QF join(" ",@qual),"\n";

      $n++;
      $total++;
    }
  }
  $dbh->do ("CLOSE $cursorname");
  print STDERR "$n sequences.\n";
}

$sq->finish();

close SF;
close QF;

print STDERR "$total sequences in $output_fname\n";

sub help {
  print <<EOF;

  Program to create input FASTA files for running a unigene assembly. At time
  of writing, the unigene assembly was performed by running 
  "unified_assembly.pl" with FASTA input files "qc-passed.seq" and 
  "qc-passed.qual" located in a directory called "basecall/" that is contained
  in the directory the unigene_assembly.pl script is being run. 

  This program functions to create these input files by querying the database.
  Note that subsequent database loading scripts expect the FASTA headers to be
  the est_id identifier, as these FASTA headers are passed along through the
  assembly process to identify the assembled ESTs.

  For input, the program expects a list of libraries fed as STDIN. This should
  be formatted as a tab-delimited file with one line per library, with the db
  library_id in the first column. This should be created using the mysql
  monitor.

  Usage: query-unigene-input.pl [-h] -H dbhost -D dbname -o organism_ids  BASENAME

  Note: Neither BASENAME.seq nor the BASENAME.qual may exist
        already. Typically, "qc-passed.seq" should be used for usage
        with unigene_assembly.pl.

        several organism_ids can be specified as a comma separated list (no spaces).

EOF
  exit 0;
}
