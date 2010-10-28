#!/usr/bin/perl

=head1 NAME

 sgndb_report_trace_processing_by_organism.pl.
 A script to get information about trace processing tables for each organism from SGN database  (version.1.0.).

=head1 SYPNOSIS
  
 sgndb_report_trace_processing_by_organism.pl [-h] -H <dbhost> -D <dbname> [-o <organism_name>] [-A] [-F]

=head1 DESCRIPTION

This script get a table from the database with the follows fields:
 
=over

=item - 

sgn.organism.organism_name
 
=item -

COUNT(sgn.library.library_id) => biological information

=item - 

COUNT(sgn.clone.clone_id) => clone names

=item -

COUNT(sgn.seqread.read_id) => sources (like chromatograms)

=item -

COUNT(sgn.est.est_id) => sequences and tags (status, flags)

=item -

COUNT(sgn.qc_report.qc_id) => sequences coordenates (for example from cleaning process)

=item -

COUNT(sgn.est_dbxref.est_dbxref_id) => crossreference entries

=item -

COUNT(public.dbxref.dbxref_id) => accessions

=over

=head1 AUTHOR

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=head1 METHODS
 
sgndb_report_trace_processing_by_organism.pl


=cut

use strict;

use Bio::SeqIO;
use CXGN::DB::InsertDBH;
use Getopt::Std;


our ($opt_D, $opt_H, $opt_o, $opt_A, $opt_F, $opt_h);
getopts("D:H:o:AFh");

if (!$opt_D && !$opt_H && !$opt_o && !$opt_A && !$opt_F && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
   help();
}

if (!$opt_D || !$opt_H) { 
      die "\nSorry, this script needs -D (database name) and -H (host) arguments.\n";
  }

our $dbh = CXGN::DB::InsertDBH::connect ({ 	dbname => $opt_D, 
						dbhost => $opt_H,
					  });

my ($data, $get_all, $flags);

if ($opt_A) {
    $get_all='LEFT';
}

my @header = ("ORGANISM","N_LIBRARIES", "N_CLONES", "N_READS", "N_ESTS", "N_QC_REPORTS", 
              "N_EST_DBXREF", "N_ACCESSIONS");
my @data;

push @data, \@header;



if ($opt_o) {
    if ($opt_F) {
        $flags="AND flags=0 and status=0 ";
    }
    my @organism = split(',', $opt_o);
    my $organism;
    foreach $organism (@organism) {
      my $query2 = "SELECT sgn.organism.organism_name, 
                           COUNT(DISTINCT sgn.library.library_id) AS lib_count, 
                           COUNT(DISTINCT sgn.clone.clone_id) AS clone_count, 
                           COUNT(DISTINCT sgn.seqread.read_id) AS read_count, 
                           COUNT(DISTINCT sgn.est.est_id) AS est_count, 
                           COUNT(DISTINCT sgn.qc_report.qc_id) AS qc_count, 
                           COUNT(DISTINCT sgn.est_dbxref.est_dbxref_id) AS est_dbxref_count, 
                           COUNT(DISTINCT public.dbxref.accession) AS accession_count 
                           FROM sgn.organism 
                           LEFT JOIN sgn.library ON sgn.organism.organism_id=sgn.library.organism_id 
                           LEFT JOIN sgn.clone ON sgn.library.library_id=sgn.clone.library_id 
                           LEFT JOIN sgn.seqread ON sgn.clone.clone_id=sgn.seqread.clone_id 
                           LEFT JOIN sgn.est ON sgn.seqread.read_id=sgn.est.read_id 
                           LEFT JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id 
                           LEFT JOIN sgn.est_dbxref ON sgn.est.est_id=sgn.est_dbxref.est_id 
                           LEFT JOIN public.dbxref ON sgn.est_dbxref.dbxref_id = public.dbxref.dbxref_id 
                           WHERE organism_name=? $flags 
                           GROUP BY sgn.organism.organism_name";
      my $sth=$dbh->prepare($query2);
      $sth->execute($organism);
      my (@trace_data)=$sth->fetchrow_array();
      push @data, \@trace_data;
      
    }
    
} else {

    if ($opt_F) {
        $flags="WHERE flags=0 AND status=0 ";
    }
    my $query = "SELECT sgn.organism.organism_name, 
                 COUNT(DISTINCT sgn.library.library_id) AS lib_count, 
                 COUNT(DISTINCT sgn.clone.clone_id) AS clone_count, 
                 COUNT(DISTINCT sgn.seqread.read_id) AS read_count, 
                 COUNT(DISTINCT sgn.est.est_id) AS est_count, 
                 COUNT(DISTINCT sgn.qc_report.qc_id) AS qc_count, 
                 COUNT(DISTINCT sgn.est_dbxref.est_dbxref_id) AS est_dbxref_count, 
                 COUNT(DISTINCT public.dbxref.accession) AS accession_count 
                 FROM sgn.organism 
                 $get_all JOIN sgn.library ON sgn.organism.organism_id=sgn.library.organism_id 
                 LEFT JOIN sgn.clone ON sgn.library.library_id=sgn.clone.library_id 
                 LEFT JOIN sgn.seqread ON sgn.clone.clone_id=sgn.seqread.clone_id 
                 LEFT JOIN sgn.est ON sgn.seqread.read_id=sgn.est.read_id 
                 LEFT JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id 
                 LEFT JOIN sgn.est_dbxref ON sgn.est.est_id=sgn.est_dbxref.est_id 
                 LEFT JOIN public.dbxref ON sgn.est_dbxref.dbxref_id = public.dbxref.dbxref_id 
                 $flags
                 GROUP BY sgn.organism.organism_name ORDER BY organism_name";
    my $sth=$dbh->prepare($query);
    $sth->execute();
    while (my ($organism, @trace_data)=$sth->fetchrow_array()) {
           $organism =~ s/ \(\w+\s+\w+\s+\w+\)//gi;
           unshift @trace_data, $organism;
           push @data, \@trace_data;
    }
}

my $data=\@data;

print_table($data);

=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0: 

    Description: 
      This script get information about trace processing tables (sgn.library, sgn.clone, sgn.seqread, sgn.est, sgn.qc_report, sgn.est_dbxref and public.dbxref) for the sequences associated to each organism in SGN database.

    Usage: 
      sgndb_report_trace_processing_by_organism.pl [-h] -H <dbhost> -D <dbname> [-o <organism_name>] [-A]

      if you want get the information in a file (and not print in the screen):

      sgndb_report_trace_processing_by_organism.pl -H <dbhost> -D <dbname> [-o <organism_name>] [-A] > file
    
    Example:
      sgndb_report_trace_processing_by_organism.pl -H localhost -D sandbox -o 'Nicotiana tabacum','Nicotiana benthamiana'

    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -A enable option        get the values for all organism (not only the organism that have libraries)
      -F enable option        get the value where flags=0 and status=0
      -o organism_names       if you want not all the organism, you can specify the organism 
                              (or organisms inside single quotes and separated by commas, example
                              'Nicotiana tabacum','Nicotiana benthamiana').
      -h print this help

EOF
exit (1);
}


=head2 print_table

  Usage: print_table ($arrayref_of_arrays);
  Desc: print a table of a array with arrayref like rows.
  Ret: none
  Args: $arrayref_of_arrays (so means, an arrayref with @rows like elements).
  Side_Effects: none
  Example: @header = ('organism', 'est_ids');
           @row1 = ('tomato', 1200);
           @row2 = ('tobacco', 320);
           @array = (\@header, \@row1, \@row2);
           print_table(\@array);

=cut


sub print_table {
    my $table=shift;
    my @table=@$table;

    my $length=0;
    my ($row, $data, $datalen);
    
    foreach $row (@table) {
	my @row=@$row;
        my $organism=shift(@row);
        my $organism_length=length($organism);
	    if ($organism_length > $length) {
		$length = $organism_length;
	    }
    }
    
    
    $length += 2;
    foreach $row (@table) {
	my @row=@$row;
        my $organism=shift(@row);
        my $organism_length=length($organism);
        while ($organism_length < $length) {
	    $organism .= " ";
            $organism_length=length($organism);
        }
        print "|  $organism|  ";
	foreach $data (@row) {
            $datalen=length($data);
	    while ($datalen < 14) {
		$data .= " ";
                $datalen=length($data);
	    }
	    print "$data|  ";   
	}
        print "\n";
    }
}
