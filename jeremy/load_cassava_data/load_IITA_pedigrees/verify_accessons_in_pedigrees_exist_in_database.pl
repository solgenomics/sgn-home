
=head1

verify_accessions_in_pedigrees_exist_in_database.pl

=head1 SYNOPSIS

    $this_script.pl -H [dbhost] -D [dbname] [-t]

=head1 COMMAND-LINE OPTIONS

 -H  host name
 -D  database name
 -i infile
 -u sgn user name
 -t  Test run . Rolling back at the end.


=head2 DESCRIPTION



=head2 AUTHORS

Jeremy Edwards (jde22@cornell.edu)
Naama Menda (nm249@cornell.edu)

June 2012

=cut


#!/usr/bin/perl
use strict;
use Getopt::Std;
use CXGN::Tools::File::Spreadsheet;
use CXGN::People::Person;

use Bio::Chado::Schema;
use CXGN::DB::InsertDBH;
use Carp qw /croak/ ;
##
##


our ($opt_H, $opt_D, $opt_i, $opt_t, $opt_u);

getopts('H:i:tD:u:');

my $dbhost = $opt_H;
my $dbname = $opt_D;
my $file = $opt_i;


my $dbh = CXGN::DB::InsertDBH->new( { dbhost=>$dbhost,
				      dbname=>$dbname,
				      dbargs => {AutoCommit => 0,
						 RaiseError => 1}
				    }
    );
my $schema= Bio::Chado::Schema->connect(  sub { $dbh->get_actual_dbh() } );
$dbh->do('SET search_path TO public');


my $accession_cvterm = $schema->resultset("Cv::Cvterm")->create_with(
    { name   => 'accession',
      cv     => 'stock type',
      db     => 'null',
      dbxref => 'accession',
    });

########################

#new spreadsheet, skip first column
my $spreadsheet=CXGN::Tools::File::Spreadsheet->new($file, 0);




my @rows = $spreadsheet->row_labels();
my @columns = $spreadsheet->column_labels();


print "accession_name\n";


my %missing_accessions;
#my %missing_pedigree_accessions;
my %pedigree_accessions;
foreach my $num (@rows ) {
 my $stock_name = $spreadsheet->value_at($num , "accession_name");
 my $pedigree_string = $spreadsheet->value_at($num , "pedigree");

 my $accession_stock = $schema->resultset("Stock::Stock")->find(
     { name => $stock_name,
     } );
 if (defined($accession_stock)) {
     print STDERR "Found: ".$stock_name."\n";
 }
 else {
     print STDERR "Not found: ".$stock_name."\n";
     $missing_accessions{$stock_name}='';

 }
 pedigree_parse($pedigree_string);

 sub pedigree_parse {
     my $pedigree = shift;

     #deal with case where pedigree is a single accession and no cross is indicated (shouldn't happen?)
     if (!($pedigree =~ m/\//)) {
	 $pedigree_accessions{$pedigree}='';
	 return;
     }


     my $max_slash = '';
     while($pedigree =~ m/(\/+)/g) {
	 if (length($1)>length($max_slash)) {
	     $max_slash=$1;
	 }
     }
     my @split_pedigree = split($max_slash,$pedigree);
     if ($split_pedigree[2]) {
	 die("Incorrect pedigree format\n");
     }
     if ($split_pedigree[0] =~ m/(\/)/g) {
	 pedigree_parse($split_pedigree[0]);
     }
     else {
	 unless ($split_pedigree[0] eq '.' || $split_pedigree[0] eq '?') {
	     $pedigree_accessions{$split_pedigree[0]}='';
	 }
     }
     if ($split_pedigree[1] =~ m/(\/)/g) {
	 pedigree_parse($split_pedigree[1]);
     }
     else {
	 unless ($split_pedigree[1] eq '.' || $split_pedigree[1] eq '?') {
	     $pedigree_accessions{$split_pedigree[1]}='';
	 }
     }
 }
}

foreach my $accession_from_pedigree (keys %pedigree_accessions) {
 my $accession_stock = $schema->resultset("Stock::Stock")->find(
     { name => $accession_from_pedigree,
     } );
 if (defined ($accession_stock)) {
     print STDERR "Found: ".$accession_from_pedigree."\n";
 }
 else {
     print STDERR "Not found: ".$accession_from_pedigree."\n";
     $missing_accessions{$accession_from_pedigree}='';
 }
}

foreach my $accession_to_print (keys %missing_accessions){
 print $accession_to_print."\n";
}



if ($@) { print "An error occured! Rolling backl!\n\n $@ \n\n "; }
else {
    #print "Transaction succeeded! Commiting ! \n\n";
    $dbh->commit();
}





