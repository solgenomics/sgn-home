#!/usr/bin/perl

=head1 NAME

 load_domainmatch_into_sgndb.pl
 A script to load the domain matches produced by an InterProScan analysis into db

=cut

=head1 SYPNOSIS

load_domainmatch_into_sgndb.pl [-h] -H <dbhost> -D <dbname> -i <input_file> [-T]
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory)

=item -i

B<input_file>           input file with 14 columns

=item -u

B<sgn_user>             sgn user that is going to load the data

=item -T

B<run as test>          run as test option

=item -h 

B<help>                  show the help

=back

=cut

=head1 DESCRIPTION

 This script load the domain matched produced by a InterProScan analysis in 
 raw format (a basic tab delimited format) with the following columns:

  + id of the input sequence (SGN-P accession).
  + crc64 (checksum) of the proteic sequence (supposed to be unique).
  + length of the sequence (in AA).
  + anaysis method launched.
  + database members entry for this match.
  + database member description for the entry.
  + start of the domain match (mandatory).
  + end of the domain match (mandatory). 
  + evalue of the match (reported by member database method) (mandatory).
  + status of the match (T: true, ?: unknown) (mandatory).
  + date of the run.
  + corresponding InterPro entry (if iprlookup requested by the user)(mandatory)
  + description of the InterPro entry.
  + GO (gene ontology) description for the InterPro entry.


=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

load_domainmatch_into_sgndb.pl


=cut

use strict;
use Getopt::Std;
use CXGN::DB::InsertDBH;
use Math::BigFloat;
use File::Basename;

our ($opt_H, $opt_D, $opt_i, $opt_u, $opt_T, $opt_h);
getopts("H:D:i:u:Th");

if (!$opt_H && !$opt_D && !$opt_i && !$opt_u && !$opt_T && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
elsif ($opt_h) {
    help();
}


## Parse the input file.

if (!$opt_H || !$opt_D) {
    die("\nArgument error: Database argument -D <database_name> and/or -H <hostname> were not supplied.\n");
}

if (!$opt_u) {
    die("Argument error: None -u <sgn_user> was supplied to the script.\n");
}

unless (-s $opt_i) {
    die("\nArgument error: The -i <input_file> does not exists or has zero size.\n");
}

open my $input_fh, '<', $opt_i || die("Sorry, I can not open the file $opt_i.\n");

my (%data_parsed, %ipr);

my $row = 0;

print STDERR "\nParsing input file.\n";

while (<$input_fh>) {
    $row++;
   
    my @data = split(/\t/, $_);
    print STDERR "...processing line $row with id=$data[0] ...\t\t\t\r";
    
    ## Checking the mandarory variables in the ipr result file

    my $id = $data[0] || die("\n\nParsing error: Line $row have not any id\n");
    my $begin_match = $data[6];
    my $end_match = $data[7]; 

    ## Use BigInt with e-value to return in always in e-value format.
    my $e_value;
    my $e_value_number = Math::BigFloat->new($data[8]);
    Math::BigFloat->accuracy(3);
    if ($e_value_number->is_nan() ) {
	$e_value = "'NA'";
    }
    else {
	$e_value = "'" . $e_value_number->bsstr() . "'";
	$e_value =~ s/E/e/;
    }
    my $status = $data[9] || die("\n\nParsing error: Line $row have not any status\n");
    my $ipr = $data[4] || die("\n\nParsing error: Line $row have not any interpro accession\n");


    ## Store the data as array reference in a hash with a key = cds_id
    
    if ($id =~ m/SGN-P(\d+)/) {
	if ( exists $data_parsed{$1} ) {
	    push @{$data_parsed{$1}}, [$id, $begin_match, $end_match, $e_value, $status, $ipr];
	}
	else {
	    $data_parsed{$1} = [ [$id, $begin_match, $end_match, $e_value, $status, $ipr] ];
	}
    }
    else {
	die("The protein ID =$id is not SGN-P accession.\n");
    }

    ## Get a list of all the different ipr_accessions using 1 as value, 
    ## it will replace by an id after the data parsing

    unless (exists $ipr{$ipr}) {
	$ipr{$ipr} = 1;
    }
}

## Now it will create a db connection

print STDERR "\n\nConnecting with the database (Database:$opt_D Hostname:$opt_H)\n";


my $dbh = CXGN::DB::InsertDBH->new(
                                     {       
					 dbname=>$opt_D ,
					 dbhost=>$opt_H ,
					 dbargs=>{ RaiseError => 1, PrintError => 1 }
				     }
                                   ) or die "Failed to connect to database ($DBI::errstr)";


## Now it will get the user_id for the sgn_username

print STDERR "\nGetting the user_id associated to the sgn_username:$opt_u...\n";

my ($user_id) = $dbh->selectrow_array("SELECT sp_person_id FROM sgn_people.sp_person WHERE username = ?", undef, $opt_u);
unless (defined $user_id) {
    die("\nDatabase coherence error: username=$opt_u does not exist into the database.\n");
}

## Get the unigene_id associated to each cds_id

print STDERR "\nGetting the unigene_id for each SGN-P accession...\n";

foreach my $cds_id (keys %data_parsed) {
    my $query1 = "SELECT unigene_id FROM sgn.cds WHERE cds_id=?";

    my ($unigene_id) = $dbh->selectrow_array($query1, undef, $cds_id);
    unless (defined $unigene_id) {
	die("\nDatabase coherence error: cds_id=$cds_id has not any unigene_id associated\n");
    }
    else {
	unshift @{$data_parsed{$cds_id}}, $unigene_id;
    }
}

## Get the domain_id associated to each ipr accession

print STDERR "\nGetting the domain_id for each ipr_accession...\n";

my ($iprdb, $iprnondb) = (0, 0);
foreach my $ipracc (sort keys %ipr) {
    my $query2 = "SELECT domain_id FROM sgn.domain WHERE domain_accession = ?";
    
    my ($domain_id) = $dbh->selectrow_array($query2, undef, $ipracc);
    unless (defined $domain_id) {
	$iprnondb++;
	print STDERR "\nThe interpro_accession=$ipracc do not exists into the sgn.domain table.\n";
    }
    else {
	$iprdb++;
	$ipr{$ipracc} = $domain_id;
    }
}
print STDERR "\nThere are:\n\t$iprnondb domain accessions that ARE NOT IN THE DB\n\tand $iprdb that EXISTS into the db\n";

## Now it will take the last domain_match_id to set to rollback


my ($last_domain_match_id) = $dbh->selectrow_array("SELECT max(domain_match_id) FROM sgn.domain_match");
my ($last_sgnmetadata_id) = $dbh->selectrow_array("SELECT max(metadata_id) FROM sgn.metadata");

unless (defined $last_domain_match_id) {
    $last_domain_match_id = 0;
}
print STDERR "\nGet the last_domain_match_id=$last_domain_match_id\n";

my ($last_run_id) = $dbh->selectrow_array("SELECT max(run_id) FROM sgn.domain_match");
my $curr_run_id = $last_run_id+1;


my ($user) = $dbh->selectrow_array("SELECT user");
if ($user eq 'web_usr') {
    die("\nInsert functions only can be used by postgres user.\n");
}

if ($opt_T) {
    print STDERR "Running the script with -T option (as test).\n";
    eval {

	## Create a new metadata_id
	my $metadata_query = "INSERT INTO sgn.metadata (create_person_id) VALUES ($user_id)";
	$dbh->do($metadata_query);
	my ($insert_sgnmetadata_id) = $dbh->selectrow_array("SELECT max(metadata_id) FROM sgn.metadata");
	if ($insert_sgnmetadata_id == $last_sgnmetadata_id) {
	    die("\nDatabase insert error: None new sgn.metadata row was inserted.\n");
	}

	foreach my $sgn_id (keys %data_parsed) {
	    my @ipr_result = @{$data_parsed{$sgn_id}};
	    my $unigene_id = shift(@ipr_result);

	    foreach my $ipr_line (@ipr_result) {
		my $i_id = $ipr_line->[0];
		my $i_begin_match = $ipr_line->[1];
		my $i_end_match = $ipr_line->[2];
		my $i_e_value = $ipr_line->[3];
		my $i_status = "'" . $ipr_line->[4] . "'";  ## Added single quotes to insert text values in the db
		my $i_domain_id = $ipr{$ipr_line->[5]};

		if (defined $i_domain_id) {
		    my $insert = "INSERT INTO sgn.domain_match ( cds_id, unigene_id, domain_id, match_begin,
                                                                 match_end, e_value, hit_status, run_id, 
                                                                 metadata_id) 
                                                        VALUES ( $sgn_id, $unigene_id, $i_domain_id, $i_begin_match, 
                                                                 $i_end_match,$i_e_value, $i_status, $curr_run_id, 
                                                                 $insert_sgnmetadata_id )";
		    $dbh->do($insert);
		}
	    }
	}
    };
    if ($@) {
	print STDERR "\n\nERROR: $@\n";
    }
    else {
	print STDERR "\n\nRunning as test did not report any error, rolling back...\n";
    }
    $dbh->rollback;
    $dbh->do("SELECT setval ('sgn.domain_match_domain_match_id_seq', $last_domain_match_id, 'true')");
    $dbh->do("SELECT setval ('sgn.metadata_metadata_id_seq', $last_sgnmetadata_id, 'true')");
}
else {

    ## Create a new metadata_id
    my $metadata_query = "INSERT INTO sgn.metadata (create_person_id) VALUES ($user_id)";
    $dbh->do($metadata_query);
    my ($insert_sgnmetadata_id) = $dbh->selectrow_array("SELECT max(metadata_id) FROM sgn.metadata");
    if ($insert_sgnmetadata_id == $last_sgnmetadata_id) {
	die("\nDatabase insert error: None new sgn.metadata row was inserted.\n");
    }

    foreach my $sgn_id (keys %data_parsed) {
	my @ipr_result = @{$data_parsed{$sgn_id}};
	my $unigene_id = shift(@ipr_result);
	
	foreach my $ipr_line (@ipr_result) {
	    my $i_id = $ipr_line->[0];
	    my $i_begin_match = $ipr_line->[1];
	    my $i_end_match = $ipr_line->[2];
	    my $i_e_value = $ipr_line->[3];
	    my $i_status = "'" . $ipr_line->[4] . "'";  ## Added single quotes to insert text values in the db
	    my $i_domain_id = $ipr{$ipr_line->[5]};
	    
	    if (defined $i_domain_id) {
		my $insert = "INSERT INTO sgn.domain_match ( cds_id, 
                                                             unigene_id, 
                                                             domain_id, 
                                                             match_begin, 
                                                             match_end,
                                                             e_value, 
                                                             hit_status, 
                                                             run_id,
                                                             metadata_id ) 
                                                    VALUES ( $sgn_id, 
                                                             $unigene_id, 
                                                             $i_domain_id, 
                                                             $i_begin_match, 
                                                             $i_end_match,
                                                             $i_e_value,
                                                             $i_status, 
                                                             $curr_run_id, 
                                                             $insert_sgnmetadata_id )";
		$dbh->do($insert);
	    }
	}
    }
    
    print STDERR "Commit?\n(yes|no, default no)> ";
    if (<STDIN> =~ m/^y(es)/i) {
	print STDERR "Committing...";
	$dbh->commit;
	print "okay.\n";
    } 
    else {
	print STDERR "Rolling back...";
	$dbh->rollback;

	print STDERR "\tsetting seq to original values...";
	$dbh->do("SELECT setval ('sgn.domain_match_domain_match_id_seq', $last_domain_match_id, 'true')");
	$dbh->do("SELECT setval ('sgn.metadata_metadata_id_seq', $last_sgnmetadata_id, 'true')");
	
	print STDERR "done.\n";
    }
}


my ($curr_domain_match_id) = $dbh->selectrow_array("SELECT max(domain_match_id) FROM sgn.domain_match");

my $matches_inserted = $curr_domain_match_id - $last_domain_match_id;
print STDERR "\nREPORT: $matches_inserted new matches have been inserted in the sgn.domain_match table.\n\n";    


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
      This script load the domain matched produced by a InterProScan analysis in 
      raw format (a basic tab delimited format) with the following columns:

	  + id of the input sequence (SGN-P accession).
	  + crc64 (checksum) of the proteic sequence (supposed to be unique).
	  + length of the sequence (in AA).
	  + anaysis method launched.
	  + database members entry for this match.
	  + database member description for the entry.
	  + start of the domain match (mandatory).
	  + end of the domain match (mandatory). 
	  + evalue of the match (reported by member database method) (mandatory).
	  + status of the match (T: true, ?: unknown) (mandatory).
	  + date of the run.
	  + corresponding InterPro entry (if iprlookup requested by the user)(mandatory)
	  + description of the InterPro entry.
	  + GO (gene ontology) description for the InterPro entry.

    Usage: 
      load_domainmatch_into_sgndb.pl [-h] -H <dbhost> -D <dbname> -i <input_file> [-T]
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -i input_file           raw file with the interpro analysis results (14 columns)(mandatory)
      -u sgn_user             sgn username that is going to load the data 
      -T run as test          run as a test option (finally it rollback)
      -h print this help

EOF
exit (1);
}


