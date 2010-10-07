#!/usr/bin/perl

=head1 NAME

set_preferred_proteins.pl
A script to decide which is the best protein prediction in the sgn database 

=head1 DESCRIPTION

This program gets all the proteins associated with a unigene (usually, there 
should be one or more ESTScan predictions and one or more 6-frame translations), 
and decides which protein is most likely the best prediction, and sets the 
preferred bit in the database. 

Currently, the best prediction is considered to be the longest prediction, 
but this should be enhanced using other parameters.

=head1 USAGE

set_preferred_proteins -H dbhost -D dbname -u unigene_build

=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory)

=item -U

B<database username>    database username (mandatory)

=item -u

B<unigene build id>     unigene build id to get the protein predicted (mandatory) 

=item -h

B<help>                 print help

=back 

=cut
  
=head1 AUTHOR(S)

Lukas Mueller <lam87@cornell.edu>

Modified by:

Aureliano Bombarely <ab782@cornell.edu>

=cut


use strict;

use Getopt::Std;
use CXGN::DB::DBICFactory;
use SGN::Schema;

use Term::ReadKey;

our ($opt_H, $opt_D, $opt_U, $opt_u, $opt_h);
getopts('H:D:U:u:h');

if (!$opt_H && !$opt_D && !$opt_U && !$opt_u && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Check the arguments

my $dbhost = $opt_H 
    || die("MANDATORY ARGUMENT ERROR: -H <dbhost> argument was not supplied to the script.\n");

my $dbname = $opt_D 
    || die("MANDATORY ARGUMENT ERROR: -D <dbname> argument was not supplied to the script.\n");

my $dbuser = $opt_U
    || die("MANDATORY ARGUMENT ERROR: -U <dbuser> argument was not supplied to the script.\n");

my $unigene_builds = $opt_u 
    || die("MANDATORY ARGUMENT ERROR: -u <unigene_build_id> argument was not supplied to the script.\n");

## First it will get the password
   
print STDERR "\nType password for database user=$dbuser:\npswd> ";

ReadMode('noecho');
my $passw = <>;
chomp($passw);
ReadMode('normal');

print STDERR "\n\n";

## Create a new db_connection


my $schema = CXGN::DB::DBICFactory->open_schema( 
                                                  'SGN::Schema',
						  search_path => ['public','sgn'],
                                                  dbconn_args => {
						                    dbhost => $dbhost,
								    dbname => $dbname,
								    dbuser => $dbuser,
								    dbpass => $passw
                                                                 }
                                                );



## Get the unigene ids, or check if they exists

print STDERR "1) GET THE UNIGENE BUILD DATA.\n";

my @unigenebuild_rows = ();

if ($unigene_builds =~ m/(c)urrent|(p)revious|(d)eprecated/i) {
    my $status = uc($1);
    
    @unigenebuild_rows = $schema->resultset('UnigeneBuild')
	                        ->search({ status => $status }, { order_by => 'unigene_build_id' })
				->all();
}
elsif ($unigene_builds =~ m/^[\d+,?]+$/) {
    @unigenebuild_rows = $schema->resultset('UnigeneBuild')
	                        ->search({ unigene_build_id => { 'IN' => $unigene_builds } }, { order_by => 'unigene_build_id' })
				->all();
}
elsif ($unigene_builds =~ m/all/i) {
    @unigenebuild_rows = $schema->resultset('UnigeneBuild')
	                        ->search( undef, { order_by => 'unigene_build_id' })
				->all();
}
else {
    print STDERR "\nARGUMENT ERROR: -u <unigene_build> has wrong format. It only can be:\n";
    die("digit or list of digits separated by comma, or the words: 'all', 'current', 'previous' or 'deprecated'\n\n");
}

my $build_c = scalar(@unigenebuild_rows);
if ( $build_c == 0) {
    die("\nARGUMENT ERROR: There are not any unigene_build associated with -u $unigene_builds.\n");
}
else {
    print STDERR "\n$build_c unigene_builds were catched for -u $unigene_builds.\n";
}

## Start the transaction

$schema->txn_begin();

## Begin the interaction with the database

my $pep_c = 0;

print STDERR "\n2) GET THE UNIGENE LIST AND CDS LIST.\n";

foreach my $unigene_build_row (@unigenebuild_rows) {

    my ($pref_c, $non_pref_c) = (0, 0);
    
    ## Create a hash to store the method count
    
    my %method_c;
    my $ubuild_prot_c = 0;
    
    ## Get the unigene build id for each column

    my $ubuild_id = $unigene_build_row->get_column('unigene_build_id');
    
    print STDERR "\t<*> Get the unigene_id for unigene_build_id=$ubuild_id\n";

    ## Get all the unigenes and print the report

    my $unigene_rs = $schema->resultset('Unigene')
	                  ->search({ unigene_build_id => $ubuild_id }, { order_by => 'unigene_id' });

    my $unigene_c = $unigene_rs->count();

    print STDERR "\t\t$unigene_c unigenes were catched from unigene_build_id=$ubuild_id\n";

    my @unigene_rows = $unigene_rs->all();

    ## Get all the cds associated with the unigene

    foreach my $unigene_row (@unigene_rows) {
    
	my $unigene_id = $unigene_row->get_column('unigene_id');

	my @cds_rows = $schema->resultset('Cd')
	                      ->search({ unigene_id => $unigene_id }, { order_by => 'cds_id' })
			      ->all();

	## Load all the cds_rows in a hash to do more easy modify the old

	my %pep_preferred;

	foreach my $cds_row (@cds_rows) {

	    $pep_c++;
	    $ubuild_prot_c++;

	    my %cds_data = $cds_row->get_columns();

	    if (exists $method_c{$cds_data{'method'}}) {
		$method_c{$cds_data{'method'}}++;
	    }
	    else {
		$method_c{$cds_data{'method'}}++;
	    }


	    ## Can happens that the cds haven't ant protein_seq (or it can be NULL). Set these always as false.
	    
	    if (defined $cds_data{'protein_seq'} && $cds_data{'protein_seq'} ne 'NULL') {

		if (exists $pep_preferred{$unigene_id}) {
		    if (length($pep_preferred{$unigene_id}->get_column('protein_seq')) > length($cds_data{'protein_seq'}) ) {
		    
			## If the length is less, it will set preferred as false

			$cds_row->set_column( preferred => 'f' ); 
		    }
		    else {

			## Set the preferred as false for the old

			$pep_preferred{$unigene_id}->set_column( preferred => 'f' );
		    
			## and store it

			$pep_preferred{$unigene_id}->update()
	 		                           ->discard_changes();

			## and overwrite with the new one

			$cds_row->set_column( preferred => 't' );
			$pep_preferred{$unigene_id} = $cds_row;
		    }
		}
		else {
		    $cds_row->set_column( preferred => 't' );
		    $pep_preferred{$unigene_id} = $cds_row;
		}
	    }
	    else {
		
		$cds_row->set_column( preferred => 'f' );
	    }
	    
	    ## Now store in the database (the change will be produced only if the data is commited).

	    $cds_row->update()
   		    ->discard_changes();
	}		
	 
	if (defined $pep_preferred{$unigene_id}) {
	    my $pref_pep = $pep_preferred{$unigene_id}->get_column('cds_id');
	    print STDERR "\t\tThe preferred protein for unigene_id=$unigene_id is cds_id=$pref_pep         \r";
	    $pref_c++;
	}
	else {
	    print STDERR "\t\tThere is not preferred protein for unigene_id=$unigene_id\r";
	    $non_pref_c++;
	}
    }

    if ($ubuild_prot_c > 0) {
	print STDERR "\n\t\tUnigene_build_id=$ubuild_id has:\n\t\t\t$pref_c preferred proteins\n\t\t\t";
	print STDERR "$non_pref_c non preferred protein\n\n";
	foreach my $method (keys %method_c) {
	    print STDERR "\t\tThere are $method_c{$method} proteins associated to the method:$method.\n";
	}
    }
    else {
	print STDERR "\n\t\tThere are not any protein associated with this unigene_build.\n\n";
    }

    print STDERR "\n";
}

print STDERR "\nThe preferred field has been set for $pep_c proteins.\nDo you want commit the results (yes|no, default no):\n";

my $answer = <>;
chomp($answer);

if ($answer =~ m/^yes$/i) {
    $schema->txn_commit();
    print STDERR "\nCOMMIT...done, the data has been commited.\n\n";
}
else {
    $schema->txn_rollback();
    print STDERR "\nROLLBACK...done.\n\n";
}

    

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
 
     This program gets all the proteins associated with a unigene (usually, there 
     should be one or more ESTScan predictions and one or more 6-frame translations), 
     and decides which protein is most likely the best prediction, and sets the 
     preferred bit in the database. 

     Currently, the best prediction is considered to be the longest prediction,  
     but this should be enhanced using other parameters.

    Usage: 
   
     set_preferred_proteins [-h] -H <dbhost> -D <dbname> -U <dbuser> -u <unigene_builds>     
    
    Example:
     
      set_preferred_proteins [-h] -H localhost -D sandbox -U postgres -u 22  

      set_preferred_proteins [-h] -H db.sgn.cornell.edu -D cxgn -U postgres -u 'current'  


    Flags:
      -H database hostname    database host (for example localhost or db.sgn.cornell.edu) (mandatory)
      -D database name        database name (for example: sandbox or cxgn) (mandatory)
      -U database username    database username (for example: postgres or web_usr) (mandatory)
      -u unigene builds       unigene build as a digit (22), list of digits (22,24), or the words 'all',
                              'current', 'previous' or 'deprecated'. (mandatory)
      -h print this help

EOF
exit (1);
}


