#!/usr/bin/perl

=head1 NAME

 experiment_dbload.pl
 A script to parse experiment file and load in a database for gem schema (version.0.1.).

=cut

=head1 SYPNOSIS

 experiment_dbload.pl [-h] -U <load_username> -H <dbhost> -D <dbname> -p <platform_file> [-T] [-X]

  To collect the data loaded report into a file:

 experiment_dbload [-h] [-X] -D <dbname> -H <dbhost> -p <platform_file> [-T] > file.log


=head1 EXAMPLE:

 perl experiment_dbload.pl -u aure -H localhost -D sandbox -p TobEA.bs
        
    
=head2 I<Flags:>

=over

=item -e

B<data load file>               data load file in bs format (mandatory).

=item -u

B<user_loader>                  user that load the data (mandatory)

=item -H

B<database_host>                database host (mandatory if you want check the relations in the database)

=item -D

B<database_name>                database name (mandatory if you want check the relations in the database)

=item -X

B<print data load file>         print a template with examples of the data load file in bs format

=item -T

B<run as a test>                run the script as test

=item -h

B<help>                         print the help  

=back

=cut

=head1 DESCRIPTION

    This script parse the experiment_dbload files and load the data into the database. Also youb can run it with -T test mode.


    Note about -T (Test mode): You can run test mode in two ways. The first using -T parameter and the second login to
                               the database as web_usr. In this mode the ids that this script will return comes from
                               the simulation of new _id_seq (get the current id_seq in an object and add using $v++). 
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

experiment_dbload.pl

=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use CXGN::GEM::Schema;
use CXGN::GEM::ExperimentalDesign;
use CXGN::GEM::Experiment;
use CXGN::GEM::Target;
use CXGN::GEM::Hybridization;
use CXGN::Biosource::Sample;
use CXGN::DB::InsertDBH;
use CXGN::Metadata::Metadbdata;

our ($opt_u, $opt_H, $opt_D, $opt_e, $opt_T, $opt_X, $opt_h);
getopts("u:H:D:e:TXh");
if (!$opt_u && !$opt_H && !$opt_D && !$opt_e && !$opt_T && !$opt_X && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
} elsif ($opt_h) {
    help();
} elsif ($opt_X) {
    print_experiment_template();
}


## Define the hash reference that is going to store the data

my $data_parsed_href = {};

## Checking the input arguments

my $loader_username = $opt_u || die("MANDATORY ARGUMENT ERROR: The -u <loader_username> argument was not supplied.\n");
$data_parsed_href->{'loader_username'} = $loader_username;

my $dbname = $opt_D || die("MANDATORY ARGUMENT ERROR: The -D <database_name> argument was not supplied.\n");
my $dbhost = $opt_H || die("MANDATORY ARGUMENT ERROR: The -H <db_hostname> argument was not supplied.\n"); 
my $experiment_file = $opt_e || die("MANDATORY ARGUMENT ERROR: The -e <experiment_dataload_file> argument was not supplied.\n");

## Connecting with the database

my $dbh =  CXGN::DB::InsertDBH->new({ dbname => $dbname, dbhost => $dbhost })->get_actual_dbh();

## The triggers need to set the search path to tsearch2 in the version of psql 8.1
my $psqlv = `psql --version`;
chomp($psqlv);

my $schema_list = 'gem,biosource,metadata,public';
if ($psqlv =~ /8\.1/) {
    $schema_list .= ',tsearch2';
}

print STDERR "\nStep 1: Connect with the database.\n";

my $schema = CXGN::GEM::Schema->connect( sub { $dbh },
                                         { on_connect_do => ["SET search_path TO $schema_list;"] },
                                        );

$data_parsed_href->{'schema'} = $schema;

## Getting the last ids for the different tables to set the database sequences values in case of rollback 
## or something wrong during the test

print STDERR "\nStep 2: Get the last ids for each table.\n";

my $all_last_ids_href = $schema->get_all_last_ids();


## Parse the sample_file and transfer the data to sample objects

print STDERR "\nStep 3: Open and parse the sample file.\n";

open my $ifh, '<', $experiment_file || die("Sorry, but I can not open the input file: $experiment_file.\n");

my $l = 0;

## The input file can store more than one element. Multiple elements will be stored as a hash
## with keys=element_name and values=element object

my ( %expdesign, %experiments, %expdesign_exp_relations, %targets, %exp_target_relations, %hyb, %target_hyb_relations, %platforms );

## Each data field will be defined by $data_type variable, also will be define element name.

my ( $dt, $expdesign_name, $experiment_name, $target_name, $target_element_name, $platform_name, $db_id );

while(<$ifh>) {
		
    $l++; ## Line counter

    ## It will do not read any line that start with #
    unless ($_ =~ m/^#/) {
	
	## First define the data_type or empty the variables defined for each type
	
	if ($_ =~ m/\*DATA_TYPE:\s+\[(\w+?)\]/) {
	    $dt = $1;
	}
	elsif ($_ =~ m/\/\//) {
	    $dt = '';
	    $expdesign_name = '';
	    $experiment_name = '';
	    $target_name = '';
	    $target_element_name = '';
	    $platform_name = '';
	    $db_id = '';
	}
	
	## Create the list of variables created to check if the data was added before

	my $target_list = join(',', keys %targets);
	my $experiment_list = join(',', keys %experiments);
	my $expdesign_list = join(',', keys %expdesign);
	
	## Parse platform and set the data in an object

	if (defined $dt) {
	    if ($dt eq 'experimental_design') {
		if ($_ =~ m/\*EXPDESIGN_NAME:\s+\[(.+?)\]/) {
		    $expdesign_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None expdesign_name data was detailed for experimental_design section.\n");
		    
		    ## Check if the expdesign exists in the database. If exists, it will take the data from the database using
		    ## CXGN::GEM::ExperimentalDesign object

		    my $expdesign_obj;
		    my ($expdesign_row) = $schema->resultset('GeExperimentalDesign')
			                         ->search({ experimental_design_name => $expdesign_name});
		    if (defined $expdesign_row) {
			$expdesign_obj = CXGN::GEM::ExperimentalDesign->new_by_name($schema, $expdesign_name);
		    }
		    else {
			$expdesign_obj = CXGN::GEM::ExperimentalDesign->new($schema);
			$expdesign_obj->set_experimental_design_name($expdesign_name);
		    }
		    
		    ## Add the expdesign object to the techtype hash. If exists this object it will overwrite it. It will made possible
		    ## modifications in the object in the same file.
		    
		    $expdesign{$expdesign_name} = $expdesign_obj;
		}
		elsif ($_ =~ m/\*EXPDESIGN_TYPE:\s+\[(.+?)\]/ ) {
		    my $expdesign_type = $1;
		    
		    ## It will set the type for the object.

		    if (defined $expdesign_type) {
			$expdesign{$expdesign_name}->set_design_type($expdesign_type);
		    }
		}
		elsif ($_ =~ m/\*EXPDESIGN_DESCRIP:\s+\[(.+?)\]/ ) {
		    my $expdesign_description = $1;
		    
		    ## It will set the description for the object.

		    if (defined $expdesign_description) {
			$expdesign{$expdesign_name}->set_description($expdesign_description);
		    }
		}
	    }
	    elsif ($dt eq 'expdesign_dbxref') {
		if ($_ =~ m/\*EXPDESIGN_NAME:\s+\[(.+?)\]/) {
		    $expdesign_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None expdesign_name data was detailed in expdesign_dbxref section.\n");
		    unless (defined $expdesign{$expdesign_name}) {		
			die("MANDATORY DATA ERROR (line $l): None expdesign_name data match with expdesign_list ($expdesign_list).\n");
		    }
		}
		elsif ($_ =~ m/\*DBNAME:\s+\[(.+?)\]/) {
		    my $dbxname = $1 ||
			die("MANDATORY DATA ERROR (line $l): None dbname data was detailed in element_dbxref section.\n");
		    
		    my ($db_row) = $schema->resultset('General::Db')
                	                  ->search( { name => $dbxname } );

		    if (defined $db_row) {
			$db_id = $db_row->get_column('db_id');
		    }
		    else {
			die("MADATORY DATA ERROR (line $l): Dbname=$dbxname do not exists in db.\n");
		    }
		} 
		elsif ($_ =~ m/\*ACCESSIONS:\s+\[(.+?)\]/) {
		    my $accessions = $1 ||
			die("MANDATORY DATA ERROR (line $l): None accessions data was detailed in element_dbxref section.\n");
		
		    my @accessions = split(/,/, $accessions);		    
		    
		    foreach my $acc (@accessions) {		    
			my ($dbxref_row) = $schema->resultset('General::Dbxref')
                		                  ->search( 
			                                    { 
					    		      accession => $acc,
							      db_id     => $db_id,
							    }  
					                  );

			if (defined $dbxref_row) {
			    my $dbxref_id = $dbxref_row->get_column('dbxref_id');
			    
			    my @expdesign_dbxref_ids = $expdesign{$expdesign_name}->get_dbxref_list();
			    
			    my $dbxrefmatch = 0;
		   
			    foreach my $prev_dbxref_id (@expdesign_dbxref_ids) {
				if ($dbxref_id == $prev_dbxref_id) {
				    $dbxrefmatch = 1;
				}
			    }
			    
			    if ($dbxrefmatch == 0) {
				
				$expdesign{$expdesign_name}->add_dbxref( $dbxref_id );
			    }
			    else {
				warn("\nSKIP WARNING: Dbxref-access=$acc exist associated to expdesign:$expdesign_name.\n");
			    }
			}
			else {
			    die("MADATORY DATA ERROR (line $l): Dbxref=$acc do not exists in db.\n");
			}
		    } 
		}
	    }
	    elsif ($dt eq 'expdesign_pub') {
		if ($_ =~ m/\*EXPDESIGN_NAME:\s+\[(.+?)\]/) {
		    $expdesign_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None expdesign_name data was detailed in expdesign_pub section.\n");
		    unless (defined $expdesign{$expdesign_name}) {		
			die("MANDATORY DATA ERROR (line $l): None expdesign_name data match with expdesign_list ($expdesign_list).\n");
		    }
		}
		elsif ($_ =~ m/\*TITLE:\s+\[(.+?)\]/) {
		    my $title = $1 ||
			die("MANDATORY DATA ERROR (line $l): None title data was detailed in expdesign_pub section.\n");
		    
		    ## To know if the publication match with previous added publications, it will compare teh title without
		    ## dots '.', spaces ' ', underlines '_' or dashes '-'

		    my $match;
		    my @pub_list = $expdesign{$expdesign_name}->get_publication_list('title');
		    foreach my $pub (@pub_list) {
			
			my $formated_title = $title;
			$formated_title =~ s/ |\.|-|_//g;

			my $formated_pub = $pub;
			$formated_pub =~ s/ |\.|-|_//g;

			if ($formated_title =~ m/$formated_pub/i) {
			    $match = $pub;
			}
		    }

		    ## If match the publication with the previous ones, it will skip add_publication function with a warn message

		    unless (defined $match) {
			$expdesign{$expdesign_name}->add_publication( { title => $title } );
		    }
		    else {
			warn("\nThe pub=$title match with a previous publication\n\t($match).\n\tIt will not be added.\n");
		    }
		}
	    }
	    elsif ($dt eq 'experiment') {
		if ($_ =~ m/\*EXPDESIGN_NAME:\s+\[(.+?)\]/) {
		    $expdesign_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None expdesign_name data was detailed in experiment section.\n");
		    
		    ## For expdesign, it will search in the hash. If it doesn't exist in the expdesign hash it will take the parsed
		    ## expdesign_name to search in the database. A new negative result will stop the process with die and a positive
		    ## will add the expdesign_object to the techtype hash

		    unless (defined $expdesign{$expdesign_name}) {		
			my $expdesign_obj;
			my ($expdesign_row) = $schema->resultset('GeExperimentalDesign')
	 		                             ->search({ experimental_design_name => $expdesign_name});
			if (defined $expdesign_row) { 
			    $expdesign_obj = CXGN::GEM::ExperimentalDesign->new_by_name($schema, $expdesign_name);
			    $expdesign{$expdesign_name} = $expdesign_obj;
			}
			else {
			    die("MANDATORY DATA ERROR (line $l):Expdesign=$expdesign_name doesn't exist db or has been parsed before.\n");
			}
		    }
		}
		elsif ($_ =~ m/\*EXPERIMENT_NAME:\s+\[(.+?)\]/) {
		    $experiment_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None experiment_name data was detailed in experiment section.\n");
		    
		    ## For the experiment, it will search if it exists into the database. If it is possitive, it will take the data
		    ## using CXGN::GEM::Experiment object. If not, it will create a new one.

		    my $experiment_obj;
		    my ($experiment_row) = $schema->resultset('GeExperiment')
			                          ->search({ experiment_name => $experiment_name});
		    if (defined $experiment_row) {
			$experiment_obj = CXGN::GEM::Experiment->new_by_name($schema, $experiment_name);
		    }
		    else {
			$experiment_obj = CXGN::GEM::Experiment->new($schema);
			$experiment_obj->set_experiment_name($experiment_name);
		    }
		    
		    ## Add the experiment object to the experiment hash. If exists this object it will overwrite it. It will made possible
		    ## modifications in the object in the same file.

		    $experiments{$experiment_name} = $experiment_obj;
		    
		    ## It will define the expdesign-experiment relation adding a key=expdesign_name value=array ref. of experiment_names
		    ## to the relation hash. The script will use this hash to set experimental_design_id in the experiment object after
		    ## store the techtype object.

		    if (defined $expdesign_exp_relations{$expdesign_name}) {
			push @{$expdesign_exp_relations{$expdesign_name}}, $experiment_name;
		    }
		    else {
			$expdesign_exp_relations{$expdesign_name} = [$experiment_name];
		    }
		}
		elsif ($_ =~ m/\*REPLICATES_NR:\s+\[(.+?)\]/) {
		    my $replicates_nr = $1;
		    
		    if (defined $replicates_nr) {
			$experiments{$experiment_name}->set_replicates_nr($replicates_nr);
		    }
		}
		elsif ($_ =~ m/\*COLOUR_NR:\s+\[(.+?)\]/) {
		    my $colour_nr = $1;
		    
		    if (defined $colour_nr) {
			$experiments{$experiment_name}->set_colour_nr($colour_nr);
		    }
		}
		elsif ($_ =~ m/\*EXPERIMENT_DESCRIP:\s+\[(.+?)\]/) {
		    my $experiment_description = $1;
		    
		    if (defined $experiment_description) {
			$experiments{$experiment_name}->set_description($experiment_description);
		    }
		}
		elsif ($_ =~ m/\*CONTACT_NAME:\s+\[(.+?)\]/ ) {
		    my $contact_name = $1;
		    
		    if (defined $contact_name) {
			
			## The contact format should be: Fist_Name{space}Last_Names
			
			my @user = split(/ /, $contact_name);
			my $first_name = shift(@user);
			my $last_name = join(' ', @user);
			    
			my $query = "SELECT sp_person_id FROM sgn_people.sp_person WHERE first_name ILIKE ? AND last_name ILIKE ?";
			my $sth = $dbh->prepare($query);
			$sth->execute($first_name, $last_name);
			my ($sp_person_id) = $sth->fetchrow_array();
			
			if (defined $sp_person_id) {
			    $experiments{$experiment_name}->set_contact_id($sp_person_id);
			}
			else {
			    warn("OPTIONAL DATA WARNING (line $l): Contact_name=$contact_name do not exists db. Skipping data.\n");
			}
		    }
		}
	    }
	    elsif ($dt eq 'experiment_dbxref') {
		if ($_ =~ m/\*EXPERIMENT_NAME:\s+\[(.+?)\]/) {
		    $experiment_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None experiment_name data was detailed in experiment_dbxref section.\n");
		    unless (defined $experiments{$experiment_name}) {		
			die("MANDATORY DATA ERROR (line $l): None experiment_name data match with experiment_list ($experiment_list).\n");
		    }
		}
		elsif ($_ =~ m/\*DBNAME:\s+\[(.+?)\]/) {
		    my $dbxname = $1 ||
			die("MANDATORY DATA ERROR (line $l): None dbname data was detailed in element_dbxref section.\n");
		    
		    my ($db_row) = $schema->resultset('General::Db')
                	                  ->search( { name => $dbxname } );

		    if (defined $db_row) {
			$db_id = $db_row->get_column('db_id');
		    }
		    else {
			die("MADATORY DATA ERROR (line $l): Dbname=$dbxname do not exists in db.\n");
		    }
		} 
		elsif ($_ =~ m/\*ACCESSIONS:\s+\[(.+?)\]/) {
		    my $accessions = $1 ||
			die("MANDATORY DATA ERROR (line $l): None accessions data was detailed in element_dbxref section.\n");
		
		    my @accessions = split(/,/, $accessions);		    
		    
		    foreach my $acc (@accessions) {		    
			my ($dbxref_row) = $schema->resultset('General::Dbxref')
                		                  ->search( 
			                                    { 
					    		      accession => $acc,
							      db_id     => $db_id,
							    }  
					                  );

			if (defined $dbxref_row) {
			    my $dbxref_id = $dbxref_row->get_column('dbxref_id');
			    
			    my @experiment_dbxref_ids = $experiments{$experiment_name}->get_dbxref_list();
			    
			    my $dbxrefmatch = 0;
		   
			    foreach my $prev_dbxref_id (@experiment_dbxref_ids) {
				if ($dbxref_id == $prev_dbxref_id) {
				    $dbxrefmatch = 1;
				}
			    }
			    
			    if ($dbxrefmatch == 0) {
				
				$experiments{$experiment_name}->add_dbxref( $dbxref_id );
			    }
			    else {
				warn("\nSKIP WARNING: Dbxref-access=$acc exist associated to experiment:$experiment_name.\n");
			    }
			}
			else {
			    die("MADATORY DATA ERROR (line $l): Dbxref=$acc do not exists in db.\n");
			}
		    } 
		}
	    }
	    elsif ($dt eq 'target') {
		if ($_ =~ m/\*EXPERIMENT_NAME:\s+\[(.+?)\]/) {
		    $experiment_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None experiment_name data was detailed in target section.\n");
		    
		    ## For experiment, it will search in the hash. If it doesn't exist in the experiment hash it will take the parsed
		    ## experiment_name to search in the database. A new negative result will stop the process with die and a positive
		    ## will add the experiment_object to the experiment hash

		    unless (defined $experiments{$experiment_name}) {		
			my $experiment_obj;
			my ($experiment_row) = $schema->resultset('GeExperiment')
	 		                              ->search({ experiment_name => $experiment_name});
			if (defined $experiment_row) { 
			    $experiment_obj = CXGN::GEM::ExperimentalDesign->new_by_name($schema, $expdesign_name);
			    $experiments{$experiment_name} = $experiment_obj;
			}
			else {
			    die("MANDATORY DATA ERROR (line $l): Exp=$experiment_name doesn't exist db or has been parsed before.\n");
			}
		    }
		}
		elsif ($_ =~ m/\*TARGET_NAME:\s+\[(.+?)\]/) {
		    $target_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None target_name data was detailed in target section.\n");
		    
		    ## For the target, it will search if it exists into the database. If it is possitive, it will take the data
		    ## using CXGN::GEM::Target object. If not, it will create a new one.

		    my $target_obj;
		    my ($target_row) = $schema->resultset('GeTarget')
			                      ->search({ target_name => $target_name});
		    if (defined $target_row) {
			$target_obj = CXGN::GEM::Target->new_by_name($schema, $target_name);
		    }
		    else {
			$target_obj = CXGN::GEM::Target->new($schema);
			$target_obj->set_target_name($target_name);
		    }
		    
		    ## Add the target object to the experiment hash. If exists this object it will overwrite it. It will made possible
		    ## modifications in the object in the same file.

		    $targets{$target_name} = $target_obj;
		    
		    ## It will define the exp-target relation adding a key=experiment_name value=array ref. of target_names
		    ## to the relation hash. The script will use this hash to set experiment_id in the target object after
		    ## store the experiment object.

		    if (defined $exp_target_relations{$experiment_name}) {
			push @{$exp_target_relations{$experiment_name}}, $target_name;
		    }
		    else {
			$exp_target_relations{$experiment_name} = [$target_name];
		    }
		}
	    }
	    elsif ($dt eq 'target_dbxref') {
		if ($_ =~ m/\*TARGET_NAME:\s+\[(.+?)\]/) {
		    $target_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None target_name data was detailed in target_dbxref section.\n");
		    unless (defined $targets{$target_name}) {		
			die("MANDATORY DATA ERROR (line $l): None target_name data match with target_list ($target_list).\n");
		    }
		}
		elsif ($_ =~ m/\*DBNAME:\s+\[(.+?)\]/) {
		    my $dbxname = $1 ||
			die("MANDATORY DATA ERROR (line $l): None dbname data was detailed in element_dbxref section.\n");
		    
		    my ($db_row) = $schema->resultset('General::Db')
                	                  ->search( { name => $dbxname } );

		    if (defined $db_row) {
			$db_id = $db_row->get_column('db_id');
		    }
		    else {
			die("MADATORY DATA ERROR (line $l): Dbname=$dbxname do not exists in db.\n");
		    }
		} 
		elsif ($_ =~ m/\*ACCESSIONS:\s+\[(.+?)\]/) {
		    my $accessions = $1 ||
			die("MANDATORY DATA ERROR (line $l): None accessions data was detailed in element_dbxref section.\n");
		
		    my @accessions = split(/,/, $accessions);		    
		    
		    foreach my $acc (@accessions) {		    
			my ($dbxref_row) = $schema->resultset('General::Dbxref')
                		                  ->search( 
			                                    { 
					    		      accession => $acc,
							      db_id     => $db_id,
							    }  
					                  );

			if (defined $dbxref_row) {
			    my $dbxref_id = $dbxref_row->get_column('dbxref_id');
			    
			    my @target_dbxref_ids = $targets{$target_name}->get_dbxref_list();
			    
			    my $dbxrefmatch = 0;
		   
			    foreach my $prev_dbxref_id (@target_dbxref_ids) {
				if ($dbxref_id == $prev_dbxref_id) {
				    $dbxrefmatch = 1;
				}
			    }
			    
			    if ($dbxrefmatch == 0) {
				
				$targets{$target_name}->add_dbxref( $dbxref_id );
			    }
			    else {
				warn("\nSKIP WARNING: Dbxref-access=$acc exist associated to target:$target_name.\n");
			    }
			}
			else {
			    die("MADATORY DATA ERROR (line $l): Dbxref=$acc do not exists in db.\n");
			}
		    } 
		}
	    }	
	    elsif ($dt eq 'target_element') {
		if ($_ =~ m/\*TARGET_NAME:\s+\[(.+?)\]/) {
		    $target_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None target_name data was detailed in target_element section.\n");
		    unless (defined $targets{$target_name}) {		
			die("MANDATORY DATA ERROR (line $l): None target_name data match with current target_list ($target_list).\n");
		    }
		}
		elsif ($_ =~ m/\*TARGET_ELEMENT_NAME:\s+\[(.+?)\]/) {
		    $target_element_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None target_element_name data was detailed in target_element section.\n");
		    
		    my %target_elements = $targets{$target_name}->get_target_elements();

		    if (defined $target_elements{$target_element_name}) {
			warn("\nSKIP WARNING: Target_element:$target_element_name exists associated to target=$target_name.\n");
		    }
		    else {
			$targets{$target_name}->add_target_element( { target_element_name => $target_element_name } );
		    }
		}
		elsif ($_ =~ m/\*SAMPLE_NAME:\s+\[(.+?)\]/) {
		    my $sample_name = $1;
		    
		    if (defined $sample_name) {
			$targets{$target_name}->edit_target_element(
			                                             $target_element_name, 
				   				     { sample_name => $sample_name }
			                                           );
		    }
		}
		elsif ($_ =~ m/\*PROTOCOL_NAME:\s+\[(.+?)\]/) {
		    my $protocol_name = $1;
		    
		    if (defined $protocol_name) {
			$targets{$target_name}->edit_target_element(
			                                                $target_element_name, 
				   				        { protocol_name => $protocol_name }
			                                              );
		    }
		}
		elsif ($_ =~ m/\*DYE:\s+\[(.+?)\]/) {
		    my $dye = $1;
		    
		    if (defined $dye) {
			$targets{$target_name}->edit_target_element(
			                                             $target_element_name, 
				   				     { dye => $dye }
			                                           );
		    }
		}
	    }
	    elsif ($dt eq 'hybridization') {

		## In the hybridization case it should create the hybridization object before and fill it because there aren't
		## hybridization name to do it
 				    
		if ($_ =~ m/\*TARGET_NAME:\s+\[(.+?)\]/) {
		    $target_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None target_name data was detailed in hybridization section.\n");
		    
		    ## For target, it will search in the hash. If it doesn't exist in the target hash it will take the parsed
		    ## target_name to search in the database. A new negative result will stop the process with die and a positive
		    ## will add the target_object to the target hash

		    unless (defined $targets{$target_name}) {		
			my $target_obj;
			my ($target_row) = $schema->resultset('GeTarget')
	 		                          ->search({ target_name => $target_name});
			if (defined $target_row) { 
			    $target_obj = CXGN::GEM::Target->new_by_name($schema, $target_name);
			    $targets{$target_name} = $target_obj;
			}
			else {
			    die("MANDATORY DATA ERROR (line $l): Target=$target_name doesn't exist db or has been parsed before.\n");
			}
		    }
		}
		elsif ($_ =~ m/\*PLATFORM_NAME:\s+\[(.+?)\]/) {
		    $platform_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None target_name data was detailed in target section.\n");
		    
		    ## For the platform, it will search if it exists into the database. If it is possitive, it will take the data
		    ## using CXGN::GEM::Platform object. If not, it will create a new one.

		    my $platform_obj;
		    my ($platform_row) = $schema->resultset('GePlatform')
			                        ->search({ platform_name => $platform_name});
		    unless (defined $platform_row) {
			die("DATA ERROR: platform_name=$platform_name does not exist into the database.\n");
		    }
		    else {
			$platform_obj = CXGN::GEM::Platform->new_by_name($schema, $platform_name);
			$platforms{$platform_name} = $platform_obj;
		    }
		}
		elsif ($_ =~ m/\*PLATFORM_BATCH:\s+\[(.+?)\]/) {
		    my $platf_batch = $1;
		    
		    if (defined $platf_batch) {
			$hyb{$target_name."+".$platform_name}->set_platform_batch($platf_batch);
		    }
		}
		elsif ($_ =~ m/\*PROTOCOL_NAME:\s+\[(.+?)\]/) {
		    my $protocol_name = $1;
		    
		    if (defined $protocol_name) {
			$hyb{$target_name."+".$platform_name}->set_protocol_by_name($protocol_name);
		    }
		}

		if ($target_name =~ m/\w+/ && $platform_name =~ m/\w+/) {
		    unless (exists $hyb{$target_name."+".$platform_name}) {

			## Search if exists previous hybridizations with this combination platform_id, 
			## target_id

			my ($hyb_row) = $schema->resultset('GeHybridization')
			                       ->search( 
			                                 { 
						           target_id   => $targets{$target_name}->get_target_id,
							   platform_id => $platforms{$platform_name}->get_platform_id,
						         } 
					               );
			if (defined $hyb_row) {
			    $hyb{$target_name."+".$platform_name} = CXGN::GEM::Hybridization->new( $schema, 
												   $hyb_row->get_column(
												                     'hybridization_id'
												                       )
				                                                                 );
			    print STDERR "\nSKIP WARNING: Hybridization=$target_name+$platform_name exist in the database.\n";
			}
			else {
			    $hyb{$target_name."+".$platform_name} = CXGN::GEM::Hybridization->new($schema);			
			
			    ## It will not set target_by_name at this point because probably it is not store
			    ## but it will add the relation to the %target_hyb relations hash

			    $hyb{$target_name."+".$platform_name}->set_platform_by_name($platform_name);			    
			}

			if (defined $target_hyb_relations{$target_name}) {
			    push @{$target_hyb_relations{$target_name}}, $target_name."+".$platform_name;
			}
			else {
			    $target_hyb_relations{$target_name} = [$target_name."+".$platform_name];
			}
		    }
		}
	    }
	}
    }
}

$data_parsed_href->{'expdesign'} = \%expdesign;
$data_parsed_href->{'experiment'} = \%experiments;
$data_parsed_href->{'expdesign_exp_rel'} = \%expdesign_exp_relations;
$data_parsed_href->{'target'} = \%targets;
$data_parsed_href->{'exp_target_rel'} = \%exp_target_relations;
$data_parsed_href->{'hyb'} = \%hyb; 
$data_parsed_href->{'target_hyb_rel'} = \%target_hyb_relations;


## The file should be parsed.
## Test mode, run as evaluation of the code

print STDERR "\nStep 4: Store (or store simulation for Test mode) the experiment data into the database.\n";

if ($opt_T) {
    print STDERR "\nRunning the TEST MODE.\n\n";
    eval {
	
	store_pipeline($data_parsed_href);	    
    };

    ## Print errors if something was wrong during the test

    if ($@) {
	print STDERR "\nTEST ERRORS:\n\n$@\n";
    }

    ## Finally, rollback, because it is a test and set the sequence values

    $schema->txn_rollback;
    $schema->set_sqlseq_values_to_original_state($all_last_ids_href);
    
} else {
    print STDERR "\nRunning the NORMAL MODE.\n\n";
   
    store_pipeline($data_parsed_href);
    
    ## Finally, commit or rollback option
    
    commit_prompt($schema, $all_last_ids_href);
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
      This script load data for experiment data into the gem schema

      In a .bs format (text format with fields delimited by [] ). To print a template file use the -X option and fill following the
      instructions. The load data file can have one or more of these data types.

    
    Usage: 
       experiment_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -p <platform_file> [-T]

      To collect the data loaded report into a file:

       experiment_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -p <platform_file> [-T] > file.log

    Example: 
      perl experiment_dbload.pl -u aure -H localhost -D sandbox -e TobEA.bs


    Flags:
      -u loader username      loader username (mandatory)
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory for check domains in the database)
      -D database name        sandbox or cxgn etc (mandatory for check domains in the database)
      -e experiment file        data load file input file (mandatory)
      -T run as test          run this script as a test
      -X create dataload file create a template for a data load file (follow the instructions to fill it)
      -h this help

EOF
exit (1);
}


################################
### GENERAL DATABASE METHODS ###
################################

=head2 commit_prompt

 Usage: commit_prompt($schema, $set_tableseq_href, $prompt_message, $after_message);
 Desc: ask if you want commit or rollback the changes in the database. If the answer (STDIN) is yes, commit the changes.
 Ret: none
 Args: $dbh (database conection object), $prompt_message and $after_message are object to print a message during prompt
       and after the answer. $initials_tableseq_href is an hash reference with keys=name of the sequence of a concrete
       table and values=current value before the process. It will be used if the option choosed is rollback.
 Side_Effects: print message
 Example: commit_prompt($schema, $initials_table_seq_href);

=cut

sub commit_prompt {
  my ($schema, $seqvalues_href, $prompt_message, $after_message) = @_;
  unless ($prompt_message) {
    $prompt_message = "Commit?\n(yes|no, default no)> ";
  }
  print STDERR $prompt_message;

  ## Ask the question... commit or rollback

  if (<STDIN> =~ m/^y(es)/i) {

      ## If is yes... commit

      print STDERR "Committing...\n\n";
      $schema->txn_commit;
      print "okay.\n";
      if ($after_message) {
	  print STDERR $after_message;
      }
  } else {
      
      ## If it is no, rollback and set the database sequences values to the initial values

      print STDERR "Rolling back...\n\n";
      $schema->set_sqlseq_values_to_original_state($seqvalues_href);
      $schema->txn_rollback;
      print STDERR "done.\n\n";
  }
}

=head2 store_pipeline 

  Usage: store_pipeline($data_parsed_href)
  Desc: store in the db the data parsed
  Ret: none
  Args: $data_parsed_href, a hash reference with the
        data parsed
  Side_Effects: die if something is wrong
  Example: store_pipeline($data_parsed_href);

=cut 

sub store_pipeline {
    my $data_parsed_href = shift ||
	die("DATA INPUT ERROR: none data_parsed hash reference was supplied to store_pipeline function.\n");

    ## First, creation of a metadata object

    my $metadbdata = CXGN::Metadata::Metadbdata->new($schema, $data_parsed_href->{'loader_username'});
    $metadbdata->store();
	
    ## Second, store the expdesign objects, experiment object, target objects and hybridization objects. 
    ## geting the primary_id for each and setting the dependencies in the next objects

    my %expdesigns = %{ $data_parsed_href->{'expdesign'} };
    my %expdesign_exp_rel = %{ $data_parsed_href->{'expdesign_exp_rel'} };
    my %experiments = %{ $data_parsed_href->{'experiment'} };

    foreach my $expdesign_name (sort keys %expdesigns ) {
	$expdesigns{$expdesign_name}->store($metadbdata);
	
	my @exp_list = @{ $expdesign_exp_rel{$expdesign_name} };
	foreach my $exp (@exp_list) {
	    $experiments{$exp}->set_experimental_design_id( $expdesigns{$expdesign_name}->get_experimental_design_id );
	}
    }

    my %exp_target_rel = %{ $data_parsed_href->{'exp_target_rel'} };
    my %targets = %{ $data_parsed_href->{'target'} };

    foreach my $exp_name (sort keys %experiments) {
	$experiments{$exp_name}->store($metadbdata);
	
	my @target_list = @{ $exp_target_rel{$exp_name} };
	foreach my $target (@target_list) {
	    $targets{$target}->set_experiment_id( $experiments{$exp_name}->get_experiment_id );
	}
    }

    my %target_hyb_rel = %{ $data_parsed_href->{'target_hyb_rel'} };
    my %hybs = %{ $data_parsed_href->{'hyb'} };

    foreach my $target_name (sort keys %targets) {
	$targets{$target_name}->store($metadbdata);
	
	if (defined $target_hyb_rel{$target_name}) {
	    my @hyb_list = @{ $target_hyb_rel{$target_name} };
	    foreach my $hyb (@hyb_list) {
		$hybs{$hyb}->set_target_id( $targets{$target_name}->get_target_id );
	    }
	}
    }

    foreach my $hyb_name (sort keys %hybs) {
	$hybs{$hyb_name}->store($metadbdata);
    }

    ## Third, print the sample data stored, to do that, it will get the different objects.

    my $schema = $data_parsed_href->{'schema'};

    print STDERR "\nStep 5: Print data loaded log.\n";

    print STDOUT "\nEXPERIMENTAL_DESIGN:\n";
    foreach my $expdesignname (sort keys %expdesigns ) {
	my $expdesign_obj = CXGN::GEM::ExperimentalDesign->new_by_name($schema, $expdesignname);

	my $expdesign_row = $expdesign_obj->get_geexpdesign_row(); 

	my %expdesign_data = $expdesign_row->get_columns();
	foreach my $col (keys %expdesign_data) {
	    my $expdesign_data = $expdesign_data{$col} || 'undef';
	    print STDOUT "\tgem.ge_experimental_design.$col\t=>\t$expdesign_data\n";
	}
	print STDOUT "\n";

	my @pub_title_list = $expdesign_obj->get_publication_list('title');
	my $expdesign_pub = join(',', @pub_title_list) || 'none';
	print STDOUT "\texpdesign publication list: $expdesign_pub\n";

	my @dbxref_id_list = $expdesign_obj->get_dbxref_list();
	my $expdesign_dbxref = join(',', @dbxref_id_list) || 'none';
	print STDOUT "\texpdesign dbxref_id list: $expdesign_dbxref\n";
    }

    print STDOUT "\nEXPERIMENT:\n";
    foreach my $expname (sort keys %experiments ) {
	my $exp_obj = CXGN::GEM::Experiment->new_by_name($schema, $expname);

	my $exp_row = $exp_obj->get_geexperiment_row(); 

	my %exp_data = $exp_row->get_columns();
	foreach my $col2 (keys %exp_data) {
	    my $exp_data = $exp_data{$col2} || 'undef';
	    print STDOUT "\tgem.ge_experiment.$col2\t=>\t$exp_data\n";
	}
	print STDOUT "\n";

	my @dbxref_id_list2 = $exp_obj->get_dbxref_list();
	my $exp_dbxref = join(',', @dbxref_id_list2) || 'none';
	print STDOUT "\texp dbxref_id list: $exp_dbxref\n";
    }

    my %target_equiv;
    print STDOUT "\nTARGET:\n";
    foreach my $targetname (sort keys %targets ) {
	my $target_obj = CXGN::GEM::Target->new_by_name($schema, $targetname);

	my $target_row = $target_obj->get_getarget_row(); 

	my %target_data = $target_row->get_columns();

	$target_equiv{$targetname} = $target_data{'target_id'};

	foreach my $col3 (keys %target_data) {
	    my $target_data = $target_data{$col3} || 'undef';
	    print STDOUT "\tgem.ge_target.$col3\t=>\t$target_data\n";
	}
	print STDOUT "\n";

	my %target_elements = $target_obj->get_target_elements();
	foreach my $target_element_name (keys %target_elements) {
	    my %target_el_data = %{ $target_elements{$target_element_name} };
	    foreach my $col4 (keys %target_el_data) {
		my $target_el_data = $target_el_data{$col4} || 'undef';
		print STDOUT "\tgem.ge_target_element.$col4\t=>\t$target_el_data\n";
	    }
	}
	print STDOUT "\n";

	my @dbxref_id_list3 = $target_obj->get_dbxref_list();
	my $target_dbxref = join(',', @dbxref_id_list3) || 'none';
	print STDOUT "\ttarget dbxref_id list: $target_dbxref\n";
    }
    
    print STDOUT "\nHYBRIDIZATIONS:\n";
    foreach my $hybname (sort keys %hybs ) {

	my @hybvar = split(/\+/, $hybname);
	my ($hyb_row) = $schema->resultset('GeHybridization')
	                       ->search({ target_id => $target_equiv{$hybvar[0]} });

	my %hyb_data = $hyb_row->get_columns();
	foreach my $col5 (keys %hyb_data) {
	    my $hyb_data = $hyb_data{$col5} || 'undef'; 
	    print STDOUT "\tgem.ge_hybridization.$col5\t=>\t$hyb_data\n";
	}
	print STDOUT "\n";
    }
    
}





=head2 print_experiment_template

  Usage: print_experiment_template()
  Desc: print a experiment_template file
  Ret: none
  Args: none
  Side_Effects: create a file and exit of the script
  Example:

=cut

sub print_experiment_template {
    my $dir = `pwd`;
    chomp($dir);
    my $template_experiment_file = $dir . '/experiment_load_file.bs';
    open my $TFH, '>', $template_experiment_file || die("Sorry, I can not open the file: $template_experiment_file.\n");

    print STDERR "PRINT GEM(Experiment) DATA LOAD FILE option...\n";
    my $info = '# Notes for the Source data load to GEM Experiment data ()
#
#	Data stingency keys:
#
#		Load the data in a flat file with different data types. The data to load will should be between [data]. 
#		 - Mandatory => means that you should have a data in this field to load in the database.
#		 - Optional  => this field could be empty
#		 - Single    => this field can have only one data
#		 - Multiple  => this field can have more than one data separated by comma
#		 
#
#               Note1: The format for contact_name should be [FIRST_NAME{space}LAST_NAMES]
#           
# NOTE FOR DATABASE CURATORS: To load all this data into the database, use the experiment_dbload.pl script
#


#####################################################################################
# EXPERIMENT_DATA ###################################################################
# FIELD	#	    	# DATA #	# DATA_STRINGENCY #	# EXAMPLE_AND_NOTES #
#####################################################################################

## To add a new Experimental Design (MANDATORY): 

*DATA_TYPE:		[]		#mandatory,single	example:[experimental_design]
*EXPDESIGN_NAME:	[]		#mandatory,single 	example:[TobEA]
*EXPDESIGN_TYPE:        []              #mandatory,single       example:[organ development comparissons]
*EXPDESIGN_DESCRIP:     []		#optional,single	example:[Expression analysis of 19 different organs]
//

## To link a dbxref to an Experimental Design (OPTIONAL):
## Note: DBNAME and ACCESSIONS have to exist into db and dbxref tables before add the relation to an experimental design

*DATA_TYPE:             []              #mandatory,single       example:[expdesign_dbxref]
*EXPDESIGN_NAME:        []              #mandatory,single       example:[TobEA]
*DBNAME:                []              #mandatory,single       example:[GEO-Platform]
*ACCESSIONS:            []              #mandatory,multiple     example:[GPL7824]
//

## To link a publication to an Experimental Design (OPTIONAL):
## Note: TITLE have to exist into pub table before add the relation to an experimental design

*DATA_TYPE:             []              #mandatory              example:[expdesign_pub]
*EXPDESIGN_NAME:        []              #mandatory,single       example:[TobEA]
*TITLE:                 []              #mandatory              example:[The example of the pub]
//

## To add a new Experiment (MANDATORY)

*DATA_TYPE:             []              #mandatory,single       example:[experiment]
*EXPDESIGN_NAME:	[]		#mandatory,single 	example:[TobEA]
*EXPERIMENT_NAME:       []              #mandatory,single       example:[TobEA seed]
*REPLICATES_NR:		[]		#mandatory,single	example:[3]
*COLOUR_NR:             []              #optional,single	example:[1]
*EXPERIMENT_DESCRIP:    []              #optional,single        example:[Expression analysis of N. tabacum seeds]
*CONTACT_NAME:          []              #optional,single,note1  example:[Aureliano Bombarely]
//

## To link a dbxref to an Experiment (OPTIONAL):
## Note: DBNAME and ACCESSIONS have to exists into db and dbxref tables before add the relation to an experiment

*DATA_TYPE:             []              #mandatory,single       example:[experiment_dbxref]
*EXPERIMENT_NAME:       []              #mandatory,single       example:[TobEA seed]
*DBNAME:                []              #mandatory,single       example:[GEO-Platform]
*ACCESSIONS:            []              #mandatory,multiple     example:[GPL7825]
//

## To associate a new target (MANDATORY)

*DATA_TYPE:             []              #mandatory              example:[target]
*EXPERIMENT_NAME:	[]		#mandatory,single	example:[TobEA seed]
*TARGET_NAME:           []              #mandatory              example:[Atlas 76.CEL]
//

## To link a dbxref to a Target (OPTIONAL):
## Note: DBNAME and ACCESSIONS have to exists into db and dbxref tables before add the relation to a target

*DATA_TYPE:             []              #mandatory,single       example:[target_dbxref]
*TARGET_NAME:           []              #mandatory,single       example:[Atlas 76.CEL]
*DBNAME:                []              #mandatory,single       example:[GEO-Platform]
*ACCESSIONS:            []              #mandatory,multiple     example:[GPL7826]
//

## To add a target element to a target (MANDATORY):

*DATA_TYPE:             []              #mandatory,single       example:[target_element]
*TARGET_NAME:           []              #mandatory,single       example:[Atlas 76.CEL]
*TARGET_ELEMENT_NAME:	[]		#mandatory,single	example:[Atlas 76.CEL]
*SAMPLE_NAME:           []              #mandatory,single       example:[ExpressionAtlas_01_seed_03]
*PROTOCOL_NAME:         []              #optional,single        example:[Microarray target synthesis for TobEA]
*DYE:                   []              #optional,single        example:[streptavidin phycoerythrin conjugate]
//  

## To associate a new hybridation (MANDATORY)

*DATA_TYPE:	        []              #mandatory,single	example:[hybridization]
*TARGET_NAME:           []              #mandatory,single       example:[Atlas 76.CEL]
*PLATFORM_NAME:         []              #mandatory,single       example:[Affymetrix TobEA]
*PLATFORM_BATCH:        []              #optional,single        example:[AtlasBatch05]
*PROTOCOL_NAME:         []              #optional,single        example:[Microarray hybridization for TobEA]
//';

    print $TFH "$info";
    print STDERR "...done (printed platform_file with the name: $template_experiment_file)\n\n";
    exit (1);
}




####
1;##
####
