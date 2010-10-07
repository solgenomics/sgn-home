
#!/usr/bin/perl
=head1 NAME

 platform_dbload.pl
 A script to parse platform file and load in a database for gem schema (version.0.1.).

=cut

=head1 SYPNOSIS

 platform_dbload.pl [-h] -U <load_username> -H <dbhost> -D <dbname> -p <platform_file> [-T] [-X]

  To collect the data loaded report into a file:

 platform_dbload [-h] [-X] -D <dbname> -H <dbhost> -p <platform_file> [-T] > file.log


=head1 EXAMPLE:

 perl platform_dbload.pl -u aure -H localhost -D sandbox -p platform-tom2.bs
        
    
=head2 I<Flags:>

=over

=item -p

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

    This script parse the platform_dbload files and load the data into the database. Also youb can run it with -T test mode.

    To load other platform associated data as template, use template_dbload.pl script.


    Note about -T (Test mode): You can run test mode in two ways. The first using -T parameter and the second login to
                               the database as web_usr. In this mode the ids that this script will return comes from
                               the simulation of new _id_seq (get the current id_seq in an object and add using $v++). 
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

platform_dbload.pl


=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use CXGN::GEM::Schema;
use CXGN::GEM::TechnologyType;
use CXGN::GEM::Platform;
use CXGN::Biosource::Sample;
use CXGN::DB::InsertDBH;
use CXGN::Metadata::Metadbdata;

our ($opt_u, $opt_H, $opt_D, $opt_p, $opt_T, $opt_X, $opt_h);
getopts("u:H:D:p:TXh");
if (!$opt_u && !$opt_H && !$opt_D && !$opt_p && !$opt_T && !$opt_X && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
} elsif ($opt_h) {
    help();
} elsif ($opt_X) {
    print_platform_template();
}

## Checking the input arguments

my $loader_username = $opt_u || die("MANDATORY ARGUMENT ERROR: The -u <loader_username> argument was not supplied.\n");
my $dbname = $opt_D || die("MANDATORY ARGUMENT ERROR: The -D <database_name> argument was not supplied.\n");
my $dbhost = $opt_H || die("MANDATORY ARGUMENT ERROR: The -H <db_hostname> argument was not supplied.\n"); 
my $platform_file = $opt_p || die("MANDATORY ARGUMENT ERROR: The -p <platform_dataload_file> argument was not supplied.\n");

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

## Getting the last ids for the different tables to set the database sequences values in case of rollback 
## or something wrong during the test

print STDERR "\nStep 2: Get the last ids for each table.\n";

my $all_last_ids_href = $schema->get_all_last_ids();


## Parse the sample_file and transfer the data to sample objects

print STDERR "\nStep 3: Open and parse the sample file.\n";

open my $ifh, '<', $platform_file || die("Sorry, but I can not open the input file: $platform_file.\n");

my $l = 0;

## The input file can store more than one sample. Multiple samples will be stored as a hash
## with keys=sample_name and values=sample object

my %techtypes;
my %platforms;
my %techtype_platform_relations;

## Each data field will be defined by $data_type variable, also will be define $sample_name and $sample_element_name.

my ($dt, $techtype_name, $platform_name, $db_id);

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
	    $techtype_name = '';
	    $platform_name = '';
	    $db_id = '';
	}
	
	## Create the list of variables created to check if the data was added before

	my $platform_list = join(',', keys %platforms);
	my $techtype_list = join(',', keys %techtypes);
	
	## Parse platform and set the data in an object

	if (defined $dt) {
	    if ($dt eq 'technology_type') {
		if ($_ =~ m/\*TECHNOLOGY_NAME:\s+\[(.+?)\]/) {
		    $techtype_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None technology_name data was detailed for technology_type section.\n");
		    

		    ## Check if the techtype exists in the database. If exists, it will take the data from the database using
		    ## CXGN::GEM::TechnologyType object

		    my $techtype_obj;
		    my ($techtype_row) = $schema->resultset('GeTechnologyType')
			                        ->search({ technology_name => $techtype_name});
		    if (defined $techtype_row) {
			$techtype_obj = CXGN::GEM::TechnologyType->new_by_name($schema, $techtype_name);
		    }
		    else {
			$techtype_obj = CXGN::GEM::TechnologyType->new($schema);
			$techtype_obj->set_technology_name($techtype_name);
		    }
		    
		    ## Add the techtype object to the techtype hash. If exists this object it will overwrite it. It will made possible
		    ## modifications in the object in the same file.
		    
		    $techtypes{$techtype_name} = $techtype_obj;
		}
		elsif ($_ =~ m/\*TECHNOLOGY_DESCRIP:\s+\[(.+?)\]/ ) {
		    my $techtype_description = $1;
		    
		    ## It will set the description for the object.

		    if (defined $techtype_description) {
			$techtypes{$techtype_name}->set_description($techtype_description);
		    }
		}
	    }
	    elsif ($dt eq 'platform') {
		if ($_ =~ m/\*TECHNOLOGY_NAME:\s+\[(.+?)\]/) {
		    $techtype_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None technology_name data was detailed in platform section.\n");
		    
		    ## For techtype, it will search in the hash. If it doesn't exist in the techtype hash it will take the parsed
		    ## techtype_name to search in the database. A new negative result will stop the process with die and a positive
		    ## will add the techtype_object to the techtype hash

		    unless (defined $techtypes{$techtype_name}) {		
			my $techtype_obj;
			my ($techtype_row) = $schema->resultset('GeTechnologyType')
			                            ->search({ technology_name => $techtype_name});
			if (defined $techtype_row) { 
			    $techtype_obj = CXGN::GEM::TechnologyType->new_by_name($schema, $techtype_name);
			    $techtypes{$techtype_name} = $techtype_obj;
			}
			else {
			    die("MANDATORY DATA ERROR (line $l): TechType=$techtype_name does not exist db or has been parsed before.\n");
			}
		    }
		}
		elsif ($_ =~ m/\*PLATFORM_NAME:\s+\[(.+?)\]/) {
		    $platform_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None platform_name data was detailed in platform section.\n");
		    
		    ## For the platform, it will search if it exists into the database. If it is possitive, it will take the data
		    ## using CXGN::GEM::Platform object. If not, it will create a new one.

		    my $platform_obj;
		    my ($platform_row) = $schema->resultset('GePlatform')
			                        ->search({ platform_name => $platform_name});
		    if (defined $platform_row) {
			$platform_obj = CXGN::GEM::Platform->new_by_name($schema, $platform_name);
		    }
		    else {
			$platform_obj = CXGN::GEM::Platform->new($schema);
			$platform_obj->set_platform_name($platform_name);
		    }
		    
		    ## Add the platform object to the platform hash. If exists this object it will overwrite it. It will made possible
		    ## modifications in the object in the same file.

		    $platforms{$platform_name} = $platform_obj;
		    
		    ## It will define the techtype-platform relation adding a key=techtype_name value=array ref. of platform_names
		    ## to the relation hash. The script will use this hash to set technology_type_id in the platform object after
		    ## store the techtype object.

		    if (defined $techtype_platform_relations{$techtype_name}) {
			push @{$techtype_platform_relations{$techtype_name}}, $platform_name;
		    }
		    else {
			$techtype_platform_relations{$techtype_name} = [$platform_name];
		    }
		}
		elsif ($_ =~ m/\*PLATFORM_DESCR:\s+\[(.+?)\]/) {
		    my $platform_description = $1;
		    
		    if (defined $platform_description) {
			$platforms{$platform_name}->set_description($platform_description);
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
			    $platforms{$platform_name}->set_contact_id($sp_person_id);
			}
			else {
			    warn("OPTIONAL DATA WARNING (line $l): Contact_name=$contact_name do not exists db. Skipping data.\n");
			}
		    }
		}
	    }
	    elsif ($dt eq 'pub') {
		if ($_ =~ m/\*PLATFORM_NAME:\s+\[(.+?)\]/) {
		    $platform_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None platform_name data was detailed in pub section.\n");
		    unless (defined $platforms{$platform_name}) {		
			die("MANDATORY DATA ERROR (line $l): None platform_name data match with curr. platform_list ($platform_list).\n");
		    }
		}
		elsif ($_ =~ m/\*TITLE:\s+\[(.+?)\]/) {
		    my $title = $1 ||
			die("MANDATORY DATA ERROR (line $l): None title data was detailed in pub section.\n");
		    
		    ## To know if the publication match with previous added publications, it will compare teh title without
		    ## dots '.', spaces ' ', underlines '_' or dashes '-'

		    my $match;
		    my @pub_list = $platforms{$platform_name}->get_publication_list('title');
		    foreach my $pub (@pub_list) {
			
			my $formated_title = $title;
			$formated_title =~ s/ //g;
			$formated_title =~ s/\.//g;
			$formated_title =~ s/-//g;
			$formated_title =~ s/_//g;

			my $formated_pub = $pub;
			$formated_pub =~ s/ //g;
			$formated_pub =~ s/\.//g;
			$formated_pub =~ s/-//g;
			$formated_pub =~ s/_//g;

			if ($formated_title =~ m/$formated_pub/i) {
			    $match = $pub;
			}
		    }

		    ## If match the publication with the previous ones, it will skip add_publication function with a warn message

		    unless (defined $match) {
			$platforms{$platform_name}->add_publication( { title => $title } );
		    }
		    else {
			warn("\nThe pub=$title match with a previous publication\n\t($match).\n\tIt will not be added.\n");
		    }
		}
	    }
	    elsif ($dt eq 'platform_design') {
		if ($_ =~ m/\*PLATFORM_NAME:\s+\[(.+?)\]/) {
		    $platform_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None platform_name data was detailed in platform_design section.\n");
		    unless (defined $platforms{$platform_name}) {		
			die("MANDATORY DATA ERROR (line $l): None platform_name data match with curr. platform_list ($platform_list).\n");
		    }
		}
		if ($_ =~ m/\*SAMPLE_NAME:\s+\[(.+?)\]/) {
		    my $sample_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None sample_name data was detailed in platform_Design section.\n");
		    
		    ## If the sample_name do not exists into the database the add_platform_design function will die
		    
		    ## Before add a new platform_design, it will check if exists associated with this platform

		    my $designmatch = 0;
		    my @sample_name_list = $platforms{$platform_name}->get_design_list('sample_name');
		    foreach my $prev_sample (@sample_name_list) {
			if ($prev_sample eq $sample_name) {
			    $designmatch = 1;
			}
		    }

		    if ($designmatch == 0) {
			$platforms{$platform_name}->add_platform_design($sample_name);
		    }
		    else {
			warn("\nSKIP WARNING: sample_name=$sample_name exist associated to platform:$platform_name.\n");
		    }
		}
	    }
	    elsif ($dt eq 'dbxref') {
		if ($_ =~ m/\*PLATFORM_NAME:\s+\[(.+?)\]/) {
		    $platform_name = $1 ||
			die("MANDATORY DATA ERROR (line $l): None platform_name data was detailed in platform_design section.\n");
		    unless (defined $platforms{$platform_name}) {		
			die("MANDATORY DATA ERROR (line $l): None platform_name data match with curr. platform_list ($platform_list).\n");
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
			    
			    my @platform_dbxref_ids = $platforms{$platform_name}->get_dbxref_list();
			    
			    my $dbxrefmatch = 0;
		   
			    foreach my $prev_dbxref_id (@platform_dbxref_ids) {
				if ($dbxref_id == $prev_dbxref_id) {
				    $dbxrefmatch = 1;
				}
			    }
			    
			    if ($dbxrefmatch == 0) {
				
				$platforms{$platform_name}->add_dbxref( $dbxref_id );
			    }
			    else {
				warn("\nSKIP WARNING: Dbxref-access=$acc exist associated to platform:$platform_name.\n");
			    }
			}
			else {
			    die("MADATORY DATA ERROR (line $l): Dbxref=$acc do not exists in db.\n");
			}
		    } 
		}
	    }
	}
    }
}

## The file should be parsed.
## Test mode, run as evaluation of the code

print STDERR "\nStep 4: Store (or store simulation for Test mode) the sample data into the database.\n";

if ($opt_T) {
    print STDERR "\nRunning the TEST MODE.\n\n";
    eval {

	## First, creation of a metadata object

	my $metadbdata = CXGN::Metadata::Metadbdata->new($schema, $loader_username);
	$metadbdata->store();
	
	## Second, store the technology objects and the platform object. After store technology_type and before store 
	## platform, it will use the %techtype_platform_relations hash to set the technology_type related with each platform

	foreach my $techtypename (sort keys %techtypes) {
	    $techtypes{$techtypename}->store($metadbdata);
	    
	    my @platf_list = @{ $techtype_platform_relations{$techtypename} };
	    foreach my $platf (@platf_list) {
		$platforms{$platf}->set_technology_type_id( $techtypes{$techtypename}->get_technology_type_id );
	    }
	}

	foreach my $platformname (sort keys %platforms) {
	    $platforms{$platformname}->store($metadbdata);
	}

	## Third, print the sample data stored

	print STDERR "\nStep 5: Print data loaded log.\n";
	print STDOUT "\nDATA NOT STORED in TEST MODE: \n";

	foreach my $techtypename (sort keys %techtypes) {
	    my $t_techtype_id = $techtypes{$techtypename}->get_technology_type_id() || 'undef';
	    my $t_techtype_name = $techtypes{$techtypename}->get_technology_name() || 'undef';
	    my $t_techtype_description = $techtypes{$techtypename}->get_description() || 'undef';
	    my $t_tt_metadata_id = $techtypes{$techtypename}->get_technology_type_metadbdata()->get_metadata_id() || 'undef';
	    print STDOUT "+ TECHNOLOGY_TYPE_DATA:\n\ttechnology_type_id =\t$t_techtype_id\n";
	    print STDOUT "\ttechnology_name =\t$t_techtype_name\n";
	    print STDOUT "\tdescription =\t$t_techtype_description\n\n";
	    print STDOUT "\tmetadata_id =\t$t_tt_metadata_id\n\n";
	}

	foreach my $platformname (sort keys %platforms) {
	    my $t_platform_id = $platforms{$platformname}->get_platform_id() || 'undef';
	    my $t_platform_name = $platforms{$platformname}->get_platform_name() || 'undef';
	    my $t_technology_type_id = $platforms{$platformname}->get_technology_type_id() || 'undef';
	    my $t_description =  $platforms{$platformname}->get_description() || 'undef';
	    my $t_contact_id = $platforms{$platformname}->get_contact_id() || 'undef';
	    my $t_p_metadata_id = $platforms{$platformname}->get_platform_metadbdata()->get_metadata_id() || 'undef';
	    print STDOUT "+ PLATFORM_DATA:\n\tplatform_id =\t$t_platform_id\n";
	    print STDOUT "\tplatform_name =\t$t_platform_name\n";
	    print STDOUT "\ttechnology_type_id =\t$t_technology_type_id\n";
	    print STDOUT "\tdescription =\t$t_description\n";
	    print STDOUT "\tcontact_id =\t$t_contact_id\n\n";
	    print STDOUT "\tmetadata_id =\t$t_p_metadata_id\n\n";
	
	    print STDOUT "  * Associated publications:\n";
	    my @pub_title_list = $platforms{$platformname}->get_publication_list('title');
	    if (scalar(@pub_title_list) == 0) {
		push @pub_title_list, 'none';
	    }
	    foreach my $pub (@pub_title_list) {
		print STDOUT "\t\tpub_title: $pub\n";
	    }

	    print STDOUT "\n";
	    print STDOUT "  * Associated platform_design:\n";
	    my @design_sample_names = $platforms{$platformname}->get_design_list('sample_name');
	    if (scalar(@design_sample_names) == 0) {
		push @design_sample_names, 'none';
	    }
	    foreach my $samplename (@design_sample_names) {
		print STDOUT "\t\tdesign (sample_name): $samplename\n";
	    }

	    print STDOUT "\n";
	    print STDOUT "  * Associated dbxref_ids:\n";
	    my @dbxref_id_list = $platforms{$platformname}->get_dbxref_list();
	    if (scalar(@dbxref_id_list) == 0) {
		push @dbxref_id_list, 'none';
	    }
	    foreach my $dbxref_id (@dbxref_id_list) {
		print STDOUT "\t\tdbxref_id: $dbxref_id\n";
	    }
	    
	    print STDOUT "\n";
	}
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

    ## TRUE RUN
    ## First, creation of a metadata object
    
    my $metadbdata = CXGN::Metadata::Metadbdata->new($schema, $loader_username);
    $metadbdata->store();       
    
    ## Second, store the technology objects and the platform object. After store technology_type and before store 
    ## platform, it will use the %techtype_platform_relations hash to set the technology_type related with each platform
    
    foreach my $techtypename (sort keys %techtypes) {
	$techtypes{$techtypename}->store($metadbdata);
	    
	    my @platf_list = @{ $techtype_platform_relations{$techtypename} };
	    foreach my $platf (@platf_list) {
		$platforms{$platf}->set_technology_type_id( $techtypes{$techtypename}->get_technology_type_id );
	    }
	}

    foreach my $platformname (sort keys %platforms) {
	    $platforms{$platformname}->store($metadbdata);
    }

    print STDERR "\nStep 5: Print the data loaded log.\n";
    print STDOUT "\nDATA PROCESSED TO BE STORED: \n";

    foreach my $techtypename (sort keys %techtypes) {
	my $t_techtype_id = $techtypes{$techtypename}->get_technology_type_id() || 'undef';
	my $t_techtype_name = $techtypes{$techtypename}->get_technology_name() || 'undef';
	my $t_techtype_description = $techtypes{$techtypename}->get_description() || 'undef';

	my $t_tt_metadata_id = $techtypes{$techtypename}->get_technology_type_metadbdata()->get_metadata_id() || 'undef';
	
	print STDOUT "+ TECHNOLOGY_TYPE_DATA:\n\ttechnology_type_id =\t$t_techtype_id\n";
	print STDOUT "\ttechnology_name =\t$t_techtype_name\n";
	print STDOUT "\tdescription =\t$t_techtype_description\n";
	print STDOUT "\tmetadata_id =\t$t_tt_metadata_id\n\n";
    }

    foreach my $platformname (sort keys %platforms) {
	my $t_platform_id = $platforms{$platformname}->get_platform_id() || 'undef';
	my $t_platform_name = $platforms{$platformname}->get_platform_name() || 'undef';
	my $t_technology_type_id = $platforms{$platformname}->get_technology_type_id() || 'undef';
	my $t_description =  $platforms{$platformname}->get_description() || 'undef';
	my $t_contact_id = $platforms{$platformname}->get_contact_id() || 'undef';
	
	my $t_p_metadata_id = $platforms{$platformname}->get_platform_metadbdata()->get_metadata_id() || 'undef';

	print STDOUT "+ PLATFORM_DATA:\n\tplatform_id =\t$t_platform_id\n";
	print STDOUT "\tplatform_name =\t$t_platform_name\n";
	print STDOUT "\ttechnology_type_id =\t$t_technology_type_id\n";
	print STDOUT "\tdescription =\t$t_description\n";
	print STDOUT "\tcontact_id =\t$t_contact_id\n";
	print STDOUT "\tmetadata_id =\t$t_p_metadata_id\n\n";
	
	print STDOUT "  * Associated publications:\n";
	my @pub_title_list = $platforms{$platformname}->get_publication_list('title');
	if (scalar(@pub_title_list) == 0) {
	    push @pub_title_list, 'none';
	}
	foreach my $pub (@pub_title_list) {
	    print STDOUT "\t\tpub_title: $pub\n";
	}
	
	print STDOUT "\n";
	print STDOUT "  * Associated platform_design:\n";
	my @design_sample_names = $platforms{$platformname}->get_design_list('sample_name');
	if (scalar(@design_sample_names) == 0) {
	    push @design_sample_names, 'none';
	}
	foreach my $samplename (@design_sample_names) {
	    print STDOUT "\t\tdesign (sample_name): $samplename\n";
	}

	print STDOUT "\n";
	print STDOUT "  * Associated dbxref_ids:\n";
	my @dbxref_id_list = $platforms{$platformname}->get_dbxref_list();
	if (scalar(@dbxref_id_list) == 0) {
	    push @dbxref_id_list, 'none';
	}
	foreach my $dbxref_id (@dbxref_id_list) {
	    print STDOUT "\t\tdbxref_id: $dbxref_id\n";
	}

	print STDOUT "\n";
    }
    
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
      This script load data for samples

      In a .bs format (text format with fields delimited by [] ). To print a template file use the -X option and fill following the
      instructions. The load data file can have one or more of these data types.

    
    Usage: 
       platform_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -p <platform_file> [-T]

      To collect the data loaded report into a file:

       platform_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -p <platform_file> [-T] > file.log

    Example: 
      perl platform_dbload.pl -u aure -H localhost -D sandbox -p TOM2.bs


    Flags:
      -u loader username      loader username (mandatory)
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory for check domains in the database)
      -D database name        sandbox or cxgn etc (mandatory for check domains in the database)
      -p platform file        data load file input file (mandatory)
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



=head2 print_platform_template

  Usage: print_platform_template()
  Desc: print a platform_template file
  Ret: none
  Args: none
  Side_Effects: create a file and exit of the script
  Example:

=cut

sub print_platform_template {
    my $dir = `pwd`;
    chomp($dir);
    my $template_platform_file = $dir . '/platform_load_file.bs';
    open my $TFH, '>', $template_platform_file || die("Sorry, I can not open the file: $template_platform_file.\n");

    print STDERR "PRINT GEM(Platform) DATA LOAD FILE option...\n";
    my $info = '# Notes for the Source data load to GEM Platform data ()
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
# NOTE FOR DATABASE CURATORS: To load all this data into the database, use the platform_dbload.pl script
#


#####################################################################################
# PLATFORM_DATA #####################################################################
# FIELD	#		# DATA #	# DATA_STRINGENCY #	# EXAMPLE_AND_NOTES #
#####################################################################################

## To add a new Technology Type:

*DATA_TYPE:		[]		#mandatory,single	example:[technology_type]
*TECHNOLOGY_NAME:	[]		#mandatory,single,fl	example:[oligo spotted microarray]
*TECHNOLOGY_DESCRIP:    []		#optional,single	example:[The oligo spotted microarrays...]
//

## To add a new Platform (MANDATORY)

*DATA_TYPE:             []              #mandatory,single       example:[platform]
*TECHNOLOGY_NAME:       []              #mandatory,single       example:[oligo spotted microarray]
*PLATFORM_NAME:		[]		#mandatory,single	example:[TOM2]
*PLATFORM_DESCR:        []              #optional,single	example:[TOM2 (Tomato Microarray) is...]
*CONTACT_NAME:          []              #optional,single,note1  example:[Aureliano Bombarely]
//


## To associate a publication to an existing platform:
## Note: title have to exists into pub table before do the association

*DATA_TYPE:             []              #mandatory              example:[pub]
*PLATFORM_NAME:		[]		#mandatory,single	example:[TOM2]
*TITLE:                 []              #mandatory              example:[The example of the pub]
//

## To add a sample link using platform_design to an existing platform:
## Note: sample_name have to exists into biosource.bs_sample table before do the association

*DATA_TYPE:             []              #mandatory,single       example:[platform_design]
*PLATFORM_NAME:		[]		#mandatory,single	example:[TOM2]
*SAMPLE_NAME:           []              #mandatory,single       example:[Tomato Unigene Build#1]
//  

## To associate external references to the platform
## Note: dbname and accessions have to exists into db and dbxref tables before do the association

*DATA_TYPE:	        []              #mandatory,single	example:[dbxref]
*PLATFORM_NAME:         []              #mandatory,single       example:[TOM2]
*DBNAME:                []              #mandatory,single       example:[GEO-Platform]
*ACCESSIONS:		[]		#mandatory,multiple     example:[GPL7824]
//';

    print $TFH "$info";
    print STDERR "...done (printed platform_file with the name: $template_platform_file)\n\n";
    exit (1);
}




####
1;##
####
