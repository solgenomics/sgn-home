#!/usr/bin/perl

=head1 NAME

 template_dbload.pl
 A script to parse template and probe data files and load into the database (version.0.1)

=cut

=head1 SYPNOSIS

 template_dbload.pl [-h] -U <load_username> -H <dbhost> -D <dbname> -i <input_file> [-T] [-X]

  To collect the data loaded report into a file:

 template_dbload [-h] [-X] -D <dbname> -H <dbhost> -i <input_file> [-T] > file.log


=head1 EXAMPLE:

 perl template_dbload.pl -u aure -H localhost -D sandbox -p tom2_templates.bs
        
    
=head2 I<Flags:>

=over

=item -i

B<input load file>              data load file in bs format (mandatory).

=item -u

B<username>                     username for the loader (mandatory)

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

  This script parse a data file, create some temp files and finally load the data into the database using COPY command.

  1) Parse bs file:

     This script parse a bs file. These files have the following structure:

     *DATA_TYPE: []
     {some fields annotated as *FIELD_NAME: []}
     //

  2) Create temp files with the data parsed 

  3) Load the data into the database and search the correspondences

   3.0) Create a metadata row

   3.1) Check template and load template data into the database.

        + gem.ge_template

   3.2) Search for dbiref and create a temp file with the new ones. Load into the database. (optional)

        + metadata.md_dbiref

   3.3) Create a template_dbiref file with template_ids and dbiref_ids. Load into the database. (optional)
   
        + gem.ge_templates_dbiref
    
   3.4) Search for dbxref and create a temp file with the new ones. Load into the database. (optional)

        + public.dbxref

   3.5) Create a template_dbxref file with template_ids and dbxref_ids. Load into the database. (optional)

        + gem.ge_templates_dbxref

   3.6) Load the probe sequence filename into the file table (optional)

        + metadata.md_file

   3.7) Create a probe file with the probe data parse, the new template ids and the files ids. (optional)
        Load into the database.

        + gem.ge_probe

   3.8) Create a spot file and load it into the database. (optional)

        + gem.ge_probe_spot

   3.9) Create a spot coordinate file and load into the database. (optional)

        + gem.ge_probe_spot_coordinates

  4) Print a report

  5) Rollback and set the sequences in the primary keys or commit the changes. 
    


=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 template_dbload.pl


=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use Cwd;

use CXGN::GEM::Schema;
use CXGN::GEM::TechnologyType;
use CXGN::GEM::Platform;
use CXGN::Biosource::Sample;
use CXGN::DB::InsertDBH;
use CXGN::Metadata::Metadbdata;

our ($opt_u, $opt_H, $opt_D, $opt_i, $opt_T, $opt_X, $opt_h);
getopts("u:H:D:i:TXh");
if (!$opt_u && !$opt_H && !$opt_D && !$opt_i && !$opt_T && !$opt_X && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
} elsif ($opt_h) {
    help();
} elsif ($opt_X) {
    print_template();
}

## Checking the input arguments

my $loader_username = $opt_u || die("MANDATORY ARGUMENT ERROR: The -u <loader_username> argument was not supplied.\n");
my $dbname = $opt_D || die("MANDATORY ARGUMENT ERROR: The -D <database_name> argument was not supplied.\n");
my $dbhost = $opt_H || die("MANDATORY ARGUMENT ERROR: The -H <db_hostname> argument was not supplied.\n"); 
my $argument_file = $opt_i || die("MANDATORY ARGUMENT ERROR: The -i <argument_file> argument was not supplied.\n");

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


## Parse the general information file

print STDERR "\nStep 3: Open and parse the argument file.\n";

my $arg_href = parse_argument_file($argument_file);

## Now it will get the platform_id for the platform_name.
## It will die if platform_name do not exists into the database


print STDERR "\nStep 4: Get platform data from DB.\n";

my $platform_name = $arg_href->{'platform_name'};

my ($platform_row) = $schema->resultset('GePlatform')
                            ->search( { platform_name => $platform_name} );

unless (defined $platform_row) {
    die("DATA ERROR: The platform_name=$platform_name does not exist into the database.\n");
}

$arg_href->{'platform_id'} = $platform_row->get_column('platform_id');


## For test option it will use eval

if ($opt_T) {
    print STDERR "\nRUNNING THE SCRIPT IN TEST MODE\n";
    eval {

	## The different functions will be group using the dataload_pipeline

	dataload_pipeline($schema, $loader_username, $all_last_ids_href, $arg_href);

    };

    ## Print errors if something was wrong during the test

    if ($@) {
        print STDERR "\nTEST ERRORS:\n\n$@\n";
    }

    ## Finally, rollback, because it is a test and set the sequence values

    $schema->txn_rollback;
    $schema->set_sqlseq_values_to_original_state($all_last_ids_href);
}
else {
    print STDERR "\nRUNNING THE SCRIPT IN LIVE MODE\n";

    dataload_pipeline($schema, $loader_username, $all_last_ids_href, $arg_href);

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
      This script load data for templates

      In a .bs format (text format with fields delimited by [] ). To print a template file use the -X option and fill following the
      instructions. The load data file can have one or more of these data types.

    
    Usage: 
       templates_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -i <templates_file> [-T]

      To collect the data loaded report into a file:

       templates_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -i <templates_file> [-T] > file.log

    Example: 
      perl templates_dbload.pl -u aure -H localhost -D sandbox -i TOM2.bs


    Flags:
      -u loader username      loader username (mandatory)
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory for check domains in the database)
      -D database name        sandbox or cxgn etc (mandatory for check domains in the database)
      -i argument file        data load file input file (mandatory)
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


########################
### DATA LOAD METHOD ###
########################

=head2 dataload_pipeline

  Usage: dataload_pipeline($schema, $loader_username, $all_last_ids_href, $args_href);
  Desc: run the different methods to load the data 
  Ret: none
  Args: $schema, a CXGN::GEM::Schema object
        $loader_username, username from sgn_people.sp_person that
        load the data
        $args_href, a hash reference with keys=script_argument and 
        value=value
  Side_Effects: Die with something wrong
                Print status messages
  Example: dataload_pipeline($schema, $loader_username, $all_last_ids_href, $args_href);

=cut

sub dataload_pipeline {
    my $schema = shift 
	|| die("Schema argument was not supplied to dataload_pipeline function in template_dbload.pl script.\n");
    my $loader_username = shift 
	|| die("Argument loader_username was not supplied to dataload_pipeline function in template_dbload.pl script.\n");
    my $all_last_ids_href = shift 
	|| die("Argument all_last_ids_href was not supplied to dataload_pipeline function in template_dbload.pl script.\n");
    my $args_href = shift 
	|| die("Argument hash reference was not supplied to dataload_pipeline function in template_dbload.pl script.\n");

    print STDERR "\nStep 5: Create a new metadata_id.\n";
    
    ## First, create a metadata and store it to get the metadata_id and add to arguments

    my $metadbdata = CXGN::Metadata::Metadbdata->new($schema, $loader_username);
    $metadbdata->store();    

    $arg_href->{'metadata_id'} = $metadbdata->get_metadata_id();
	
    ## Second parse TEMPLATE file

    print STDERR "\nStep 6: Parse the template_file ($arg_href->{template_file}) and create the temp_template_dbload file.\n";

    my $template_dbload_file = create_template_dbload_file($arg_href);

    print STDERR "\nStep 7: Load the $template_dbload_file file into the gem.ge_template table.\n";
    
    ## Use of pg_putline to use the copy command as stdin (it will let copy from a machine different from the db host)

    $schema->storage()
	   ->dbh()
	   ->do("COPY gem.ge_template (template_name, template_type, platform_id, metadata_id) FROM STDIN");

    open my $temp_io, '<', $template_dbload_file;
    while (<$temp_io>) {
	$schema->storage()->dbh()->pg_putline($_);
    }
    close $temp_io;

    $schema->storage()->dbh()->pg_endcopy();

    my $curr_last_template = $schema->resultset('GeTemplate')
	                            ->search( undef, { order_by => 'template_id DESC' })
				    ->first()
				    ->get_column('template_id');

    my $template_count = $curr_last_template - $all_last_ids_href->{'gem.ge_template_template_id_seq'};

    print STDERR "\t\t$template_count template rows has been loaded into the database.\n";

    print STDERR "\nStep 8: Get the template_id for the data loaded for platform_id=$arg_href->{'platform_id'}\n";

    my $c = 0;
    my %templates_equiv;
    my @template_rows = $schema->resultset('GeTemplate')
	                       ->search( { platform_id => $arg_href->{'platform_id'} });
    
    foreach my $template_row (@template_rows) {
	$c++;
	my $temp_name = $template_row->get_column('template_name');
	my $temp_id = $template_row->get_column('template_id');
	
	$templates_equiv{$temp_name} = $temp_id;
	print STDERR "\t\tGetting template_ids ($c) ( $temp_name => $temp_id )                                \r";
    }
    print STDERR "\n\t\t\t...Done\n";

    $arg_href->{'templates_equiv'} = \%templates_equiv;

    ## Third parse iref_mapping if it is supplied

    my $dbiref_dbload_file;
 
    print STDERR "\nStep 9: Parse the iref_mapping file and create the temp_dbiref_dbload file.\n";

    if (defined $arg_href->{'iref_mapping_file'}) {
	
	$dbiref_dbload_file = create_dbiref_dbload_file($arg_href);
    }
    else {
	print STDERR "\t\tNONE IREF_MAPPING_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    print STDERR "\nStep 10: Load the dbiref_dbload_file file into the gem.ge_template_dbiref table.\n";

    if (defined $dbiref_dbload_file) {

	$schema->storage()
	       ->dbh()
	       ->do("COPY gem.ge_template_dbiref (template_id, dbiref_id, metadata_id) FROM STDIN");

	open my $temp_dbiref_io, '<', $dbiref_dbload_file;
	while (<$temp_dbiref_io>) {
	    $schema->storage()->dbh()->pg_putline($_);
	}
	close $temp_dbiref_io;

	$schema->storage()->dbh()->pg_endcopy();

	my $curr_last_dbiref = $schema->resultset('GeTemplateDbiref')
	                              ->search( undef, { order_by => 'template_dbiref_id DESC' })
				      ->first()
				      ->get_column('template_dbiref_id');

	my $dbiref_count = $curr_last_dbiref - $all_last_ids_href->{'gem.ge_template_dbiref_template_dbiref_id_seq'};
	    
	print STDERR "\t\t$dbiref_count template rows has been loaded into the database.\n";	
    }
    else {
	print STDERR "\t\tNONE IREF_MAPPING_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    ## Forth parse xref_mapping if it is supplied

    my $dbxref_dbload_file;
 
    print STDERR "\nStep 11: Parse the xref_mapping file and create the temp_dbxref_dbload file.\n";

    if (defined $arg_href->{'xref_mapping_file'}) {
	$dbxref_dbload_file = create_dbxref_dbload_file($arg_href);
    }
    else {
	print STDERR "\t\tNONE XREF_MAPPING_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    print STDERR "\nStep 12: Load the dbxref_dbload_file file into the gem.ge_template_dbxref table.\n";

    if (defined $dbxref_dbload_file) {
	
	$schema->storage()
	       ->dbh()
	       ->do("COPY gem.ge_template_dbxref (template_id, dbxref_id, metadata_id) FROM STDIN");

	open my $temp_dbxref_io, '<', $dbxref_dbload_file;
	while (<$temp_dbxref_io>) {
	    $schema->storage()->dbh()->pg_putline($_);
	}
	close $temp_dbxref_io;

	$schema->storage()->dbh()->pg_endcopy();

	my $curr_last_dbxref = $schema->resultset('GeTemplateDbxref')
	                              ->search( undef, { order_by => 'template_dbxref_id DESC' })
				      ->first()
				      ->get_column('template_dbxref_id');

	my $dbxref_count = $curr_last_dbxref - $all_last_ids_href->{'gem.ge_template_dbxref_template_dbxref_id_seq'};
	
	print STDERR "\t\t$dbxref_count template rows has been loaded into the database.\n";	
    }
    else {
	print STDERR "\t\tNONE XREF_MAPPING_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    ## Sixth load the file data into the database (it isnot necessary to parse a file because all the data needed should be 
    ## in the arg_href variable).

    print STDERR "\nStep 13: Associate file to probe. Loading data step.\n";
    my $file_id;

    if (defined $arg_href->{'basename'}) {
	my ($file_row) = $schema->resultset('MdFiles')
	                        ->search( 
	                                  { 
					    basename => $arg_href->{'basename'},
					    dirname  => $arg_href->{'dirname'},
					  }
				        );

	if (defined $file_row) {
	    $file_id = $file_row->get_column('file_id');
	    print STDERR "\t\tFile with file_id=$file_id was found in db for basename=$arg_href->{'basename'}. Using this file_id.\n";
	}
	else {

	    ## Collect the file data into a hash ref.
	    my @file_fields = ('basename', 'dirname', 'filetype', 'alt_filename', 'comment', 
			       'md5checksum', 'urlsource', 'urlsource_md5checksum', 'metadata_id');

	    my %file_args;
	    foreach my $file_field (@file_fields) {
		$file_args{$file_field} = $arg_href->{$file_field};
	    }
	    
	    $file_row = $schema->resultset('MdFiles')
		               ->new(\%file_args)
			       ->insert()
			       ->discard_changes();
	    $file_id = $file_row->get_column('file_id');
	    print STDERR "\t\tInserted a new file row with basename=$arg_href->{'basename'} (file_id=$file_id).\n";
	}
	$arg_href->{'file_id'} = $file_id;
    }
    else {
	print STDERR "\t\tNONE BASENAME was supplied to this script. SKIPPING STEP.\n";
    }
		
    ## Seventh, parse and load the probe data (if exists).
    
    my $probe_dbload_file;
 
    print STDERR "\nStep 14: Parse the probe_file and create the temp_probe_dbload file.\n";

    if (defined $arg_href->{'probe_file'}) {
	    
	$probe_dbload_file = create_probe_dbload_file($arg_href);
    }
    else {
	print STDERR "\t\tNONE PROBE_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    print STDERR "\nStep 15: Load the probe_dbload_file file into the gem.ge_probe table.\n";

    if (defined $probe_dbload_file) {

	$schema->storage()
	       ->dbh()
	       ->do("COPY gem.ge_probe (platform_id, probe_name, probe_type, sequence_file_id, template_id, template_start,
                                        template_end, metadata_id ) FROM STDIN");

	open my $temp_probe_io, '<', $probe_dbload_file;
	while (<$temp_probe_io>) {
	    $schema->storage()->dbh()->pg_putline($_);
	}
	close $temp_probe_io;

	$schema->storage()->dbh()->pg_endcopy();

	my $curr_last_probe = $schema->resultset('GeProbe')
	                             ->search( undef, { order_by => 'probe_id DESC' })
				     ->first()
				     ->get_column('probe_id');

	my $probe_count = $curr_last_probe - $all_last_ids_href->{'gem.ge_probe_probe_id_seq'};
	    
	print STDERR "\t\t$probe_count probe rows has been loaded into the database.\n";	
    }
    else {
	print STDERR "\t\tNONE PROBE_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    print STDERR "\nStep 16: Get the probe_id for the data loaded for platform_id=$arg_href->{'platform_id'}\n";

    if (defined $probe_dbload_file) {

	my $d = 0;
	my %probes_equiv;
	my @probes_rows = $schema->resultset('GeProbe')
	                         ->search( { platform_id => $arg_href->{'platform_id'} });
	    
	foreach my $probe_row (@probes_rows) {
	    $d++;
	    my $prob_name = $probe_row->get_column('probe_name');
	    my $prob_id = $probe_row->get_column('probe_id');
	    
	    $probes_equiv{$prob_name} = $prob_id;
	    print STDERR "\t\tGetting probes_ids ($d) ( $prob_name => $prob_id )                                \r";
	}
	print STDERR "\n\t\t\t...Done\n";
	
	$arg_href->{'probes_equiv'} = \%probes_equiv;
    }
    else {
	print STDERR "\t\tNONE PROBE_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    ## Eigth, parse and load the spot data (if exists), in two steps, load the spots names and the coordinates
	
    my $spot_dbload_file;
 
    print STDERR "\nStep 17: Parse the spot_file and create the temp_spot_dbload file.\n";

    if (defined $arg_href->{'spot_file'}) {
	    
	$spot_dbload_file = create_spot_dbload_file($arg_href);
    }
    else {
	print STDERR "\t\tNONE SPOT_FILE was supplied to this script. SKIPPING STEP.\n";
    }


    print STDERR "\nStep 18: Load the spot_dbload_file file into the gem.ge_probe_spot table.\n";

    if (defined $spot_dbload_file) {

	$schema->storage()
	       ->dbh()
	       ->do("COPY gem.ge_probe_spot (probe_id, spot_name, metadata_id ) FROM STDIN");

	open my $temp_spot_io, '<', $spot_dbload_file;
	while (<$temp_spot_io>) {
	    $schema->storage()->dbh()->pg_putline($_);
	}
	close $temp_spot_io;

	$schema->storage()->dbh()->pg_endcopy();

	my $curr_last_spot = $schema->resultset('GeProbeSpot')
	                            ->search( undef, { order_by => 'probe_spot_id DESC' })
				    ->first()
				    ->get_column('probe_spot_id');

	my $spot_count = $curr_last_spot - $all_last_ids_href->{'gem.ge_probe_spot_probe_spot_id_seq'};
	    
	print STDERR "\t\t$spot_count probe spot rows has been loaded into the database.\n";	
    }
    else {
	print STDERR "\t\tNONE SPOT_FILE was supplied to this script. SKIPPING STEP.\n";
    }
    
    print STDERR "\nStep 19: Get the probe_spot_id for the data loaded for platform_id=$arg_href->{'platform_id'}\n";

    if (defined $spot_dbload_file) {

	my $e = 0;
	my %spot_equiv;
	my %probes_equiv = %{ $arg_href->{'probes_equiv'} };
	foreach my $probe_name (keys %probes_equiv) {
	    my $prob_id = $probes_equiv{$probe_name};
	    my @spot_rows = $schema->resultset('GeProbeSpot')
	                           ->search( { probe_id => $prob_id } );
	    
	    foreach my $spot_row (@spot_rows) {
		$e++;
		my $spot_name = $spot_row->get_column('spot_name');
		my $spot_id = $spot_row->get_column('probe_spot_id');
	    
		$spot_equiv{$spot_name} = $spot_id;
		print STDERR "\t\tGetting probe_spot_ids ($e) ( $spot_name => $spot_id )                                \r";
	    }
	}
	print STDERR "\n\t\t\t...Done\n";
	
	$arg_href->{'spots_equiv'} = \%spot_equiv;
    }
    else {
	print STDERR "\t\tNONE PROBE_FILE was supplied to this script. SKIPPING STEP.\n";
    }

    ## create_spot_coord_dbload_file will use the same file than create_spot_dbload, but it will depend of 
    ## spot_coordinates keys in the arg_href.

	
    my $spot_coord_dbload_file;
 
    print STDERR "\nStep 20: Parse the spot_file and create the temp_spot_coord_dbload file.\n";

    if (defined $arg_href->{'spot_coordinates'}) {
	
	$spot_coord_dbload_file = create_spot_coord_dbload_file($arg_href);
    }
    else {
	print STDERR "\t\tNONE SPOT_COORDINATES was supplied to this script. SKIPPING STEP.\n";
    }


    print STDERR "\nStep 21: Load the spot_coord_dbload_file file into the gem.ge_probe_spot_coordinate table.\n";

    if (defined $spot_coord_dbload_file) {

	$schema->storage()
	       ->dbh()
	       ->do("COPY gem.ge_probe_spot_coordinate (probe_spot_id, coordinate_type, coordinate_value, metadata_id ) 
                     FROM STDIN");

	open my $temp_spot_coord_io, '<', $spot_coord_dbload_file;
	while (<$temp_spot_coord_io>) {
	    $schema->storage()->dbh()->pg_putline($_);
	}
	close $temp_spot_coord_io;

	$schema->storage()->dbh()->pg_endcopy();

	my $curr_last_spotcoord = $schema->resultset('GeProbeSpotCoordinate')
	                                 ->search( undef, { order_by => 'probe_spot_coordinate_id DESC' })
					 ->first()
					 ->get_column('probe_spot_coordinate_id');

	my $spotcoord_ct = $curr_last_spotcoord - $all_last_ids_href->{'gem.ge_probe_spot_coordinate_probe_spot_coordinate_id_seq'};
	
	print STDERR "\t\t$spotcoord_ct probe spot coordinate rows has been loaded into the database.\n";	
    }
    else {
	print STDERR "\t\tNONE SPOT_COORDINATES was supplied to this script. SKIPPING STEP.\n";
    }
}



#######################
### PARSING METHODS ###
#######################

=head2 parse_argument_file

  Usage: my $arg_href = parse_argument_file($file);
  Desc: parse the argument file and return a hash
        with the arguments
  Ret: A hash ref. with keys=argument_name and value=value
  Args: $file, a file name to open
  Side_Effects: Check if the files exists, if not die
  Example: my $arg_href = parse_argument_file($file);

=cut

sub parse_argument_file {
    my $file = shift || die("Schema argument was not supplied to parse_argument_file function in template_dbload.pl script.\n");

    ## Define the argument hash

    my $arg_href = {};

    ## Define the keys and the types.

    my @fields = ( 'PLATFORM_NAME', 'TEMPLATE_FILE', 'TEMPLATE_TYPE', 'IREF_MAPPING_FILE', 'IPATH',
		   'XREF_MAPPING_FILE', 'XDB_NAME', 'BASENAME', 'DIRNAME', 'FILETYPE', 'ALT_FILENAME',
		   'COMMENT', 'MD5CHECKSUM', 'URLSOURCE', 'URLSOURCE_MD5CHECKSUM', 'PROBE_FILE', 'SPOT_FILE', 'SPOT_COORDINATES' );

    ## Open the file

    open my $ifh, '<', $file || die("Sorry, I can not open the file: $file.\n");

    ## Define the variables

    my $l = 0;

    while (<$ifh>) {
	 $l++; ## Line counter

	 ## It will do not read any line that start with #
	 unless ($_ =~ m/^#/) {
	     foreach my $field (@fields) {
		 $arg_href = get_variable($field, $arg_href, $l, $_);
	     }
	 }
    }
    close $ifh;

    ## Now check the dependencies

    unless (defined $arg_href->{'platform_name'}) {
	die("MANDATORY ARGUMENT ERROR: Platform_name was not supplied from argument_file.\n");
    }
    unless (defined $arg_href->{'template_file'} && defined $arg_href->{'template_type'}) {
	die("MANDATORY ARGUMENT ERROR: Template_file and template_type were not supplied from argument_file.\n");
    }
    if (defined $arg_href->{'iref_mapping_file'} || defined $arg_href->{'ipath'}) {
	unless (defined $arg_href->{'iref_mapping_file'} && defined $arg_href->{'ipath'}) {
	    die("MANDATORY ARGUMENT ERROR: iref_mapping_file or ipath argument were not supplied to load dbiref.\n");
	}
    }
    if (defined $arg_href->{'xref_mapping_file'} || defined $arg_href->{'xdb_name'}) {
	unless (defined $arg_href->{'xref_mapping_file'} && defined $arg_href->{'xdb_name'}) {
	    die("MANDATORY ARGUMENT ERROR: xref_mapping_file or xdb_name argument were not supplied to load dbxref.\n");
	}
    }
 
    return $arg_href;
}

=head2 get_variable

  Usage: my $arg_href = get_variable($variable, $arg_href, $l, $line)
  Desc: check if exists a variable
  Ret: The variable parsed
  Args: $variable to get,
        $arg_href, the hash reference of the argument array
        $l, line counter
        $line, the line to parse
  Side_Effects: Die if the variable do not exists
  Example: my $arg_href = get_variable($variable, $arg_href, $l, $line)

=cut

sub get_variable {
    my $variable = shift;
    my $arg_href = shift;
    my $l = shift;
    my $line = shift;

    my $form_variable = lc($variable);

    if ($variable =~ m/_FILE/) {

	if ($_ =~ m/\*$variable:\s+\[(.+)\]/) {
	    my $data1 = $1 ||
		die("MANDATORY DATA ERROR (line $l): None $variable data was supplied.\n");
	    	    
	    unless (-s $data1) {
		die("FILE ERROR: The $variable=$data1 does not exist or has zero size.\n");
	    }
	    else {
		unless (defined $arg_href->{$form_variable}) {
		    $arg_href->{$form_variable} = $data1;
		}
		else {
		    die("DUPLICATION ERROR (line $l): There are more than one $form_variable arg. It should be unique per load.\n ");
		}
	    }	
	}
    }
    else {
	if ($_ =~ m/\*$variable:\s+\[(.+)\]/) {
	    my $data2 = $1 ||
		die("MANDATORY DATA ERROR (line $l): None $variable data was supplied.\n");	    
	    
	    unless ( defined $arg_href->{$form_variable} ) {
		$arg_href->{$form_variable} = $data2;
	    }
	    else {
		die("DUPLICATION ERROR (line $l): There are more than one $form_variable arg. It should be unique per load.\n ");
	    }
	}
    }
    return $arg_href;
}

=head2 create_template_dbload_file

  Usage: my $template_dbload_file = create_template_dbload_file(arg_href);
  Desc: Parse the template file and create a template file with the data
  Ret: $template_dbload_file, a file name
  Args: $arg_href, a hash reference with keys=arg and value=values
  Side_Effects: Die if something is wrong
  Example: my $template_dbload_file = create_template_dbload_file($arg_href);

=cut

sub create_template_dbload_file {
    my $arg_href = shift 
	|| die("DATA ERROR: None argument hash reference was supplied to create_template_dbload_file.\n");

    ## Get the variables

    my $template_type = $arg_href->{'template_type'} 
       || die("DATA ERROR: None template_type was supplied to create_template_dbload_file.\n");
    my $platform_id = $arg_href->{'platform_id'}
       || die("DATA ERROR: None platform_id was supplied to create_template_dbload_file.\n");
    my $metadata_id = $arg_href->{'metadata_id'} 
       || die("DATA ERROR: None metadata_id was supplied to create_template_dbload_file.\n");

    ## template redundancy count

    my $template_c = 0;

    ## It will use as path `pwd`

    my $path = getcwd();

    open my $in_tfh, '<', $arg_href->{'template_file'} 
        || die("Sorry but I can not open the file:" . $arg_href->{'template_file'} . " ($!)\n");

    my $temp_template = $path ."/temp_template_dbload.tab";

    open my $out_tfh, '>', $temp_template 
	|| die("Sorry but I can not open the file:$temp_template ($!)\n");

    ## The template file should have only a column

    my $l = 0;
    while (<$in_tfh>) {
	$l++;
	print STDERR "\t\tParsing line:$l from file $arg_href->{'template_file'}\r";
	chomp($_);
	my @data = split(/\t/, $_);
	my $template_name = $data[0];
	
	## Filter spaces in the begining and in the end

	$template_name =~ s/^\s+//g;
	$template_name =~ s/\s+$//g;

	## Check if exists in the database as: template_name + platform_id

	my ($template_row) = $schema->resultset('GeTemplate')
	                            ->search( { template_name => $template_name, platform_id => $platform_id });
	
	unless (defined $template_row) {
	    print $out_tfh "$template_name\t$template_type\t$platform_id\t$metadata_id\n";
	}
	else {
	    $template_c++;
	}
    }
    print STDERR "\n\t\t\tDone... file parsed.\n\n";
    if ($template_c > 0) {
	print STDERR "\t\tThere are $template_c templates that exist into the database. These templates have been skipped to create the dbload file.\n\n";
    }
    close($in_tfh);
    
    return $temp_template;
}

=head2 create_dbiref_dbload_file

  Usage: my $dbiref_dbload_file = create_dbiref_dbload_file(arg_href);
  Desc: Parse the dbiref file and create a dbiref dbload file with the data
  Ret: $dbiref_dbload_file, a file name
  Args: $arg_href, a hash reference with keys=arg and value=values
  Side_Effects: Store a new dbipath and dbiref if don't exist
                Die if something is wrong
  Example: my $dbiref_dbload_file = create_dbiref_dbload_file($arg_href);

=cut

sub create_dbiref_dbload_file {
    my $arg_href = shift 
	|| die("DATA ERROR: None argument hash reference was supplied to create_dbiref_dbload_file.\n");

    ## Get the variables

    my $ipath = $arg_href->{'ipath'} 
       || die("DATA ERROR: None ipath was supplied to create_dbiref_dbload_file.\n");
    my $metadata_id = $arg_href->{'metadata_id'} 
       || die("DATA ERROR: None metadata_id was supplied to create_dbiref_dbload_file.\n");
    my $templates_equiv_href = $arg_href->{'templates_equiv'}
       || die("DATA ERROR: None templates_equiv was supplied to create_dbiref_dbload_file.\n");

    ## Redundancy count

    my $rel_dbiref_c = 0;

    my @ipath = split(/\./, $ipath);

    ## Get the dbipath_id for $ipath, if do not exists it will add a new one

    my ($dbipath_row) = $schema->resultset('MdDbipath')
	                       ->search(
	                                 { 
				           schema_name => $ipath[0], 
					   table_name  => $ipath[1],
					   column_name => $ipath[2],
				         }
			               );
    my $dbipath_id;
    
    if (defined $dbipath_row) {
	$dbipath_id = $dbipath_row->get_column('dbipath_id');
    }
    else {
	$dbipath_row = $schema->resultset('MdDbipath')
	                      ->new(
			             { 
				       schema_name => $ipath[0], 
				       table_name  => $ipath[1],
				       column_name => $ipath[2],
				       metadata_id => $arg_href->{'metadata_id'}
				     }
			           );
	$dbipath_row->insert()
	            ->discard_changes();
	$dbipath_id = $dbipath_row->get_column('dbipath_id');
    }

    ## It will use as path `pwd`

    my $path = getcwd();

    open my $in_ifh, '<', $arg_href->{'iref_mapping_file'} 
        || die("Sorry but I can not open the file:" . $arg_href->{'iref_mapping_file'} . " ($!)\n");

    my $temp_dbiref = $path ."/temp_dbiref_dbload.tab";

    open my $out_ifh, '>', $temp_dbiref 
	|| die("Sorry but I can not open the file:$temp_dbiref ($!)\n");

    ## The iref_mapping_file should have at least two columns (-f1 template_name and -f2 iaccession)

    my $l = 0;
    while (<$in_ifh>) {
	$l++;
	print STDERR "\t\tParsing line:$l from file $arg_href->{'iref_mapping_file'}\r";
	chomp($_);
	my @data = split(/\t/, $_);
	my $template_name = $data[0];
	my $iaccession = $data[1];
	
	## Check if exists dbiref in the database as: accession + dbipath_id

	my ($dbiref_row) = $schema->resultset('MdDbiref')
	                          ->search( { iref_accession => $iaccession, dbipath_id => $dbipath_id });
	
	my $dbiref_id;
	if (defined $dbiref_row) {
	    $dbiref_id = $dbiref_row->get_column('dbiref_id');
	}
	else {
	    $dbiref_row = $schema->resultset('MdDbiref')
	                         ->new(
			                { 
 				          iref_accession => $iaccession, 
				          dbipath_id     => $dbipath_id,
 				          metadata_id    => $arg_href->{'metadata_id'}
				        }
			              );
	    $dbiref_row->insert()
	               ->discard_changes();
	    $dbiref_id = $dbiref_row->get_column('dbiref_id');
	}

	## Now it will search the template_id based in the template name and it will print them into the output file

	my $template_id = $templates_equiv_href->{$template_name};
	
	unless (defined $template_id) {
	    print STDERR "\nSKIPPING template_name=$template_name, it have not template_id.\n";
	}
	else {

	    ## Finally it will check if exists the relation for a concrete template_id and dbiref_id

	    my ($ge_template_dbiref_row) = $schema->resultset('GeTemplateDbiref')
		                                  ->search({ template_id => $template_id, dbiref_id => $dbiref_id });

	    unless (defined $ge_template_dbiref_row) {
		print $out_ifh "$template_id\t$dbiref_id\t$arg_href->{'metadata_id'}\n";
	    }
	    else {
		$rel_dbiref_c++;
	    }
	}
	    
    }
    print STDERR "\n\t\t\tDone... file parsed.\n\n";
    if ($rel_dbiref_c > 0) {
	print STDERR "\t\tThere are $rel_dbiref_c template-dbiref relations that exist in the db, these template_id-dbiref_id have been skipped to create the dbload.\n\n";
    }

    close($in_ifh);
    
    return $temp_dbiref;
}


=head2 create_dbxref_dbload_file

  Usage: my $dbxref_dbload_file = create_dbxref_dbload_file(arg_href);
  Desc: Parse the dbxref file and create a dbxref dbload file with the data
  Ret: $dbxref_dbload_file, a file name
  Args: $arg_href, a hash reference with keys=arg and value=values
  Side_Effects: Store a new db and dbxref if don't exist
                Die if something is wrong
  Example: my $dbxref_dbload_file = create_dbxref_dbload_file($arg_href);

=cut

sub create_dbxref_dbload_file {
    my $arg_href = shift 
	|| die("DATA ERROR: None argument hash reference was supplied to create_dbxref_dbload_file.\n");

    ## Get the variables

    my $xdb_name = $arg_href->{'xdb_name'} 
       || die("DATA ERROR: None xdb_name was supplied to create_dbxref_dbload_file.\n");
    my $metadata_id = $arg_href->{'metadata_id'} 
       || die("DATA ERROR: None metadata_id was supplied to create_dbxref_dbload_file.\n");
    my $templates_equiv_href = $arg_href->{'templates_equiv'}
       || die("DATA ERROR: None templates_equiv was supplied to create_dbxref_dbload_file.\n");

    ## Count the redundacy

    my $rel_dbxref_c = 0;

    ## Get the db_id for $xdb_name, if do not exists it will die

    my $db_id;

    my ($db_row) = $schema->resultset('General::Db')
	                       ->search(
	                                 { 
				           name => $xdb_name, 
				         }
			               );
    
    
    
    if (defined $db_row) {
	$db_id = $db_row->get_column('db_id');
    }
    else {
	die("DATA ERROR: xdb_name=$xdb_name does not exist in db table.\n");
    }

    ## It will use as path

    my $path = getcwd();

    open my $in_xfh, '<', $arg_href->{'xref_mapping_file'} 
        || die("Sorry but I can not open the file:" . $arg_href->{'xref_mapping_file'} . " ($!)\n");

    my $temp_dbxref = $path ."/temp_dbxref_dbload.tab";

    open my $out_xfh, '>', $temp_dbxref 
	|| die("Sorry but I can not open the file:$temp_dbxref ($!)\n");

    ## The xref_mapping_file should have at least two columns (-f1 template_name and -f2 accession)

    my $l = 0;
    while (<$in_xfh>) {
	$l++;
	print STDERR "\t\tParsing line:$l from file $arg_href->{'xref_mapping_file'}\r";
	chomp($_);
	my @data = split(/\t/, $_);
	my $template_name = $data[0];
	my $accession = $data[1];
	
	## Check if exists dbxref in the database as: accession + db_id

	my ($dbxref_row) = $schema->resultset('General::Dbxref')
	                          ->search( { accession => $accession, db_id => $db_id });
	
	## If do not exists, it will added the new ones

	my $dbxref_id;
	if (defined $dbxref_row) {
	    $dbxref_id = $dbxref_row->get_column('dbxref_id');
	}
	else {
	    $dbxref_row = $schema->resultset('General::Dbxref')
	                         ->new(
			                { 
 				          accession   => $accession, 
				          db_id       => $db_id,
				        }
			              );
	    $dbxref_row->insert()
	               ->discard_changes();
	    $dbxref_id = $dbxref_row->get_column('dbxref_id');
	}

	## Now it will search the template_id based in the template name and it will print them into the output file

	my $template_id = $templates_equiv_href->{$template_name};
	
	unless (defined $template_id) {
	    print STDERR "\nSKIPPING template_name=$template_name, it have not template_id.\n";
	}
	else {

	    ## Finally it will check if exists the relation for a concrete template_id and dbxref_id

	    my ($ge_template_dbxref_row) = $schema->resultset('GeTemplateDbxref')
		                                  ->search({ template_id => $template_id, dbxref_id => $dbxref_id });
	    
	    unless (defined $ge_template_dbxref_row) {
		print $out_xfh "$template_id\t$dbxref_id\t$arg_href->{'metadata_id'}\n";
	    }
	    else {
		$rel_dbxref_c++;
	    }	    
	}
	    
    }
    print STDERR "\n\t\t\tDone... file parsed.\n\n";
    if ($rel_dbxref_c > 0) {
	print STDERR "\t\tThere are $rel_dbxref_c template-dbxref relations that exist in the db, these template_id-dbxref_id have been skipped to create the dbload files.\n\n";
    }
    close($in_xfh);
    
    return $temp_dbxref;
}

=head2 create_probe_dbload_file

  Usage: my $probe_dbload_file = create_probe_dbload_file(arg_href);
  Desc: Parse the probe file and create a probe dbload file with the data
  Ret: $probe_dbload_file, a file name
  Args: $arg_href, a hash reference with keys=arg and value=values
  Side_Effects: Die if something is wrong
  Example: my $probe_dbload_file = create_probe_dbload_file($arg_href);

=cut

sub create_probe_dbload_file {
    my $arg_href = shift 
	|| die("DATA ERROR: None argument hash reference was supplied to create_probe_dbload_file.\n");

    ## Get the variables

    my $probe_file = $arg_href->{'probe_file'} 
       || die("DATA ERROR: None probe was supplied to create_probe_dbload_file.\n");
    my $platform_id = $arg_href->{'platform_id'}
       || die("DATA ERROR: None platform_id was supplied to create_probe_dbload_file.\n");
    my $metadata_id = $arg_href->{'metadata_id'} 
       || die("DATA ERROR: None metadata_id was supplied to create_probe_dbload_file.\n");

    ## Redundancy count

    my $probe_c = 0;

    ## Templates and file_id are conditionals so it will not die if they are no present.
    ## file_id also will be added to the database load file, so by default it will be \N (NULL VALUE)

    my $templates_equiv_href = $arg_href->{'templates_equiv'};
    my $file_id = $arg_href->{'file_id'} || '\N';

    ## It will use as path

    my $path = getcwd();

    open my $in_pfh, '<', $probe_file 
        || die("Sorry but I can not open the file:" . $probe_file . " ($!)\n");

    my $temp_probe = $path ."/temp_probe_dbload.tab";

    open my $out_pfh, '>', $temp_probe 
	|| die("Sorry but I can not open the file:$temp_probe ($!)\n");

    ## The template file should have only a column

    my $l = 0;
    while (<$in_pfh>) {
	$l++;
	print STDERR "\t\tParsing line:$l from file $probe_file\r";
	chomp($_);
	my @data = split(/\t/, $_);
	
	## Asign data to the variables using \N as default (NULL VALUE to load data into the database)

	my $probe_name = $data[0];
	my $probe_type = $data[1];
	
	my $template_id = '\N';
	if (defined $templates_equiv_href && defined $data[2]) {
	    $template_id = $templates_equiv_href->{$data[2]} || '\N';
	}
	else {
	    $template_id = '\N';
	}
	
	my $temp_start = $data[3] || '\N';
	my $temp_end = $data[4] || '\N';

	## Filter spaces in the begining and in the end

	$probe_name =~ s/^\s+//g;
	$probe_name =~ s/\s+$//g;

	## Check if exists in the database as: probe_name + platform_id.
	## If don't exists, it will add to the file to load data into the database

	my ($probe_row) = $schema->resultset('GeProbe')
	                         ->search( { probe_name => $probe_name, platform_id => $platform_id });
	
	unless (defined $probe_row) {
	    print $out_pfh "$platform_id\t$probe_name\t$probe_type\t$file_id\t$template_id\t$temp_start\t$temp_end\t$metadata_id\n";
	}
	else {
	    $probe_c++;
	}
    }
    print STDERR "\n\t\t\tDone... file parsed.\n\n";
    if ($probe_c > 0) {
	print STDERR "\t\tThere are $probe_c probes that exist into the database. These probes have been skipped to create the dbload file.\n\n";
    }

    close($in_pfh);

    return $temp_probe;
}

=head2 create_spot_dbload_file

  Usage: my $spot_dbload_file = create_spot_dbload_file(arg_href);
  Desc: Parse the spot file and create a spot dbload file with the data
  Ret: $spot_dbload_file, a file name
  Args: $arg_href, a hash reference with keys=arg and value=values
  Side_Effects: Die if something is wrong
  Example: my $spot_dbload_file = create_spot_dbload_file($arg_href);

=cut

sub create_spot_dbload_file {
    my $arg_href = shift 
	|| die("DATA ERROR: None argument hash reference was supplied to create_spot_dbload_file.\n");

    ## Get the variables

    my $spot_file = $arg_href->{'spot_file'} 
       || die("DATA ERROR: None spot file was supplied to create_spot_dbload_file.\n");
    my $probes_equiv_href = $arg_href->{'probes_equiv'}
       || die("DATA ERROR: None probe_equiv hash reference was supplied to create_spot_dbload_file.\n");
    my $metadata_id = $arg_href->{'metadata_id'} 
       || die("DATA ERROR: None metadata_id was supplied to create_probe_dbload_file.\n");

    ## Count the redundancy

    my $spot_rc = 0;

    my $path = getcwd();

    open my $in_spfh, '<', $spot_file 
        || die("Sorry but I can not open the file:" . $spot_file . " ($!)\n");

    my $temp_spot = $path ."/temp_spot_dbload.tab";

    open my $out_spfh, '>', $temp_spot 
	|| die("Sorry but I can not open the file:$temp_spot ($!)\n");

    my $l = 0;
    while (<$in_spfh>) {
	$l++;
	print STDERR "\t\tParsing line:$l from file $spot_file\r";
	chomp($_);
	my @data = split(/\t/, $_);
	
	## Asign data to the variables using \N as default (NULL VALUE to load data into the database)

	my $probe_name = $data[0];
	my $spot_name = $data[1];
	
	my $probe_id = $probes_equiv_href->{$probe_name}; 
	
	## Filter spaces in the begining and in the end

	$spot_name =~ s/^\s+//g;
	$spot_name =~ s/\s+$//g;

	## Check if exists in the database as: spot_name + probe_id.
	## If don't exists, it will add to the file to load data into the database

	if (defined $probe_id) {
	    my ($spot_row) = $schema->resultset('GeProbeSpot')
	                            ->search( { spot_name => $spot_name, probe_id => $probe_id });
	
	    unless (defined $spot_row) {
		print $out_spfh "$probe_id\t$spot_name\t$metadata_id\n";
	    }
	    else {
		$spot_rc++;
	    }
	}
	else {
	    print STDERR "DATA WARNING (line $l): Probe_name=$probe_name has not been parsed in previous files. SKIP SPOT.\n";
	}
    }
    print STDERR "\n\t\t\tDone... file parsed.\n\n";
    if ($spot_rc > 0) {
	print STDERR "\t\tThere are $spot_rc spots that exist into the database. These spots have been skipped to create the dbload file.\n\n";
    }
    close($in_spfh);

    return $temp_spot;
}

=head2 create_spot_coord_dbload_file

  Usage: my $spot_coord_dbload_file = create_spot_coord_dbload_file(arg_href);
  Desc: Parse the spot file and create a spot coordinate dbload file with the data
  Ret: $spot_coord_dbload_file, a file name
  Args: $arg_href, a hash reference with keys=arg and value=values
  Side_Effects: Die if something is wrong
  Example: my $spot_coord_dbload_file = create_spot_coord_dbload_file($arg_href);

=cut

sub create_spot_coord_dbload_file {
    my $arg_href = shift 
	|| die("DATA ERROR: None argument hash reference was supplied to create_spot_coord_dbload_file.\n");

    ## Get the variables

    my $spot_coord = $arg_href->{'spot_coordinates'} 
       || die("DATA ERROR: None spot coordinates was supplied to create_spot_coord_dbload_file.\n");
    my $spot_file = $arg_href->{'spot_file'}
       || die("DATA ERROR: None spot file was supplied to create_spot_coord_dbload_file.\n");
    my $spot_equiv_href = $arg_href->{'spots_equiv'}
       || die("DATA ERROR: None spots_equiv hash reference was supplied to create_spot_coord_dbload_file.\n");
    my $metadata_id = $arg_href->{'metadata_id'} 
       || die("DATA ERROR: None metadata_id was supplied to create_probe_dbload_file.\n");

    ## Count the redundancy

    my $coord_rc = 0;


    my @coord_types = split(/,/, $spot_coord);

    my $path = getcwd();

    open my $in_scfh, '<', $spot_file 
        || die("Sorry but I can not open the file:" . $spot_file . " ($!)\n");

    my $temp_spot_coord = $path ."/temp_spot_coord_dbload.tab";

    open my $out_scfh, '>', $temp_spot_coord 
	|| die("Sorry but I can not open the file:$temp_spot_coord ($!)\n");

    ## The template file should have only a column

    my $l = 0;
    while (<$in_scfh>) {
	$l++;
	print STDERR "\t\tParsing line:$l from file $spot_file\r";
	chomp($_);
	my @data = split(/\t/, $_);
	
	## Asign data to the variables using \N as default (NULL VALUE to load data into the database)

	my $probe_name = shift(@data);
	my $spot_name = shift(@data);

	$spot_name =~ s/^\s+//g;
	$spot_name =~ s/\s+$//g;

	my $probe_spot_id = $spot_equiv_href->{$spot_name};
	
	my $a = 0;
	foreach my $coord_type (@coord_types) {
	    my $coord_value = $data[$a];
	    my ($spot_coord_row) = $schema->resultset('GeProbeSpotCoordinate')
	                                  ->search( { probe_spot_id => $probe_spot_id, coordinate_type => $coord_type });
	
	    unless (defined $spot_coord_row) {
		print $out_scfh "$probe_spot_id\t$coord_type\t$coord_value\t$metadata_id\n";
	    }
	    else {
		$coord_rc++;
	    }
	    $a++;
	}
    }
    print STDERR "\n\t\t\tDone... file parsed.\n\n";
    if ($coord_rc > 0) {
	print STDERR "\t\tThere are $coord_rc spots_coordinates that exist into the database. These coords. have been skipped to create the dbload file.\n\n";
    }
    close($in_scfh);

    return $temp_spot_coord;
}




=head2 print_template

  Usage: print_template()
  Desc: print a template file
  Ret: none
  Args: none
  Side_Effects: create a file and exit of the script
  Example:

=cut

sub print_template {
    my $dir = `pwd`;
    chomp($dir);
    my $template_file = $dir . '/argument_template_dbload_file.bs';
    open my $TFH, '>', $template_file || die("Sorry, I can not open the file: $template_file.\n");

    print STDERR "PRINT ARGUMENT FIELD FILE option...\n";
    my $info = '# Notes for the Source data load to GEM Template data ()
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


#####################################################################################
# TEMPLATE_DATA #####################################################################
# FIELD	#		# DATA #	# DATA_STRINGENCY #	# EXAMPLE_AND_NOTES #
#####################################################################################

## For templates it will use a template file with one column (-f1: template_name)

*PLATFORM_NAME:	        []		#mandatory,single	example:[GeneChip Tomato Genome Array]
*TEMPLATE_FILE:         []		#mandatory,single	example:[TomGenChip_templates.tab]
*TEMPLATE_TYPE:         []              #mandatory,single       example:[Consensus sequence]

## To add a new internal reference (Dbiref) it will use a mapping file with two columns (-f1: template_name, -f2:iref_accession)

*IREF_MAPPING_FILE:     []		#optional,single	example:[TomGenChip_templates_iref.map]
*IPATH:                 []		#optional,single	example:[sgn.unigene.unigene_id]

## To associate external references (Dbxref) it will use a mapping file with two columns (-f1: template_name, -f2:dbxref_accession)

*XREF_MAPPING_FILE:     []              #optional,single       example:[TomGenChip_templates_xref.map]
*XDB_NAME:		[]		#optional,single       example:[DB:GenBank_accesss]

## To add a file or url (all the probes will have the same file_id)

*BASENAME:              []              #opptional,single       example:[GeneChipTomato_probes.fasta]
*DIRNAME:               []              #optional*,single      example:[expression/platforms/affy/]
*FILETYPE:		[]		#optional*,single      example:[fasta]
*ALT_FILENAME:          []              #optional,single        example:[AffyTomato_probes.fasta]
*COMMENT:               []              #optional,single        example:[This is a file with...]
*MD5CHECKSUM:           []              #optional,single        example:[ddcdf85aaa8431a2e05d322e8cdc6bda]
*URLSOURCE:             []              #optional*,single      example:[ftp://ftp.ncbi.nih.gov/pub/geo/DATA/MINiML/by_platform/GPL4741/]
*URLSOURCE_MD5CHECKSUM: []              #optional,single        example:[fd2d8bb0501a38deecaf559c7165f6d1]

## To add a probe it will use a probe file with two (-f1 probe_name, -f2 probe_type), three (-f3 template_name) 
## or five (-f4 template_start, -f5 template_end) depending if you want add only the probes, link them with 
## templates or also map them using coordinates

*PROBE_FILE:            []		#optional,single	example:[TomGenChip_probes.tab]

## To add a probe spot it will use a mapping spot file with a minimun of two columns (-f1 probe_name, -f2 spot_name)
## The rest of the column will be coordinates (from 0 to as many as you need)

*SPOT_FILE:             []              #optional,single       example:[TomGenChip_spots.tab]
*SPOT_COORDINATES:      []              #optional,multiple     example:[X,Y,Interrogation Position]

';

    print $TFH "$info";
    print STDERR "...done (printed template_file with the name: $template_file)\n\n";
    exit (1);
}


###
1;#
###
