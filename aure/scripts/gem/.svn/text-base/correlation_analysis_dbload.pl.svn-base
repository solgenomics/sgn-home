
#!/usr/bin/perl
=head1 NAME

 correlation_analysis_dbload.pl
 A script to parse correlation files and load in a database for gem schema (version.0.1.).

=cut

=head1 SYPNOSIS

  correlation_analysis_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -a <argument_file> [-T]

  To collect the data loaded report into a file:

  correlation_analysis_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -a <argument_file> [-T] > file.log

=head1 EXAMPLE:

 perl expression_dbload.pl -u aure -H localhost -D sandbox -a CorrelationArgumentFile.bs
          
=head2 I<Flags:>

=over

=item -a

B<data load file>               argument file in bs format (mandatory).

=item -u

B<user_loader>                  user that load the data (mandatory)

=item -H

B<database_host>                database host (mandatory if you want check the relations in the database)

=item -D

B<database_name>                database name (mandatory if you want check the relations in the database)

=item -X

B<print data load file>         print an argument file with examples of the data load file in bs format

=item -T

B<run as a test>                run the script as test

=item -h

B<help>                         print the help  

=back

=cut

=head1 DESCRIPTION

 This script copy correlation analysis data into:

    + gem.ge_experiment_analysis_group
    + gem.ge_experiment_analysis_member
    + gem.ge_correlation_analysis
    + gem.ge_correlation_member

 The correlation analysis files are composed by millions of lines. To a better way to work with them, this script will
 create subsets of 1,000,000 lines, and after that load each one using SQL copy command.
 
 
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

use Cwd;
use File::Basename;
use Getopt::Std;
use CXGN::GEM::Schema;
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
    print_argument_file();
}

## Checking the input arguments

my $loader_username = $opt_u || die("MANDATORY ARGUMENT ERROR: The -u <loader_username> argument was not supplied.\n");
my $dbname = $opt_D || die("MANDATORY ARGUMENT ERROR: The -D <database_name> argument was not supplied.\n");
my $dbhost = $opt_H || die("MANDATORY ARGUMENT ERROR: The -H <db_hostname> argument was not supplied.\n"); 
my $argument_file = $opt_e || die("MANDATORY ARGUMENT ERROR: The -e <argument_dataload_file> argument was not supplied.\n");

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


## It will store the different parse arguments in the arg_href hash reference.
## The first argument to store will be $schema (to use dbic) and $loader_username (to create a metadata object)

my $args_href = {
                  'schema'          => $schema,
		  'loader_username' => $loader_username,
		  'exp_group_name'  => [],

                };

## Now it will parse the argument file.

print STDERR "\nStep 3: Open and parse the argument file.\n";

$args_href = parse_arguments($argument_file, $args_href);

## Test mode, run as evaluation of the code

print STDERR "\nStep 4: Store (or store simulation for Test mode) the expression data into the database.\n";

if ($opt_T) {
    print STDERR "\nRunning the TEST MODE.\n\n";
    eval {

	## it will use this pipeline to get, create and copy the expression data

	store_pipeline($args_href);

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

    ## it will use this pipeline to get, create and copy the expression data

    store_pipeline($args_href);

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
      
      This script copy correlation analysis data into:

       + gem.ge_experiment_analysis_group
       + gem.ge_experiment_analysis_member
       + gem.ge_correlation_analysis
       + gem.ge_correlation_member

      The correlation analysis files are composed by millions of lines. To a better way to work with them, this script will
      create subsets of 1,000,000 lines, and after that load each one using SQL copy command.
    
    Usage: 

       correlation_analysis_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -a <argument_file> [-T]

     To collect the data loaded report into a file:

       correlation_analysis_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -a <argument_file> [-T] > file.log
      

    Example: 
      
        perl expression_dbload.pl -u aure -H localhost -D sandbox -a CorrelationArgumentFile.bs

    Flags:
      -u loader username      loader username (mandatory)
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory for check domains in the database)
      -D database name        sandbox or cxgn etc (mandatory for check domains in the database)
      -e argument file        data load file input file (mandatory)
      -T run as test          run this script as a test
      -X create dataload file create an argument file for a data load file (follow the instructions to fill it)
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

#######################
### PARSING METHODS ###
#######################

=head2 parse_arguments

  Usage: my $arg_href = parse_arguments($file, $arg_href);
  Desc: parse the argument from the argument file
  Ret: $arg_href, a hash reference with the arguments.
       keys=argument_type and value=hash reference with
            keys=argument_name and value=value

       exceptions are: schema, that supply an schema object
                     : loader_username, that supply a scalar   
  Args: $args_href, a hash reference with arguments
        $file, an argument file
  Side_Effects: die if something it is wrong
  Example: my $arg_href = parse_arguments($file, $arg_href);

=cut

sub parse_arguments {
    my $argument_file = shift ||
	die("None argument file was supplied to the parse_argument function.\n");
    my $args_href = shift ||
	die("None argument hash reference was supplied to the parse_argument_fuction.\n");

    open my $ifh, '<', $argument_file || die("Sorry, but I can not open the input file: $argument_file.\n");

    my $l = 0;

    ## It will check if the file have the right fields, DATA_TYPE is mandatory, so it will check the presence with dt variable       
    my ($exp_group_name, $dt);

    while(<$ifh>) {
		
	$l++; ## Line counter

	## It will do not read any line that start with #
	unless ($_ =~ m/^#/) {
	
	    ## First define the data_type or empty the variables defined for each type
	
	    if ($_ =~ m/\*DATA_TYPE:\s+\[(\w+)\]/) {
		$dt = $1;
	    }
	    elsif ($_ =~ m/\/\//) {
		$dt = '';
		$exp_group_name = '';
	    }

	    if (defined $dt) {
		if ($dt eq 'experiment_group_analysis') {
		    if ($_ =~ m/\*EXPERIMENT_GROUP_NAME:\s+\[(.+)\]/) {
			$exp_group_name = $1 ||
			    die("MANDATORY DATA ERROR (line $l): None experiment_group_name data was detailed.\n");
			
			## It will check if exists in the database later, because if it do not exists it will create a new group

			my $analysis_group_href = {};
			push @{$args_href->{'exp_group_name'}}, $exp_group_name;
			$args_href->{$exp_group_name} = $analysis_group_href;
		    }
		    elsif ($_ =~ m/\*EXPERIMENT_GROUP_DESC:\s+\[(.+)\]/ ) {
			my $exp_group_description = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None experiment_group_description was supplied");
			
			$args_href->{$exp_group_name}->{'experiment_group_description'} = $exp_group_description;
		    }
		    elsif ($_ =~ m/\*EXPERIMENT_NAME_LIST:\s+\[(.+)\]/ ) {
			my $exp_name_list = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None experiment_name_list was supplied");
			
			$args_href->{$exp_group_name}->{'experiment_name_list'} = $exp_name_list;
		    }
		}
	
		elsif ($dt eq 'correlation_analysis') {
		    if ($_ =~ m/\*EXPERIMENT_GROUP_NAME:\s+\[(.+)\]/) {
			$exp_group_name = $1 ||
			    die("MANDATORY DATA ERROR (line $l): None experiment_group_name data was detailed.\n");
			
			## It will check if exists in the database later, because if it do not exists it will create a new group

			unless (defined $args_href->{$exp_group_name} ) {
			    die("DATA ERROR (line $l): experiment_group_name has not been parsed before.\n");
			}
		    }
		    elsif ($_ =~ m/\*PLATFORM_NAME:\s+\[(.+)\]/ ) {
			my $platf_name = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None platform_name was supplied");
			
			my ($platf_row) = $args_href->{'schema'}
			                            ->resultset('GePlatform')
					            ->search({ platform_name => $platf_name });
			if (defined $platf_row) { 
			    $args_href->{$exp_group_name}->{'platform_id'} = $platf_row->get_column('platform_id');
			}
			else {
			    die("DATA ERROR: platform_name=$platf_name doesn't exist in gem.ge_platform table.\n");
			}
		    }
		    elsif ($_ =~ m/\*DATASET_NAME:\s+\[(.+)\]/ ) {
			my $dataset_name = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None dataset_name was supplied");
			
			my ($sample_row) = $args_href->{'schema'}
			                            ->resultset('BsSample')
					            ->search({ sample_name => $dataset_name });
			if (defined $sample_row) { 
			    $args_href->{$exp_group_name}->{'sample_id'} = $sample_row->get_column('sample_id');
			}
			else {
			    die("DATA ERROR: dataset_name=$dataset_name doesn't exist in biosource.bs_samplem table.\n");
			}
		    }
		    elsif ($_ =~ m/\*METHODOLOGY:\s+\[(.+)\]/ ) {
			my $methodology = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None methodology was supplied");
			
			$args_href->{$exp_group_name}->{'methodology'} = $methodology;
		    }
		    elsif ($_ =~ m/\*DESCRIPTION:\s+\[(.+)\]/ ) {
			my $description = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None experiment_name_list was supplied");
			
			$args_href->{$exp_group_name}->{'description'} = $description;
		    }
		    elsif ($_ =~ m/\*CORR_ANALYSIS_FILES:\s+\[(.+)\]/ ) {
			my $corr_analysis_files = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None correlation_analysis_file were supplied");
			
			my @files = split(/,/, $corr_analysis_files);
			
			foreach my $file (@files) {
			    unless (-s $file) {
				die("DATA ERROR (line $l): file=$file do not exists or has zero size.\n");
			    }
			}
			$args_href->{$exp_group_name}->{'correlation_analysis_file_list'} = \@files;
		    }
		    elsif ($_ =~ m/\*CORR_TYPE:\s+\[(.+)\]/ ) {
			my $corr_type = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None correlation_type was supplied");
			
			$args_href->{$exp_group_name}->{'correlation_type'} = $corr_type;
		    }
		    elsif ($_ =~ m/\*CORR_VALUE_FILTER:\s+\[(.+)\]/ ) {
			my $corr_value_filter = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None correlation_value_filter was supplied");
			
			$args_href->{$exp_group_name}->{'correlation_value_filter'} = $corr_value_filter;
		    }
		}		
	    }
	}
    }

    unless (defined $dt) {
	die("WRONG INPUT ARGUMENT FILE: The argument file have not any DATA_TYPE field. Please check that it have the right format.\n");
    }

    return $args_href;
}

=head2 find_or_store_experiment_group

  Usage: my %exp_group_id = find_or_store_experiment_group($args_href);
  Desc: find the experiment_group_id for the experiments from the arguments
        hash reference.
        If they does not exist in the db, it will add them
  Ret: A hash with keys=experiment_group_name
                   value=experiment_group_id
  Args: $arg_href, a hash reference with the arguments.
       keys=argument_type and value=hash reference with
            keys=argument_name and value=value

       exceptions are: schema, that supply an schema object
                     : loader_username, that supply a scalar
                       metadata_id, a scalar
  Side_Effects: die if something it is wrong
  Example: my %exp_group_id = find_or_store_experiment_group($args_href);

=cut

sub find_or_store_experiment_group {
    my $args_href = shift ||
	die("None argument hash reference was supplied to the find_or_store_experiment_group function.\n");

    my %exp_group_id;

    my @exp_group_names = @{ $args_href->{'exp_group_name'} };

    foreach my $exp_group_name (@exp_group_names) {

	my ($exp_analysis_group_row) = $args_href->{'schema'}
	                                         ->resultset('GeExperimentAnalysisGroup')
					         ->search( { group_name => $exp_group_name } );
	
	my $exp_analysis_group_id;

	if (defined $exp_analysis_group_row) {
	    $exp_analysis_group_id = $exp_analysis_group_row->get_column('experiment_analysis_group_id');
	    
	}
	else {

	    print STDERR "\t\tMESSAGE: experiment_group_name=$exp_group_name doesn't exist in gem.ge_experiment_analysis_grouup table.\n";
	    print STDERR "\t\tINSERTING a new experiment_analysis_group in the database.\n";

	    my $exp_analysis_group_row = $args_href->{'schema'}
	                                           ->resultset('GeExperimentAnalysisGroup')
			  		           ->new( { group_name => $exp_group_name } );

	    my $group_descript = $args_href->{$exp_group_name}->{'experiment_group_description'};

	    if (defined $group_descript) {
		$exp_analysis_group_row->set_column( group_description => $group_descript);
	    }

	    my $metadata_id = $args_href->{'metadata_id'} ||
		die("DATA ARGUMENT ERROR: metadata_id was not supplied in the argument hash to find_or_store_experiment_group.\n");
	    
	    $exp_analysis_group_row->set_column( metadata_id => $metadata_id);

	    $exp_analysis_group_row->insert()
		                   ->discard_changes();

	    $exp_analysis_group_id = $exp_analysis_group_row->get_column('experiment_analysis_group_id');
	    
	    my $group_members = $args_href->{$exp_group_name}->{'experiment_name_list'} ||
		die("DATA ARGUMENT ERROR: exper_name_list was not supplied in the argument hash to find_or_store_experiment_group.\n");

	    my @exp_names = split(/,/, $group_members);

	    foreach my $exp_name (@exp_names) {
		my ($exp_row) = $args_href->{'schema'}
		                          ->resultset('GeExperiment')
					  ->search({ experiment_name => $exp_name });
		
		if (defined $exp_row) {
		    my $exp_analysis_member_row = $args_href->{'schema'}
	                                           ->resultset('GeExperimentAnalysisMember')
			  		           ->new( 
			                                  { 
						            experiment_analysis_group_id => $exp_analysis_group_id,
							    experiment_id                => $exp_row->get_column('experiment_id'),
							  } 
						        );
		    $exp_analysis_member_row->insert();
		}
		else {
		    die("DATA ERROR: experiment_name=$exp_name doesn't exist in gem.ge_experiment table.\n");
		}
	    }
	}
	$exp_group_id{$exp_group_name} = $exp_analysis_group_id;
    }
    return %exp_group_id;
}


=head2 create_correlation_member_dbload
    
  Usage: my @corr_analysis_member_dbload = create_correlation_member_dbload($args_href);
  Desc: parse the correlation files and create correlation_member_dbload files
  Ret: a list of file names
  Args: $ags_href, arguments_hash reference
  Side_Effects: print status messages
  Example: my @corr_analysis_member_dbload = create_correlation_member_dbload($args_href);

=cut

sub create_correlation_member_dbload {
    my $args_href = shift ||
	die("None argument hash reference was supplied to the create_template_dbload fuction.\n");

    ## Create the array to store the file names

    my @corr_dbload_files;

    ## Get the data stored in the argument hash

    my @exp_group_names = @{ $args_href->{'exp_group_name'} };
    my %corr_analysis = %{ $args_href->{'correlation_analysis_id'} };
    
    my $metadata_id = $args_href->{'metadata_id'} ||
	die("DATA ARGUMENT ERROR: metadata_id was not supplied in the argument hash to create_correlation_member_dbload.\n");

    foreach my $exp_group_name (@exp_group_names) {

	my @files = @{ $args_href->{$exp_group_name}->{'correlation_analysis_file_list'} };

	my $corr_type = $args_href->{$exp_group_name}->{'correlation_type'} || '\N';
	my $corr_value_filter = $args_href->{$exp_group_name}->{'correlation_value_filter'};
	my $platf_id = $args_href->{$exp_group_name}->{'platform_id'} ||
	    die("DATA ARGUMENT ERROR: platform_id was not supplied in the argument hash to create_correlation_member_dbload.\n");;
	my $sample_id = $args_href->{$exp_group_name}->{'sample_id'} || '\N';
	my $corr_id = $corr_analysis{$exp_group_name};

	## get all the templates associated to this platform

	print STDERR "\t\tGetting template_ids for templates associated to platform_id=$platf_id.\n";
	my @templ_rows = $args_href->{'schema'}
	                           ->resultset('GeTemplate')
				   ->search({ platform_id => $platf_id });

	my %templ_id;
	foreach my $templ_row (@templ_rows) {
	    $templ_id{$templ_row->get_column('template_name')} = $templ_row->get_column('template_id');
	}
	
	## First, if exist filter see if the value is in the filter range or use by default.

	my $low_val = -$corr_value_filter || -10;
	my $high_val = $corr_value_filter || 10; 

	print STDERR "\n\tFilter value used: LowCutoff=$low_val and HighCutOff=$high_val\n";

	foreach my $file (@files) {
	    
	    my $dir = getcwd;
	    my $basename = basename($file);

	    open my $ifh, '<', $file || die("Sorry, I can not open the file:$file.\n");

	    my ($l, $a, $n, $t, $z) = (1, 1, 1, 0, 1000000); 
	    
	    my ($ofh, $filename);

	    while (<$ifh>) {
		chomp($_);

		unless ($_ =~ m/^#/) {

		    $t++;

		    ## Filter the end of the line

		    $_ =~ s/\s*$//;

		    ## Open a new fh when the counter is 1

		    if ($t == 1 && $a == 1) {
			$filename = $dir . "/temp_" . $basename . "_dbload_" . $n . ".tab" ;
			
			print STDERR "\n\t$filename has been created.\n";
			
			open $ofh, '>', $filename || die("Sorry, I can not open the file:$filename.\n");
			
			push @corr_dbload_files, $filename;
		    }
	     		    		    
		    my ($template_name_a, $template_name_b, $corrval) = split(/\t/, $_);
		    
		    unless (defined $template_name_a && defined $template_name_b && defined $corrval) {
			print STDERR "\nWARNING: line:$l has not three fields (#template_name_a\t#template_name_b\t#corr_value)\n";
			print STDERR "\t[skipping $_]\n";
		    }
		    else {

			## Print in the out file and increase the counter

			my $template_id_a = $templ_id{$template_name_a};
			my $template_id_b = $templ_id{$template_name_b};

			if (defined $template_id_a && defined $template_id_b && defined $corrval) {
			    unless ($low_val < $corrval && $high_val > $corrval) {
				print $ofh "$corr_id\t$template_id_a\t$template_id_b\t$corrval\t$corr_type\t$sample_id\t$metadata_id\n";
				$a++;
			    }
			}
		    }
		    
		    ## Print the status message
		    
		    my $catch = $a-1;
		    print STDERR "\t\tProcessing line:$l in file:$file ( values catched= $catch )         \r";
		    
                    ## When the counter is equal to the limit, reset that to 1 (it will create a new one)

		    if ($a == $z) {
			close $ofh;
			$a = 1;
			$t = 0;
			$n++;			
		    }
		    $l++;
		}
	    }
	    print STDERR "\n";
	}
    }
    return @corr_dbload_files;
}


=head2 store_pipeline

  Usage: store_pipeline($args_href)
  Desc: copy the data parsed using $dbh->do functions
        print status messages
  Ret: none
  Args: $ags_href, arguments_hash reference
  Side_Effects: print status messages
  Example:

=cut 

sub store_pipeline {
    my $args_href = shift ||
	die("None argument hash reference was supplied to the store_pipeline fuction.\n");

    my $metadbdata = CXGN::Metadata::Metadbdata->new( $args_href->{'schema'}, $args_href->{'loader_username'} );
    $metadbdata->store();
	
    ## Second add the metadata_id to the arg_href
	
    $args_href->{'metadata_id'} = $metadbdata->get_metadata_id();

    ## Third get the experiment_analysis_group_ids

    print STDERR "\t4.1) Getting experiment_analysis_group_id ...\n";
    my %exp_group_id = find_or_store_experiment_group($args_href);
    print STDERR "\t\t... done\n\n";

    my %corr_analysis;
    print STDERR "\t4.2) Adding new correlation_analysis to gem.ge_correlation_analysis table.\n";
    foreach my $exp_group_name (sort keys %exp_group_id) {
	my $exp_group_id = $exp_group_id{$exp_group_name};

	my $corr_analysis_row = $args_href->{'schema'}
	                                  ->resultset('GeCorrelationAnalysis')
			  		  ->new( 
			                          { 
						    experiment_analysis_group_id => $exp_group_id,
						    metadata_id                  => $args_href->{'metadata_id'},
						  } 
					       );
	my $methodology = $args_href->{$exp_group_name}->{'methodology'};

	if (defined $methodology) {
	    $corr_analysis_row->set_column( methodology => $methodology );
	}
	
	my $description = $args_href->{$exp_group_name}->{'description'};

	if (defined $description) {
	    $corr_analysis_row->set_column( description => $description );
	}

	$corr_analysis_row->insert()
	                  ->discard_changes();

	$corr_analysis{$exp_group_name} = $corr_analysis_row->get_column('correlation_analysis_id');
    }
    $args_href->{'correlation_analysis_id'} = \%corr_analysis;
    print STDERR "\t\t... done.\n\n";

    print STDERR "\t4.3) Parsing and creating the correlation_member_dbload files.\n";
    my @corr_member_dbload = create_correlation_member_dbload($args_href);

    my $corr_dbload_n = scalar(@corr_member_dbload);
    print STDERR "\t\t... done ($corr_dbload_n files were parsed)\n\n";

    ## Forth load the files into the database.

    my $last_corr_member_id = $all_last_ids_href->{'gem.ge_correlation_member_correlation_member_id_seq'};
 
    print STDERR "\t4.4) Copying correlation_member_dbload files (if there are any).\n\n";
    foreach my $corr_member_file (@corr_member_dbload) {
	$schema->storage()
	       ->dbh()
	       ->do("COPY gem.ge_correlation_member ( correlation_analysis_id,
                                                      template_a_id, 
                                                      template_b_id, 
                                                      correlation_value, 
                                                      correlation_type,  
                                                      dataset_id,  
                                                      metadata_id) 
                         FROM '$corr_member_file'");

	my $curr_corr_member_id = $schema->resultset('GeCorrelationMember')
		                         ->search(undef, { order_by => "correlation_member_id DESC"})
		      	                 ->first()
				         ->get_column('correlation_member_id');

	my $corr_load_n = $curr_corr_member_id - $last_corr_member_id;
	$last_corr_member_id = $curr_corr_member_id;

	print STDERR "\t\t$corr_load_n expression values from\n\t$corr_member_file file\n";
	print STDERR "\t\t\thas been loaded into gem.ge_correlation_member table.\n";
    }
    print STDERR "\n";
}




=head2 print_argument_file

  Usage: print_argument_file()
  Desc: print a argument file
  Ret: none
  Args: none
  Side_Effects: create a file and exit of the script
  Example:

=cut

sub print_argument_file {
    my $dir = `pwd`;
    chomp($dir);
    my $argument_file = $dir . '/correlation_analysis_argument_loadfile.bs';
    open my $TFH, '>', $argument_file || die("Sorry, I can not open the file: $argument_file.\n");

    print STDERR "PRINT GEM (correlation analysis) DATA LOAD FILE option...\n";
    my $info = '# Notes for the Source data load to GEM Expression data ()
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
#               Note1: An empty TEMPL_EXPRESS_FILE in the "experiment_expression" will take the 
#                      all the template_id and target_id associated to a experiment and will calculate
#                      this values. 
#           
# NOTE FOR DATABASE CURATORS: To load all this data into the database, use the platform_dbload.pl script
#


#####################################################################################
# EXPRESSION_DATA ###################################################################
# FIELD	#		# DATA #	# DATA_STRINGENCY #	# EXAMPLE_AND_NOTES #
#####################################################################################

## Definition of the experiment group used in the analysis

*DATA_TYPE:		[]		#mandatory,single	example:[experiment_group_analysis]
*EXPERIMENT_GROUP_NAME: []              #mandatory,single       example:[]
*EXPERIMENT_GROUP_DESC: []              #optional,single        example:[] 
*EXPERIMENT_NAME_LIST:  []              #mandatory,multiple     example:[]
//

## Definition of the correlation analysis arguments

*DATA_TYPE:		[]		#mandatory,single	example:[correlation_analysis]
*EXPERIMENT_GROUP_NAME: []              #mandatory,single       example:[]
*PLATFORM_NAME:         []              #mandatory,single       example:[]
*METHODOLOGY:           []              #optional,single        example:[]
*DESCRIPTION:           []              #optional,single        example:[]
*CORR_ANALYSIS_FILES:   []              #mandatory,multiple     example:[]     
*CORR_TYPE:             []              #mandatory,single       example:[]
*CORR_VALUE_FILTER:     []              #optional,single        example:[]
//

';

    print $TFH "$info";
    print STDERR "...done (printed argument_file with the name: $argument_file)\n\n";
    exit (1);
}




####
1;##
####
