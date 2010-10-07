#!/usr/bin/perl

=head1 NAME

 expression_dbload.pl
 A script to parse expression files and load in a database for gem schema (version.0.1.).

=cut

=head1 SYPNOSIS

  expression_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -e <argument_file> [-T]

  To collect the data loaded report into a file:

  expression_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -e <argument_file> [-T] > file.log

=head1 EXAMPLE:

 perl expression_dbload.pl -u aure -H localhost -D sandbox -e TobEA_Expression.bs
          
=head2 I<Flags:>

=over

=item -e

B<data load file>               argument data load file in bs format (mandatory).

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

 This script load data for expression

      In a .bs format (text format with fields delimited by [] ). To print a template file use the -X option and fill following the
      instructions. The load data file can have one or more of these data types. 

      There are four different ways to use this script:

      1) Load probe_expression data (associate to a target_element)
    
      2) Load template_expression data (associate to an hybridization)

      3-A) Load experiment_expression data (associate to an experiment)

      3-B) Calculate mean, median, sd and cv for all the hybridization of an experiment, create files and
           load them.

      1, 2 and 3 can be used at the same time, it depends of the argument in the -e <argument_file>     
 
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
use CXGN::GEM::Platform;
use CXGN::GEM::Hybridization;
use CXGN::GEM::Target;
use CXGN::GEM::Experiment;
use CXGN::Biosource::Sample;
use CXGN::DB::InsertDBH;
use CXGN::Metadata::Metadbdata;
use Statistics::Descriptive;

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
    print STDERR "\n\n\tRunning in TEST MODE... executing rollback.\n";
    $schema->set_sqlseq_values_to_original_state($all_last_ids_href);
    print STDERR "\n\tRunning in TEST MODE... executing set the last ids.\n\n";

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
      This script load data for expression

      In a .bs format (text format with fields delimited by [] ). To print a template file use the -X option and fill following the
      instructions. The load data file can have one or more of these data types.

      There are four different ways to use this script:

      1) Load probe_expression data (associate to a target_element)
    
      2) Load template_expression data (associate to an hybridization)

      3-A) Load experiment_expression data (associate to an experiment)

      3-B) Calculate mean, median, sd and cv for all the hybridization of an experiment, create files and
           load them.

      1, 2 and 3 can be used at the same time, it depends of the argument in the -e <argument_file> 

    
    Usage: 
       expression_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -e <argument_file> [-T]

      To collect the data loaded report into a file:

       expression_dbload [-h] [-X] -u <loader_username> -D <dbname> -H <dbhost> -e <argument_file> [-T] > file.log

    Example: 
      perl expression_dbload.pl -u aure -H localhost -D sandbox -p expression_argument.bs


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
    
    my ($target_el_name, $hyb_name, $experiment_name, $dt);

    my (%prob_expr, %temp_expr, %exp_expr);

    $args_href->{'probe_expression'} = \%prob_expr;
    $args_href->{'template_expression'} = \%temp_expr;
    $args_href->{'experiment_expression'} = \%exp_expr;

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
		$target_el_name = '';
		$hyb_name = '';
		$experiment_name = '';
	    }
	    
	    if (defined $dt) {
		if ($dt eq 'probe_expression') {
		    if ($_ =~ m/\*TARGET_ELEMENT_NAME:\s+\[(.+)\]/) {
			$target_el_name = $1 ||
			    die("MANDATORY DATA ERROR (line $l): None target_element_name data was detailed.\n");
			
			## Now it will create a new hash reference to store the data associated with this target_element
			## First, get the target_element_id or die

			my ($target_el_row) = $args_href->{'schema'}
		                                       ->resultset('GeTargetElement')
			                               ->search( { target_element_name => $target_el_name } );

			if (defined $target_el_row) {
			    my $target_el_id = $target_el_row->get_column('target_element_id');
			    my $target_href = { 'target_element_id'  => $target_el_id };
			    $prob_expr{$target_el_name} = $target_href;
			}
			else {
			    die("DATA ERROR (line $l):target_element_name=$target_el_name does not exist gem.ge_target_element table.\n");
			}
		    }
		    elsif ($_ =~ m/\*PROBE_EXPRESS_FILE:\s+\[(.+)\]/ ) {
			my $prob_express_file = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None probe_express_file was supplied");
			
			unless (-s $prob_express_file) {
			    die("DATA ERROR (line $l): probe_expression_file:$prob_express_file doesn't exist or has zero size.\n");
			}
			else {
			    $prob_expr{$target_el_name}->{'probe_expression_file'} = $prob_express_file;
			}
		    }
		    elsif ($_ =~ m/\*PROBE_EXPRESS_HEADERS:\s+\[(.+)\]/ ) {
			my $prob_express_head = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None probe_express_headers was supplied");
			
			$prob_expr{$target_el_name}->{'probe_expression_header'} = $prob_express_head;
		    }
		    elsif ($_ =~ m/\*DATASET_NAME\s+\[(.+)\]/ ) {
			my $dataset_name = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None dataset_name was supplied");
			
			my ($sample_row) = $args_href->{'schema'}
		                                     ->resultset('BsSample')
				       		     ->search( { sample_name => $dataset_name } );

			if (defined $sample_row) {
			    my $sample_id = $sample_row->get_column('sample_id');
			    $prob_expr{$target_el_name}->{'sample_id'} = $sample_id;
			}
			else {
			    die("DATA ERROR (line $l): dataset_name=$dataset_name does not exist in biosource.bs_sample table.\n");
			}
		    }
		    elsif ($_ =~ m/\*SIGNAL_TYPE:\s+\[(.+)\]/ ) {
			my $signal_type = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None signal_type was supplied");
			
			$prob_expr{$target_el_name}->{'signal_type'} = $signal_type;
		    }
		    elsif ($_ =~ m/\*BACKGROUND_TYPE:\s+\[(.+)\]/ ) {
			my $backg_type = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None background_type was supplied");
			
			$prob_expr{$target_el_name}->{'background_type'} = $backg_type;
		    }
		}
		elsif ($dt eq 'template_expression') {
		    if ($_ =~ m/\*HIBRIDIZATION:\s+\[(.+)\+(.+)\]/) {
			my $target_name = $1 ||
			    die("MANDATORY DATA ERROR (line $l): None target_name data was detailed.\n");
			my $platform_name = $2 ||
			    die("MANDATORY DATA ERROR (line $l): None platform_name data was detailed.\n");
		    
			## Hybridization name is a combination of target_name and platform_name bind for a + sign

			$hyb_name = $target_name . "+" . $platform_name;
		    
			my ($target_id, $platform_id);

			## Now it will create a new hash reference to store the data associated with this hybridization
			## First, get the target_id and platform_id or die

			my ($target_row) = $args_href->{'schema'}
		                                     ->resultset('GeTarget')
			                             ->search( { target_name => $target_name } );

			if (defined $target_row) {
			    $target_id = $target_row->get_column('target_id');
			}
			else {
			    die("DATA ERROR (line $l): target_name=$target_name does not exist in gem.ge_target table.\n");
			}
			
			my ($platform_row) = $args_href->{'schema'}
		                                       ->resultset('GePlatform')
			                               ->search( { platform_name => $platform_name } );

			if (defined $platform_row) {
			    $platform_id = $platform_row->get_column('platform_id');
			}
			else {
			    die("DATA ERROR (line $l): platform_name=$target_name does not exist in gem.ge_platform table.\n");
			}
			
			my ($hyb_row) = $args_href->{'schema'}
		                                  ->resultset('GeHybridization')
			                          ->search( { target_id => $target_id, platform_id => $platform_id } );

			if (defined $hyb_row) {
			    my $hyb_id = $target_row->get_column('target_id');
			    my $hyb_href = { 'hybridization_id'  => $hyb_id };

			    $temp_expr{$hyb_name} = $hyb_href;
			}
			else {
			    die("DATA ERROR (line $l):It doesnot exist hyb_row with target_id=$target_id & platform_id=$platform_id.\n");
			}
		    }
		    elsif ($_ =~ m/\*TEMPL_EXPRESS_FILE:\s+\[(.+)\]/ ) {
			my $temp_express_file = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None template_express_file was supplied");

			unless (-s $temp_express_file) {
			    die("DATA ERROR (line $l): template_expression_file:$temp_express_file doesn't exist or has zero size.\n");
			}
			else {
			    $temp_expr{$hyb_name}->{'template_expression_file'} = $temp_express_file;
			}
		    }
		    elsif ($_ =~ m/\*TEMPL_EXPRESS_HEADERS:\s+\[(.+)\]/ ) {
			my $temp_express_head = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None template_express_headers was supplied");

			$temp_expr{$hyb_name}->{'template_expression_headers'} = $temp_express_head;
		    }
		    elsif ($_ =~ m/\*DATASET_NAME\s+\[(.+)\]/ ) {
			my $dataset_name = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None dataset_name was supplied");
			
			my ($sample_row) = $args_href->{'schema'}
		                                     ->resultset('BsSample')
				  	  	     ->search( { sample_name => $dataset_name } );

			if (defined $sample_row) {
			    my $sample_id = $sample_row->get_column('sample_id');
			    $temp_expr{$hyb_name}->{'sample_id'} = $sample_id;
			}
			else {
			    die("DATA ERROR (line $l): dataset_name=$dataset_name does not exist in biosource.bs_sample table.\n");
			}
		    }
		    elsif ($_ =~ m/\*SIGNAL_TYPE:\s+\[(.+)\]/ ) {
			my $signal_type = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None signal_type was supplied");
			
			$temp_expr{$hyb_name}->{'signal_type'} = $signal_type;
		    }
		    elsif ($_ =~ m/\*STAT_VALUE_TYPE:\s+\[(.+)\]/ ) {
			my $stat_type = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None stats_value_type was supplied");
			
			$temp_expr{$hyb_name}->{'stat_value_type'} = $stat_type;
		    }
		}
		elsif ($dt eq 'experiment_expression') {
		    if ($_ =~ m/\*EXPERIMENT_NAME:\s+\[(.+)\]/) {
			$experiment_name = $1 ||
			    die("MANDATORY DATA ERROR (line $l): None experiment_name data was detailed.\n");

			## Now it will create a new hash reference to store the data associated with this target_element
			## First, get the target_element_id or die

			my ($experiment_row) = $args_href->{'schema'}
		                                         ->resultset('GeExperiment')
			                                 ->search( { experiment_name => $experiment_name } );

			if (defined $experiment_row) {
			    my $exp_id = $experiment_row->get_column('experiment_id');
			    my $exp_href = { 'experiment_id'  => $exp_id };
			    $exp_expr{$experiment_name} = $exp_href;
			}
			else {
			    die("DATA ERROR (line $l): experiment_name=$experiment_name does not exist in gem.ge_experiment table.\n");
			}
		    }
		    elsif ($_ =~ m/\*PLATFORM_NAME:\s+\[(.+)\]/) {
			my $platform_name = $1 ||
			    die("MANDATORY DATA ERROR (line $l): None platform_name data was detailed.\n");
		    
			## Now it will create a new hash reference to store the data associated with this target_element
			## First, get the target_element_id or die

			my ($platform_row) = $args_href->{'schema'}
		                                       ->resultset('GePlatform')
			                               ->search( { platform_name => $platform_name } );

			if (defined $platform_row) {
			    my $platf_id = $platform_row->get_column('platform_id');
			    $exp_expr{$experiment_name}->{'platform_id'} = $platf_id;
			}
			else {
			    die("DATA ERROR (line $l): experiment_name=$experiment_name does not exist in gem.ge_experiment table.\n");
			}
		    }
		    elsif ($_ =~ m/\*TEMPL_EXPRESS_FILE:\s+\[(.+)\]/ ) {
			my $temp_express_file = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None template_express_file was supplied");
			
			unless (-s $temp_express_file) {
			    die("DATA ERROR (line $l): template_expression_file:$temp_express_file doesn't exist or has zero size.\n");
			}
			else {
			    $exp_expr{$experiment_name}->{'template_expression_file'} = $temp_express_file;
			}
		    }
		    elsif ($_ =~ m/\*TEMPL_EXPRESS_HEADERS:\s+\[(.+)\]/ ) {
			my $temp_express_head = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None template_express_headers was supplied");

			$exp_expr{$experiment_name}->{'template_expression_header'} = $temp_express_head;
		    }
		    elsif ($_ =~ m/\*DATASET_NAME\s+\[(.+)\]/ ) {
			my $dataset_name = $1 ||
			    die("MANDATORY DATA ERROR (line $1): None dataset_name was supplied");
			
			my ($sample_row) = $args_href->{'schema'}
		                                     ->resultset('BsSample')
						     ->search( { sample_name => $dataset_name } );

			if (defined $sample_row) {
			    my $sample_id = $sample_row->get_column('sample_id');
			    $exp_expr{$experiment_name}->{'sample_id'} = $sample_id;
			}
			else {
			    die("DATA ERROR (line $l): dataset_name=$dataset_name does not exist in biosource.bs_sample table.\n");
			}
		    }
		}
	    }
	}
    }
    return $args_href;
}

=head2 create_probe_dbload

  Usage: my @probe_filedata = create_probe_dbload($arg_href);
  Desc: parse the probe_files contained in the argument hash reference 
  Ret: a list of file names to load into the database
  Args: $arg_href, a hash reference with the arguments.
       keys=argument_type and value=hash reference with
            keys=argument_name and value=value

       exceptions are: schema, that supply an schema object
                     : loader_username, that supply a scalar
                       metadata_id, a scalar
  Side_Effects: die if something it is wrong
  Example: my @probe_filedata = create_probe_dbload($arg_href);

=cut

sub create_probe_dbload {
    my $args_href = shift ||
	die("None argument hash reference was supplied to the create_probe_dbload function.\n");

    ## Create the array to store the file names

    my @probe_dbload_files;

    ## Get the data stored in the argument hash reference for probe_expression and get all the file names

    my %prob_express = %{ $args_href->{'probe_expression'}};  

    foreach my $targetel_name (sort keys %prob_express) {

	my $file = $prob_express{$targetel_name}->{'probe_expression_file'};
	
	## Create the out file

	my $path = `pwd`;
	my $basename = basename($file);

	chomp($path);

	my $outfile = $path . "/" . $basename . ".probe_dbload.tab";

	open my $ipfh, '<', $file || die("Sorry, I can not open the probe_expression_file=$file.\n");
	
	open my $opfh, '>', $outfile || die("Sorry, I can not open the probe_expression_file output=$outfile.\n");

	my $header_line = $prob_express{$targetel_name}->{'probe_expression_header'};
	my @headers = split(/,/, $header_line);

	my $l = 0;

	## Parse the file

	while (<$ipfh>) {
	    $l++;
	    chomp($_);
	    unless ($_ =~ m/^#/) {
		my @data = split(/\t/, $_);
	    
		if (scalar(@data) != scalar(@headers)) {
		    die("DATA PROBE ERROR: File=$file has different number of column than headers in the line=$l.\n");
		}

		print STDERR "\t\tProcessing file:$file line $l (probe_name=$data[0])                 \r";

		## To put all the data in a hash to get the correspondences

		my %datacol;
		my $a = 0;
		foreach my $head (@headers) {
		    $datacol{$head} = $data[$a];
		    $a++;
		}

		## Now it will replace probe_name by probe_id

		my ($probe_row) = $args_href->{'schema'}
	                                    ->resultset('GeProbe')
				   	    ->search( { probe_name => $datacol{'probe_name'} } );
	    
		unless (defined $probe_row) {
		    warn("WARNING: probe_name=$datacol{'probe_name'} does not exists into the table: gem.ge_probe. SKIPING DATA.\n")
		}
		else {
		    my $probe_id = $probe_row->get_column('probe_id');

		    my $error_head = "PROBE EXPRESS PARSING ERROR (file:$file, line:$l)";

		    ## Die if there aren't some of the not null variables

		    my $target_e_id = $prob_express{$targetel_name}->{'target_element_id'} || 
			die("$error_head: None target_element_id was supplied to parse_probe_files.\n");
		    
		    my $signal = $datacol{'signal'} || 
			die("$error_head: None signal was supplied to parse_probe_files.\n");

		    my $signal_type = $prob_express{$targetel_name}->{'signal_type'} ||
			die("$error_head: None signal_type was supplied to parse_probe_files.\n");

		    my $backg = $datacol{'background'} || 
			die("$error_head: None background was supplied to parse_probe_files.\n");

		    my $backg_type = $prob_express{$targetel_name}->{'background_type'} ||
			die("$error_head: None background_type was supplied to parse_probe_files.\n");;
		    
		    my $flags = $datacol{'flags'} || '\N';
		      
		    my $data_id = $prob_express{$targetel_name}->{'sample_id'} || '\N';
		    
		    my $metadata_id = $args_href->{'metadata_id'} ||
			die("$error_head: None metadata_id was supplied to parse_probe_files.\n");
		    
		    ## Print the output in a file

		    print $opfh "$probe_id\t$target_e_id\t$signal\t$signal_type\t$backg\t$backg_type\t$flags\t$data_id\t$metadata_id\n";
		}		
	    }
	}
	close $ipfh;
	push @probe_dbload_files, $outfile;
    }
    return @probe_dbload_files;
}

=head2 create_template_dbload

  Usage: my @template_filedata = create_template_dbload($arg_href);
  Desc: parse the template_files contained in the argument hash reference 
  Ret: a list of file names to load into the database
  Args: $arg_href, a hash reference with the arguments.
       keys=argument_type and value=hash reference with
            keys=argument_name and value=value

       exceptions are: schema, that supply an schema object
                     : loader_username, that supply a scalar
                       metadata_id, a scalar
  Side_Effects: die if something it is wrong
  Example: my @template_filedata = create_template_dbload($arg_href);

=cut

sub create_template_dbload {
    my $args_href = shift ||
	die("None argument hash reference was supplied to the create_template_dbload fuction.\n");

    ## Create the array to store the file names

    my @template_dbload_files;

    ## Get the data stored in the argument hash reference for probe_expression and get all the file names

    my %temp_express = %{ $args_href->{'template_expression'}};  

    foreach my $hyb_name (sort keys %temp_express) {

	my $file = $temp_express{$hyb_name}->{'template_expression_file'};
	
	## Create the out file

	my $path = `pwd`;
	my $basename = basename($file);

	chomp($path);

	my $outfile = $path . "/" . $basename . ".template_dbload.tab";

	open my $itfh, '<', $file || die("Sorry, I can not open the template_expression_file=$file.\n");
	
	open my $otfh, '>', $outfile || die("Sorry, I can not open the template_expression_file output=$outfile.\n");

	my $header_line = $temp_express{$hyb_name}->{'template_expression_headers'} || 
	    die("DATA ERROR: None template_expression_headers were supplied.\n");

	my @headers = split(/,/, $header_line);

	my $l = 0;

	## Parse the file

	while (<$itfh>) {
	    $l++;
	    chomp($_);
	    unless ($_ =~ m/^#/) {
		my @data = split(/\t/, $_);
	    		
		if (scalar(@data) != scalar(@headers)) {
		    die("DATA TEMPLATE ERROR: File=$file has different number of column than headers in the line=$l.\n");
		}

		## To put all the data in a hash to get the correspondences

		my %datacol;
		my $a = 0;
		foreach my $head (@headers) {
		    $datacol{$head} = $data[$a];
		    if ($data[$a] =~ m/^\d+$/ && $data[$a] == 0) {
			$datacol{$head} = '0.000';
		    }
		    $a++;
		}

		print STDERR "\t\tProcessing file:$file line $l (template_name=$data[0])                 \r";

		## Now it will replace probe_name by probe_id

		my ($template_row) = $args_href->{'schema'}
	                                       ->resultset('GeTemplate')
				       	       ->search( { template_name => $datacol{'template_name'} } );
	    
		unless (defined $template_row) {
		    warn("WARNING: template_name=$datacol{'template_name'} doesn't exist into the table: gem.ge_template.SKIPING DATA.\n")
		}
		else {
		    my $template_id = $template_row->get_column('template_id');

		    my $error_head = "TEMPLATE EXPRESS PARSING ERROR (file:$file, line:$l)";

		    ## Die if there aren't some of the not null variables

		    my $hyb_id = $temp_express{$hyb_name}->{'hybridization_id'} || 
			die("$error_head: None hybridization_id was supplied to parse_template_files.\n");
		    
		    my $signal = $datacol{'signal'} || 
			die("$error_head: None template_signal was supplied to parse_template_files.\n");

		    my $signal_type = $temp_express{$hyb_name}->{'signal_type'} ||
			die("$error_head: None template_signal_type was supplied to parse_template_files.\n");

		    my $stat = $datacol{'statistical_value'} || 
			die("$error_head: None statistical_value was supplied to parse_template_files.\n");

		    my $stat_type = $temp_express{$hyb_name}->{'stat_value_type'} ||
			die("$error_head: None stat_value_type was supplied to parse_template_files.\n");;
		    
		    my $flags = $datacol{'flags'} || '\N'; 

		    my $dataset_id = $temp_express{$hyb_name}->{'sample_id'} || '\N';
		    
		    my $metadata_id = $args_href->{'metadata_id'} ||
			die("$error_head: None metadata_id was supplied to parse_template_files.\n");
		    
		    ## Print the output in a file

		    print $otfh "$hyb_id\t$template_id\t$signal\t$signal_type\t$stat\t$stat_type\t$flags\t$dataset_id\t$metadata_id\n";
		}		
	    }
	}
	print "\n\n";
	close $itfh;
	push @template_dbload_files, $outfile;
    }
    return @template_dbload_files;
}

=head2 create_experiment_dbload

  Usage: my @experiment_filedata = create_experiment_dbload($arg_href);
  Desc: parse the experiment_files contained in the argument hash reference.
        if these files are omited, it will take the hybridization associated
        with theses experiments and it will calculate the different data for
        each template id
  Ret: a list of file names to load into the database
  Args: $arg_href, a hash reference with the arguments.
       keys=argument_type and value=hash reference with
            keys=argument_name and value=value

       exceptions are: schema, that supply an schema object
                     : loader_username, that supply a scalar
                       metadata_id, a scalar
  Side_Effects: die if something it is wrong
  Example: my @experiment_filedata = create_experiment_dbload($arg_href);

=cut

sub create_experiment_dbload {
    my $args_href = shift ||
	die("None argument hash reference was supplied to the create_template_dbload fuction.\n");

    ## Create the array to store the file names

    my @experiment_dbload_files;

    ## Get the data stored in the argument hash reference for probe_expression and get all the file names

    my %exp_express = %{ $args_href->{'experiment_expression'}};  

    foreach my $exp_name (sort keys %exp_express) {

	## Create the out file
	
	my $path = `pwd`;
	chomp($path);
	my $basename = $exp_name;
	$basename =~ s/\s+/_/g;
	$basename =~ s/\\/_/g;
	$basename =~ s/\//_/g;

	my $outfile = $path . "/" . $basename . "_template_dbload.tab";
	

	open my $otfh, '>', $outfile || die("Sorry, I can not open the template_expression_file output=$outfile.\n");

	my $file = $exp_express{$exp_name}->{'experiment_expression_file'};

	if (defined $file) {
  
	    open my $itfh, '<', $file || die("Sorry, I can not open the template_expression_file=$file.\n");
	
	    my $header_line = $exp_express{$exp_name}->{'template_expression_header'};
	    my @headers = split(/,/, $header_line);
	
	    my $l = 0;

	    ## Parse the file

	    while (<$itfh>) {
		$l++;
		chomp($_);
		unless ($_ =~ m/^#/) {
		    my @data = split(/\t/, $_);
	    
		    if (scalar(@data) != scalar(@headers)) {
			die("DATA TEMPLATE ERROR: File=$file has different number of column than headers in the line=$l.\n");
		    }

		    print STDERR "\t\tProcessing file:$file line $l (template_name=$data[0])                 \r";

		    ## To put all the data in a hash to get the correspondences

		    my %datacol;
		    my $a = 0;
		    foreach my $head (@headers) {
			$datacol{$head} = $data[$a];
			$a++;
		    }

		    ## Now it will replace probe_name by probe_id

		    my ($template_row) = $args_href->{'schema'}
	                                           ->resultset('GeTemplate')
				         	   ->search( { template_name => $datacol{'template_name'} } );
	    
		    unless (defined $template_row) {
			warn("WARNING: template_name=$datacol{template_name} doesn't exist in the table: gem.ge_template.SKIPING DATA.\n")
		    }
		    else {
			my $template_id = $template_row->get_column('template_id');

			my $error_head = "TEMPLATE EXPRESS PARSING ERROR (file:$file, line:$l)";

			## Die if there aren't some of the not null variables

			my $exp_id = $exp_express{$exp_name}->{'experiment_id'} || 
			    die("$error_head: None experiment_id was supplied to create_experiment_dbload.\n");

			my $repl_n = $datacol{'replicate_n'} || '0';
			
			my $mean = $datacol{'mean'} || '\N';
		 
			my $median = $datacol{'median'} || '\N';
			
			my $sd = $datacol{'sd'} || '\N';
			
			my $cv = $datacol{'cv'} || '\N';
		    		    
			my $dataset_id = $exp_express{$exp_name}->{'sample_id'} || '\N';
			
			my $metadata_id = $args_href->{'metadata_id'} ||
			    die("$error_head: None metadata_id was supplied to parse_template_files.\n");
		    
			## Print the output in a file

			print $otfh "$exp_id\t$template_id\t$repl_n\t$mean\t$median\t$sd\t$cv\t$dataset_id\t$metadata_id\n";
		    }		
		}
	    }
	    close $itfh;
	}
	else {

	    ## If don't exist any experiment file, it will create these files using the data from gem.ge_template_expression

	    print STDERR "\t\tSCRIPT MESSAGE: There aren't any file associated with exp_name=$exp_name.\n";
	    print STDERR "\t\tRETRIEVING DATA FROM DATABASE to CALCULATE experiment expression data.\n\n";

	    ## First, get the list of the target_ids associated with a experiment and the platform_id prased before

	    my $exp_id = $exp_express{$exp_name}->{'experiment_id'} || 
			    die("DATA ERROR ($exp_name): None experiment_id was supplied to create_experiment_dbload.\n");

	    my $platf_id = $exp_express{$exp_name}->{'platform_id'} || 
			    die("DATA ERROR ($exp_name): None platform_id was supplied to create_experiment_dbload.\n");

	    my $exp = CXGN::GEM::Experiment->new($args_href->{'schema'}, $exp_id);

	    ## Define the hash to store the expression data

	    my %exptemp;

	    my @target_list = $exp->get_target_list();
	    foreach my $target (@target_list) {
		my $target_id = $target->get_target_id();
		
		## Using target_id and platform_id, it will search the hybrdization_id
		
		my ($hyb_row) = $args_href->{'schema'}
		                        ->resultset('GeHybridization')
       				        ->search({ target_id => $target_id, platform_id => $platf_id });

		if (defined $hyb_row) {

		    my $hyb_id = $hyb_row->get_column('hybridization_id');

		    print STDERR "\t\tRetrieving data for hybridization_id=$hyb_id (associated to experiment=$exp_name).\r";
		    
		    ## With hyb_id it will get all the temnplate_expression rows from gem.ge_template_expression table

		    my @temp_rows = $args_href->{'schema'}
		                              ->resultset('GeTemplateExpression')
					      ->search({ hybridization_id => $hyb_id });

		    foreach my $temp_row (@temp_rows) {

			my %tempdata = $temp_row->get_columns();
			
			## It will store the data as a array reference of signal values associated with each template_id

			if (exists $exptemp{$tempdata{'template_id'}} ) {
			    push @{$exptemp{$tempdata{'template_id'} } }, $tempdata{'template_signal'};
			}
			else {
			    $exptemp{$tempdata{'template_id'} } = [$tempdata{'template_signal'}];
			}
		    }

		    print STDERR "\n\t\t\t..Done.\n";
		}
	    }

	    ## When all the target associated to the experiment will be gotten from the database, it will calculate
	    ## the stats values (median, mean, sd and cv)

	    print STDERR "\n";

	    ## First take the global data used in the table

	    my $dataset_id = $exp_express{$exp_name}->{'sample_id'} || '\N';
			
	    my $metadata_id = $args_href->{'metadata_id'} ||
		die("DATA ERROR ($exp_name): None metadata_id was supplied to parse_template_files.\n");

	    foreach my $temp_id (sort {$a <=> $b} keys %exptemp) {
		
		print STDERR "\t\tProcessing STATS for template_id=$temp_id for experiment=$exp_name.\r";
		
		my ($mean, $median, $std, $cv, $repl_n);
		my $template_values_aref = $exptemp{$temp_id};

		if (defined $template_values_aref) {
		    my @values = @{$template_values_aref};
		    $repl_n = scalar(@values);
		    
		    if ($repl_n > 1) {
			my $stats = Statistics::Descriptive::Full->new();
			$stats->add_data(@values);
			$mean = $stats->mean();
			$median = $stats->median();
			$std = $stats->standard_deviation();
			if ($mean > 0) {
			    $cv = $std / $mean;
			} else {
			    $cv = '\N';
			}
		    } elsif ($repl_n == 1) {
			$mean = $values[0];
			$median = $values[0];
			$std = '\N';
			$cv = '\N';
		    }
		} else {
		    $repl_n = 0;
		    $mean = '\N';
		    $median = '\N';
		    $std = '\N';
		    $cv = '\N';
		}
		
		print $otfh "$exp_id\t$temp_id\t$repl_n\t$mean\t$median\t$std\t$cv\t$dataset_id\t$metadata_id\n";
	    }
	    print STDERR "\n\t\t\t...Done\n";
	}
	push @experiment_dbload_files, $outfile;
    }
    return @experiment_dbload_files;
}

=head2 store_pipeline

  Usage: store_pipeline($args_href)
  Desc: copy the data parsed using $dbh->do functions
        use the create_xxxxx_functions to parse the files
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

    ## Third parse and create the dbload files (they need metadata_id)

    print STDERR "\t4.1) Parsing probe_expression_files ...\n";
    my @probe_expression_dbload = create_probe_dbload($args_href);
	
    my $prob_dbload_n = scalar(@probe_expression_dbload);
    print STDERR "\t\t... done ($prob_dbload_n files were parsed)\n\n";

    print STDERR "\t4.2) Parsing template_expression_files ...\n";
    my @template_expression_dbload = create_template_dbload($args_href);

    my $templ_dbload_n = scalar(@template_expression_dbload);
    print STDERR "\t\t... done ($templ_dbload_n files were parsed)\n\n";

    ## Forth load the files into the database.

    my $last_probe_expression_id = $all_last_ids_href->{'gem.ge_probe_expression_probe_expression_id_seq'};
    my $last_template_expression_id = $all_last_ids_href->{'gem.ge_template_expression_template_expression_id_seq'};

    print STDERR "\t4.3) Copying probe_expression_dbload files (if there are any).\n\n";
    foreach my $probe_dbload (@probe_expression_dbload) {
	$schema->storage()
	       ->dbh()
	       ->do("COPY gem.ge_probe_expression ( target_element_id,
                                                    probe_id, 
                                                    signal, 
                                                    signal_type, 
                                                    background, 
                                                    background_type, 
                                                    flag, 
                                                    dataset_id, 
                                                    metadata_id) 
                         FROM STDIN");

	open my $temp_probe_io, '<', $probe_dbload;
	while (<$temp_probe_io>) {
	    $schema->storage()->dbh()->pg_putline($_);
	}
	close $temp_probe_io;

	$schema->storage()->dbh()->pg_endcopy();


	my $curr_prob_expr_id = $schema->resultset('GeProbeExpression')
		                       ->search(undef, { order_by => "probe_expression_id DESC"})
		      	               ->first()
				       ->get_column('probe_expression_id');

	my $prob_load_n = $curr_prob_expr_id - $last_probe_expression_id;
	$last_probe_expression_id = $curr_prob_expr_id;

	print STDERR "\t\t$prob_load_n expression values from $probe_dbload file\n";
	print STDERR "\t\t\thas been loaded into gem.ge_probe_expression table.\n";
    }
    print STDERR "\n";
    
    print STDERR "\t4.4) Copying template_expression_dbload files (if there are any).\n\n";
    foreach my $template_dbload (@template_expression_dbload) {
	  $schema->storage()
	         ->dbh()
	         ->do("COPY gem.ge_template_expression ( hybridization_id,
                                                         template_id, 
                                                         template_signal, 
                                                         template_signal_type, 
                                                         statistical_value, 
                                                         statistical_value_type, 
                                                         flag, 
                                                         dataset_id, 
                                                         metadata_id ) 
                            FROM STDIN");
	  
	  open my $temp_template_io, '<', $template_dbload;
	  while (<$temp_template_io>) {
	      $schema->storage()->dbh()->pg_putline($_);
	  }
	  close $temp_template_io;

	  $schema->storage()->dbh()->pg_endcopy();


	  my $curr_temp_expr_id = $schema->resultset('GeTemplateExpression')
		                         ->search(undef, { order_by => "template_expression_id DESC"})
		       	                 ->first()
				         ->get_column('template_expression_id');

	  my $temp_load_n = $curr_temp_expr_id - $last_template_expression_id;
	  $last_template_expression_id = $curr_temp_expr_id;
	    
	  print STDERR "\t\t$temp_load_n expression values from $template_dbload file\n";
	  print STDERR "\t\t\thas been loaded into gem.ge_template_expression table.\n";
    }
    print STDERR "\n";

    ## After store the possible template_expression data it will work over experiment expression data
    
    my $last_experiment_expression_id = $all_last_ids_href->{'gem.ge_expression_by_experiment_expression_by_experiment_id_seq'};

    print STDERR "\t4.5) Parsing experiment_expression_files ...\n";
    
    my @experiment_expression_dbload = create_experiment_dbload($args_href);
    
    my $exper_dbload_n = scalar(@experiment_expression_dbload);
    print STDERR "\t\t... done ($exper_dbload_n files were parsed)\n\n";

    print STDERR "\t4.6) Copying template_expression_dbload files (if there are any).\n\n";
    foreach my $experiment_dbload (@experiment_expression_dbload) {
	$schema->storage()
       	       ->dbh()
	       ->do("COPY gem.ge_expression_by_experiment ( experiment_id, 
                                                            template_id, 
                                                            replicates_used, 
                                                            mean, 
                                                            median, 
                                                            standard_desviation, 
                                                            coefficient_of_variance, 
                                                            dataset_id, 
                                                            metadata_id ) 
                     FROM STDIN");

	open my $temp_exp_io, '<', $experiment_dbload;
	while (<$temp_exp_io>) {
	    $schema->storage()->dbh()->pg_putline($_);
	}
	close $temp_exp_io;

	$schema->storage()->dbh()->pg_endcopy();


	my $curr_expe_expr_id = $schema->resultset('GeExpressionByExperiment')
		                       ->search(undef, { order_by => "expression_by_experiment_id DESC"})
				       ->first()
			               ->get_column('expression_by_experiment_id');

	my $expe_load_n = $curr_expe_expr_id - $last_experiment_expression_id;
	$last_experiment_expression_id = $curr_expe_expr_id;
	
	print STDERR "\t\t$expe_load_n expression values from $experiment_dbload file\n";
	print STDERR "\t\t\thas been loaded into gem.ge_expression_by_experiment table.\n";
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
    my $argument_file = $dir . '/expression_argument_loadfile.bs';
    open my $TFH, '>', $argument_file || die("Sorry, I can not open the file: $argument_file.\n");

    print STDERR "PRINT GEM(Expression) DATA LOAD FILE option...\n";
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

## To use with a probe expression file

*DATA_TYPE:		[]		#mandatory,single	example:[probe_expression]
*TARGET_ELEMENT_NAME:   []              #mandatory,single       example:[]
*PROBE_EXPRESS_FILE:    []		#mandatory,multiple	example:[]
*PROBE_EXPRESS_HEADERS: []              #mandatory,multiple     example:[]
*DATASET_NAME:          []              #mandatory,single       example:[]
*SIGNAL_TYPE:           []              #mandatory,single       example:[]
*BACKGROUND_TYPE:       []              #mandatory,single       example:[]
//

## To use with a template expression file

*DATA_TYPE:		[]		#mandatory,single	example:[template_expression]
*HIBRIDIZATION:         []              #mandatory,single       example:[]     
*TEMPL_EXPRESS_FILE:    []		#mandatory,multiple	example:[]
*TEMPL_EXPRESS_HEADERS: []              #mandatory,multiple     example:[]
*DATASET_NAME:          []              #mandatory,single       example:[]
*SIGNAL_TYPE:           []              #mandatory,single       example:[]
*STAT_VALUE_TYPE:       []              #mandatory,single       example:[]
//

## To use with a experiment expression file

*DATA_TYPE:		[]		#mandatory,single	example:[experiment_expression]
*EXPERIMENT:            []              #mandatory,single       example:[]
*PLATFORM:              []              #mandatory,single       example:[]     
*TEMPL_EXPRESS_FILE:    []		#optional,single,note1	example:[]
*TEMPL_EXPRESS_HEADERS: []              #optional,single        example:[]
*DATASET_NAME:          []              #optional,single        example:[]
//

';

    print $TFH "$info";
    print STDERR "...done (printed argument_file with the name: $argument_file)\n\n";
    exit (1);
}




####
1;##
####
