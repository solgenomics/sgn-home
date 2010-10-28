#!/usr/bin/perl

=head1 NAME

 load_blastresults_into_sgndb.pl
 A script to load blast results for unigenes into the sgn database (version.0.1.).

=cut

=head1 SYPNOSIS

 load_blastresults_into_sgndb [-h] -H <dbhost> -D <dbname> -U <dbuser> -b <blast_result_file_m8> -a <defline_file> -d <external_dbname> 
                              [-T] [-e <e_value_filter>] [-s <score_filter>] [-i <identity_filer>] [-l <length_filter>] [-n <hits_filter>]

=head1 EXAMPLE:

 load_blastresults_into_sgndb [-h] -H localhost -D sandbox -U postgres -b potato.blastx.gb.m8 -a gb_defline.gz -d genbank -e 1e-20 [-T]
        
    
=head2 I<Flags:>

=over

=item -H

B<database_host>                database host (mandatory)

=item -D

B<database_name>                database name (mandatory)

=item -U

B<database_user>                database user (mandatory)

=item -b

B<blast_results_file>           blast results file in m8 format (mandatory)

=item -a

B<defline_file>                 defline file (it can be compressed) (mandatory)
 
=item -d

B<external_dbname>              name of the database used in the blast (example: genbank) (mandatory)

=item -e

B<e_value_filter>               filter by e_value, only will be selected hits with e_value < e_value_cutoff

=item -s

B<hit_score_filter>             filter by hit_score, only will be selected hits with hit_score > hit_score_cutoff

=item -i

B<identity_filter>              filter by identity, only will be selected hits identity_score > identity_cutoff

=item -l

B<alignment_length_filter>      filter by alignment length, only will be selected hits alignment length > align_length_cutoff

=item -n

B<number_of_hits_filter>        filter by number of hits, only will be selected this number of hits (50 by default)

=item -T

B<run as a test>                run the script as test

=item -h

B<help>                         print the help  

=back

=cut

=head1 DESCRIPTION

    This script parse blast result file in m8 format associated with an unigene build and 
    the defline file. Load the blast results into sgn.blast_annotation and sgn.blast_hits tables 
    and the deflines into the sgn.blast_defline table.

    Note about -T: Enable the test mode. It run the script as a datatabase transaction and the 
                   code as eval{} function, rollback the results at the end of the script.

    Note about deflines: A defline file can be produced using the following commands:

    fastacmd -D 1 -d my_blast_db_formated | sed -n '/^>/{s/[[:space:]]+/\t/; p;}' | sort | gzip -c > defline.gz
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

load_blastresults_into_sgndb.pl


=cut

use strict;
use warnings;

use Cwd;
use File::Basename;
use Getopt::Std;
use Term::ReadKey;

use CXGN::DB::Connection;


our ($opt_H, $opt_D, $opt_U, $opt_b, $opt_a, $opt_d, $opt_e, $opt_s, $opt_i, $opt_l, $opt_n, $opt_T, $opt_h);
getopts("H:D:U:b:a:d:e:s:i:l:n:Th");
if (!$opt_H && !$opt_D && !$opt_U && !$opt_b && !$opt_a && !$opt_d && !$opt_e && !$opt_s && !$opt_i 
    && !$opt_l && !$opt_n && !$opt_T && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
} elsif ($opt_h) {
    help();
} 

## Checking the input arguments

my $dbname = $opt_D || die("MANDATORY ARGUMENT ERROR: The -D <database_name> argument was not supplied.\n");
my $dbhost = $opt_H || die("MANDATORY ARGUMENT ERROR: The -H <database_hostname> argument was not supplied.\n"); 
my $dbuser = $opt_U || die("MANDATORY ARGUMENT ERROR: The -U <database_user> argument was not supplied.\n");
my $blastf = $opt_b || die("MANDATORY ARGUMENT ERROR: The -b <blast_result_file> argument was not supplied.\n");
my $deflif = $opt_a || die("MANDATORY ARGUMENT ERROR: The -a <defline_file> argument was not supplied.\n");
my $bla_db = $opt_d || die("MANDATORY ARGUMENT ERROR: The -d <blast_database_name> argument was not supplied.\n");
my $hits_cutoff = $opt_n || 50;

## Start the process

my $start_time = time();

print STDERR "\n\n*** STARTING BLAST RESULT LOADING PROCESS.\n";

## Connecting with the database

print STDERR "\n1) CONNECTING WITH THE DATABASE.\n";

## First, get the password as prompt

print STDERR "\n\tType password for database user=$dbuser:\n\tpswd> ";

ReadMode('noecho');
my $passw = <>;
chomp($passw);
ReadMode('normal');

print STDERR "\n\n";

## Create a new db_connection

my $dbh = CXGN::DB::Connection->new(
                                     {
					 dbhost => $dbhost,
					 dbname => $dbname,
					 dbuser => $dbuser,
					 dbpass => $passw
				     }
                                   );

## This always will run as a transaction

$dbh->{AutoCommit} = 0;  # enable transactions, if possible
$dbh->{RaiseError} = 1;

## Check if the -d supplied exists in sgn.blast tables

print STDERR "\n2) CHECKING blast_dbname.\n";

my $t = time() - $start_time;
my $f_time = secs_to_complex($t);

my $blast_db_query = "SELECT blast_target_id FROM sgn.blast_targets WHERE db_name ILIKE ?";
my $f_blast_db = "%".$bla_db."%";

my ($blast_target_id) = $dbh->selectrow_array($blast_db_query, undef, $f_blast_db);

unless (defined $blast_target_id) {
    die("ARGUMENT ERROR: -d $bla_db does not match with any name in sgn.blast_target database table.\n");
}
else {
    print STDERR "\n\tBlast database = $bla_db, exists in the db with blast_target_id=$blast_target_id.\n";
    print STDERR "\t(op. time = $f_time).\n";
}

## Before parsing the defline and add all the lines, it will select from the blast file the subject_ids 
## to select the data from the defline file

open my $b_fh, '<', $blastf || die("OPEN FILE ERROR: $blastf could not be openned (system error: $!).\n");

## It will store the filtered blast file 

print STDERR "\n3) PARSING AND FILTERING BLAST FILE.\n\n";

my $path = getcwd();
my $blast_basename = basename($blastf);
my $filter_blastf = $path . '/blastload_filter_' . $blast_basename;

open my $fb_fh, '>', $filter_blastf || die("OPEN FILE ERROR: $filter_blastf could not be openned (system error: $!).\n");

## Also it will keep the object_ids and the subject_ids in a hashes

my %obj_ids;
my %subj_ids;
my %blast_results;

my $l = 0;
while(<$b_fh>) {
    chomp($_);

    $l++;
    print STDERR "\tParsing line:$l for file $blastf              \r";
    
    my @data = split(/\t/, $_);

    my $obj_id = $data[0];
    my $subj_id = $data[1];
  
    if (exists $obj_ids{$obj_id}) {
	$obj_ids{$obj_id}++;
    }
    else {
	$obj_ids{$obj_id}++;
    }

    ## Now it will print the blast data into the filter file according the
    ## filter cutoffs

    if ($obj_ids{$obj_id} <= $hits_cutoff) {
	my $select = 1;
	if (defined $opt_e) {
	    my $e_value_cut = Math::BigFloat->new($opt_e);
	    my $e_value = Math::BigFloat->new($data[10]);

	    if ($e_value->bmcp($e_value_cut) > 0) {
		$select = 0;		
	    }
	}
	if (defined $opt_s) {
	    if ($data[11] < $opt_s) {
		$select = 0;		
	    }
	}
	if (defined $opt_i) {
	    if ($data[2] < $opt_i) {
		$select = 0;		
	    }
	}
	if (defined $opt_l) {
	    if ($data[3] < $opt_l) {
		$select = 0;		
	    }
	}
	
	## It it has passed the filter, will be selected and printed

	if($select == 1) {

	    ## Add the subj_id to the hash to be selected in the defline parsing

	    if (exists $subj_ids{$subj_id}) {
		$subj_ids{$subj_id}++;
	    }
	    else {
		$subj_ids{$subj_id}++;
	    }
	    
	    if (exists $blast_results{$obj_id}) {
		push @{$blast_results{$obj_id}}, \@data;
	    }
	    else {
		$blast_results{$obj_id} = [\@data];
	    }

	    print $fb_fh "$_\n";
	}
    }
}

print STDERR "\n\n\t$blastf blast file has been parsed and filtered.\n";
$t = time() - ($start_time + $t);
$f_time = secs_to_complex($t);
print STDERR "\t(op. time = $f_time).\n";

## Now it will process the deflines

## First, check if the defline is unzip

print STDERR "\n4) UNCOMPRESSING AND PARSING $deflif file\n\n";

if ($deflif =~ m/\.gz$/) {  ## It is compress with .gz

    my $uncompress = $deflif;
    $uncompress =~ s/\.gz$//;
    
    my $g = system("gunzip -c $deflif > $uncompress");
    unless ($g == 0) {
	die("SYSTEM ERROR: $deflif file could not be uncompress using gunzip (system error $?).\n");
    }

    $t = time() - ($start_time + $t);
    print STDERR "\t$deflif file has been uncompressed.\n";
    $f_time = secs_to_complex($t);
    print STDERR "\t(op. time = $f_time).\n";

    $deflif = $uncompress;
}    

## With the defline uncompressed and the list of subject_id it will prepare the
## file to load the data into the database

open my $def_fh, '<', $deflif ||
    die("OPEN FILE ERROR: $deflif defline file could not be openned (system error: $!).\n");

my $defline_basename = basename($deflif);
my $load_deflif = $path . '/blastload_' . $defline_basename;

open my $def_out_fh, '>', $load_deflif || die("OPEN FILE ERROR: $load_deflif could not be openned (system error: $!).\n");

my $d = 0;
while(<$def_fh>) {
    chomp($_);
    
    $d++;
    print STDERR "\tParsing line:$d for file $deflif              \r";

    my @defdata = split(/\s/, $_);
    my $id = shift(@defdata);
    $id =~ s/>//;

    my $defline = join(' ', @defdata);
        
    ## If the blast database is arabidopsis it will remove lcl

    if ($blast_basename =~ m/[ath|arabidopsis]/i) {
	$id =~ s/lcl\|//;
	$defline =~ s/\s*\|//;
	$defline =~ s/Symbols:\s*\|\s*//;
    }

    if (exists $subj_ids{$id}) {
	print $def_out_fh "$blast_target_id\t$id\t$defline\n";
	$subj_ids{$id} = $defline;
    }
}

print STDERR "\n\n\t$deflif defline file has been parsed and filtered.\n";
$t = time() - ($start_time + $t);
$f_time = secs_to_complex($t);
print STDERR "\t(op. time = $f_time).\n";

## Now it will start the work where the database will be modify.
## This script have a test option -T that will be runned as eval {};
## process. Also all the process is run in a transaction so, if 
## the changes are not committed, none will be changed except the
## database sequences. They will be set to original state for
## test or rollback options

## Get last id for all the tables that this scripts modify.

my ($l_defline_id) = 
    $dbh->selectrow_array("SELECT defline_id FROM sgn.blast_defline ORDER BY defline_id DESC LIMIT 1");
my ($l_annot_id) = 
    $dbh->selectrow_array("SELECT blast_annotation_id FROM sgn.blast_annotations ORDER BY blast_annotation_id DESC LIMIT 1");
my ($l_hit_id) = 
    $dbh->selectrow_array("SELECT blast_hit_id FROM sgn.blast_hits ORDER BY blast_hit_id DESC LIMIT 1");

if ($opt_T) {

    print STDERR "\n\n*** TEST MODE (Operations will be done under eval function).\n\n";

    eval {
	
	store_pipeline(\%obj_ids, \%subj_ids, \%blast_results);
    };

    if ($@) {
        print STDERR "\nTEST ERRORS:\n\n$@\n";
    }

    ## Finally, rollback, because it is a test and set the sequence values

    print STDERR "\nDATABASE TRANSACTION ROLLBACK...";
    $dbh->rollback();
    print STDERR "done.\n";
    

    ## And set the sequence values

    print STDERR "\nSETTING database sequences to original value...";
    $dbh->do("SELECT setval ('sgn.blast_hits_blast_hit_id_seq', $l_hit_id, true)");
    $dbh->do("SELECT setval ('sgn.blast_annotations_blast_annotation_id_seq', $l_annot_id, true)");
    $dbh->do("SELECT setval ('sgn.blast_defline_defline_id_seq', $l_defline_id, true)");
    print STDERR "done.\n\n\n";

}
else {

    print STDERR "\n\n*** LIVE MODE (Changes can be commited or rollbacked at the end of the script).\n\n";

    store_pipeline(\%obj_ids, \%subj_ids, \%blast_results);

    print STDERR "\n\n\tDo you want commit the results (yes/no, default no)?\n";
    print STDERR "\tcommit> ";

    my $ans = <>;
    chomp($ans);

    if ($ans =~ m/^y/i) {

	print STDERR "\nDATABASE TRANSACTION COMMIT...";
	$dbh->commit();
	print STDERR "done.\n";	
    }
    else {

	print STDERR "\nDATABASE TRANSACTION ROLLBACK...";
	$dbh->rollback();
	print STDERR "done.\n";

	print STDERR "\nSETTING database sequences to original value...";
	$dbh->do("SELECT setval ('sgn.blast_hits_blast_hit_id_seq', $l_hit_id, true)");
	$dbh->do("SELECT setval ('sgn.blast_annotations_blast_annotation_id_seq', $l_annot_id, true)");
	$dbh->do("SELECT setval ('sgn.blast_defline_defline_id_seq', $l_defline_id, true)");
	print STDERR "done.\n\n\n";
    }

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
  
      This script parse blast result file in m8 format associated with an unigene build and 
      the defline file. Load the blast results into sgn.blast_annotation and sgn.blast_hits tables 
      and the deflines into the sgn.blast_defline table.

      Note about -T: Enable the test mode. It run the script as a datatabase transaction and the 
                     code as eval{} function, rollback the results at the end of the script.

      Note about deflines: A defline file can be produced using the following commands:

      fastacmd -D 1 -d my_blast_db_formated | sed -n '/^>/{s/[[:space:]]+/\t/; p;}' | sort | gzip -c > defline.gz     

    
    Usage: 
  
      load_blastresults_into_sgndb [-h] -H <dbhost> -D <dbname> -U <dbuser> -b <blast_result_file_m8> -a <defline_file> 
                                   -d <external_dbname> [-e <e_value_filter>] [-s <score_filter>] [-i <identity_filer>] 
				   [-l <length_filter>] [-n <hits_filter>] [-T]

    Example: 

      load_blastresults_into_sgndb [-h] -H localhost -D sandbox -U postgres -b potato.blastx.gb.m8 -a gb_defline.gz 
                                        -d genbank -e 1e-20 [-T]   

    Flags:

      -H <database_host>                database host (mandatory)
      -D <database_name>                database name (mandatory)
      -U <database_user>                database user (mandatory)
      -b <blast_results_file>           blast results file in m8 format (mandatory)
      -a <defline_file>                 defline file (it can be compressed) (mandatory)
      -d <external_dbname>              name of the database used in the blast (example: genbank) (mandatory)
      -e <e_value_filter>               filter by e_value, only will be selected hits with e_value < e_value_cutoff
      -s <hit_score_filter>             filter by hit_score, only will be selected hits with hit_score > hit_score_cutoff
      -i <identity_filter>              filter by identity, only will be selected hits identity_score > identity_cutoff
      -l <alignment_length_filter>      filter by alignment length, only will be selected hits alignment length > align_length_cutoff
      -n <number_of_hits_filter>        filter by number of hits, only will be selected this number of hits (50 by default)
      -T <run_as_a_test>                run the script as test
      -h <help>                         print the help       

EOF
exit (1);
}



=head2 store_pipeline

  Usage: store_pipeline(\%obj_ids, \%subj_ids, \%blast_results);
  Desc: function to store data into sgn.blast_defline, 
        sgn.blast_annotation and sgn.blast_hits
  Ret: none
  Args: %obj_ids, a hash with the blast object id as keys and hits
        count as value
        %subj_ids, a hash with the blast subject id as keys and defline
        as value
        %blast_results, a hash with object_id as keys and array reference
        with array reference with blast lines as values
  Side_Effects: print status messages
  Example: store_pipeline();

=cut

sub store_pipeline {

    my $obj_ids_href = shift;
    my $subj_ids_href = shift;
    my $blast_results_href = shift;

    my %obj_ids = %{$obj_ids_href};
    my %subj_ids = %{$subj_ids_href};
    my %blast_results = %{$blast_results_href};


    ### First operation, delete the old data for this unigene build (if exists)
    ### and this target_blast_id

    print STDERR "5) DELETING OLD BLAST DATA.\n";

    ## It will get the first unigene_id and get the unigene_build_id for it

    my @unigene_ids = keys %obj_ids;
    my $first_unigene_id = $unigene_ids[0];
    
    $first_unigene_id =~ s/SGN-U//;

    my ($unigene_build_id) = $dbh->selectrow_array( 
	                                             "SELECT unigene_build_id FROM sgn.unigene WHERE unigene_id = ?", 
	                                             undef,
	                                             $first_unigene_id
	                                          );

    ## Count the number of annotations for this unigene build and this target_blast_id

    my ($annot_count) = $dbh->selectrow_array( "SELECT COUNT(*), blast_target_id FROM sgn.blast_annotations 
                                                JOIN sgn.unigene ON(sgn.blast_annotations.apply_id=sgn.unigene.unigene_id) 
                                                WHERE unigene_build_id = ? AND blast_target_id = ?
                                                GROUP BY sgn.blast_annotations.blast_target_id ",
	                                        undef, 
	                                        $unigene_build_id, $blast_target_id );

    unless (defined $annot_count) {
	$annot_count = 0;
    }

    my $t = time() - ($start_time + $t);
    my $f_time = secs_to_complex($t);
    print STDERR "\n\tThere are $annot_count unigene with annotations for $bla_db database.\n";
    print STDERR "\t(op. time = $f_time).\n";


    ## Now it will delete all the sgn.blast_annotations associated with this build_id and this 
    ## target_blast_id
    ## Also it will delete all the entries associated with sgn.blast_hits table.

    if ($annot_count > 0) {

	my $delete = "DELETE FROM sgn.blast_annotations 
                      USING sgn.unigene WHERE sgn.blast_annotations.apply_id=sgn.unigene.unigene_id 
                      AND unigene_build_id=$unigene_build_id AND blast_target_id=$blast_target_id ";

	my $rows_deleted = $dbh->do($delete) 
	    || die($dbh->errstr);
		
	print STDERR "\n\t$rows_deleted rows have been deleted from sgn.blast_annotations and sgn.blast_hits tables.\n";
	$t = time() - ($start_time + $t);
	$f_time = secs_to_complex($t);
	print STDERR "\t(op. time = $f_time).\n";
    }
    else {
	print STDERR "\n\tThere are not annotations associated to the unigene_build_id=$unigene_build_id ($bla_db).\n";
    }

    ## For sgn.blast_defline it will check if exists or not each defline
    ## If exists it will update the description
    ## If does not exists, it will add
	
    print STDERR "\n6) UPDATING/ADDING DEFLINES.\n\n";
	
    my ($ins, $upd) = (0, 0);

    my %defline_ids = ();

    foreach my $subj_id (keys %subj_ids) {
	    
	my ($defline_id) = $dbh->selectrow_array("SELECT defline_id FROM sgn.blast_defline WHERE target_db_id = ?", undef, $subj_id);

	if (defined $defline_id) {
	    
	    ## It will UPDATE
	    
	    my $rows_updated = $dbh->do( "UPDATE sgn.blast_defline SET defline=? WHERE defline_id=?",
					  undef,
				          $subj_ids{$subj_id}, $defline_id );
		
	    print STDERR "\tUPDATED target_db_id=$subj_id ($rows_updated rows have been updated, $upd in total).            \r";
	    
	    ## and add to the defline_id hash

	    $defline_ids{$subj_id} = $defline_id;

	    $upd += $rows_updated;
	}
	else {
	    
	    ## It will INSERT

	    my $row_inserted = $dbh->do( "INSERT INTO sgn.blast_defline (blast_target_id, target_db_id, defline) 
                                          VALUES (?,?,?)",
					  undef,
					  $blast_target_id, $subj_id, $subj_ids{$subj_id} );

	    ## postgres 8.1 have not RETURNING command so it will get the new defline_id by simple SELECT command

	    ($defline_id) = $dbh->selectrow_array( "SELECT defline_id FROM sgn.blast_defline WHERE target_db_id = ?", 
						    undef, 
						    $subj_id 
		                                 );

	    $defline_ids{$subj_id} = $defline_id;

	    print STDERR "\tINSERTED target_db_id=$subj_id ($row_inserted rows have been inserted, $ins in total).         \r";
	    
	    $ins += $row_inserted;
	}
	
    }

    print STDERR "\n\n\t$ins deflined have been inserted and $upd deflines have been updated.\n";	
    $t = time() - ($start_time + $t);
    $f_time = secs_to_complex($t);
    print STDERR "\t(op. time = $f_time).\n";


    ## After update the blast deflines it will store sgn.blast_annotations
    ## It will store: apply_id (unigene_id), apply_type (15), blast_target_id, n_hits (stored in obj_ids hash)
    ## hits stored (scalar@{$blast_results{$obj_id}}), last_updated ($start_time)

    ## It will take all the unigene id for this unigene build, if do not exists hits it will store n_hits=0

    print STDERR "\n7) ADDING BLAST ANNOTATIONS.\n\n";

    my @u_ids = ();
    my $query = "SELECT unigene_id FROM sgn.unigene WHERE unigene_build_id=?";
    my $sth = $dbh->prepare($query);
    $sth->execute($unigene_build_id);

    while(my ($u_id) = $sth->fetchrow_array() ) {
	print STDERR "\tGetting unigene_id=$u_id for unigene_build_id=$unigene_build_id            \r";
	push @u_ids, $u_id;
    }
	  
    my $unig_count = scalar(@u_ids);
    print STDERR "\n\n\t$unig_count unigenes have been selected.\n\n";

    ## After the insertion it will retrieve the blast_annotation_id by select and store in blast_annot hash

    my %blast_annot = ();
    
    my $bl_annot_ins = 0;

    foreach my $unigene_id (sort @u_ids) {

	if (exists $blast_results{'SGN-U' . $unigene_id}) {
		
	    my @blast_results = @{$blast_results{'SGN-U' . $unigene_id} };

	    my $hits_stored = scalar(@blast_results);
	    my $n_hits = $obj_ids{'SGN-U' . $unigene_id};

	    my ($rows_inserted) = $dbh->do("INSERT INTO sgn.blast_annotations 
                                            (apply_id, apply_type, blast_target_id, n_hits, hits_stored, last_updated) 
                                            VALUES (?,?,?,?,?,?)", 
		                            undef, 
		                            $unigene_id, 15, $blast_target_id, $n_hits, $hits_stored, $start_time);
		
	    my ($blast_annot_id) = $dbh->selectrow_array("SELECT blast_annotation_id FROM sgn.blast_annotations 
                                                          WHERE apply_id=? AND blast_target_id=? AND last_updated=?", 
		                                          undef, 
		                                          $unigene_id, $blast_target_id, $start_time );

	    $blast_annot{$unigene_id} = $blast_annot_id;
	    
	    print STDERR "\tINSERTED apply_id=$unigene_id ($rows_inserted rows have been inserted, $bl_annot_ins in total).  \r";
	    
	    $bl_annot_ins += $rows_inserted;
	}
	else {
	    
	    my ($rows_inserted) = $dbh->do("INSERT INTO sgn.blast_annotations 
                                            (apply_id, apply_type, blast_target_id, n_hits, hits_stored, last_updated) 
                                            VALUES (?,?,?,?,?,?)", 
		                            undef, 
		                            $unigene_id, 15, $blast_target_id, 0, 0, $start_time);
		
	    my ($blast_annot_id) = $dbh->selectrow_array("SELECT blast_annotation_id FROM sgn.blast_annotations 
                                                          WHERE apply_id=? AND blast_target_id=? AND last_updated=?", 
		                                          undef, 
		                                          $unigene_id, $blast_target_id, $start_time );

	    $blast_annot{$unigene_id} = $blast_annot_id;

	    print STDERR "\tINSERTED apply_id=$unigene_id ($rows_inserted rows have been inserted, $bl_annot_ins in total).  \r";
	    
	    $bl_annot_ins += $rows_inserted;
	}
    }
    
    print STDERR "\n\n\t$bl_annot_ins blast_annotations have been inserted in sgn.blast_annotations table.\n";	
    $t = time() - ($start_time + $t);
    $f_time = secs_to_complex($t);
    print STDERR "\t(op. time = $f_time).\n";

    ## Finally it will add the sgn.blast_hits data

    print STDERR "\n8) ADDING BLAST HITS.\n\n";

    my $bl_hits_ins = 0;

    foreach my $unig_id (sort keys %blast_annot) {

	my $blast_annotation_id = $blast_annot{$unig_id};
	
	if (defined $blast_results{'SGN-U' . $unig_id}) {
				 
	    my @blast_results = @{$blast_results{'SGN-U' . $unig_id} };

	    foreach my $blast_line_aref (@blast_results) {
		    
		my $target_db_id = $blast_line_aref->[1];
		my $evalue = $blast_line_aref->[10];
		my $hitscore = $blast_line_aref->[11];
		my $identity = $blast_line_aref->[2];
		my $apply_start = $blast_line_aref->[6];
		my $apply_end = $blast_line_aref->[7];
		    
		my $defline_id = $defline_ids{$target_db_id};

		my ($rows_inserted) = $dbh->do("INSERT INTO sgn.blast_hits 
                                                (blast_annotation_id, target_db_id, evalue, score, 
                                                identity_percentage, apply_start, apply_end, defline_id) 
                                                VALUES (?,?,?,?,?,?,?,?)", 
	                                        undef, 
                                                $blast_annotation_id, $target_db_id, $evalue, $hitscore, $identity, 
		                                $apply_start, $apply_end, $defline_id );

		print STDERR "\tINSERTED blast_annotation_id=$blast_annotation_id ";
		print STDERR "($rows_inserted rows have been inserted, $bl_hits_ins in total).  \r";
	    
		$bl_hits_ins += $rows_inserted;
	    }
	}
    }

    print STDERR "\n\n\t$bl_hits_ins blast_hits have been inserted in sgn.blast_hits table.\n";	
    $t = time() - ($start_time + $t);
    $f_time = secs_to_complex($t);
    print STDERR "\t(op. time = $f_time).\n";
}


=head2 secs_to_complex

  Usage: my @time = secs_to_complex($secs);
  Desc: function to convert time from secs to days, hours, min and seconds
  Ret: an array with days, hours, min and seconds
  Args: $time in seconds
  Side_Effects: print status messages
  Example: my ($day, $hour, $min, $secs) = secs_to_complex($time);

=cut



sub secs_to_complex {
    my $time = shift;

    ## It will convert the time to days and get the parts from there
    ## 1 day = 24x60x60 = 86400 secs

    my ($days, $hours, $min, $secs) = (0, 0, 0, 0);

    ## Convert to days

    my $t;

    if ($time > 86400) {
	$t = $time / 86400;
	$t =~ m/^(\d+)\.?/;
	$days += $1;	
	$time = $t - $days;
    
	## Convert the rest to hours

	$t = $time * 24;
	$t =~ m/^(\d+)\.?/;
	$hours += $1;
	$time = $t - $hours;

	## Convert the rest to min
	
	$t = $time * 60;
	$t =~ m/^(\d+)\.?/;
	$min += $1;	
	$time = $t - $min;
    
	## Convert the rest to secs
	$t = $time * 60;
	$t =~ m/^(\d+)\.?/;
	$secs += $1;
    }
    elsif ($time > 3600) {
	$t = $time / 3600;
	$t =~ m/^(\d+)\.?/;
	$hours += $1;
	$time = $t - $hours;

	## Convert the rest to min
	
	$t = $time * 60;
	$t =~ m/^(\d+)\.?/;
	$min += $1;	
	$time = $t - $min;
    
	## Convert the rest to secs
	$t = $time * 60;
	$t =~ m/^(\d+)\.?/;
	$secs += $1;
    }
    elsif ($time > 60) {

	## Convert the rest to min
	
	$t = $time / 60;
	$t =~ m/^(\d+)\.?/;
	$min += $1;	
	$time = $t - $min;
    
	## Convert the rest to secs
	$t = $time * 60;
	$t =~ m/^(\d+)\.?/;
	$secs += $1;
    }
    else {
	$secs += $time;
    }


    my $f_time = $days .'D:'. $hours .'H:'. $min .'M:'. $secs . 'S';
    return $f_time;
}






####
1;##
####
