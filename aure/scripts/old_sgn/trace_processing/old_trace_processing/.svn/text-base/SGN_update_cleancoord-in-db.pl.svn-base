=head1 NAME

 SGN_update_cleancoord-in-db.pl.
 A script to update the cleaning coordenates and cleaning tags in the SGN databases.(version.1.0.).

=cut

=head1 SYPNOSIS

SGN_update_cleancoord-in-db.pl -H <dbhost> -D <dbname> -c <cleaning_coordenates_file> -q <chimera_screening_file> -o organism_name [-r <restore file>]
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>      for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>          sandbox or cxgn etc (mandatory)
      
=item -c

B<cleaning_coord. file>   cleaning coordenates file

=item -q

B<chimera screen file>	  chimera screen file

=item -r

B<restore file>           file with old clean data
           
=item -h 

B<help>                   show the help

=back

=cut

=head1 DESCRIPTION

This script update the cleaning coordenates (sgn.qc_report.hqi_start, sgn.qc_report.hqi_legth) and the cleaning tags (sgn.est.status and sgn.est.flags). Also store the old values in a file or restore the old values from these files.
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

SGN_update_cleancoord-in-db.pl


=cut

use strict;
use Getopt::Std;
use File::Basename;
use CXGN::DB::InsertDBH;

our ($opt_H, $opt_D, $opt_c, $opt_q, $opt_o, $opt_r, $opt_h);
getopts("H:D:c:q:o:r:h");
if (!$opt_H && !$opt_D && !$opt_c && !$opt_q && !$opt_o && !$opt_r && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

our ($cleanfile, $chimerafile, $organism_name, $dbh) = validate_the_input();

our $outdir;
if ($opt_c && !$opt_q) {
    $outdir=File::Basename::dirname($cleanfile)."/update_cleaning";
} elsif ($opt_q && !$opt_c) {
    $outdir=File::Basename::dirname($chimerafile)."/update_chimera_tags";
} elsif ($opt_r) {
    $outdir=File::Basename::dirname($opt_r);
} else {
    $outdir=File::Basename::dirname($cleanfile)."/update_cleaning_and_chimera";
}
mkdir "$outdir", 0775 or warn "Cannot make $outdir directory: $!";

if ($opt_r) {
    my ($old_register_file, $count_entries_register_file, $count_updates)=restore_old_data();
    print "Cleaning data restore from the file: $opt_r\n";
    print "\nOld cleaning data for $opt_o store in the file:\n\t$old_register_file";
    print "\nOld file register has $count_entries_register_file entries.\n";
    print "\nThe restore function of the script has updated $count_updates.\n";
    commit_prompt($dbh);
} else {
    my ($old_temp_table, $old_register_file, $count_entries_register_file)=store_old_data();
    print "\nOld cleaning data for $opt_o store in the file:\n\t$old_register_file";
    print "\nOld file register has $count_entries_register_file entries.\n";
    my ($temp_clean_coord, $temp_chimera_screen)=load_temp_data();
    my $master_table_file=create_master_table($temp_clean_coord, $temp_chimera_screen);
    my $count_updates=update_clean_data($master_table_file, $old_temp_table);
    print "Updated $count_updates entries.\n";
    commit_prompt($dbh);
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
      This script update the cleaning coordenates (sgn.qc_report.hqi_start, sgn.qc_report.hqi_legth) 
      and the cleaning tags (sgn.est.status and sgn.est.flags). Also store the old values in a file, or restore the old values using these files.      
    
    Usage: 
      SGN_update_cleancoord-in-db.pl -H <dbhost> -D <dbname> -c <cleaning_coordeantes_file> -q <chimera_screen_file> 
      -o <organism_name> [-r <restore file>]

    Example:
      SGN_update_cleancoord-in-db.pl -H localhost -D sandbox -c /home/aure/workdir/GenBank_sequences/
      Nicotiana_tabacum/update_clean/nt_global_seq.clean_coord.tab -q /home/aure/workdir/GenBank_sequences/
      Nicotiana_tabacum/update_clean/nt_global_Seq.chimera_screen.tab -o 'Nicotiana tabacum'
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -c cleaning coord file  input tab file with the cleaning coordenates and cleaning tags 
      -q chimera screen file  input tab file with the results of the chimera screening
      -o organism_name        organism name to save the old cleaning coordenates (mandatory)
      -r restore file         tab file with old cleaning values
      -h this help

EOF
exit (1);
}

=head2 validate_the_input

  Usage: my @input_var=validate_the_input();
  Desc: this subroutine check if exist the mandatory input variables and count the entries in the input files. 
  Ret: four scalars. $cleanfile (cleaning coordenates file name), $chimerafile (chimera id file name), 
       $organism_name and $dbh (database conection object)
  Args: none
  Side_Effects: die the process if the input data is not right
  Example: my ($cleanfile, $chimerafile, $organism_name, $dbh)=validate_the_input();

=cut
 
sub validate_the_input {
  print "Check the variables...\n";
  if ($opt_h) {
    help ();
  }

  unless ($opt_c || $opt_q || $opt_r) {
      die ("Required argument -c <cleaning_coord.file> or/and -q <chimera_screen_file> or <restore cleandata file> were not supplied.\n");
  }
  unless ($opt_o) {
      die ("Required argument -o <organism_name> was not supplied.\n");
  }
  if ($opt_c && $opt_r) {
      die ("The option -c (update clean data) is incompatible with the -r (restore) option.\n");
  } elsif ($opt_q && $opt_r) {
      die ("The option -q (update chimera tags) is incompatible with the -r (restore) option.\n");
  }

  if ($opt_c) {
      my $clean_count_columns=count_columns($opt_c);
      if ($clean_count_columns != 5) {
	  die "Sorry the cleaning coordenates file ($opt_c) have $clean_count_columns.\nThese must have 5 columns\n";
      } else {
	  print "\tColumns number in the $opt_c file ... Ok.\n";
      }
  }
  if ($opt_q) {
      my $chimera_count_columns=count_columns($opt_q);
      if ($chimera_count_columns != 3) {
	  die "Sorry the chimera screening file ($opt_q) have $chimera_count_columns.\nThese must have 3 columns\n";
      } else {
	  print "\tColumns number in the $opt_q file ... Ok.\n";
      }
  }
  if ($opt_r) {
      my $restore_count_columns=count_columns($opt_r);
      if ($restore_count_columns !=6) {
          die "Sorry the restore cleandata file ($opt_r) have $restore_count_columns.\nThese must have 6 columns\n";
      } else {
	  print "\tColumns number in the $opt_r file ... Ok.\n";
      }
  }
  

  print "\nCount process...\n";
  if ($opt_c) {
      my $clean_count_id = `cut -f1 $opt_c | sort -u | wc -l`;
      chomp ($clean_count_id);
      print "\tId count for $opt_c file:\t$clean_count_id\n";
  }
  if ($opt_q) {
      my $chimera_count_id = `cut -f1 $opt_q | sort -u | wc -l`; 
      chomp ($chimera_count_id);
      print "\tId count for $opt_q file:\t$chimera_count_id\n";
  }

  if (!$opt_D || !$opt_H) { 
      die "\nneed -D (database name) and -H (host) options";
  }

  my $dbh = CXGN::DB::InsertDBH::connect({	dbname=>$opt_D ,
						dbhost=>$opt_H ,
					  	dbargs=>{ RaiseError => 1, PrintError => 1 }
					})
	or die "Failed to connect to database ($DBI::errstr)";

  print "\nChecking organism...\n";
  my $query ="SELECT sgn.organism.organism_id, COUNT(sgn.est.est_id), COUNT(sgn.qc_report.qc_id) FROM sgn.organism 
              JOIN sgn.library ON sgn.organism.organism_id=sgn.library.organism_id 
              JOIN sgn.clone ON sgn.library.library_id=sgn.clone.library_id 
              JOIN sgn.seqread ON sgn.clone.clone_id=sgn.seqread.clone_id 
              JOIN sgn.est ON sgn.seqread.read_id=sgn.est.read_id 
              JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id
              WHERE organism_name = ? GROUP BY sgn.organism.organism_id";
  my $sth=$dbh->prepare($query);
  $sth->execute($opt_o);
  my ($organism_id, $count_est, $count_qc)=$sth->fetchrow_array();
  if (!$organism_id) {
      die "Sorry, the organism $opt_o does not exists in the sgn.organism table in $opt_D database.\n";
  } elsif ($count_est eq 0) {
      die "Sorry, the organism $opt_o have not associated any est_id.\n";
  } else {
      print "Organism: $opt_o \n\t=> organism_id: $organism_id\n\t=> est_id associated: $count_est";
      print "\n\t=> qc_id associated: $count_qc\n";
  }
 
  return ($opt_c, $opt_q, $opt_o, $dbh);
}

=head2 count_columns

  Usage: my $columns_number=count_columns($file);
  Desc: this subroutine counts the number of columns in a file and returns the number
  Ret: A scalar, the number of columns of the file
  Args: $file (the file)
  Side_Effects: none
  Example: my $columns_number=count_columns($file);

=cut

sub count_columns {
  my $filename=shift;
  my @columns;
  my ($line, $max_value);
  open my $input, '<', "$filename" or die "Sorry, I can not open the input file: $filename.\n";
      while (<$input>) {
         chomp($line = $_);
	 @columns = split(/\t/,$line);
         my $columns_number=@columns;
         if ($columns_number > $max_value) {
	     $max_value=$columns_number;
	 }
     }
  
  close $input;
  return $max_value;
}

=head2 commit_prompt

 Usage: commit_prompt($dbh);
 Desc: ask if you want commit or rollback the changes in the database. If the answer (STDIN) is yes, commit the changes.
 Ret: none
 Args: $dbh (database conection object)
 Side_Effects: print message
 Example: commit_prompt($dbh);

=cut

sub commit_prompt {
  my ($dbh, $prompt_message, $after_message) = @_;
  unless ($prompt_message) {
    $prompt_message = "Commit?\n(yes|no, default no)> ";
  }
  print $prompt_message;
  if (<STDIN> =~ m/^y(es)/i) {
    print "Committing...";
    $dbh->commit;
    print "okay.\n";
    if ($after_message) {
      print $after_message;
    }
  } else {
    print "Rolling back...";
    $dbh->rollback;
    print "done.\n";
  }
};

#This subroutine create the file to copy the data before load in the database.
 
=head2 create_dbloadfile

  Usage: my $newfilename=create_dbloadfile($sufix);
  Desc: make a new file and change the permissions to (-rw-rw-rw-).
  Ret: the name of the new file.
  Args: $sufix (sufix for the filename, for example if find 'est' the filename will be 'est.tab'.
  Side_Effects: none
  Example: my $estload_name=create_dbloadfile('est');

=cut
 
sub create_dbloadfile {
    my $sufix = shift;
   
    my $outfile_name = $outdir ."/". $sufix . ".tab";
    print "Create $outfile_name file.\n";
    open my $outfile, '>', "$outfile_name" || die "Cannot create the file $outfile_name\n";
    chmod 0666, "$outfile_name";
    return $outfile_name;
}

=head2 store_old_data

  Usage: my ($old_temp_table, $old_store_file, $count_store_file)=store_old_data();
  Desc: this subroutine store the old cleaning data in a file (with the date)
  Ret: three scalars, the $old_temp_table (table name for the selected data), the $old_store_file (the name of the 
       file where the old data will be stored) and $count_store_file (the number of stored data)
  Args: none
  Side_Effects: none
  Example: my ($old_temptable, $old_storefile, $old_count_storefile)=store_old_data();

=cut

sub store_old_data {
    my @time=localtime();
    my $year=$time[5];
    $year += 1900;
    my $month=$time[4];
    $month += 1;
    my $time=$time[2]."-".$time[1]."-".$time[0]."_".$year."_".$month."_".$time[3];
    my $filename="OLD_cleandata_".$time."_".$opt_o.".tab";
    $filename =~ s/ /_/gi;
    my $temp_table="old_cleandata_".$time."_".$opt_o."_temp";
    $temp_table =~ s/\W//gi;

    my $store_file=create_dbloadfile($filename);

    my $query="SELECT sgn.est.est_id, sgn.est.status, sgn.est.flags, sgn.qc_report.qc_id, sgn.qc_report.hqi_start, sgn.qc_report.hqi_length INTO TEMP TABLE $temp_table FROM sgn.est JOIN sgn.qc_report ON sgn.qc_report.est_id=sgn.est.est_id JOIN sgn.seqread ON sgn.est.read_id=sgn.seqread.read_id JOIN sgn.clone ON sgn.seqread.clone_id=sgn.clone.clone_id JOIN sgn.library ON sgn.clone.library_id=sgn.library.library_id JOIN sgn.organism ON sgn.library.organism_id=sgn.organism.organism_id WHERE organism_name = ?";
    my $sth=$dbh->prepare($query);
    $sth->execute($opt_o);
    $dbh->do("COPY $temp_table TO '$store_file'");
    
    my $count_store_file=`cut -f1 $store_file | sort -u | wc -l`;
    chomp ($count_store_file);
    
    return ($temp_table, $store_file, $count_store_file);
}
    
=head2 restore_old_data

  Usage: my ($old_storefile, $count_oldstorefile)=restore_old_data();
  Desc: restore the old clean data from a old store clean data file
  Ret: Two scalars, the $old_storefile (the store clean data file) and the $count_oldstorefile 
       (entries counts for the file)
  Args: none
  Side_Effects: store the data in a file (like the store_old_data subroutine)
  Example: my ($old_storefile, $count_oldstorefile)=restore_old_data();

=cut

sub restore_old_data {
    $dbh->do("CREATE TEMP TABLE restore_cleandata (est_id int, status int, flags int, qc_id int, hqi_start int, hqi_length int)");
    $dbh->do("COPY restore_cleandata FROM '$opt_r'");
  
    my ($old_temp_table, $old_storefile, $count_oldstorefile)=store_old_data();
    print "Updating sgn.est.status and sgn.est.flags...\t";
    
    $dbh->do("UPDATE sgn.est SET 
              status = 
             (SELECT restore_cleandata.status FROM restore_cleandata WHERE restore_cleandata.est_id=sgn.est.est_id), 
              flags = 
             (SELECT restore_cleandata.flags FROM restore_cleandata WHERE restore_cleandata.est_id=sgn.est.est_id) 
              FROM restore_cleandata WHERE restore_cleandata.est_id=sgn.est.est_id");
    print "...updated\nUpdating sgn.qc_report.hqi_start and sgn.qc_report.hqi_length...\t";
    $dbh->do("UPDATE sgn.qc_report SET 
              hqi_start = 
            (SELECT restore_cleandata.hqi_start FROM restore_cleandata WHERE restore_cleandata.qc_id=sgn.qc_report.qc_id),
              hqi_length = 
            (SELECT restore_cleandata.hqi_length FROM restore_cleandata WHERE restore_cleandata.qc_id=sgn.qc_report.qc_id)
              FROM restore_cleandata WHERE restore_cleandata.qc_id=sgn.qc_report.qc_id");
    print "...updated\n";

    $dbh->do("SELECT sgn.est.est_id, sgn.qc_report.qc_id, sgn.est.status-$old_temp_table.status AS diff_status, sgn.est.flags-$old_temp_table.flags AS diff_flags, sgn.qc_report.hqi_start-$old_temp_table.hqi_start AS diff_hqi_start, sgn.qc_report.hqi_length-$old_temp_table.hqi_length AS diff_hqi_length INTO TEMP TABLE new_restoreupdates FROM sgn.est JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id JOIN $old_temp_table ON sgn.est.est_id=$old_temp_table.est_id");    
    
    my $query3="SELECT COUNT(*) FROM new_restoreupdates WHERE diff_status != 0 OR diff_flags != 0 OR diff_hqi_start != 0 OR diff_hqi_length != 0";
    my $sth=$dbh->prepare($query3);
    $sth->execute();
    my ($count_entries) = $sth->fetchrow_array();


    return ($old_storefile, $count_oldstorefile, $count_entries);
}

=head2 load_temp_data

  Usage: my @load_temp_tables=load_temp_data();
  Desc: load the input file data in the database in temp tables. 
  Ret: two scalar, the temp table names ($temp_clean_coord, $temp_chimera_screen)
  Args: none (use the $opt)
  Side_Effects: none
  Example: my ($temp_clean_coord, $temp_chimera_screen)=load_temp_data();

=cut
 
sub load_temp_data { 
    print "Load files in the temp tables...\n\n";
    my $temp_clean_coord;
    if ($opt_c) {
	$temp_clean_coord = "clean_coord";
	$dbh->do("CREATE TEMP TABLE $temp_clean_coord (preest_id varchar(250), est_id int, hqi_start int, hqi_length int, clean_code varchar(250), n_percentage real)");
	$dbh->do("COPY $temp_clean_coord (preest_id, hqi_start, hqi_length, clean_code, n_percentage) FROM '$opt_c'");
	$dbh->do("UPDATE $temp_clean_coord SET preest_id = ltrim(preest_id, 'SGN-E')");
        $dbh->do("UPDATE $temp_clean_coord SET est_id = to_number(preest_id, '99999999')");
	$dbh->do("ALTER TABLE $temp_clean_coord DROP COLUMN preest_id");
	print "\tLoad the clean file in the temp tables\n";
    }
    my $temp_chimera_screen;
    if ($opt_q) {
	$temp_chimera_screen = "chimera_screen";
	$dbh->do("CREATE TEMP TABLE $temp_chimera_screen (preest_id varchar(250), est_id int, ath_b_q5 varchar(250), ath_b_q3 varchar(250))");
	print "CREATE TEMP TABLE $temp_chimera_screen\n";
	$dbh->do("COPY $temp_chimera_screen (preest_id, ath_b_q5, ath_b_q3) FROM '$opt_q'");
	$dbh->do("UPDATE $temp_chimera_screen SET preest_id = ltrim(preest_id, 'SGN-E')");
        $dbh->do("UPDATE $temp_chimera_screen SET est_id = to_number(preest_id, '99999999')");
	$dbh->do("ALTER TABLE $temp_chimera_screen DROP COLUMN preest_id");
        print "\tLoad the chimera file in the temp tables\n";
    }
    return ($temp_clean_coord, $temp_chimera_screen);
}

=head2 create_master_table

  Usage: my $master_table_file=create_master_table($temp_clean_coord, $temp_chimerascreen);
  Desc: make the master table join the different temp tables, add new fields like flags and status, update them using
        the clean codes and copy in a file.
  Ret: A scalar, the master table file name $master_table_file.
  Args: the temp table names ($temp_clean_coord, $temp_chimerascreen)
  Side_Effects: none
  Example: my $load_master_table = create_master_table($temp_clean_coord, $temp_chimera_screen);

=cut
 
sub create_master_table {
    my $temp_clean_coord = shift;
    my $temp_chimera_screen = shift;

    my ($query);

    if ($opt_c && $opt_q) {
        $query = "SELECT $temp_clean_coord.est_id, hqi_start, hqi_length, clean_code, n_percentage, ath_b_q5, ath_b_q3 INTO TEMP TABLE load_master_table FROM $temp_clean_coord LEFT JOIN $temp_chimera_screen ON $temp_clean_coord.est_id=$temp_chimera_screen.est_id";	
    } elsif ($opt_c && !$opt_q) {
        $query = "SELECT est_id, hqi_start, hqi_length, clean_code, n_percentage INTO TEMP TABLE load_master_table FROM $temp_clean_coord";	
    } elsif (!$opt_c && $opt_q) {
	$query = "SELECT est_id, ath_b_q5, ath_b_q3 INTO TEMP TABLE load_master_table FROM $temp_chimera_screen";
    }
       
    print "\nThe constructor for the load_master_table is ...\n\n$query\n\n";
    $dbh->do("$query");
   
    print "Adding flags and status to load_master_table ...\t";
    $dbh->do("ALTER TABLE load_master_table ADD COLUMN flags int");
    $dbh->do("ALTER TABLE load_master_table ADD COLUMN status int");

    if ($opt_c) {
	$dbh->do("UPDATE load_master_table SET flags = 8, status = 1 WHERE clean_code = 'seq_remove_in_phred_trimming'");
	$dbh->do("UPDATE load_master_table SET flags = 4, status = 1 WHERE clean_code LIKE '%short%'");
	$dbh->do("UPDATE load_master_table SET flags = 20, status = 1 WHERE clean_code LIKE '%UniVec'");
	$dbh->do("UPDATE load_master_table SET flags = 20, status = 1 WHERE clean_code LIKE '%ColiBank95.lib'");
	$dbh->do("UPDATE load_master_table SET flags = 20, status = 1 WHERE clean_code LIKE '%dust'");
	$dbh->do("UPDATE load_master_table SET flags = 8, status = 1 WHERE clean_code LIKE '%low_qual'");
        $dbh->do("UPDATE load_master_table SET flags = 20, status = 1 WHERE clean_code LIKE '%trim%'");
        $dbh->do("UPDATE load_master_table SET flags = 20, status = 1 WHERE clean_code LIKE '%contam%'");
    }
    if ($opt_q) {
    $dbh->do("UPDATE load_master_table SET flags = 80, status = 1 WHERE ath_b_q5 IS NOT NULL AND ath_b_q3 IS NOT NULL");
    }

    $dbh->do("UPDATE load_master_table SET flags = 0, status = 0 WHERE flags IS NULL");
    print "... added.\n\n";

    print "Adding the qc_id to load_master_table ...\t";
    $dbh->do("ALTER TABLE load_master_table ADD COLUMN qc_id int");
    $dbh->do("UPDATE load_master_table SET qc_id = (SELECT sgn.qc_report.qc_id FROM sgn.qc_report WHERE sgn.qc_report.est_id=load_master_table.est_id)");
    print "... added.\n\n";    

 
    my $master_table_file=create_dbloadfile('master_table');
    $dbh->do("COPY load_master_table TO '$master_table_file'");
    return $master_table_file;
}


=head2 update_clean_data

  Usage: my $clean_count_updated=update_clean_data($master_table_file, $old_temp_table);
  Desc: update the clean data (flags, status, hqi_start and hqi_length) in sgn.est and sgn.qc_report
  Ret: A scalar, the count of data updated.
  Args: the master_table_file name
  Side_Effects: it do not update the sequences that have not qc_id
  Example: my $clean_count_updated=update_clean_data($master_table_file, $old_temp_table);

=cut
 
sub update_clean_data {
    my $master_table_file=shift;
    my $old_temp_table=shift;

    my $query2="SELECT COUNT(*) FROM load_master_table WHERE qc_id IS NULL";
    my $sth=$dbh->prepare($query2);
    $sth->execute();
    my ($count_qcnulls) = $sth->fetchrow_array();

    if ($count_qcnulls > 0) {
	print "There are $count_qcnulls est_id that have not associated qc_id.They will not be updated.\n";
    } else {
	print "All the est_id entries have associated qc_id values.\n";
    }
    print "Updating sgn.est.status and sgn.est.flags...\t";
    $dbh->do("UPDATE sgn.est SET 
              status = 
             (SELECT load_master_table.status FROM load_master_table WHERE load_master_table.est_id=sgn.est.est_id), 
              flags = 
             (SELECT load_master_table.flags FROM load_master_table WHERE load_master_table.est_id=sgn.est.est_id)
              FROM load_master_table
              WHERE sgn.est.est_id=load_master_table.est_id AND load_master_table.qc_id IS NOT NULL");
    print "...updated\n";
    if ($opt_c) {
	print "Updating sgn.qc_report.hqi_start and sgn.qc_report.hqi_length...\t";
    $dbh->do("UPDATE sgn.qc_report SET 
              hqi_start = 
            (SELECT load_master_table.hqi_start FROM load_master_table WHERE load_master_table.qc_id=sgn.qc_report.qc_id),
              hqi_length = 
            (SELECT load_master_table.hqi_length FROM load_master_table WHERE load_master_table.qc_id=sgn.qc_report.qc_id)
              FROM load_master_table
              WHERE sgn.qc_report.qc_id=load_master_table.qc_id");
	print "...updated\n";
    }
    
    
    $dbh->do("SELECT sgn.est.est_id, sgn.qc_report.qc_id, sgn.est.status-$old_temp_table.status AS diff_status, sgn.est.flags-$old_temp_table.flags AS diff_flags, sgn.qc_report.hqi_start-$old_temp_table.hqi_start AS diff_hqi_start, sgn.qc_report.hqi_length-$old_temp_table.hqi_length AS diff_hqi_length INTO TEMP TABLE new_updates FROM sgn.est JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id JOIN $old_temp_table ON sgn.est.est_id=$old_temp_table.est_id");    
    
    my $query3="SELECT COUNT(*) FROM new_updates WHERE diff_status != 0 OR diff_flags != 0 OR diff_hqi_start != 0 OR diff_hqi_length != 0";
    $sth=$dbh->prepare($query3);
    $sth->execute();
    my ($count_entries) = $sth->fetchrow_array();

    $dbh->do("COPY load_master_table TO '$master_table_file'");
    return $count_entries;
}
