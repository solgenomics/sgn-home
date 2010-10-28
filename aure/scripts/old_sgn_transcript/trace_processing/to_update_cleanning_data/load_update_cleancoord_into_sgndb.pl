#!/usr/bin/perl

=head1 NAME

 load_update_cleancoord_into_sgndb.pl.
 A script to update the cleaning coordenates and cleaning tags in the SGN databases.(version.2.0.).

=cut

=head1 SYPNOSIS

load_update_cleancoord_into_sgndb.pl -H <dbhost> -D <dbname> [-c <cleaning_coordenates_file>] [-q <chimera_screening_file>] [-u <chimera_screen_from_unigene_assebly_file>] [-m <manually_censored_file>] [-e <clean_code_equivalences_file>] -o organism_name [-r <restore file>]
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>      for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>          sandbox or cxgn etc (mandatory)
      
=item -c

B<cleaning_coord. file>   cleaning coordenates file

=item -q

B<ath chimera screen>	  chimera screen file from chimera screening using Arabidopsis thaliana genes as dataset

=item -u

B<unigene chimera screen> chimera screen file from the chimera discovered during the unigene assembly (one column)

=item -m

B<manually censored>      id manually censored file (one column)

=item -r

B<restore file>           file with old clean data
           
=item -h 

B<help>                   show the help

=back

=cut

=head1 DESCRIPTION

This script update the cleaning coordenates (sgn.qc_report.hqi_start, sgn.qc_report.hqi_legth) and the cleaning tags (sgn.est.status and sgn.est.flags). Also store the old values in a file or restore the old values from these files.

This script change the clean_code before try to match with the default clean codes equivalences, so means that UniVec change to univec and High Expected Error to high_expected_error. This means that by default a clean code like 'ColiBank95.lib' will match with 'coli' and will have flags=32, but something like 'rRNA_contamination' will match with 'rrna' and 'contmination', in alphabetical order. The flags and status, will be the last match (in this case flag=64). For better equivalences you can load them using -e <equivalences_file> parameter.      

 * cleaning_tags and its conversion into SGN flags (these values are the conversion to decimal of the hexadecimal tags from EST.pm module):
        flags = 2    => <Possibly chimeric (vector parsing triggered)>                  => 'trim' 
	                                                                                => 'vector'
											=> 'adaptor'
	flags = 4    => <Insert too short>                                              => 'short'
	flags = 8    => <High expected error (low base calling quality values overall)> => 'low_qual' 
                                                                                        => 'high_expected_error' 
                                                                                        => 'seq_remove_in_phred'
        flags = 16   => <Low complexity>                                                => 'dust'
                                                                                        => 'nseg'
                                                                                        => 'trf'
	flags = 32   => <E.coli or cloning host contamination>                          => 'coli'
	                                                                                => 'contamination'
											=> 'host'
        flags = 64   => <rRNA contamination>                                            => 'rrna'
	flags = 128  => <Possibly chimeric (arabidopsis screen)>                        => see below
        flags = 256  => <Possibly chimeric (internal screen during unigene assembly)>   => see below
        flags = 512  => <Manually censored (reason may not be recorded)>                => see below 
        flags = 1024 => <plastid contamination>                                         => 'plastid';

 * Also: 
        flags = 128 => <Possibly chimeric (arabidopsis screen)> IF id is into the chimera screen file
        flags = 256 => <Possibly chimeric (internal screen during unigene assembly)> IF appears into chimeric unigene file
        flags = 512 => <Manually censored (reason may not be recorded)> IF appears into manually censored file     

 * All the tags will be status = 0;
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

load_update_cleancoord_into_sgndb.pl


=cut

use strict;
use warnings;

use Getopt::Std;
use File::Basename;
use CXGN::DB::InsertDBH;

our ($opt_H, $opt_D, $opt_c, $opt_q, $opt_o, $opt_r, $opt_u, $opt_m, $opt_e, $opt_h);
getopts("H:D:c:q:o:r:u:m:e:h");
if (!$opt_H && !$opt_D && !$opt_c && !$opt_q && !$opt_o && !$opt_r && !$opt_e && !$opt_u && !$opt_m && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

our ($dbh) = validate_the_input();

our $outdir;
if ($opt_c && !$opt_q) {
    $outdir=File::Basename::dirname($opt_c)."/update_cleaning";
} elsif ($opt_q && !$opt_c) {
    $outdir=File::Basename::dirname($opt_q)."/update_chimera_tags";
} elsif ($opt_u && !$opt_c) {
    $outdir=File::Basename::dirname($opt_q)."/update_chimera_tags_from_unigene_assembly";
} elsif ($opt_m && !$opt_c) {
    $outdir=File::Basename::dirname($opt_q)."/update_manually_censored_tags";
} elsif ($opt_r) {
    $outdir=File::Basename::dirname($opt_r);
} else {
    $outdir=File::Basename::dirname($opt_c)."/update_cleaning_and_chimera";
}
mkdir "$outdir", 0775 or warn "Cannot make $outdir directory: $!";
print "\n";
if ($opt_r) {
    my ($old_register_file, $count_entries_register_file, $count_updates_ref)=restore_old_data();
    print "Cleaning data restore from the file: $opt_r\n";
    print "\nOld cleaning data for $opt_o store in the file:\n\t$old_register_file";
    print "\nOld file register has $count_entries_register_file entries.\n";
    my @count_updates=@$count_updates_ref;
    print "\n---------------------------------\n";
    print "REPORT: Restore function.";
    print "\n---------------------------------\n";
    print "Est: $count_updates[0]\n";
    print "status updated in the sgn.est table: $count_updates[1]\n";
    print "flags updated in the sgn.est table: $count_updates[2]\n";
    print "hqi_start updated in the sgn.qc_report table: $count_updates[3]\n";
    print "hqi_length updated in the sgn.qc_report table: $count_updates[4]\n";
    print "---------------------------------\n\n";
    my $query = "SELECT flags, status, count(est_id) FROM sgn.est 
                                                     JOIN sgn.seqread USING(read_id) 
                                                     JOIN sgn.clone USING(clone_id) 
                                                     JOIN sgn.library USING(library_id) 
                                                     JOIN sgn.organism USING(organism_id) 
                                                     WHERE organism_name = ? 
                                                     GROUP BY flags, status";
    my $sth = $dbh->prepare($query);
    $sth->execute($opt_o);
    while (my @data = $sth->fetchrow_array()) {
	print "|\tflags:\t$data[0]\t|\tstatus:\t$data[1]\t|\tcount(est_id):\t$data[2]\t|\t\n";
    }
    print "---------------------------------\n\n";
    commit_prompt($dbh);
} else {
    print "STORING THE OLD CLEAN DATA for:\n\t- database:$opt_D\n\t- organism:$opt_o\n";
    my ($old_temp_table, $old_register_file, $count_entries_register_file)=store_old_data();
    print "\n\tOld cleaning data for $opt_o store in the file:\n\t$old_register_file";
    print "\n\tOld file register has $count_entries_register_file entries.\n\n";
    print "LOADING DATA INTO TEMP TABLE...\n";
    load_temp_data();
    my $master_table_file=create_master_table();
    print "UPDATING (as transaction) SGN TABLES (sgn.est and sgn.qc_report)\n";
    my @count_updates=update_clean_data($master_table_file, $old_temp_table);
    print "\n---------------------------------\n";
    print "REPORT: Update function.";
    print "\n---------------------------------\n";
    print "Est: $count_updates[0]\n";
    print "status updated in the sgn.est table: $count_updates[1]\n";
    print "flags updated in the sgn.est table: $count_updates[2]\n";
    print "hqi_start updated in the sgn.qc_report table: $count_updates[3]\n";
    print "hqi_length updated in the sgn.qc_report table: $count_updates[4]\n";
    print "\nnew qc entries added to sgn.qc_report table: $count_updates[5]\n";
    print "\nsequences with start position before the old start: $count_updates[6]\n";
    print "sequences with the start position after the old start: $count_updates[7]\n";
    print "sequences shorter than the old ones: $count_updates[8]\n";
    print "sequences longer than the old ones: $count_updates[9]\n";
    print "---------------------------------\n\n";
    my $query = "SELECT flags, status, count(est_id) FROM sgn.est 
                                                     JOIN sgn.seqread USING(read_id) 
                                                     JOIN sgn.clone USING(clone_id) 
                                                     JOIN sgn.library USING(library_id) 
                                                     JOIN sgn.organism USING(organism_id) 
                                                     WHERE organism_name = ? 
                                                     GROUP BY flags, status";
    my $sth = $dbh->prepare($query);
    $sth->execute($opt_o);
    while (my @data = $sth->fetchrow_array()) {
	print "|\tflags:\t$data[0]\t|\tstatus:\t$data[1]\t|\tcount(est_id):\t$data[2]\t|\t\n";
    }
    print "---------------------------------\n\n";
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

      This script change the clean_code before try to match with the default clean codes equivalences, so means that UniVec change to univec and High Expected Error to high_expected_error. This means that by default a clean code like 'ColiBank95.lib' will match with 'coli' and will have flags=32, but something like 'rRNA_contamination' will match with 'rrna' and 'contmination', in alphabetical order. The flags and status, will be the last match (in this case flag=64). For better equivalences you can load them using -e <equivalences_file> parameter.      

      * cleaning_tags and its conversion into SGN flags:
        flags = 2    => <Possibly chimeric (vector parsing triggered)>                  => 'trim' 
	                                                                                => 'vector'
											=> 'adaptor'
	flags = 4    => <Insert too short>                                              => 'short'
	flags = 8    => <High expected error (low base calling quality values overall)> => 'low_qual' 
                                                                                        => 'high_expected_error' 
                                                                                        => 'seq_remove_in_phred'
											=> 'removed'
        flags = 16   => <Low complexity>                                                => 'dust'
	                                                                                => 'nseg'
											=> 'trf'
	flags = 32   => <E.coli or cloning host contamination>                          => 'coli'
	                                                                                => 'contamination'
											=> 'host'
        flags = 64   => <rRNA contamination>                                            => 'rrna'
	flags = 128  => <Possibly chimeric (arabidopsis screen)>                        => see below
        flags = 256  => <Possibly chimeric (internal screen during unigene assembly)>   => see below
        flags = 512  => <Manually censored (reason may not be recorded)>                => see below 
        flags = 1024 => <plastid contamination>                                         => 'plastid';

      * Also: 
        flags = 128 => <Possibly chimeric (arabidopsis screen)> IF id is into the chimera screen file
        flags = 256 => <Possibly chimeric (internal screen during unigene assembly)> IF appears into chimeric unigene file
        flags = 512 => <Manually censored (reason may not be recorded)> IF appears into manually censored file     

      * All the tags will be status = 0;

    Usage: 
      load_update_cleancoord_into_sgndb.pl -H <dbhost> -D <dbname> -c <cleaning_coordeantes_file> -q <chimera_screen_file> 
      -o <organism_name> [-r <restore file>]

    Example:
      load_update_cleancoord_into_sgndb.pl -H localhost -D sandbox -c /home/aure/workdir/GenBank_sequences/
      Nicotiana_tabacum/update_clean/nt_global_seq.clean_coord.tab -q /home/aure/workdir/GenBank_sequences/
      Nicotiana_tabacum/update_clean/nt_global_Seq.chimera_screen.tab -o 'Nicotiana tabacum'
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -c cleaning coord file  input tab file with the cleaning coordenates and cleaning tags. There are two options:
                               - seqclean format (*.cln) with 7 columns (id|N%|start(1)|end|seq_lenght|clean_tag|trim_tag)
			       - tab format (*.tab) with 4 columns (id|start(0)|length|clean_tag
      -q at chimera screen    input tab file with the results of the chimera screening in arabidopsis (three column)
      -u unig chimera screen  input tab file with the results of the chimera screening in unigene assembly (one column)
      -m manually censored    input tab file with a list of manually censored ids (one column)
      -o organism_name        organism name to save the old cleaning coordenates (mandatory)
      -r restore file         tab file with old cleaning values
      -e equivalences file    tab file with three columns (clean_code|flags|status)
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
  print "\nCHECKING VARIABLES...\n";
  if ($opt_h) {
    help ();
  }

  unless ($opt_c || $opt_q || $opt_u || $opt_m || $opt_r) {
      die ("Required argument -c <cleaning_coord_file> or/and -q <chimera_screen_file> or/and -u <chimera_screen_in_unigene_assembly> or/and -m <manually_censored_list> or <restore cleandata file> were not supplied.\n");
  }
  unless ($opt_o) {
      die ("Required argument -o <organism_name> was not supplied.\n");
  }
  
  if ($opt_c && $opt_r) {
      die ("The option -c (update clean data) is incompatible with the -r (restore) option.\n");
  } elsif ($opt_q && $opt_r) {
      die ("The option -q (update chimera tags from datasets screening) is incompatible with the -r (restore) option.\n");
  } elsif ($opt_u && $opt_r) { 
      die ("The option -u (update chimera tags from unigene assembly) is incompatible with the -r (restore) option.\n");
  } elsif ($opt_m && $opt_r) {
      die ("The option -m (manually censored list) is incompatible with the -r (restore) option.\n");
  }

  if ($opt_c) {
      my @parse_opt_c=split(/\./, $opt_c);
      my $fileC_extension=pop(@parse_opt_c);
      unless ($fileC_extension eq 'tab' || 'cln') {
	  die "Sorry, the -c <cleaning_coord_file> should have .tab (4 columns) or .cln (seqclean format) extensions\n";
      }
      my $clean_count_columns=count_columns($opt_c);      
      if ($clean_count_columns != 4 && $fileC_extension eq 'tab') {
	  print "Sorry the cleaning coordenates file ($opt_c) with have $clean_count_columns.\n";
          die "The -c files with .tab extension must have 4 columns\n";
      } elsif ($clean_count_columns != 7 && $fileC_extension eq 'cln') {
	  print "Sorry the cleaning coordenates file ($opt_c) with have $clean_count_columns.\n";
	  die "The -c files with .cln extension must have 7 columns\n";
      } else {
	  print "\tColumns number for -c <clean_coordenates_file> ... Ok.\n";
      }
  }
  if ($opt_q) {
      my $chimera_count_columns=count_columns($opt_q);
      if ($chimera_count_columns != 3) {
	  die "Sorry the chimera screening file ($opt_q) have $chimera_count_columns.\nThese must have 3 columns\n";
      } else {
	  print "\tColumns number in the -q <screen_chimera_file_using_arabidopsis_dataset> ... Ok.\n";
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
  

  print "\n\tColumn count process...\n";
  if ($opt_c) {
      my $clean_count_id = `cut -f1 $opt_c | sort -u | wc -l`;
      chomp ($clean_count_id);
      print "\tId count for -c <clean_coordenates_file>:\t$clean_count_id\n";
  }
  if ($opt_q) {
      my $chimera_count_id = `cut -f1 $opt_q | sort -u | wc -l`; 
      chomp ($chimera_count_id);
      print "\tId count for -q <screen_chimera_file_using_arabidopsis_dataset>:\t$chimera_count_id\n";
  }

  if (!$opt_D || !$opt_H) { 
      die "\nneed -D (database name) and -H (host) options";
  }
  print "\n";
  my $dbh = CXGN::DB::InsertDBH::connect({	dbname=>$opt_D ,
						dbhost=>$opt_H ,
					  	dbargs=>{ RaiseError => 1, PrintError => 1 }
					})
	or die "Failed to connect to database ($DBI::errstr)";

  print "\nCHECKING ORGANISM...\n";
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
 
  return ($dbh);
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
  my $max_value=0;
  open my $input, '<', "$filename" or die "Sorry, I can not open the input file: $filename.\n";
      while (<$input>) {
         chomp($_);
	 @columns = split(/\t/,$_);
         my $columns_number=scalar(@columns);
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
    $sufix =~ s/\.tab$//;
    $sufix =~ s/_$//;
    my $outfile_name = $outdir ."/". $sufix . ".tab";
    #print "Create $outfile_name file.\n";
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
    my $name=$opt_o;
    $name =~ s/\(.+\)//;
    my $filename="OLD_cleandata_".$time."_".$name.".tab";
    $filename =~ s/ /_/gi;
    my $temp_table="old_cleandata_".$time."_".$name."_temp";
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

  Usage: my ($old_storefile, $count_oldstorefile, $count_update_arrayref)=restore_old_data();
  Desc: restore the old clean data from a old store clean data file
  Ret: Three scalars, the $old_storefile (the store clean data file), the $count_oldstorefile 
       (entries counts for the file) and $count_update_arrayref, the array reference of the different 
       updates counts 
  Args: none
  Side_Effects: store the data in a file (like the store_old_data subroutine)
  Example: my ($old_storefile, $count_oldstorefile, $count_update_arrayref)=restore_old_data();

=cut

sub restore_old_data {
    $dbh->do("CREATE TEMP TABLE restore_cleandata (est_id int, status int, flags int, qc_id int, 
              hqi_start int, hqi_length int)");
    $dbh->do("COPY restore_cleandata FROM '$opt_r'");
    my $count_table_q="SELECT COUNT(*) FROM restore_cleandata";
    my $sth=$dbh->prepare($count_table_q);
    $sth->execute();
    my ($count_entries) = $sth->fetchrow_array();
    print "There are $count_entries entries in the restore_cleandata table.\n\n";
  
    my ($old_temp_table, $old_storefile, $count_oldstorefile)=store_old_data();

    my ($c_est, $update_status_count, $update_flags_count, $update_hqistart_count, $update_hqilength_count)=(0,0,0,0,0);
    my ($update_est, $update_qc);

    my $query="SELECT est_id, status, flags, qc_id, hqi_start, hqi_length FROM restore_cleandata";
    $sth=$dbh->prepare($query);
    $sth->execute();
    while (my ($est_id, $status, $flags, $qc_id, $hqi_start, $hqi_length)=$sth->fetchrow_array()) {
        my $predata_est="SELECT status, flags FROM sgn.est WHERE est_id = ?";
        my $sth2=$dbh->prepare($predata_est);
        $sth2->execute($est_id);
        my ($prestatus, $preflags)=$sth2->fetchrow_array();
        $c_est++;
        print "\n=> Updating est: SGN-E$est_id ($c_est of $count_entries entries).\n";
        if ($prestatus eq $status) {
	    print "\tFor SGN-E$est_id, old status ($prestatus) is the same than new ($status)\n";
	} else {
            $update_est="UPDATE sgn.est SET status = ? WHERE est_id = ?";
            $sth2=$dbh->prepare($update_est);
            $sth2->execute($status, $est_id);
            print "\tUpdate SGN-E$est_id, from status: $prestatus to status: $status\n";
            $update_status_count++;
	}
        if ($preflags eq $flags) {
            print "\tFor SGN-E$est_id, old flags ($preflags) is the same than new ($flags)\n";
	} else {
            $update_est="UPDATE sgn.est SET flags = ? WHERE est_id = ?";
            $sth2=$dbh->prepare($update_est);
            $sth2->execute($flags, $est_id);
            print "\tUpdate SGN-E$est_id, from flags: $preflags to flags: $flags\n";
            $update_flags_count++;
	}

        my $predata_qc="SELECT hqi_start, hqi_length FROM sgn.qc_report WHERE qc_id = ?";
        $sth2=$dbh->prepare($predata_qc);
        $sth2->execute($qc_id);
        my ($prehqi_start, $prehqi_length)=$sth2->fetchrow_array();
        if ($prehqi_start eq $hqi_start) {
	    print "\tFor SGN-E$est_id, old hqi_start ($prehqi_start) is the same than new ($hqi_start)\n";
	} else {
            $update_qc="UPDATE sgn.qc_report SET hqi_start = ? WHERE qc_id = ?";
            $sth2=$dbh->prepare($update_qc);
            $sth2->execute($hqi_start, $qc_id);
            print "\tUpdate SGN-E$est_id, from hqi_start: $prehqi_start to hqi_start: $hqi_start\n";
            $update_hqistart_count++;
	}
        if ($prehqi_length eq $hqi_length) {
            print "\tFor SGN-E$est_id, old hqi_length ($prehqi_length) is the same than new ($hqi_length)\n";
	} else {
            $update_qc="UPDATE sgn.qc_report SET hqi_length = ? WHERE qc_id = ?";
            $sth2=$dbh->prepare($update_qc);
            $sth2->execute($hqi_length, $qc_id);
            print "\tUpdate SGN-E$est_id, from hqi_length: $prehqi_length to hqi_length: $hqi_length\n";
            $update_hqilength_count++;
	}
    }
   
    my @countupdates=($c_est, $update_status_count, $update_flags_count, $update_hqistart_count, $update_hqilength_count);


    return ($old_storefile, $count_oldstorefile, \@countupdates);
}

=head2 load_temp_data

  Usage: load_temp_data();
  Desc: load the input file data in the database in temp tables. 
  Ret: two scalar, the temp table names ($temp_clean_coord, $temp_chimera_screen)
  Args: none (use our variables and the temp table always have the same name load_master_table)
  Side_Effects: none
  Example: load_temp_data();

=cut
 
sub load_temp_data { 
  
  $dbh->do("CREATE TEMP TABLE load_master_table (est_id int, hqi_start int, hqi_length int, clean_code varchar(250))");
  print "\tCREATE TEMP TABLE load_master_table\n";


  if ($opt_c) {
      my $c=0;
      my $n_c=`cut -f1 $opt_c | wc -l`;
      chomp($n_c);
      my @fileC_parsename=split(/\./, $opt_c);
      my $fileC_extension=pop(@fileC_parsename);
      if ($fileC_extension eq 'tab') {
	  open my $fileC_fh, '<', $opt_c || die "I can not open the file -c $opt_c.\n";
	  while (<$fileC_fh>) {
	      chomp($_);
	      my @variables=split(/\t/, $_);
	      my $id=$variables[0];
	      $id =~ s/SGN-E//;
	      my $start=$variables[1];
              my $length=$variables[2];
	      my $clean_tag = lc $variables[3] || 'without_clean_code';
	      $clean_tag =~ s/^\s+//g;
	      $clean_tag =~ s/\s+$//g;
	      $clean_tag =~ s/\s+/_/g;
              if ($start == 0 && $length == 0) {
		  $length=1;
	      }

	      my $query="INSERT INTO load_master_table VALUES (?,?,?,?)";
	      my $sth=$dbh->prepare($query);
	      $sth->insert($id, $start, $length, $clean_tag);
	      $c++;
	      print "\tLoading line $c of $n_c for -c <clean_coordenates_file>\r";
	  }
	  close $fileC_fh;
      } elsif ($fileC_extension eq 'cln') {	  
          open my $fileC_fh, '<', $opt_c || die "I can not open the file -c $opt_c.\n";
	  while (<$fileC_fh>) {
	      chomp($_);
	      my @variables=split(/\t/, $_);
	      my $id=$variables[0];
	      $id =~ s/SGN-E//;
	      my $n_percent=$variables[1];
	      my $begin=$variables[2];
              my $end=$variables[3];
	      my $raw_length=$variables[4];
	      my $clean_tag = lc $variables[5] || 'without_clean_code';
	      $clean_tag =~ s/^\s+//g;
	      $clean_tag =~ s/\s+$//g;
	      $clean_tag =~ s/\s+/_/g;
	      my $trim_tag=$variables[6];
	      if ($begin  == 0 && $end > 0) {
		  print "The id:$id has start=$begin and end=$end. It looks like other format (0/something).\n";
		  print "Do you want continue? (yes|no) (0 value will change to -1)\n";
		  if (<STDIN> =~ m/n(o)/i) {
		      die "Stop process by probably wrong format.\n";
		  }
	      } elsif ($begin == 0 && $end == 0) {
		  $begin=1;
		  $end=1;
	      }
	      my $start=$begin-1;
	      my $length=$end-$begin+1;

	      my $query="INSERT INTO load_master_table VALUES (?,?,?,?)";
	      my $sth=$dbh->prepare($query);
	      $sth->execute($id, $start, $length, $clean_tag);
	      $c++;
	      print "\tLoading line $c of $n_c -c <clean_coordenates_file>\r";
	  }
	  close $fileC_fh;
      }
      print "\n\n";
  }
			  
  my ($n, $up, $in);
  if ($opt_q) {
      ($n,$up,$in)=load_single_entries_in_temp_table($opt_q, 'ath_chimera_screen_positive');     
      print "\n\tLoaded the chimera file in the load_master_table ($n entries, $up updates and $in inserts).\n";
  }
			   
  if ($opt_u) {
      ($n,$up,$in)=load_single_entries_in_temp_table($opt_u,'unigene_assembly_chimera_screen_positive');
      print "\n\tLoaded the chimera from unigene file in the load_master_table ($n entries, $up updates and $in inserts).\n";
  }

  if ($opt_m) {
      ($n,$up,$in)=load_single_entries_in_temp_table($opt_m,'manually_censored');
      print "\n\tLoaded the manually censored file in the load_master_table ($n entries, $up updates and $in inserts).\n";
  }
}

=head2 create_master_table

  Usage: my $master_table_file=create_master_table();
  Desc: make the master table join the different temp tables, add new fields like flags and status, update them using
        the clean codes and copy in a file.
  Ret: A scalar, the master table file name $master_table_file.
  Args: none (use our variables)
  Side_Effects: none
  Example: my $load_master_table = create_master_table();

=cut
 
sub create_master_table {

  my $equiv_file=$opt_e;    
  my %clean_code_clusters;
  my ($equiv_flags_href, $equiv_status_href)=flags_status_equivalences($equiv_file);  
  my %equiv_flags=%$equiv_flags_href;
  my %equiv_status=%$equiv_status_href;
  my @flags_codes= sort keys %equiv_flags;
  my @status_codes= sort keys %equiv_status;


  print "\n\tAdding flags and status to load_master_table ...\t";
  $dbh->do("ALTER TABLE load_master_table ADD COLUMN flags int");
  $dbh->do("ALTER TABLE load_master_table ADD COLUMN status int");
  
  my $query="SELECT COUNT(est_id), clean_code FROM load_master_table GROUP BY clean_code ORDER BY clean_code";
  my $sth=$dbh->prepare($query);
  $sth->execute();
  print "\n\n\tClean code cluster and est_id count:\n";
  while (my ($est_id_count,$clean_code_cluster)=$sth->fetchrow_array) {
      $clean_code_clusters{$clean_code_cluster}=$est_id_count;
      print "\t\t$est_id_count est_ids with clean code:\t$clean_code_cluster\n";
  }
  print "\n";
  my @clean_code_list = keys %clean_code_clusters;

  foreach my $clean_code (@clean_code_list) {
      my ($all_flags,$all_status)=(0,0);
      if ($clean_code ne 'without_clean_code') {
	  my @codes=split(/;/, $clean_code);
	  foreach my $precode (@codes) {
	      my $code=lc $precode;
	      my ($flag,$status)=(-1,-1);
	      foreach my $flag_kw (@flags_codes) {
		  if ($code =~ m/$flag_kw/) {
		      $flag=$equiv_flags{$flag_kw};
		  }
	      }
	      if ($flag < 0) {
		  print "\nSorry, the clean code: $precode (match_form:$code) is not into the equivalences status\n";
		  $flag=2048;
	      }
	      foreach my $status_kw (@status_codes) {
		  
		  if ($code =~ m/$status_kw/) {
		      $status=$equiv_status{$status_kw};
		  }
	      }
	      if ($status < 0) {
		  print "\nSorry, the clean code: $precode (match_form:$code) is not into the equivalences status\n";
		  $status=2048;
	      }

	      $all_flags += $flag;
	      $all_status += $status;
	  }
      }
      print "\tThe clean_code $clean_code will have flags:$all_flags and status:$all_status\n";
      my $update="UPDATE load_master_table SET flags=?, status=? WHERE clean_code=?";
      $sth=$dbh->prepare($update);
      $sth->execute($all_flags, $all_status, $clean_code);
      

  }

  print "\n\tAdding the qc_id to load_master_table ...\t";
  $dbh->do("ALTER TABLE load_master_table ADD COLUMN qc_id int");
  $dbh->do("UPDATE load_master_table SET qc_id = (SELECT sgn.qc_report.qc_id FROM sgn.qc_report WHERE sgn.qc_report.est_id=load_master_table.est_id)");
  print "... added.\n\n";    

 
  my $master_table_file=create_dbloadfile('master_table');
 
  $dbh->do("COPY load_master_table TO '$master_table_file'");
  return $master_table_file;
}


=head2 update_clean_data

  Usage: my $clean_count_updated_arrayref=update_clean_data($master_table_file, $old_temp_table);
  Desc: update the clean data (flags, status, hqi_start and hqi_length) in sgn.est and sgn.qc_report
  Ret: A scalar, the array reference of the count of data updated.
  Args: the master_table_file name
  Side_Effects: it do not update the sequences that have not qc_id
  Example: my $clean_count_updated_arrayref=update_clean_data($master_table_file, $old_temp_table);

=cut
 
sub update_clean_data {
    my $master_table_file=shift;
    my $old_temp_table=shift;

    my $count1="SELECT COUNT(*) FROM load_master_table";
    my $sth=$dbh->prepare($count1);
    $sth->execute();
    my ($count_entries) = $sth->fetchrow_array();

    my $count2="SELECT COUNT(*) FROM load_master_table WHERE qc_id IS NULL";
    $sth=$dbh->prepare($count2);
    $sth->execute();
    my ($count_qcnulls) = $sth->fetchrow_array();

    my $opt_add='no';

    if ($count_qcnulls > 0) {
	print "There are $count_qcnulls est_id (of $count_entries) that have not associated qc_id.\n";
        print "Do you want create a new qc entry for these est_id\'s (yes|no, default no):\n";
        $opt_add=<STDIN>;
        chomp($opt_add);
    } else {
	print "All the est_id entries ($count_entries) have associated qc_id values.\n";
    }
    my ($c_est, $update_status_count, $update_flags_count, $update_hqistart_count, $update_hqilength_count, 
        $new_qc_entries_added, $start_before, $start_after, $shorter_seq, $longer_seq)=(0,0,0,0,0,0,0,0,0,0);
    my ($update_est, $update_qc);

    if (!$opt_c) {
	my $query = "SELECT est_id, status, flags FROM load_master_table";
        $sth->$dbh->prepare($query);
        $sth->execute();
        while (my ($est_id, $status, $flags)=$sth->fetchrow_array()) {
	    my $predata_est="SELECT status, flags FROM sgn.est WHERE est_id = ?";
            my $sth2=$dbh->prepare($predata_est);
            $sth2->execute($est_id);
            my ($prestatus, $preflags)=$sth2->fetchrow_array();
            $c_est++;
            print "\n=> Updating est: SGN-E$est_id ($c_est of $count_entries entries).\n";
            if ($prestatus eq $status) {
	        print "\tFor SGN-E$est_id, old status ($prestatus) is the same than new ($status)\n";
	    } else {
                $update_est="UPDATE sgn.est SET status = ? WHERE est_id = ?";
                $sth2=$dbh->prepare($update_est);
                $sth2->execute($status, $est_id);
                print "\tUpdate SGN-E$est_id, from status: $prestatus to status: $status\n";
                $update_status_count++;
	    }
            if ($preflags eq $flags) {
                print "\tFor SGN-E$est_id, old flags ($preflags) is the same than new ($flags)\n";
	    } else {
                $update_est="UPDATE sgn.est SET flags = ? WHERE est_id = ?";
                $sth2=$dbh->prepare($update_est);
                $sth2->execute($flags, $est_id);
                print "\tUpdate SGN-E$est_id, from flags: $preflags to flags: $flags\n";
                $update_flags_count++;
	    }
         } 
    } else {
	my $query="SELECT est_id, status, flags, hqi_start, hqi_length, qc_id FROM load_master_table";
        $sth=$dbh->prepare($query);
        $sth->execute();
        while (my ($est_id, $status, $flags, $hqi_start, $hqi_length, $qc_id)=$sth->fetchrow_array()) {  
            my $predata_est="SELECT status, flags FROM sgn.est WHERE est_id = ?";
            my $sth2=$dbh->prepare($predata_est);
            $sth2->execute($est_id);
            my ($prestatus, $preflags)=$sth2->fetchrow_array();
            $c_est++;
            print "\n=> Updating est: SGN-E$est_id ($c_est of $count_entries entries).\n";
            if ($prestatus eq $status) {
	        print "\tFor SGN-E$est_id, old status ($prestatus) is the same than new ($status)\n";
	    } else {
                $update_est="UPDATE sgn.est SET status = ? WHERE est_id = ?";
                $sth2=$dbh->prepare($update_est);
                $sth2->execute($status, $est_id);
                print "\tUpdate SGN-E$est_id, from status: $prestatus to status: $status\n";
                $update_status_count++;
	    }
                if ($preflags eq $flags) {
                print "\tFor SGN-E$est_id, old flags ($preflags) is the same than new ($flags)\n";
	    } else {
                $update_est="UPDATE sgn.est SET flags = ? WHERE est_id = ?";
                $sth2=$dbh->prepare($update_est);
                $sth2->execute($flags, $est_id);
                print "\tUpdate SGN-E$est_id, from flags: $preflags to flags: $flags\n";
                $update_flags_count++;
	    }

            if ($qc_id) {
                my $predata_qc="SELECT hqi_start, hqi_length FROM sgn.qc_report WHERE qc_id = ?";
                $sth2=$dbh->prepare($predata_qc);
                $sth2->execute($qc_id);
                my ($prehqi_start, $prehqi_length)=$sth2->fetchrow_array();
                if ($prehqi_start eq $hqi_start) {
	            print "\tFor SGN-E$est_id, old hqi_start ($prehqi_start) is the same than new ($hqi_start)\n";
	        } else {
                    $update_qc="UPDATE sgn.qc_report SET hqi_start = ? WHERE qc_id = ?";
                    $sth2=$dbh->prepare($update_qc);
                    $sth2->execute($hqi_start, $qc_id);
                    print "\tUpdate SGN-E$est_id, from hqi_start: $prehqi_start to hqi_start: $hqi_start\n";
                    if ($prehqi_start > $hqi_start) {
			$start_before++;
		    } elsif ($prehqi_start < $hqi_start) {
			$start_after++;
		    }
                    $update_hqistart_count++;
	        }
                if ($prehqi_length eq $hqi_length) {
                    print "\tFor SGN-E$est_id, old hqi_length ($hqi_length) is the same than new ($hqi_length)\n";
	        } else {
                    $update_qc="UPDATE sgn.qc_report SET hqi_length = ? WHERE qc_id = ?";
                    $sth2=$dbh->prepare($update_qc);
                    $sth2->execute($hqi_length, $qc_id);
                    print "\tUpdate SGN-E$est_id, from hqi_length: $prehqi_length to hqi_length: $hqi_length\n";
                    if ($prehqi_length > $hqi_length) {
			$shorter_seq++;
		    } elsif ($prehqi_length < $hqi_length) {
			$longer_seq++;
		    }
                    $update_hqilength_count++;
	        }
	    } elsif (!$qc_id) {
		if ($opt_add =~ m/^y(es)/i) {
                    my $previous_check="SELECT qc_id FROM sgn.qc_report WHERE est_id = ?";
                    $sth2=$dbh->prepare($previous_check);
                    $sth2->execute($est_id);
                    my ($old_qc_id) = $sth2->fetchrow_array();
                    if ($old_qc_id) {
			print "Sorry, I can\'t add a new entry fo the est SGN-E$est_id.";
                        print " There is one with qc_id=$old_qc_id.\n";
		    } else {
			my $insert = "INSERT INTO sgn.qc_report (est_id, hqi_start, hqi_length) VALUES (?, ?, ?)";
			$sth2=$dbh->prepare($insert);
			$sth2->execute($est_id, $hqi_start, $hqi_length);
			my $new_qcid_q="SELECT qc_id FROM sgn.qc_report WHERE est_id = ?";
			$sth2=$dbh->prepare($new_qcid_q);
			$sth2->execute($est_id);
			my ($new_qc_id) = $sth2->fetchrow_array();
			print "\tAdded a new qc report entry for the est SGN-E$est_id with hqi_start=$hqi_start, ";
			print "hqi_length=$hqi_length and qc_id=$new_qc_id.\n";
			$new_qc_entries_added++;
		    }
		} else {
		    print "\tThis est (SGN-E$est_id) haven\'t qc_id associated\n";
		}
	    }
        }
    }
    $dbh->do("COPY load_master_table TO '$master_table_file'");
    my @countupdates=($c_est, $update_status_count, $update_flags_count, $update_hqistart_count, 
                      $update_hqilength_count, $new_qc_entries_added,$start_before, $start_after, 
                      $shorter_seq, $longer_seq);
    return @countupdates;
}

=head2 load_entries_in_temp_table

  Usage: my ($n_entries, $n_updated, $n_inserted)=load_single_entries_in_temp_table($input_file, $tag, $temp_clean_coord);
  Desc: Parse the first column of the inpu file, check if exists into temp table and add (or update) entries with tag
  Ret: Three scalars, the count for entries parsed from input file, update count and insert count
  Args: $input_file, $tag (tag for the field), $temp_clean_coord
  Side_Effects: None
  Example:my ($n, $up, $in)=load_single_entries_in_temp_table($input_file, $tag, $temp_clean_coord);

=cut

sub load_single_entries_in_temp_table {
  my $inputfile=shift;
  my $tag=shift;
  
  my ($n,$up,$in)=(0,0,0);
  my $n_f=`cut -f1 $inputfile | wc -l`;
  chomp($n_f);

  open my $file_fh, '<', $inputfile || die "I can not open the file: $inputfile\n";
  while (<$file_fh>) {
      chomp($_);
      my @variables=split(/\t/, $_);
      my $id=$variables[0];
      $id =~ s/SGN-E//;
      my $clean_code;

      my $check="SELECT est_id, clean_code FROM load_master_table WHERE est_id=?";
      my $sth=$dbh->prepare($check);
      $sth->execute($id);
      my ($check_id, $check_clean_code)=$sth->fetchrow_array();
      if ($check_id) {
	  if ($check_clean_code ne 'without_clean_code') {
	      if ($check_clean_code =~ m/\w+/) {
		  $clean_code=$check_clean_code.';'.$tag;
	      }
	      else {
		  $clean_code = $tag;
	      }
	  } else {
	      $clean_code=$tag;
	  }
	  my $update="UPDATE load_master_table SET clean_code=? WHERE est_id=?";
	  $sth=$dbh->prepare($update);
	  $sth->execute($clean_code,$id);
	  $up++;
      } else {
          my $insert="INSERT INTO load_master_table (est_id, clean_code) VALUES (?,?)";
          $sth=$dbh->prepare($insert);
          $sth->execute($id,$tag);
          $in++;
      }
      $n++;
      print "\tLoading line $n of $n_f <$tag file>\r";          
  }
  close $file_fh;
  return ($n, $up, $in);
}

=head2 flags_status_equivalences

Usage: my ($equiv_flags_href, $equiv_status_href)=flags_status_equivalences($equiv_file);
Desc: this subroutine parse a file with the equivalences of clean code, if there are not a file, use the default values
Ret: Two hash references with the equivalences (key='keyword', value='flag or status value')
Args: $equiv_file file name
Side_Effects: none
Example: my ($equiv_flags_href, $equiv_status_href)=flags_status_equivalences($equiv_file);

=cut

sub flags_status_equivalences {
  my $equiv_file=shift;
  my (%equiv_flags,%equiv_status);

  if ($equiv_file) {
      open my $equiv_fh, '<', $equiv_file || die "I can not open the equiv file:$equiv_file\n";
      while (<$equiv_fh>) {
	  chomp($_);
	  my @input = split(/\t/, $_);
          my $tag = lc $input[0];
	  $tag =~ s/^\s+//g;
	  $tag =~ s/\s+$//g;
	  $tag =~ s/\s+/_/g;
	  my $flag=$input[1];
	  my $status=$input[2];

	  $equiv_flags{$tag}=$flag;
	  $equiv_status{$tag}=$status;
      }
      close $equiv_fh;
  } else {                      #DEFAULT VALUES (the match rules will be lowercase, with _ instead [\s,.;])
      # <Possibly chimeric (vector parsing triggered)>, the array for vector keywords will be:
      my @vector_keywords=('trim','vector','adaptor','univec');
      foreach my $vector_kw (@vector_keywords) {
         $equiv_flags{$vector_kw}=2;
         $equiv_status{$vector_kw}=0;
      }
      # <Insert too short>
      my @short_keywords=('short');
      foreach my $short_kw (@short_keywords) {
	  $equiv_flags{$short_kw}=4;
	  $equiv_status{$short_kw}=0;
      }
      # <High expected error (low base calling quality values overall)>
      my @high_expected_error_keywords=('low_qual','high_expected_error','seq_remove_in_phred','removed');
      foreach my $high_expected_error_kw (@high_expected_error_keywords) {
	  $equiv_flags{$high_expected_error_kw}=8;
	  $equiv_status{$high_expected_error_kw}=0;
      }
      # <Low complexity>
      my @low_complexity_keywords=('dust','low_complexity','nseg','trf');
      foreach my $low_compl_kw (@low_complexity_keywords) {
          $equiv_flags{$low_compl_kw}=16;
          $equiv_status{$low_compl_kw}=0;
      }
      # <E.coli or cloning host contamination>  
      my @cloning_host_contamination_keywords=('contam','coli','host');
      foreach my $host_contam_kw (@cloning_host_contamination_keywords) {
          $equiv_flags{$host_contam_kw}=32;
          $equiv_status{$host_contam_kw}=0;
      }
      # <rRNA contamination>  
      my @rRNA_contam_keywords=('rrna');
      foreach my $rrna_contam_kw (@rRNA_contam_keywords) {
	  $equiv_flags{$rrna_contam_kw}=64;
	  $equiv_status{$rrna_contam_kw}=0;
      }
      # <plastid contamination>
      my @plastid_contam_keywords=('plastid');
      foreach my $plastid_contam_kw (@plastid_contam_keywords) {
	  $equiv_flags{$plastid_contam_kw}=1024;
	  $equiv_status{$plastid_contam_kw}=0;
      }
  }
  # the follow tags are from this script and always should have the same tag;
  $equiv_flags{'ath_chimera_screen_positive'}=128;
  $equiv_status{'ath_chimera_screen_positive'}=0;
  $equiv_flags{'unigene_assembly_chimera_screen_positive'}=256;
  $equiv_status{'unigene_assembly_chimera_screen_positive'}=0;
  $equiv_flags{'manually_censored'}=512;
  $equiv_status{'manually_censored'}=0;

  return (\%equiv_flags, \%equiv_status);
}

