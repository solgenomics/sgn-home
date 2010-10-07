#!/usr/bin/perl

=head1 NAME

 load_estseq_into_sgndb.pl.
 A script to load sequences in the SGN databases.(version.2.0.).

 Need a new version of this script, with the follow things:
 * Use of hash to manipulate data (or not, depends of the sequences)
 * Use of check tables.
 * Order the accessions to load the data in the same order
 * Clean de code (always) 

=cut

=head1 SYPNOSIS

load_estseq_into_sgndb.pl -H <dbhost> -D <dbname> -f <P|G|F> -c <cleancoord_file> -q <chimerascreen_file> -s <sequencedata_file> -l <libraries_file> [-o <output_dir>][-L] -i <index_file> 
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory)
      
=item -f 

B<type of file source>  (P for P.hred, G for G.enBank or F for F.asta) (mandatory)
      
=item -c 

B<cleaning coord file>   cleaning coordenates file produced by the SGN_seqclean.pl script (usually with the "db_fasta_clean_coord.tab" name) (seqclean format (*.cln) with 7 columns (id|N%|start(1)|end|seq_lenght|clean_tag|trim_tag))
      
=item -q 

B<chimera seq file>      chimera sequence file produced by the SGN_chimera_screen.pl script (usually with the "chimera_screen_analyze.tab" name)
      
=item -s 

B<sequences table>       the file with the sequences (remember the option -F) (mandatory)

=item -l 

B<libraries table>       (This could be originated by the SGN_PI_GBseq.pl script for the GenBank files or by the SGN_library_form.pl) (mandatory)
      
=item -o 

B<output folder>         (with the files used to load the datas in the database)

=item -i

B<index_file>            index file that order the row in the tables before load them (first columns accessions, and integer for the rest). It can have from 2 to 6 columns (to order clones, reads, ests, qcs, est_dbxref and dbxref)

=item -L

B<load by library>       divide the master table by differents libraries and load one by one (option when have big files)      

=item -h 

B<help>                  show the help

=back

=cut

=head1 DESCRIPTION

This script load sequences and their associated datas in the database in the tables sgn.library, sgn.clone, sgn.seqread, sgn.est and sgn.qc_report. If the source of the sequences file is GenBank (-f G tag) add GenBank accesions to the public.dbxref and sgn.est_dbxref tables. It is necesary the use of the complete path to describe the location of the files.
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

load_estseq_into_sgndb.pl


=cut

use strict;
use Getopt::Std;
use CXGN::DB::InsertDBH;
use File::Basename;

our ($opt_H, $opt_D, $opt_f, $opt_c, $opt_q, $opt_s, $opt_l, $opt_o, $opt_i, $opt_L, $opt_h);
getopts("H:D:f:c:q:s:l:i:o:Lh");
if (!$opt_H && !$opt_D && !$opt_f && !$opt_c && !$opt_q && !$opt_s && !$opt_l && !$opt_o && !$opt_i && !$opt_L && 
    !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my ($cleanfile, $chimerafile, $seqfile, $libraryfile, $indexfile, $dbname, $dbhost)  = validate_the_input();

our $dbh = CXGN::DB::InsertDBH::connect({	dbname=>$dbname ,
						dbhost=>$dbhost ,
					  	dbargs=>{ RaiseError => 1, PrintError => 1 }
					})
	or die "Failed to connect to database ($DBI::errstr)";

check_organism();

my $dbload_dir;
if (!$opt_o) {
    my $path=`pwd`;
    chomp ($path);
    $dbload_dir = $path . "/dbload_files";
}else{
    $dbload_dir = $opt_o;
}

mkdir "$dbload_dir", 0755 or warn "Cannot make $dbload_dir directory: $!";
print_separator('Loading temp data');
my ($temp_seq_data, $temp_library, $temp_clean_coord, $temp_chimera_screen, $temp_indextable) = load_temp_data();
print_separator('Creating mastertable');
my $global_master_table_file = create_dbloadfile ('global_mastertable');

my  ($global_master_table_file, $global_load_master_table) =
   create_master_table($temp_seq_data, $temp_library, $temp_clean_coord, $temp_chimera_screen, $global_master_table_file);

print_separator('Loading libraries data');
my $lib_add=load_library($temp_library, $global_load_master_table, $global_master_table_file);
my $lib_loadmaster;

my ($report_by_libraries);
my $clon_add=0;
my $read_add=0;
my $est_add=0;
my $qcreport_add=0;
my $dbxref_add=0;
my $est_dbxref_add=0;

if ($opt_L) {
    print "\n\n>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<";
    print "\n>>>>>>>>>>>OPTION -L ENABLE<<<<<<<<<<<<";
    print "\n>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<";
    print "\n         Loading data by library";
    print "\n>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<\n\n";
    
    my $total_lib_n= `cut -f3 $opt_l | sort -u | wc -l`;
    chomp($total_lib_n);
    my $lib_count=1;
    
    my @lib_master_tables = split_master_table($temp_library, $global_load_master_table);
    
    my @header = ("LIBRARY", "ENTRIES", "N_CLONES_ADDED", "N_READS_ADDED", "N_ESTS_ADDED", "N_QC_REP_ADDED");
    if ($opt_f =~ m/^G/i) {
        push @header, "N_EST_DBXREF";
        push @header, "N_ACCESSIONS";
    }
    my @report_table;
    push @report_table, \@header;

    foreach $lib_loadmaster (@lib_master_tables) {
        print "\n=> Loading data for $lib_loadmaster... ($lib_count / $total_lib_n)\n";
        my $lib_master_table_file=create_dbloadfile($lib_loadmaster);
        my @entries_added = sequence_of_loads($lib_loadmaster, $lib_master_table_file);
        my @data_entries=@entries_added;
        push @report_table, \@data_entries;

        my $mastertable_name=shift(@entries_added);
        my $mastertable_entries=shift(@entries_added);

        my $libclon_add=shift(@entries_added);
        $clon_add += $libclon_add;
        my $libread_add=shift(@entries_added);
        $read_add += $libread_add;
        my $libest_add=shift(@entries_added);
        $est_add += $libest_add;
        my $libqcreport_add=shift(@entries_added);
        $qcreport_add += $libqcreport_add;
        if ($opt_f =~ m/^G/i) {
            my $libdbxref_add=shift(@entries_added);
            $dbxref_add += $libdbxref_add;
            my $libest_dbxref_add=shift(@entries_added);
            $est_dbxref_add += $libest_dbxref_add;
        }
        $lib_count++;
    }
    $report_by_libraries=\@report_table;
} else {
    my @entries_added=sequence_of_loads($global_load_master_table, $global_master_table_file, $temp_indextable);
    my $mastertable_name=shift(@entries_added);
    my $mastertable_entries=shift(@entries_added);

    $clon_add=shift(@entries_added);
    $read_add=shift(@entries_added);
    $est_add=shift(@entries_added);
    $qcreport_add=shift(@entries_added);
    if ($opt_f =~ m/^G/i) {
        $dbxref_add=shift(@entries_added);
        $est_dbxref_add=shift(@entries_added);
    }
}

if ($opt_L) {
    print_separator('REPORT BY LIBRARIES:');
    print_table($report_by_libraries);
}


print_separator('GLOBAL REPORT:');
print "\tLibraries:\t$lib_add\n\tClones:\t\t$clon_add\n\tReads:\t\t$read_add\n\tESTs(seqs):\t$est_add\n\tQC_reports:\t$qcreport_add\n";
if ($opt_f =~ m/^G/i) {
    print "\tDBxref:\t\t$dbxref_add\n\tEst_Dbxref:\t$est_dbxref_add\n";
}

print_separator();
commit_prompt ($dbh);

#help subroutine

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
      Load sequences and their associated datas in the database in the tables sgn.library, sgn.clone, sgn.seqread, sgn.est and sgn.qc_report. If the source of the sequences file is GenBank (-f G tag) add GenBank accesions to the public.dbxref and sgn.est_dbxref tables. Use the complete path to describe the location of the files.
    
    Usage: 
      load_estseq_into_sgndb.pl -H <dbhost> -D <dbname> -f <P|G|F> -c <cleancoord_file> -q <chimerascreen_file> -s <sequencedata_file> -l <libraries_file> [-o <output_dir>] [-L] 
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -f type of file source  (P for P.hred, G for G.enBank or F for F.asta) (mandatory)
      -c cleanning coordenates file produced by the SGN_seqclean.pl script (usually with the "db_fasta_clean_coord.tab" name)
         (seqclean format (*.cln) with 7 columns (id|N%|start(1)|end|seq_lenght|clean_tag|trim_tag))
      -q chimera sequences file produced by the SGN_chimera_screen.pl script (usually with the "chimera_screen_analyze.tab" name)
      -s sequences table (remember the option -F) (mandatory)
      -l libraries table (This could be originated by the SGN_PI_GBseq.pl script for the GenBank files or by the SGN_library_form.pl) (mandatory)
      -o output folder (with the files used to load the datas in the database)
      -i use index file.      you can use a index file with accessions in the first column and a integer that represent a index to order the accession for load different data in different tables (for example clone). To order clone use 2nd column, read, 3rd column, est 4th column, qc 5th column, est_dbxref 6th column and dbxref 7th column. (this option only can be used with -f G and without -L) 
      -L split the master table in submaster tables by library and load them (option recomended when have big files)
      -h this help

EOF
exit (1);
}

#This subroutine check two things. If exists the input files and how many rows (members) have the diferent files, so means, how many id's and libraries have the seq file, how many id's have the clean file and the chimera file and how many libraries have the library file. Finally check if id's of seq files are the same than the clean file and if the number of libraries in the seq file are the same than in the library file.

=head2 validate_the_input

  Usage: my @input_var=validate_the_input();
  Desc: this subroutine check if exist the mandatory input variables and count the entries in the input files. 
  Ret: six scalars. $cleanfile (cleaning coordenates file name), $chimerafile (chimera id file name), 
       $seqfile (file name of the file with all the sequences and sequence datas), $libraryfile (file
       with the libraries datas), $dbname and $dbhost.
  Args: none
  Side_Effects: die the process if the input data is not right
  Example: my ($cleanfile, $chimerafile, $seqfile, $libraryfile, $dbname, $dbhost)=validate_the_input();

=cut
 
sub validate_the_input {
  print "Check the variables...\n";
  if ($opt_h) {
    help ();
  }

  unless ($opt_f) {
    die ("Required argument -f <type of sequence source> was not supplied. (P for P.hred, G for G.enBank or F for F.asta)\n");
  }
    print "\ttag -f\t$opt_f\t...Ok\n";
  unless ($opt_s) {
    die ("Required argument -s <sequence table file> was not supplied.\n");
  }
    print "\ttag -s\t$opt_s\t...Ok\n";
  unless ($opt_l) {
    die ("Required argument -l <libraries table file> was not supplied.\n");
  }
    print "\ttag -l\t$opt_l\t...Ok\n";

  my $seqfile_count_columns=count_columns($opt_s);
  if ( ($opt_f =~ m/^G/i && $seqfile_count_columns != 25) || 
       ($opt_f =~ m/^P/i && $seqfile_count_columns != 6) || 
       ($opt_f =~ m/^F/i && $seqfile_count_columns != 5) ) {
	  die "Sorry the input sequence file ($opt_s) have $seqfile_count_columns.\nThese must have 24 for GenBank origin files, 6 for Phred origin files and 5 for Free Fasta origin files.\n";
      } else {
	  print "\tColumns number in the $opt_s file ... Ok.\n";
      }
   
  my $libfile_count_columns=count_columns($opt_l);
  

  if ($libfile_count_columns != 12) {
      die "Sorry the library file ($opt_l) have $libfile_count_columns.\nIt must have 12 columns.\n";
  } else {
      print "\tColumns number in the $opt_l file ... Ok.\n";
  }
  
  if ($opt_i) {
      unless ($opt_f =~ m/^G/i) {
	  die "Sorry, -i <index_file> only can be used with -f G option\n";
      }
      if ($opt_L) {
	  die "Sorry, -i <index_table> can not be used with -L option\n";
      }
      my $index_count_column=count_columns($opt_i);
      if ($index_count_column == 2) {
	  print "Index file with 2 columns will order clone\n";
      } elsif ($index_count_column == 3) {
	  print "Index file with 3 columns will order clone and reads\n";
      } elsif ($index_count_column == 4) {
	  print "Index file with 4 columns will order clone, reads and ests\n";
      } elsif ($index_count_column == 5) {
	  print "Index file with 5 columns will order clone, reads, ests and qc_report\n";
      } elsif ($index_count_column == 6) {
	  print "Index file with 6 columns will order clone, reads, ests, qc_report and est_dbxref\n";
      } elsif ($index_count_column == 7) {
	  print "Index file with 7 columns will order clone, reads, ests, qc_report, est_dbxref and dbxref\n";
      } else {
	  die "the -i <index_file> has a wrong number of column ($index_count_column)\n";
      }
  }

  print "\nCount process...\n";
  my ($id_count_seqfile,  $lib_count_seqfile);
  if ($opt_f =~ m/^G/i) {
      $id_count_seqfile = `cut -f2 $opt_s | sort -u | wc -l`;
      $lib_count_seqfile = `cut -f13 $opt_s | sort -u | wc -l`; 
  } elsif ($opt_f =~ m/^P/gi) {
      $id_count_seqfile = `cut -f1 $opt_s | sort -u | wc -l`;
      $lib_count_seqfile = `cut -f5 $opt_s | sort -u | wc -l`;
  } elsif ($opt_f =~ m/^F/gi) {
      $id_count_seqfile = `cut -f1 $opt_s | sort -u | wc -l`;
      $lib_count_seqfile = `cut -f4 $opt_s | sort -u | wc -l`;
  }
    chomp($id_count_seqfile);
    print "\tId_number in $opt_s file:\t$id_count_seqfile\n";
    chomp($lib_count_seqfile);
    print "\tLibraries_number in $opt_s file:\t$lib_count_seqfile\n";
  my $id_count_cleanfile;
  if ($opt_c) {
      $id_count_cleanfile = `cut -f1 $opt_c | sort -u | wc -l`;
      chomp($id_count_cleanfile);
      print "\tId_number in $opt_c file:\t$id_count_cleanfile\n";
  } else {
      print "\nYou haven\'t specified -c option\n";
      possible_abort_process();
  }
  if ($opt_q) {
      my $id_count_chimerafile = `cut -f1 $opt_q | sort -u | wc -l`;
      chomp($id_count_chimerafile);
      print "\tId_number in $opt_q file:\t$id_count_chimerafile\n";
  } else {
      print "\nYou haven\'t specified -q option\n";
      possible_abort_process();
  }
    my $id_count_libraryfile = `cut -f3 $opt_l | sort -u | wc -l`;
    print "\tId_number in $opt_l file:\t$id_count_libraryfile";

  if ($opt_c) {
      if ($id_count_seqfile != $id_count_cleanfile) { 
	  die ("\nthe number of entries in $opt_s file ($id_count_seqfile entries) is different than the number of entries in $opt_c file ($id_count_cleanfile entries)\n");
      }
  }
  if ($lib_count_seqfile != $id_count_libraryfile) {
	die ("\nthe number of entries in $opt_s file ($lib_count_seqfile libraries entries) is different than the number of entries in $opt_l file ($id_count_libraryfile entries)\n");
  }

  if (!$opt_D || !$opt_H) { 
      die "\nneed -D (database name) and -H (host) options";
  }

  return ($opt_c, $opt_q, $opt_s, $opt_l, $opt_i, $opt_D, $opt_H);
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

=head2 check_organism

  Usage: check_organism();
  Desc: check if there are more than organism in the seq file or if this organism in not in the database
  Ret: none
  Args: none
  Side_Effects: die the process
  Example: check_organism();

=cut

sub check_organism {
  my $organism_count;
  if ($opt_f =~ m/^G/i) {
	$organism_count =`cut -f9 $opt_s | sort -u | wc -l`;
  } elsif ($opt_f =~ m/^P/i) {
	$organism_count =`cut -f6 $opt_s | sort -u | wc -l`;
  } elsif ($opt_f =~ m/^F/i) {
	$organism_count =`cut -f5 $opt_s | sort -u | wc -l`;
  }
  chomp($organism_count);
  if ($organism_count > 1) {
	die "Sorry, there are more than one organism in the $opt_s file.\n";
  } else {
	my $organism;
        if ($opt_f =~ m/^G/i) {
	    $organism =`cut -f9 $opt_s | sort -u`;
        } elsif ($opt_f =~ m/^P/i) {
	    $organism =`cut -f6 $opt_s | sort -u`;
        } elsif ($opt_f =~ m/^F/i) {
	    $organism =`cut -f5 $opt_s | sort -u`;
        }
        chomp($organism);
	my $query1="SELECT organism_id FROM sgn.organism WHERE organism_name = ?";
        my $sth=$dbh->prepare($query1);
	$sth->execute($organism);
	my ($organism_id) = $sth->fetchrow_array();
	if (!$organism_id) {
	    die "Sorry, the organism $organism is not in the sgn.organism table in $opt_D database.\n";
        } else {
	    print "Organism in the table $opt_s...Ok (The organism count is $organism_count).\n";
	    print "                            ...Ok ($organism -> organism_id: $organism_id).\n";
	}
  }
}


#This subroutine is for commit or rollback all the datas load in the database.

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

#There are some optionals tables. This subroutine lets abort the process if the option isn't right. For example, if you don't specify the chimera file, could be for two reasons. You haven't or you forgot it. This script tell that you don't specify it and ask if you want abort all the process.

=head2 possible_abort_process

  Usage: possible_abort_process();
  Desc: give to user the possibility of abort the process with a STDIN 
  Ret: none
  Args: none
  Side_Effects: none (die the process)
  Example: if (!$var) {
            print "There are not $var variables.\n";
            possible_abort_process();
	}

=cut
 
sub possible_abort_process {
  my $alert_message = "Do you want abort the process?\n(yes|no, default no)> ";
  print $alert_message;
    
  if (<STDIN> =~ m/^y(es)/i) {
      die "...abort the process\n";
  }
}

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
    if (!$opt_o) {
        my $path=`pwd`;
        chomp ($path);
        $dbload_dir = $path . "/dbload_files";
    }else{
        $dbload_dir = $opt_o;
    }
    my $outfile_name = $dbload_dir ."/". $sufix . ".tab";
    print "Create $outfile_name file.\n";
    open my $outfile, '>', "$outfile_name" || die "Cannot create the file $outfile_name\n";
    chmod 0666, "$outfile_name";
    return $outfile_name;
}

#This subroutine gets the last value for a sequence table...is used to check how many new data have added in a table

=head2 get_last_number

  Usage: my $last_number=get_last_number($sequence_number_id);
  Desc: get the number of the last entrie in a database table.
  Ret: a scalar, the last_number
  Args: $dbseq (the name of the sequence that control the counts in the table).
  Side_Effects: none.
  Example: my $preload_last_lib_id = get_last_number('sgn.library_library_id_seq');

=cut
 
sub get_last_number {
    my $dbseq = shift;
    my $query1 = "SELECT last_value FROM $dbseq";
    my $sth=$dbh->prepare($query1);
    $sth->execute();
    my ($last_count) = $sth->fetchrow_array();
    return $last_count;
}

# Print a separator line

=head2 print_separator

  Usage: print_separator($message);
  Desc: print in the screen a separator with a message.
  Ret: none
  Args: $var (a message).
  Side_Effects: none
  Example: print_separator('this is a test');

=cut
 
sub print_separator {
    my $var=shift;
    print "\n-------------------------------------------------------------\n";
    print "$var\n";
    print "-------------------------------------------------------------\n";
}

# This subroutine load all the data table in the database for a futher process.

=head2 load_temp_data

  Usage: my @load_temp_tables=load_temp_data();
  Desc: load the input file data in the database in temp tables. 
  Ret: four scalar, the temp table names ($temp_seq_data, $temp_library, $temp_clean_coord, $temp_chimera_screen)
  Args: none (use the $opt)
  Side_Effects: if the seq data is from phred or fasta sequences, add the library name.
  Example: my ($temp_seq_data, $temp_library, $temp_clean_coord, $temp_chimera_screen)=load_temp_data();

=cut
 
sub load_temp_data { 
    print "Load files in the temp tables...\n\n";
    my $temp_clean_coord;
    if ($opt_c) {
	$temp_clean_coord = "clean_coord";
	$dbh->do("CREATE TEMP TABLE $temp_clean_coord (preclone_name varchar(250), clone_name varchar(250), n_percentage real, hqi_start int, hqi_length int, length int, clean_code varchar(250), trim_code text)");
	$dbh->do("COPY $temp_clean_coord (preclone_name, n_percentage, hqi_start, hqi_length, length, clean_code, trim_code) FROM '$opt_c'");
	$dbh->do("UPDATE $temp_clean_coord SET clone_name = rtrim(preclone_name, ' ')"); 
	$dbh->do("ALTER TABLE $temp_clean_coord DROP COLUMN preclone_name");
	print "\tLoad the clean file in the temp tables\n";
    }
    my $temp_chimera_screen;
    if ($opt_q) {
	$temp_chimera_screen = "chimera_screen";
	$dbh->do("CREATE TEMP TABLE $temp_chimera_screen (preclone_name varchar(250), clone_name varchar(250), ath_b_q5 varchar(250), ath_b_q3 varchar(250))");
	$dbh->do("COPY $temp_chimera_screen (preclone_name, ath_b_q5, ath_b_q3) FROM '$opt_q'");
	$dbh->do("UPDATE $temp_chimera_screen SET clone_name = rtrim(preclone_name, ' ')"); 
        $dbh->do("ALTER TABLE $temp_chimera_screen DROP COLUMN preclone_name");
        print "\tLoad the chimera file in the temp tables\n";
    }
    my $temp_seq_data = "sequences_data";
    if ($opt_f =~ m/^G/i){
	 $dbh->do("CREATE TEMP TABLE $temp_seq_data (locus varchar(50), accession varchar(50), version varchar(50), primary_id varchar(250), description text, seq text, qscore text, mol_type varchar(50), organism varchar(250), cell_line varchar(250), cell_type varchar(250), clone_name varchar(250), clone_lib varchar(250), cultivar varchar(250), dbxref varchar(250), dev_stage varchar(250), ecotype varchar(250), environmental_sample varchar(250), lab_host varchar(250),note text, PCR_primers text, plasmid varchar(250), tissue_lib varchar(250), tissue_type varchar(250), author text, UNIQUE (accession, version))");
    } 
    elsif ($opt_f =~ m/^P/i) {
	$dbh->do("CREATE TEMP TABLE $temp_seq_data (clone_name varchar(250), seq text, qscore text, call_positions text, library_shortname varchar(16), organism varchar(250))");
    }
    elsif ($opt_f =~ m/^F/i) {
	$dbh->do("CREATE TEMP TABLE $temp_seq_data (clone_name varchar(250), seq text, qscore text, library_shortname varchar(16), organism varchar(250))");
    }

    $dbh->do("COPY $temp_seq_data FROM '$opt_s'");
    if ($opt_f =~ m/^G/i) {
	$dbh->do("UPDATE $temp_seq_data SET clone_lib = 'virtual mRNA library for genbank sequences from '||organism WHERE clone_lib LIKE '%without library tag'");
    }
    print "\tLoad the sequence data file in the temp tables.\n";

    my $temp_indextable;
    if ($opt_i) {
        my $index_count_column=count_columns($opt_i);
	$temp_indextable='temp_index_table_c'.$index_count_column;	
        if ($index_count_column == 2) {
	    $dbh->do("CREATE TEMP TABLE $temp_indextable (accession varchar(250), clone_order int)");
        } elsif ($index_count_column == 3) {
	    $dbh->do("CREATE TEMP TABLE $temp_indextable (accession varchar(250), clone_order int, read_order int)");
        } elsif ($index_count_column == 4) {
     	    $dbh->do("CREATE TEMP TABLE $temp_indextable (accession varchar(250), clone_order int, read_order int, 
                      est_order int)");
        } elsif ($index_count_column == 5) {
	    $dbh->do("CREATE TEMP TABLE $temp_indextable (accession varchar(250), clone_order int, read_order int, 
                      est_order int, qc_order int)");
        } elsif ($index_count_column == 6) {
            $dbh->do("CREATE TEMP TABLE $temp_indextable (accession varchar(250), clone_order int, read_order int, 
                      est_order int, qc_order int, est_dbxref_order int)");
        } elsif ($index_count_column == 7) {
	    $dbh->do("CREATE TEMP TABLE $temp_indextable (accession varchar(250), clone_order int, read_order int, 
                      est_order int, qc_order int, est_dbxref_order int, dbxref_order int)");
        } else {
	    die "the -i <index_file> has a wrong number of column ($index_count_column)\n";
        }
        $dbh->do("COPY $temp_indextable FROM '$opt_i'");
	print "\tLoad the index file into the $temp_indextable table\n";
    }

    my $temp_library = "library_data";
    $dbh->do("CREATE TEMP TABLE $temp_library (type int, clone_lib varchar(250), library_shortname varchar(250), authors text, organism varchar(250), cultivar varchar(250), tissue_type varchar(250), dev_stage varchar(250), lab_host varchar(250), plasmid varchar(250), note text, acc_count int)");
    $dbh->do("COPY $temp_library FROM '$opt_l'");
    $dbh->do("ALTER TABLE $temp_library ADD COLUMN organism_id int");
    $dbh->do("UPDATE $temp_library SET organism_id = (SELECT sgn.organism.organism_id FROM sgn.organism WHERE sgn.organism.organism_name=$temp_library.organism)");
    print "\tLoad the library file in the temp tables.\n";
  
    return ($temp_seq_data, $temp_library, $temp_clean_coord, $temp_chimera_screen, $temp_indextable);
}

# This subroutine create the master table. This is the table that this script use to load the diferent fields in the diferents table. When load a field data in a table, the new id's are update in this table in a new fields... so means that if you load clone_name in the sgn.clone table. This scripts add new table to the master table, the clone_id, and update it.

=head2 create_master_table

  Usage: my $master_table_file=create_master_table($temp_seq_data, $temp_library, $temp_clean_coord, $temp_chimerascreen);
  Desc: make the master table join the different temp tables, add new fields like flags and status, update them using
        the clean codes and copy in a file.
  Ret: A scalar, the master table file name $master_table_file.
  Args: the temp table names ($temp_seq_data, $temp_library, $temp_clean_coord, $temp_chimerascreen)
  Side_Effects: none
  Example: my $load_master_table = create_master_table($temp_seq_data, $temp_library, $temp_clean_coord, 
           $temp_chimera_screen);

=cut
 
sub create_master_table {
    my $temp_seq_data = shift;
    my $temp_library = shift;
    my $temp_clean_coord = shift;
    my $temp_chimera_screen = shift;
    my $master_table_file=shift;

    my $load_master_table='global_master_table';
    my ($prequery1, $prequery2, $prequery3, $prequery4, $prequery5, $prequery6, $prequery7);

    $prequery1 = "SELECT $temp_seq_data.*,";
    if ($opt_c) {
	$prequery2 = "hqi_start, hqi_length, clean_code, n_percentage,";
        if ($opt_f =~ m/^G/i) {
	    $prequery5 = "LEFT JOIN $temp_clean_coord ON $temp_seq_data.accession=$temp_clean_coord.clone_name";
        } elsif ($opt_f =~ m/^P|^F/i) {
	     $prequery5 = "LEFT JOIN $temp_clean_coord ON $temp_seq_data.clone_name=$temp_clean_coord.clone_name";
	}
    }
    if ($opt_q) {
	$prequery3 = "ath_b_q5, ath_b_q3,";
	if ($opt_f =~ m/^G/i) {
	    $prequery6 = "LEFT JOIN $temp_chimera_screen ON $temp_seq_data.accession=$temp_chimera_screen.clone_name";
	} elsif ($opt_f =~ m/^P|^F/i) {
	    $prequery6 ="LEFT JOIN $temp_chimera_screen ON $temp_seq_data.clone_name =$temp_chimera_screen.clone_name";
	}
    }
    $prequery4 = "sgn.clone.clone_id, sgn.seqread.read_id, sgn.est.est_id INTO TEMP TABLE $load_master_table FROM $temp_seq_data";
    if ($opt_f =~ m/^G/i) {
	$prequery7 = "LEFT JOIN public.dbxref ON $temp_seq_data.accession||$temp_seq_data.version=public.dbxref.accession||public.dbxref.version LEFT JOIN sgn.est_dbxref ON public.dbxref.dbxref_id=sgn.est_dbxref.dbxref_id LEFT JOIN sgn.est ON sgn.est_dbxref.est_id=sgn.est.est_id LEFT JOIN sgn.seqread ON sgn.est.read_id=sgn.seqread.read_id LEFT JOIN sgn.clone ON sgn.seqread.clone_id=sgn.clone.clone_id WHERE mol_type = 'EST' OR mol_type = 'mRNA'";
    }
    elsif ($opt_f =~ m/^P|^F/i) {
	$prequery7 = "LEFT JOIN sgn.clone ON $temp_seq_data.clone_name||$temp_seq_data.library_shortname=sgn.clone.clone_name||(SELECT sgn.library.library_shortname FROM sgn.library WHERE sgn.library.library_id=sgn.clone.library_id) LEFT JOIN sgn.seqread ON sgn.clone.clone_id=sgn.seqread.clone_id LEFT JOIN sgn.est ON sgn.seqread.read_id=sgn.est.read_id";
    }
   
    my $query = "$prequery1 $prequery2 $prequery3 $prequery4 $prequery5 $prequery6 $prequery7";
    print "\nThe constructor for the $load_master_table is ...\n\n$query\n\n";
    $dbh->do("$query");
    print "Adding qc_id to $load_master_table ...\t";
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN qc_id int");
    $dbh->do("UPDATE $load_master_table SET qc_id = (SELECT qc_id FROM sgn.qc_report WHERE sgn.qc_report.est_id=$load_master_table.est_id)");
    print "... added.\n\n";

    if ($opt_f =~ m/^G/i) {
	$dbh->do("ALTER TABLE $load_master_table ADD COLUMN est_dbxref_id int");
	$dbh->do("ALTER TABLE $load_master_table ADD COLUMN dbxref_id int");
	$dbh->do("UPDATE $load_master_table SET est_dbxref_id = (SELECT sgn.est_dbxref.est_dbxref_id FROM sgn.est_dbxref WHERE sgn.est_dbxref.est_id=$load_master_table.est_id AND sgn.est_dbxref.dbxref_id=(SELECT public.dbxref.dbxref_id FROM public.dbxref WHERE public.dbxref.accession=$load_master_table.accession AND public.dbxref.version=$load_master_table.version AND public.dbxref.db_id=(SELECT db_id FROM public.db WHERE name= 'DB:GenBank_Accession')))");
	$dbh->do("UPDATE $load_master_table SET dbxref_id = (SELECT dbxref_id FROM sgn.est_dbxref WHERE sgn.est_dbxref.est_id=$load_master_table.est_id AND sgn.est_dbxref.dbxref_id=(SELECT public.dbxref.dbxref_id FROM public.dbxref WHERE public.dbxref.accession=$load_master_table.accession AND public.dbxref.version=$load_master_table.version AND public.dbxref.db_id=(SELECT db_id FROM public.db WHERE name= 'DB:GenBank_Accession')))");
    }

    print "Adding flags and status to $load_master_table ...\t";
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN flags int");
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN status int");

    if ($opt_c) {
	$dbh->do("UPDATE $load_master_table SET flags = 8, status = 8 WHERE seq = '' AND clean_code = 'seq_remove_in_phred_trimming'");
        $dbh->do("UPDATE $load_master_table SET flags = 8, status = 0 WHERE clean_code ILIKE '%low_quality%'");
        $dbh->do("UPDATE $load_master_table SET flags = 8, status = 0 WHERE clean_code ILIKE '%high expected error%'");
	$dbh->do("UPDATE $load_master_table SET flags = 4, status = 0 WHERE clean_code ILIKE '%short%'");
	$dbh->do("UPDATE $load_master_table SET flags = 2, status = 0 WHERE clean_code ILIKE '%UniVec'");
	$dbh->do("UPDATE $load_master_table SET flags = 32, status = 0 WHERE clean_code ILIKE '%ColiBank95.lib'");
	$dbh->do("UPDATE $load_master_table SET flags = 32, status = 0 WHERE clean_code ILIKE '%dust'");
	$dbh->do("UPDATE $load_master_table SET flags = 8, status = 0 WHERE clean_code ILIKE '%low_qual'");
        $dbh->do("UPDATE $load_master_table SET flags = 2, status = 0 WHERE clean_code ILIKE '%trim%'");
        $dbh->do("UPDATE $load_master_table SET flags = 32, status = 0 WHERE clean_code ILIKE '%contam%'");
        $dbh->do("UPDATE $load_master_table SET flags = 1024, status = 0 WHERE clean_code ILIKE '%plastid_contam%'");
	$dbh->do("UPDATE $load_master_table SET flags = 1024, status = 0 WHERE clean_code ILIKE '%solanaceae-plastid-mito-noncoding%'");
    }
    if ($opt_q) {
    $dbh->do("UPDATE $load_master_table SET flags = 128, status = 0 WHERE ath_b_q5 IS NOT NULL AND ath_b_q3 IS NOT NULL");
    }

    $dbh->do("UPDATE $load_master_table SET flags = 0, status = 0 WHERE flags IS NULL");
    print "... added.\n\n";

    my $load_master_table_idx1=$load_master_table."_idx1";
    my $load_master_table_idx2=$load_master_table."_idx2";
    if ($opt_f =~ m/^G/i) {
	$dbh->do("CREATE INDEX $load_master_table_idx1 ON $load_master_table (accession)");
	$dbh->do("CREATE INDEX $load_master_table_idx2 ON $load_master_table (clone_name)");
	$dbh->do("UPDATE $load_master_table SET clone_name = accession WHERE clone_name LIKE 'NULL'");
	$dbh->do("UPDATE $load_master_table SET clone_name = accession WHERE clone_name IS NULL");
    }
    print "Adding organism_id to $load_master_table ...\t";
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN organism_id int");
    $dbh->do("UPDATE $load_master_table SET organism_id = (SELECT organism_id FROM sgn.organism WHERE sgn.organism.organism_name=$load_master_table.organism)");
    print "... added.\n\n";

    $dbh->do("COPY $load_master_table TO '$master_table_file'");
    print "Copy the basic master table\n";
    return ($master_table_file, $load_master_table);
}

# This subroutine counts the entries for a field. 

=head2 count_new_entries

  Usage: my $new_seq_count=count_new_entries($field, $constraint_null);
  Desc: count the number of new entries for a concrete field.
  Ret: a scalar, the number of new entries ($new_seq).
  Args: $field (the name of the field of the master table where you want count the new entries)
  Side_Effects: none
  Example: my $new_clones = count_new_entries('clone_name', 'clone_id');

=cut
 
sub count_new_entries {
    my $field=shift;
    my $null_constraint=shift;
    my $load_master_table=shift;
    my $query = "SELECT COUNT(DISTINCT $field) FROM $load_master_table WHERE $null_constraint IS NULL";
    my $sth = $dbh->prepare($query);
    $sth->execute();
    my ($new_seq) = $sth->fetchrow_array();
    return $new_seq;
}

#This subroutine get the library data, check how many new libraries have (use organism_id + library_name) and add them in a file. Use this file to load the data in the sgn.library by COPY SQL command. Also check how many new libraries have added to the sgn.library. If not is the same than in the load_library file die the process.

=head2 load_library

  Usage: my $libraries_added=load_library($temp_library_table);
  Desc: select the library data from the $temp_library_table, copy in a file and load in the database in the 
        sgn.library table. Check if it have loaded the same number of entries in the table that are in the file. Search
        the new library_id and add to the master table.
  Ret: a scalar, the number of libraries added ($libraries_added);
  Args: $temp_library_table (libraries temp table)
  Side_Effects: die if the check is wrong
  Example: my $libraries_added=load_library($temp_library_table);

=cut
 
sub load_library {
    my $temp_library=shift;
    my $load_master_table=shift;
    my $master_table_file=shift;
    
    my $id_count_libraryfile =  `cut -f3 $opt_l | sort -u | wc -l`;
    $dbh->do("SELECT $temp_library.type, $temp_library.clone_lib, $temp_library.library_shortname, $temp_library.authors, $temp_library.organism_id, $temp_library.cultivar, $temp_library.tissue_type, $temp_library.dev_stage, $temp_library.lab_host, $temp_library.plasmid, $temp_library.note, sgn.library.library_id INTO TEMP TABLE preload_library FROM $temp_library LEFT JOIN sgn.library ON $temp_library.clone_lib||$temp_library.organism_id = sgn.library.library_name||sgn.library.organism_id");
    my $query1 = "SELECT COUNT(*) FROM preload_library WHERE library_id IS NULL";
    my $sth = $dbh->prepare($query1);
    $sth->execute();
    my ($new_lib_n) = $sth->fetchrow_array();
    my $query2 = "SELECT COUNT(sgn.library.library_id) FROM sgn.library JOIN $temp_library ON sgn.library.library_shortname = $temp_library.library_shortname";
    $sth->execute();
    my ($new_lib_shn) = $sth->fetchrow_array();
    print "There are $new_lib_n new libraries in the library file checked by library_name and $new_lib_shn new libraries in the library_file checked by library_shortname.\n";

    my $libraries_added;
    if ($new_lib_n > 0 && $new_lib_shn > 0) {
	$dbh->do("SELECT type, clone_lib, library_shortname, authors, organism_id, cultivar, tissue_type, dev_stage, lab_host, plasmid, note INTO TEMP TABLE load_library FROM preload_library WHERE library_id IS NULL");
  
	my $library_preload = create_dbloadfile('library');
	$dbh->do("COPY load_library TO '$library_preload'");
	my $library_load_file_count = `cut -f3 $library_preload | sort -u | wc -l`;
        chomp ($library_load_file_count);
	print "There are $library_load_file_count libraries in the library_preload_file\n";
	my $preload_last_lib_id = get_last_number('sgn.library_library_id_seq');
	$dbh->do("COPY sgn.library (type, library_name, library_shortname, authors, organism_id, cultivar, tissue, development_stage, cloning_host, vector, comments) FROM '$library_preload'");
	my $postload_last_lib_id = get_last_number('sgn.library_library_id_seq');
	$libraries_added = $postload_last_lib_id - $preload_last_lib_id; 
	if ( $libraries_added != $library_load_file_count) {
	    die "Sorry, the number of libraries added to sgn.library table ($libraries_added) is not the same than the libraries of the $library_preload file.\nCheck it before rerun this script.\n\n";
	} else {
	    print "Added $libraries_added to sgn.library table.\n\n";
        }
    } else {
	$libraries_added=0;
    }
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN library_id int");
	
    if ($opt_f =~ m/^G/i) {
	$dbh->do("UPDATE $load_master_table SET library_id = (SELECT library_id FROM sgn.library WHERE sgn.library.library_name||sgn.library.organism_id = $load_master_table.clone_lib||$load_master_table.organism_id)");
    } elsif ($opt_f =~ m/^P/i) {
	$dbh->do("UPDATE $load_master_table SET library_id = (SELECT library_id FROM sgn.library WHERE sgn.library.library_shortname||sgn.library.organism_id = $load_master_table.library_shortname||$load_master_table.organism_id)");
    } elsif ($opt_f =~ m/^F/i) {
	$dbh->do("UPDATE $load_master_table SET library_id = (SELECT library_id FROM sgn.library WHERE sgn.library.library_shortname||sgn.library.organism_id = $load_master_table.library_shortname||$load_master_table.organism_id)");
    }
	$dbh->do("COPY $load_master_table TO '$master_table_file'");
	print "Copy the sgn.library.library_id to mastertable\n";
        return $libraries_added;
    
}

# This script load the clone_name in the sgn.clone table and check if this number is the same than the number of new clone name of the sequence data file (in the mastertable).

=head2 load_clone

  Usage: my $clones_added=load_clone();
  Desc: load new clones in the database, check if it have added the same number of entries that can find in the 
        master load table. If it is not the same, die. Remember that perhaps the number the clone not is the same
        that the number of accessions (some accessions could have the same clone_name). Search the new clone_id and
        add to master table.
  Ret: a scalar, the number of clones added to sgn.clone table
  Args: none
  Side_Effects: create load tables, add the new clone_id to the master table
  Example: my $clones_added=load_clone();

=cut
 
sub load_clone {
    my $load_master_table=shift;
    my $master_table_file=shift;
    my $temp_indextable=shift;
    my $index_counts = 0;
    if ($opt_i) {
	my @split_table=split(/_/, $temp_indextable);
	$index_counts=$split_table[3];
    }

    my $new_clones = count_new_entries('clone_name', 'clone_id', $load_master_table);
    print "There are $new_clones new clones associated to $opt_s file\n";
    my $load_clone="load_clone_".$load_master_table;
    if ($index_counts > 1) {
	$dbh->do("SELECT library_id, clone_name INTO TEMP TABLE $load_clone FROM $load_master_table 
                  LEFT JOIN $temp_indextable ON $load_master_table.accession=$temp_indextable.accession  
                  WHERE $load_master_table.clone_id IS NULL GROUP BY library_id, clone_name 
                  ORDER BY $temp_indextable.clone_order ASC");
    } else {
        $dbh->do("SELECT library_id, clone_name INTO TEMP TABLE $load_clone FROM $load_master_table 
                  WHERE $load_master_table.clone_id IS NULL GROUP BY library_id, clone_name");
    }

    my $clone="clone_".File::Basename::basename($master_table_file);
    my $clone_preload=create_dbloadfile($clone);
    $dbh->do("COPY $load_clone TO '$clone_preload'");
    my $clone_n_preload = `cut -f2 $clone_preload | sort -u | wc -l`;
    chomp($clone_n_preload);
    print "There are $clone_n_preload UNIQUES clones in the $clone_preload file.\n";
    if ( $new_clones != $clone_n_preload) {
            die "Sorry, the number of new clones in the $opt_s file ($new_clones) is not the same than the number of clones fo the preload clone file ($clone_n_preload).\n";
    }
    my $preload_last_clone_id = get_last_number('sgn.clone_clone_id_seq');
    $dbh->do("COPY sgn.clone (library_id, clone_name) FROM '$clone_preload'");
    my $postload_last_clone_id = get_last_number('sgn.clone_clone_id_seq');
    my $clones_added = $postload_last_clone_id - $preload_last_clone_id;
    if ($clones_added != $clone_n_preload) {
        my $clone_n_preload_non_u=`cut -f2 $clone_preload | sort | wc -l`;
        chomp($clone_n_preload_non_u);
        if ($clone_n_preload != $clone_n_preload_non_u) {
            my $diff_clones=$clone_n_preload_non_u-$clone_n_preload;
	    print "The number of clone names uniques ($clone_n_preload) is not the same than the number of clone names ($clone_n_preload_non_u) in the file.\nThere are $diff_clones repeat clone names in the file: $clone_preload (Possibly because there are more than one clone with the same name but different library)\n";
            possible_abort_process();
        } else { 
            die "Sorry, the number of clones added to sgn.clone table ($clones_added) is not the same than the new clones in the $opt_s file.\n Check it before rerun this script.\n\n";
	}
    } else {
	print "Added $clones_added to sgn.clone table\n\n";
    }
    my $load_master_table_idx3=$load_master_table."_idx3";
    $dbh->do("CREATE INDEX $load_master_table_idx3 ON $load_master_table (library_id)");
    my $sub_clone="subSGN_clone_".$load_master_table;
    $dbh->do("SELECT clone_id, clone_name, library_id INTO TEMP TABLE $sub_clone FROM sgn.clone WHERE clone_id BETWEEN $preload_last_clone_id AND $postload_last_clone_id");
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN new_clone_id int");
    $dbh->do("UPDATE $load_master_table SET new_clone_id = (SELECT clone_id FROM $sub_clone WHERE $sub_clone.clone_name||$sub_clone.library_id =$load_master_table.clone_name||$load_master_table.library_id) WHERE clone_id IS NULL");
    $dbh->do("UPDATE $load_master_table SET new_clone_id = clone_id WHERE clone_id IS NOT NULL");
    $dbh->do("COPY $load_master_table TO '$master_table_file'");
    print "Copy the sgn.clone.clone_id to $load_master_table\n";
    return $clones_added;
}

#This subroutine load the read data in the sgn.seqread table. If you have chromatograms (and choose the -f P tag) this subroutine ask about the trace location to add to the table. After that, check if have inserted the same number of entries that the preload file have.

=head2 load_read

  Usage: my $reads_added=load_read($load_master_table, $master_table_file);
  Desc: For -f G (GenBank) select clone_id (added in the load_clone()) and accession (like trace_name) in a new temp
        table. For -f P select clone_id and clone_name (like trace_name) and ask about the location pf the chromatograms.
        Copy this tables to a file and use this file to load the data in the sgn.seqread. Check that the number the new 
        entries is the same that the number of entries in the preload file. Search the new read_id and add to the master
        table.
  Ret: a scalar, the number of clones added. 
  Args: $load_master_table and $master_table_file
  Side_Effects: if the option is P, ask for the trace_location using STDIN. die the process if the check is wrong. 
  Example: my $reads_added=load_read($load_master_table, $master_table_file);

=cut
 
sub load_read {
    my $load_master_table=shift;
    my $master_table_file=shift;
    my $temp_indextable=shift;
    my $index_counts = 0;
    if ($opt_i) {
	my @split_table=split(/_/, $temp_indextable);
	$index_counts=$split_table[3];
    }
    my $trace_location;

    my $new_read = count_new_entries('clone_name', 'read_id', $load_master_table);
    print "There are $new_read new reads associated to $opt_s file\n";
    
    my $load_seqread="load_seqread_".$load_master_table;
    if ($opt_f =~ m/^G/i ) {
	if ($index_counts > 2) {
	    $dbh->do("SELECT new_clone_id, accession INTO TEMP TABLE $load_seqread FROM $load_master_table 
                      LEFT JOIN $temp_indextable ON $load_master_table.accession=$temp_indextable.accession
                      WHERE $load_master_table.read_id IS NULL ORDER BY $temp_indextable.read_order ASC");
	} else {
	    $dbh->do("SELECT new_clone_id, accession INTO TEMP TABLE $load_seqread FROM $load_master_table 
                      WHERE $load_master_table.read_id IS NULL");
	}
    }
    elsif ($opt_f =~ m/^P/i ) {
	$dbh->do("SELECT new_clone_id, clone_name INTO TEMP TABLE $load_seqread FROM $load_master_table WHERE $load_master_table.read_id IS NULL");
	print "Do you want specify the trace location for the chromatograms? (yes|no, default no)\n";
	if (<STDIN> =~ m/^y(es)/i) {
	    print "Write the complete path for the directory with the chromatograms>";
	    $trace_location = <STDIN>;
	    chomp($trace_location);
	    $dbh->do("ALTER TABLE $load_seqread ADD COLUMN trace_location text");
	    $dbh->do("UPDATE $load_seqread SET trace_location = '$trace_location'");
	}
    }
    elsif ($opt_f =~ m/^F/i) {
	$dbh->do("SELECT new_clone_id, clone_name INTO TEMP TABLE $load_seqread FROM $load_master_table WHERE $load_master_table.read_id IS NULL");
    }
    my $seqread="seqread_".File::Basename::basename($master_table_file);
    my $preload_seqread=create_dbloadfile($seqread);
    $dbh->do("COPY $load_seqread TO '$preload_seqread'");
    my $read_n_preload = `cut -f2 $preload_seqread | sort -u | wc -l`;
    chomp($read_n_preload);
    print "There are $read_n_preload trace_names in the $preload_seqread file.\n";
    my $preload_last_read_id = get_last_number('sgn.seqread_read_id_seq');
    if ($trace_location) {
	$dbh->do("COPY sgn.seqread (clone_id, trace_name, trace_location) FROM '$preload_seqread'");
    } else {
	$dbh->do("COPY sgn.seqread (clone_id, trace_name) FROM '$preload_seqread'");
    }
    my $postload_last_read_id = get_last_number('sgn.seqread_read_id_seq');
    my $read_added = $postload_last_read_id - $preload_last_read_id;
    if ( $read_added != $read_n_preload) {
	die "Sorry, the number of reads added to sgn.seqread table ($read_added) is not the same than the reads in the $opt_s file.\n Check it before rerun this script.\n\n";
    } else {
	print "Added $read_added reads to sgn.seqread table\n\n";
    }
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN new_read_id int");
    
    my $sub_read="subSGN_read_".$load_master_table;
    $dbh->do("SELECT * INTO TEMP TABLE $sub_read FROM sgn.seqread WHERE read_id BETWEEN $preload_last_read_id AND $postload_last_read_id");
    if ($opt_f =~ m/^G/i) {
	$dbh->do("UPDATE $load_master_table SET new_read_id = (SELECT read_id FROM $sub_read WHERE $sub_read.trace_name=$load_master_table.accession) WHERE read_id IS NULL");
    }
    elsif ($opt_f =~ m/^P|^F/i) {
	$dbh->do("UPDATE $load_master_table SET new_read_id = (SELECT read_id FROM $sub_read WHERE $sub_read.trace_name=$load_master_table.clone_name) WHERE read_id IS NULL");
    }
    $dbh->do("UPDATE $load_master_table SET new_read_id = read_id WHERE read_id IS NOT NULL");
    $dbh->do("COPY $load_master_table TO '$master_table_file'");
    print "Copy the sgn.seqread.read_id to $load_master_table\n";
    return $read_added;
}

# This subroutine load the est's data in the sgn.est table. If the sequences come from chromatograms (and they have been processed by Phred, -f P tag) load the call positions too. After that check if it have loaded the same number of entries that the preload file have.

=head2 load_est

  Usage: my $ests_added=load_est($load_master_table, $master_table_file);
  Desc: select read_id (added to master table in the load_read() subroutine), seq, qscore, status and flags from
        master table into temp table. If the opt f is P also select call_postions. Copy it to a file and load the data
        in the sgn.est table. Check if the number of the entries added is the same that the number of entries of the
        file. Search the new est_id and add to master table.   
  Ret: a scalar, the number of ests added.
  Args: $load_master_table and $master_table_file
  Side_Effects: die the process if the check is wrong.
  Example: my $ests_added=load_est($load_master_table, $master_table_file);

=cut
 
sub load_est {
    my $load_master_table=shift;
    my $master_table_file=shift;
    my $temp_indextable=shift;
    
    my $index_counts = 0;
    if ($opt_i) {
	my @split_table=split(/_/, $temp_indextable);
	$index_counts=$split_table[3];
    }

    my $new_est = count_new_entries('clone_name', 'est_id', $load_master_table);
    print "There are $new_est new sequences entries associated to $opt_s file\n";   
 
    my $load_est="load_est_".$load_master_table;
    if ($opt_f =~ m/^G|^F/i) {
	if ($index_counts > 3) {
	    $dbh->do("SELECT new_read_id, seq, qscore, status, flags INTO TEMP TABLE $load_est FROM $load_master_table 
                      LEFT JOIN $temp_indextable ON $load_master_table.accession=$temp_indextable.accession
                      WHERE $load_master_table.est_id IS NULL ORDER BY $temp_indextable.est_order ASC");
	} else {
	    $dbh->do("SELECT new_read_id, seq, qscore, status, flags INTO TEMP TABLE $load_est FROM $load_master_table 
                      WHERE $load_master_table.est_id IS NULL");
	}
    }
    elsif ($opt_f =~ m/^P/i) {
	$dbh->do("SELECT new_read_id, seq, qscore, call_positions, status, flags INTO TEMP TABLE $load_est 
                  FROM $load_master_table WHERE $load_master_table.est_id IS NULL ");
    }
    my $est="est_".File::Basename::basename($master_table_file);
    my $preload_est=create_dbloadfile($est);
    $dbh->do("COPY $load_est TO '$preload_est'");
    my $est_n_preload = `cut -f1 $preload_est | sort -u | wc -l`;
    chomp($est_n_preload);
    print "There are $est_n_preload ESTs in the $preload_est file.\n";
    my $preload_last_est_id=get_last_number('sgn.est_est_id_seq');
    if ($opt_f =~ m/^G|^F/i) {
	$dbh->do("COPY sgn.est (read_id, seq, qscore, status, flags) FROM '$preload_est'");
    }
    elsif ($opt_f =~ m/^P/i) {
	$dbh->do("COPY sgn.est (read_id, seq, qscore, call_positions, status, flags) FROM '$preload_est'");
    }
    my $postload_last_est_id=get_last_number('sgn.est_est_id_seq');
    my $est_added = $postload_last_est_id - $preload_last_est_id;
    if ( $est_added != $est_n_preload) {
	die "Sorry, the number of ests added to sgn.est table ($est_added) is not the same than the est in the $opt_s, $opt_c and $opt_q files.\n Check them before rerun this script.\n\n";
    } else {
	print "Added $est_added ESTs to sgn.est table\n\n";
    }
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN new_est_id int");
    
    my $sub_est="subSGN_est_".$load_master_table;
    $dbh->do("SELECT * INTO TEMP TABLE $sub_est FROM sgn.est WHERE est_id BETWEEN $preload_last_est_id AND $postload_last_est_id");
    $dbh->do("UPDATE $load_master_table SET new_est_id = (SELECT est_id FROM $sub_est WHERE $sub_est.read_id=$load_master_table.new_read_id) WHERE est_id IS NULL");
    $dbh->do("UPDATE $load_master_table SET new_est_id = est_id WHERE est_id IS NOT NULL");
    $dbh->do("COPY $load_master_table TO '$master_table_file'");
    print "Copy the sgn.est.est_id to $load_master_table\n";
    return $est_added;
}

#This subroutine load the qc report in the sgn.qc_report table. After that check if have inserted the same numbre of entries that have the preload file. 

=head2 load_qc_report

  Usage: my $qcreport_added=load_qc_report($load_master_table, $master_table_file);
  Desc: select est_id (added to master table in the load_est), hqi_start and hqi_length from master table into temp table.
        Copy it in a file and use it to load the data in the sgn.qc_report table. Check if the number of entries added
        is the same than the entries of the preload file. Finally add the new qc_id to master table.
  Ret: a scalar, the number of entries added to qc_report table.
  Args: $load_master_table and $master_table_file
  Side_Effects: die if the check is wrong
  Example: my $qcreport_added=load_qc_report($load_master_table, $master_table_file);

=cut
 
sub load_qc_report {
    my $load_master_table=shift;
    my $master_table_file=shift;
    my $temp_indextable=shift;
    
    my $index_counts = 0;
    if ($opt_i) {
	my @split_table=split(/_/, $temp_indextable);
	$index_counts=$split_table[3];
    }
    my $new_qcreport = count_new_entries('clone_name', 'qc_id', $load_master_table);
    print "There are $new_qcreport new qc reports associated to $opt_s file\n";
    
    my $load_qc_report="load_qcreport".$load_master_table;
    if ($index_counts > 4) {
        $dbh->do("SELECT new_est_id, hqi_start, hqi_length INTO TEMP TABLE $load_qc_report FROM $load_master_table 
                  LEFT JOIN $load_master_table.accession=$temp_indextable.accession
                  WHERE $load_master_table.qc_id IS NULL ORDER BY $temp_indextable.qc_order ASC;");
    } else {
	$dbh->do("SELECT new_est_id, hqi_start, hqi_length INTO TEMP TABLE $load_qc_report FROM $load_master_table 
                  WHERE $load_master_table.qc_id IS NULL;");
    }
    my $qc_report="qc_report_".File::Basename::basename($master_table_file);
    my $preload_qcreport=create_dbloadfile($qc_report);
    $dbh->do("COPY $load_qc_report TO '$preload_qcreport'");
    my $qc_report_n_preload=`cut -f1 $preload_qcreport | sort -u | wc -l`;
    chomp($qc_report_n_preload);
    print "There are $qc_report_n_preload reports in the $preload_qcreport file.\n";
    my $preload_last_qc_id=get_last_number('sgn.qc_report_qc_id_seq');
    $dbh->do("COPY sgn.qc_report (est_id, hqi_start, hqi_length) FROM '$preload_qcreport'");
    my $postload_last_qc_id=get_last_number('sgn.qc_report_qc_id_seq');
    my $qcreport_added = $postload_last_qc_id - $preload_last_qc_id;
    if ( $qcreport_added != $qc_report_n_preload) {
	die "Sorry, the number of qc report added to sgn.qc_report table ($qcreport_added) is not the same than the qc report in the $opt_c file.\n Check it before rerun this script.\n\n";
    } else {
	print "Added $qcreport_added QC reports to sgn.qc_report table\n\n";
    }
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN new_qc_id int");
    my $sub_qcreport="subSGN_qcreport_".$load_master_table;
    $dbh->do("SELECT * INTO TEMP TABLE $sub_qcreport FROM sgn.qc_report WHERE qc_id BETWEEN $preload_last_qc_id AND $postload_last_qc_id");
    $dbh->do("UPDATE $load_master_table SET new_qc_id = (SELECT qc_id FROM $sub_qcreport WHERE $sub_qcreport.est_id=$load_master_table.new_est_id) WHERE qc_id IS NULL");
    $dbh->do("UPDATE $load_master_table SET new_qc_id = qc_id WHERE qc_id IS NOT NULL");
    $dbh->do("COPY $load_master_table TO '$master_table_file'");
    print "Copy the sgn.qc_report.qc_id to $load_master_table\n";
    return $qcreport_added;
}

# This subroutine load the GenBank accessions in the public.dbxref table. After that check if have inserted the same number of accessions that the preload file have. 

=head2 load_dbxref

  Usage: my $dbxref_added=load_dbxref($load_master_table, $master_table_file);
  Desc: select accession and version from the master table into temp table. Copy it to a file and load the data into the 
        public.dbxref table using this file. Check that the number of new entries is the same that the number of entries 
        of the file. Add the new dbxref_id to master table. 
  Ret: a scalar, the number of new entries added to public.dbxref table.
  Args: $load_master_table and $master_table_file
  Side_Effects: die if the check is wrong.
  Example: my $accessions_added=load_dbxref($load_master_table, $master_table_file);

=cut
 
sub load_dbxref {
    my $load_master_table=shift;
    my $master_table_file=shift;
    my $temp_indextable=shift;
    
    my $index_counts = 0;
    if ($opt_i) {
	my @split_table=split(/_/, $temp_indextable);
	$index_counts=$split_table[3];
    }
    my $new_accessions = count_new_entries('accession', 'dbxref_id', $load_master_table);
    print "There are $new_accessions new accessions associated to $opt_s file\n";
    
    my $load_dbxref="load_dbxref_".$load_master_table;
    if ($index_counts > 5) {
	$dbh->do("SELECT (SELECT db_id FROM public.db WHERE name = 'DB:GenBank_Accession'), accession, version, organism 
                  INTO TEMP TABLE $load_dbxref FROM $load_master_table 
                  LEFT JOIN $temp_indextable ON $load_master_table.accession=$temp_indextable.accession
                  WHERE $load_master_table.dbxref_id IS NULL ORDER BY $temp_indextable.dbxref_order ASC");
    } else {
	$dbh->do("SELECT (SELECT db_id FROM public.db WHERE name = 'DB:GenBank_Accession'), accession, version, organism 
                  INTO TEMP TABLE $load_dbxref FROM $load_master_table WHERE $load_master_table.dbxref_id IS NULL");
    }
    $dbh->do("ALTER TABLE $load_dbxref ADD COLUMN description text");
    $dbh->do("UPDATE $load_dbxref SET description = 'This is the GenBank accession for a sequence from '||organism");
    $dbh->do("ALTER TABLE $load_dbxref DROP COLUMN organism");
    
    my $dbxref="dbxref_".File::Basename::basename($master_table_file);
    my $preload_dbxref=create_dbloadfile($dbxref);

    $dbh->do("COPY $load_dbxref TO '$preload_dbxref'");
    my $dbxref_n_preload=`cut -f2 $preload_dbxref | sort -u | wc -l`;
    chomp($dbxref_n_preload);
    print "There are $dbxref_n_preload GenBank accessions in the $preload_dbxref file.\n";
    my $preload_last_dbxref_id=get_last_number('dbxref_dbxref_id_seq');
    $dbh->do("COPY public.dbxref (db_id, accession, version, description) FROM '$preload_dbxref'");
    my $postload_last_dbxref_id=get_last_number('dbxref_dbxref_id_seq');
    my $dbxref_added = $postload_last_dbxref_id - $preload_last_dbxref_id;
    if ( $dbxref_added != $dbxref_n_preload) {
	die "Sorry, the number of accessions added to public.dbxref table ($dbxref_added) is not the same than the accessions in the $opt_s file.\n Check it before rerun this script.\n\n";
    } else {
	print "Added $dbxref_added accessions to public.dbxref table\n\n";
    }
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN new_dbxref_id int");
    my $sub_dbxref="subSGN_dbxref_".$load_master_table;
    $dbh->do("SELECT * INTO TEMP TABLE $sub_dbxref FROM public.dbxref WHERE dbxref_id BETWEEN $preload_last_dbxref_id AND $postload_last_dbxref_id");
    $dbh->do("UPDATE $load_master_table SET new_dbxref_id = (SELECT dbxref_id FROM $sub_dbxref WHERE $sub_dbxref.accession||$sub_dbxref.version=$load_master_table.accession||$load_master_table.version) WHERE dbxref_id IS NULL");
    $dbh->do("UPDATE $load_master_table SET new_dbxref_id = dbxref_id WHERE dbxref_id IS NOT NULL");
    $dbh->do("COPY $load_master_table TO '$master_table_file'");
    print "Copy the public.dbxref.dbxref_id to $load_master_table.\n";
    return $dbxref_added;
}

# This subroutine load the crossreferences in the sgn.est_dbxref table. It let to map SGN est id's and GenBank accessions. Like other load subroutines check if have inserted the same nunmber of entries that entries have the preload file. 

=head2 load_est_dbxref

  Usage: my $est_dbxref_added=load_est_dbxref($load_master_table, $master_table_file);
  Desc: select dbxref_id and est_id from the master table into temp table and copy it to a file. Use this file to load
        the datas into the sgn.est_dbxref table. Check that the number of new entries is the same that the entries
        of the file.
  Ret: a scalar, the number of new entries.
  Args: $load_master_table (master table name) and $master_table_file (master_table_file);
  Side_Effects: die if the check is wrong.
  Example: my $new_est_dbxref_added=load_est_dbxref($load_master_table, $master_table_file);

=cut
 
sub load_est_dbxref {
    my $load_master_table=shift;
    my $master_table_file=shift;
    my $temp_indextable=shift;
    my $index_counts = 0;
    if ($opt_i) {
	my @split_table=split(/_/, $temp_indextable);
	$index_counts=$split_table[3];
    }

    my $new_accessions = count_new_entries('accession', 'est_dbxref_id', $load_master_table);
    print "There are $new_accessions new accessions without entries in sgn.est_dbxref table associated to $opt_s file\n";
    
    my $load_est_dbxref="load_est_dbxref_".$load_master_table;
    if ($index_counts > 6) {
        $dbh->do("SELECT new_est_id, new_dbxref_id INTO TEMP TABLE $load_est_dbxref FROM $load_master_table
                  LEFT JOIN $temp_indextable ON $load_master_table.accession=$temp_indextable.accession
                  WHERE $load_master_table.est_dbxref_id IS NULL ORDER BY $temp_indextable.est_dbxref_order ASC");
    } else {
	$dbh->do("SELECT new_est_id, new_dbxref_id INTO TEMP TABLE $load_est_dbxref FROM $load_master_table 
                  WHERE $load_master_table.est_dbxref_id IS NULL");
    }

    my $est_dbxref="est_dbxref_".File::Basename::basename($master_table_file);
    my $preload_est_dbxref=create_dbloadfile($est_dbxref);
    $dbh->do("COPY $load_est_dbxref TO '$preload_est_dbxref'");
    my $est_dbxref_n_preload=`cut -f1 $preload_est_dbxref | sort -u | wc -l`;
    chomp($est_dbxref_n_preload);
    print "There are $est_dbxref_n_preload entries in the $preload_est_dbxref file.\n";
    my $preload_last_est_dbxref_id=get_last_number('sgn.est_dbxref_est_dbxref_id_seq');
    $dbh->do("COPY sgn.est_dbxref (est_id, dbxref_id) FROM '$preload_est_dbxref'");
    my $postload_last_est_dbxref_id=get_last_number('sgn.est_dbxref_est_dbxref_id_seq');
    my $est_dbxref_added = $postload_last_est_dbxref_id - $preload_last_est_dbxref_id;
    if ( $est_dbxref_added != $est_dbxref_n_preload) {
	die "Sorry, the number of entries added to sgn.est_dbxref table ($est_dbxref_added) is not the same than the entries in the $preload_est_dbxref file.\n Check it before rerun this script.\n\n";
    } else {
	print "Added $est_dbxref_added accessions to sgn.est_dbxref table\n\n";
    }
    my $sub_est_dbxref="subSGN_est_dbxref_".$load_master_table;
    $dbh->do("ALTER TABLE $load_master_table ADD COLUMN new_est_dbxref_id int");
    $dbh->do("SELECT * INTO TEMP TABLE $sub_est_dbxref FROM sgn.est_dbxref WHERE est_dbxref_id BETWEEN $preload_last_est_dbxref_id AND $postload_last_est_dbxref_id");
    $dbh->do("UPDATE $load_master_table SET new_est_dbxref_id = (SELECT est_dbxref_id FROM $sub_est_dbxref WHERE $sub_est_dbxref.est_id=$load_master_table.new_est_id) WHERE est_dbxref_id IS NULL");
    $dbh->do("UPDATE $load_master_table SET new_est_dbxref_id = est_dbxref_id WHERE est_dbxref_id IS NOT NULL");
    $dbh->do("COPY $load_master_table TO '$master_table_file'");
    print "Copy the public.dbxref.dbxref_id to $load_master_table.\n";
    return $est_dbxref_added;
}


=head2 sequence_of_loads

  Usage: my @entries_added=sequence_of_loads($load_master_table, $master_table_file);
  Desc: this subrotine coordinates the different load subroutines and get the number of entries added.
  Ret: An array with the name of the mastertable, the entries and the different new entries added to the database
  Args: $load_master_table (the name of the master table in the database) and $master_table_file (the name of the
        file where the master table is copying)
  Side_Effects: none
  Example: my @entries_added=sequence_of_loads($load_master_table, $master_table_file);

=cut

sub sequence_of_loads {
    my $load_master_table=shift;
    my $master_table_file=shift;
    my $temp_indextable=shift;
    my @entries_added;


    my $tablename=$load_master_table;
    $tablename =~ s/_mastertable$//gi;
    push @entries_added, $tablename;

    my $querycount="SELECT COUNT(*) FROM $load_master_table";
    my $sth=$dbh->prepare($querycount);
    $sth->execute();
    my ($entries_count) = $sth->fetchrow_array();
    push @entries_added, $entries_count;

    print_separator('Loading clone data');
    my $clon_add=load_clone($load_master_table, $master_table_file, $temp_indextable);
    push @entries_added, $clon_add;

    print_separator('Loading read data');
    my $read_add=load_read($load_master_table, $master_table_file, $temp_indextable);
    push @entries_added, $read_add;

    print_separator('Loading ESTs (sequences) data');
    my $est_add=load_est($load_master_table, $master_table_file, $temp_indextable);
    push @entries_added, $est_add;

    print_separator('Loading clean coordenates in QC report');
    my $qcrep_add=load_qc_report($load_master_table, $master_table_file, $temp_indextable);
    push @entries_added, $qcrep_add;

    my ($dbxref_add, $est_dbxref_add);
    if ($opt_f =~ m/^G/i) {
        print_separator('Loading accessions');
        $dbxref_add=load_dbxref($load_master_table, $master_table_file, $temp_indextable);
        push @entries_added, $dbxref_add;

        print_separator('Loading crossreferences between accessions and ESTs');
        $est_dbxref_add=load_est_dbxref($load_master_table, $master_table_file, $temp_indextable);
        push @entries_added, $est_dbxref_add;
    }
    return @entries_added;
}


=head2 split_master_table

  Usage: my @submastertables_by_library=split_master_table($temp_library, $global_master_table);
  Desc: this subroutine splits the global master table in anothers, smaller, group by library.
  Ret: An array with the submaster tables names
  Args: $temp_library (name of the temp table for the libraries) and $global_master_table (the name of the global 
        master table with all the entries)
  Side_Effects: none
  Example: my @mastertable_libs=split_master_table($temp_library, $global_master_table);

=cut

sub split_master_table {
    my $temp_library=shift;
    my $global_master_table=shift;
    my (@lib_master_table, @list_lib_sname);
    my ($lib, $lib_mastertable);
   
    print "Splicing the $global_master_table in subtables by libraries...\n\n";
    my $query = "SELECT library_shortname FROM $temp_library ORDER BY library_shortname";
    my $sth=$dbh->prepare($query);
    $sth->execute();
    while (my ($lib_shortname) = $sth->fetchrow_array()) {
        push @list_lib_sname, $lib_shortname;
    }
   
    print "@list_lib_sname\n";
    if ($opt_f =~ m/^G/i) {
       
        foreach $lib (@list_lib_sname) {
	    $lib_mastertable="lib_".lc($lib)."_mastertable";
            $lib_mastertable =~ s/\W+/_/gi;
            $dbh->do("SELECT $global_master_table.* INTO TEMP TABLE $lib_mastertable FROM $global_master_table 
                        WHERE $global_master_table.clone_lib = 
                        (SELECT clone_lib FROM $temp_library WHERE library_shortname = '$lib')");
            my $querycount="SELECT COUNT(*) FROM $lib_mastertable";
            $sth=$dbh->prepare($querycount);
            $sth->execute();
            my ($entries_count) = $sth->fetchrow_array();
            print "Library shortname:$lib\n";
            print "\tlibrary master table: $lib_mastertable ===> $entries_count entries\n";
  
            push @lib_master_table, $lib_mastertable;
	}
    } elsif ($opt_f =~ m/^P|^F/i) {
        foreach $lib (@list_lib_sname) {
	    $lib_mastertable="lib_".lc($lib)."_mastertable";
            $lib_mastertable =~ s/\W+/_/gi;
            $dbh->do("SELECT $global_master_table.* INTO TEMP TABLE $lib_mastertable FROM $global_master_table 
                        WHERE $global_master_table.library_shortname = '$lib'");
            my $querycount="SELECT COUNT(*) FROM $lib_mastertable";
            $sth=$dbh->prepare($querycount);
            $sth->execute();
            my ($entries_count) = $sth->fetchrow_array();
            print "Library shortname:$lib\n";
            print "\tlibrary master table: $lib_mastertable ===> $entries_count entries\n";

            push @lib_master_table, $lib_mastertable;
	}
    }
    return @lib_master_table;
}



=head2 print_table

  Usage: print_table ($arrayref_of_arrays);
  Desc: print a table of a array with arrayref like rows.
  Ret: none
  Args: $arrayref_of_arrays (so means, an arrayref with @rows like elements).
  Side_Effects: none
  Example: @header = ('organism', 'est_ids');
           @row1 = ('tomato', 1200);
           @row2 = ('tobacco', 320);
           @array = (\@header, \@row1, \@row2);
           print_table(\@array);

=cut

sub print_table {
    my $table=shift;
    my @table=@$table;

    my $length=0;
    my ($row, $row_ar, $data, $datalen);
    
    foreach $row_ar (@table) {
	my @row=@$row_ar;
        my $name=$row[0];
        my $name_length=length($name);
	    if ($name_length > $length) {
		$length = $name_length;
	    }
    }
    
    
    $length += 2;
    my $count=0;
    foreach $row (@table) {
        my @row=@$row;
        my $name=shift(@row);
	
        my $name_length=length($name);
        while ($name_length < $length) {
	    $name .= " ";
            $name_length=length($name);
        }
        print "|  $name|  ";
	foreach $data (@row) {
            $datalen=length($data);
	    while ($datalen < 16) {
		$data .= " ";
                $datalen=length($data);
	    }
	    print "$data|  ";   
	}
        print "\n";
        if ($count == 0) {
		my $line="-";
                $line =~ s/./-/gi;
                my $line_length=length($line);
		while ($line_length < $length) {
         	    $line .= "-";
            	    $line_length=length($line);
                }
		print "+--$line+";
                my ($line2, $linelen);
                foreach $line2 (@row) {
		    $line2 =~ s/./-/gi;
                    $linelen=length($line2);
	            while ($linelen < 16) {
		        $line2 .= "-";
                        $linelen=length($line2);
		    }
                print "--$line2+";
	    }
	    print "\n";
        }
        $count++;
    }
}
