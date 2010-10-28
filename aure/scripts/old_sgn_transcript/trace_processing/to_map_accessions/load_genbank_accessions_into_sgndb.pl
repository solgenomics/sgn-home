#!/usr/bin/perl

=head1 NAME

 load_genbank_accession_into_sgndb.pl.
 A script to load the accessions in the SGN databases.(version.2.0.).

=cut

=head1 SYPNOSIS

load_genbank_accession_into_sgndb.pl -H <dbhost> -D <dbname> -i <input_file> -o <organism_name> [-a <database_accession_name>] [-P]
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory)
      
=item -i

B<input file>           input tab file with the accessions, versions and est_ids (this is a output file of the script SGN_map_gb_accession_in_db.pl) (mandatory)

=item -o

B<organism_name>	Organism name for the accessions (mandatory)

=item -a

B<database_accession>   Database accession name (by default 'DB:GenBank_Accession')

=item -P

B<Print_load_data>      Enable the option to write in two files the load data (load_dbxref.tab and load_est_dbxref.tab) (In these files, the first column is new, new=0 the accession or est_id+dbxref_id is into the table, with new=1 is new).
           
=item -h 

B<help>                  show the help

=back

=cut

=head1 DESCRIPTION

Sometimes you can find sequences in the database but not a GenBank accessions associated them (perhaps the sequences were incorporated to the database before they was submit to GenBank). This script load the accessions and the versions of a GenBank sequences in the database in the tables public.dbxref and sgn.est_dbxref tables. Previously you need map the accessions in the database. To do it you can use the script SGN_map_gb_accession_in_db.pl. It is necesary the use of the complete path to describe the location of the files.
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

load_genbank_accession_into_sgndb.pl


=cut

use strict;
use Getopt::Std;
use File::Basename;
use CXGN::DB::InsertDBH;

our ($opt_H, $opt_D, $opt_i, $opt_o, $opt_a, $opt_P, $opt_h);
getopts("H:D:i:o:a:Ph");
if (!$opt_H && !$opt_D && !$opt_i && !$opt_o && !$opt_a && !$opt_P && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my ($inputfile, $organism_name, $dbname, $dbhost) = validate_the_input();

our $dbh = CXGN::DB::InsertDBH->new({	dbname=>$dbname ,
						dbhost=>$dbhost ,
					  	dbargs=>{ RaiseError => 1, PrintError => 1 }
					})
	or die "Failed to connect to database ($DBI::errstr)";

my $dbload_dirpath=File::Basename::dirname($opt_i);
our $dbload_dir=$dbload_dirpath."/accession_db_load_files";

my $query="SELECT organism_id FROM sgn.organism WHERE organism_name = ?";
my $sth=$dbh->prepare($query);
$sth->execute($organism_name);
my ($organism_id) = $sth->fetchrow_array();
print "CHECK1: $organism_id\n";
if (!$organism_id) {
   die ("Sorry, the organism $opt_o is not in the database (sgn.organism table).\n");
}

my $accession_dbname=$opt_a||'DB:GenBank_Accession';
my $query2="SELECT db_id FROM public.db WHERE name = ?";
$sth=$dbh->prepare($query2);
$sth->execute($accession_dbname);
my ($db_id) = $sth->fetchrow_array();
print "CHECK2: $db_id\n";
if (!$db_id) {
    die "Sorry, the accession database name ($accession_dbname) is not in the database (public.db).\n";
}
my $description="This is a $accession_dbname for a sequence from $organism_name";

print_separator('Parsing the input file');
my @input=parse_input($inputfile);
my ($element_aref, $outdbxref, $outestdbxref);
my ($count, $new_dbxref_entries, $new_est_dbxref_entries)=(0,0,0);
print_separator('Inserting the variables into the database');
if ($opt_P) {
    print "Enable the option -P: Print the variables into output files (first column -f1 is new, 0=not 1=yes)\n" ;
    my $dbload_dirpath=File::Basename::dirname($opt_i);
    my $dbload_dir=$dbload_dirpath."/accession_db_load_files";
    mkdir "$dbload_dir", 0755 or warn "Cannot make $dbload_dir directory: $!";
    my $load_dbxref_file=$dbload_dir."/load_dbxref.tab";
    my $load_est_dbxref_file=$dbload_dir."/load_est_dbxref.tab";
    open $outdbxref, '>', $load_dbxref_file || die "I can not open the output file: $load_dbxref_file.\n";
    open $outestdbxref, '>', $load_est_dbxref_file || die "I can not open the output file: $load_est_dbxref_file.\n";
}
my $total=`cut -f1 $opt_i | wc -l`;
chomp($total);

my $last_number_dbxref_beforeinsert=get_last_number($dbh, 'dbxref_dbxref_id_seq');
my $last_number_estdbxref_beforeinsert=get_last_number($dbh, 'sgn.est_dbxref_est_dbxref_id_seq');

foreach $element_aref (@input) {
   $count++;
   my @element=@$element_aref;
   my $accession=$element[0];
   my $version=$element[1];
   my $est_id=$element[2];
   print "\nInserting element ($count of $total) (accession:$accession, version:$version and est_id:$est_id\n";
   my ($dbxref_id, $dbxref_new)=insert_new_dbxref($dbh, $db_id, $accession, $version, $description);
   my ($est_dbxref_id, $est_dbxref_new)=insert_new_est_dbxref($dbh, $est_id, $dbxref_id);
   if ($opt_P) {
       print $outdbxref "$dbxref_new\t$db_id\t$accession\t$version\t$description\n";
       print $outestdbxref "$est_dbxref_new\t$est_id\t$dbxref_id\n";
   }
   $new_dbxref_entries += $dbxref_new;
   $new_est_dbxref_entries += $est_dbxref_new;
}

print_separator('REPORT:');
print "Rows into the input file: $total\n";
print "Rows processed: $count\n";
print "New dbxref_id:$new_dbxref_entries\n";
print "New est_dbxref_id:$new_est_dbxref_entries\n";

my $last_number_dbxref_afterinsert=get_last_number($dbh, 'dbxref_dbxref_id_seq');
my $last_number_estdbxref_afterinsert=get_last_number($dbh, 'sgn.est_dbxref_est_dbxref_id_seq');
my $delta_dbxref=$last_number_dbxref_afterinsert-$last_number_dbxref_beforeinsert;
my $delta_estdbxref=$last_number_estdbxref_afterinsert-$last_number_estdbxref_beforeinsert;

if ( $delta_dbxref ne $new_dbxref_entries) {
   die "The number of new entries in the public.dbxref table ($delta_dbxref) IS NOT the same that the number of counts of new inserts ($new_dbxref_entries)\n";
} elsif ($delta_estdbxref ne $new_est_dbxref_entries) {
   die "The number of new entries in the sgn.est_dbxref table ($delta_estdbxref) IS NOT the same that the number of counts of new inserts ($new_est_dbxref_entries)\n";
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
      Sometimes you can find sequences in the database but not a GenBank accessions associated them (perhaps the sequences were incorporated to the database before they was submit to GenBank). This script load the accessions and the versions of a GenBank sequences in the database in the tables public.dbxref and sgn.est_dbxref tables. Previously you need map the accessions in the database. To do it you can use the script SGN_map_gb_accession_in_db.pl. It is necesary the use of the complete path to describe the location of the files.
    
    Usage: 
      load_genbank_accessions_into_sgndb.pl -H <dbhost> -D <dbname> -i <input_file> -o <organism_name> [-a <accession_dbname>] [-P]

    Example:
      load_genbank_accessions_into_sgndb.pl -H localhost -D sandbox -i /home/aure/workdir/GenBank_sequences/
      Nicotiana_tabacum/tab_files/db_mapping_nicotiana_tabacum_db_gbseqdata.tab.tab -o 'Nicotiana tabacum'
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -i input file           input tab file with the accessions, versions and est_ids (this is a output file 
                              of the script SGN_map_gb_accession_in_db.pl) (mandatory)
      -o organism name        organism name of the sequences (with single quotes) (mandatory)
      -a accession database   accesssion database name (by default DB:GenBank_Accession)
      -P enable print         enable print into two output files the insert data (first column is new, that means new=0 for dbxref or est_dbxref that exists into the table and new=1 for new entries)
      -h this help

EOF
exit (1);
}

=head2 validate_the_input

  Usage: my @input_var=validate_the_input();
  Desc: this subroutine check if exist the mandatory input variables and count the entries in the input files. 
  Ret: six scalars. $cleanfile (cleaning coordenates file name), $chimerafile (chimera id file name), 
       $seqfile (file name of the file with all the sequences and sequence datas), $libraryfile (file
       with the libraries datas), $dbname and $dbhost.
  Args: none
  Side_Effects: die the process if the input data is not right
  Example: my ($input_file, $dbname, $dbhost)=validate_the_input();

=cut
 
sub validate_the_input {
  print "Check the variables...\n";
  if ($opt_h) {
    help ();
  }

  unless ($opt_i) {
    die ("Required argument -i <input file with accessions, versions and est_id> was not supplied.\n");
  }
  unless ($opt_o) {
    die ("Required argument -o <organism name> was not supplied.\n");
  }
    
  my $inputfile_count_columns=count_columns($opt_i);
  if ($inputfile_count_columns != 3) {
	  die "Sorry the input sequence file ($opt_i) have $inputfile_count_columns.\nThese must have 3 columns\n";
      } else {
	  print "\tColumns number in the $opt_i file ... Ok.\n";
      }
    
  print "\nCount process...\n";

  my $accessions_count_inputfile = `cut -f1 $opt_i | sort -u | wc -l`;
  chomp ($accessions_count_inputfile);
  my $est_id_count_inputfile = `cut -f3 $opt_i | sort -u | wc -l`; 
  chomp ($est_id_count_inputfile);
  my $input_file=File::Basename::basename($opt_i);
  print "\t> Accessions count for $input_file file:\t$accessions_count_inputfile\n";
  print "\t> Est_ids_count for $input_file file:\t$est_id_count_inputfile\n";

  if (!$opt_D || !$opt_H) { 
      die "\nneed -D (database name) and -H (host) options";
  }

  return ($opt_i, $opt_o, $opt_D, $opt_H);
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
    my $dbh=shift;
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

=head2 parse_input

  Usage: my @input=parse_input($input_file)
  Desc: This subroutine parse the input file and store the variables into an array (rows) of arrays (elements)
  Ret: An array of array references.
  Args: A scalar, the input file.
  Side_Effects: Die if parse a file with a column number different of 3
  Example: my @input=parse_input($input_file)

=cut

sub parse_input {
  my $input_file=shift;
  my @columns;
  my $n=0;
  open my $input, '<', "$input_file" or die "Sorry, I can not open the input file: $input_file.\n";
      while (<$input>) {
         chomp($_);
         $n++;
	 my @variables = split(/\t/,$_);
         my $column_n=scalar(@variables);
         if ($column_n != 3) {
	     die "The row $n has not 3 columns\n";
	 }
         my $variables_aref=\@variables;
	 push @columns, $variables_aref;
      }
  
  close $input;
  return @columns;
}

=head2 insert_new_dbxref

  Usage: my ($dbxref_id, $new)=insert_new_dbxref($db, $db_id, $accession, $version, $description)
  Desc: Insert new accessions in the public.dbxref table, get the new dbxref_id. If find old dbxref get the old and return $new object equal to 0.
  Ret: Two scalar, $dbxref_id and $new (1 if it is new and 0 if it is old)
  Args: Five scalars, $db (database conection object), and the input variables ($db_id, $accession, $version and $description).
  Side_Effects: None.
  Example: my ($dbxref_id, $new)=insert_new_dbxref($db, $db_id, $accession, $version, $description);

=cut

sub insert_new_dbxref {
    my $dbh=shift;
    my $db_id=shift;
    my $accession=shift;
    my $version=shift;
    my $description=shift;

    my $check = "SELECT dbxref_id FROM public.dbxref WHERE db_id=? AND accession=? AND version=?";
    my $sth=$dbh->prepare($check);
    $sth->execute($db_id, $accession, $version);
    my ($dbxref_id) = $sth->fetchrow_array();
    my $new;

    if (!$dbxref_id) {
	my $insert="INSERT INTO public.dbxref (db_id, accession, version, description) VALUES (?,?,?,?)";
        $sth=$dbh->prepare($insert);
        $sth->execute($db_id, $accession, $version, $description);
        print "INSERT INTO public.dbxref (db_id, accession, version, description) VALUES ($db_id, $accession, $version, $description)\n";
        my $query = "SELECT dbxref_id FROM public.dbxref WHERE db_id=? AND accession=? AND version=?";
        $sth=$dbh->prepare($query);
        $sth->execute($db_id, $accession, $version);
        ($dbxref_id) = $sth->fetchrow_array();
        print "RETURN dbxref_id=$dbxref_id\n";
        $new=1;
    } else {
	print "There is a dbxref_id associated with the values db_id=$db_id, accession=$accession, version=$version\n";
        $new=0;
    }
    return ($dbxref_id, $new);
}

=head2 insert_new_est_dbxref

  Usage: my ($est_dbxref_id, $new)=insert_new_est_dbxref($db, $est_id, $dbxref_id);
  Desc: Insert new accessions into the sgn.est_dbxref table, get the new est_dbxref_id. If find old dbxref get the old and return $new object equal to 0.
  Ret: Two scalar, $est_dbxref_id and $new (1 if it is new and 0 if it is old)
  Args: Three scalars, $db (database conection object), and the input variables ($est_id and $dbxref_id).
  Side_Effects: None.
  Example: my ($est_dbxref_id, $new)=insert_new_est_dbxref($db, $est_id, $dbxref_id);

=cut

sub insert_new_est_dbxref {
  my $dbh=shift;
  my $est_id=shift;
  my $dbxref_id=shift;
  my $new;

  my $first_check="SELECT est_id FROM sgn.est WHERE est_id=?";
  my $sth=$dbh->prepare($first_check);
  $sth->execute($est_id);
  my ($check_est_id) = $sth->fetchrow_array();
  if (!$check_est_id) {
      print "WARNING LOAD:The est_id=$est_id do not exist into the table sgn.est. It will not load into the est_dbxref.\n";
  } else {
      my $second_check="SELECT est_dbxref_id FROM sgn.est_dbxref WHERE est_id=? AND dbxref_id=?";
      my $sth2=$dbh->prepare($second_check);
      $sth2->execute($est_id, $dbxref_id);
      my ($est_dbxref_id) = $sth2->fetchrow_array();
      
      if (!$est_dbxref_id) {
	  my $insert="INSERT INTO sgn.est_dbxref (est_id, dbxref_id) VALUES (?,?)";
	  my $sth3=$dbh->prepare($insert);
          $sth3->execute($est_id, $dbxref_id);
	  print "INSERT INTO sgn.est_dbxref (est_id, dbxref_id) VALUES ($est_id, $dbxref_id)\n";

          my $sth2=$dbh->prepare($second_check);
          $sth2->execute($est_id, $dbxref_id);
          ($est_dbxref_id) = $sth2->fetchrow_array();
          print "RETURN: est_dbxref_id:$est_dbxref_id\n";
          $new=1;
      } else {
	  print "There is a est_dbxref_id associated with the values est_id=$est_id and dbxref_id=$dbxref_id\n";
	  $new=0;
      }
      return ($est_dbxref_id, $new);
  }
}
