#!/usr/bin/perl

=head1 NAME

sgnpl_create_library_file.pl.
A script to add library information to a tab file (with 3 or 4 columns) (version.2.0.).

=head1 SYPNOSIS
  
sgnpl_create_library_file.pl <input tab file> -H <dbhost> -D <dbname>

=head2 I<tags:>

=over 

=item -H 

database host

=item -D

database name
      
=item -h 

print this help

=back

=head1 DESCRIPTION

This script make a library tab file. Ask to user some information about the library using STDIN, check if exits in the database and write a tab file with them. Also add a new library_shortname to the sequence tab file.


=head1 AUTHOR

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=head1 METHODS
 
sgnpl_create_library_file.pl


=cut

use strict;
use CXGN::DB::InsertDBH;
use File::Basename;
use Getopt::Std;


my $file=shift(@ARGV);
if ($file =~ m/-h/) {
	help()
}

our ($opt_H, $opt_D, $opt_f);
getopts("H:D:f:");

if (!$opt_H && !$opt_D && !$opt_f && !$file) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my ($dbhost, $dbname, $filesource, $filename)=check_input($file);

my $dbh = CXGN::DB::InsertDBH::connect ({ 	dbname => $dbname, 
						dbhost => $dbhost,
					  });




my $preoutlibrary = File::Basename::basename($filename);
my $pathlibrary = File::Basename::dirname($filename);
$preoutlibrary =~ s/.tab$//g;
my $outlibrary = $pathlibrary."/library_".$preoutlibrary.".tab";
open my $outlib_handle, '>', "$outlibrary"
	or die "Can\'t open $outlibrary for output:$!";
chmod 0666, "$outlibrary";

my $outtab = $filename;
$outtab =~ s/.tab$//g;
$outtab .= ".with_lib.tab";
open my $outtab_handle, '>', "$outtab"
        or die "Can\'t open $outtab for output:$!";
chmod 0666, "$outtab";

my $temp_table =File::Basename::basename($filename);
$temp_table =~ s/\W/_/gi;

my ($type, $type_n, $lib_field);
my @library;

while ($type ne 'cDNA library' || $type ne 'virtual library') {
	$type = &get_library_field('type (cDNA library|virtual library)', 'none');
	if ($type eq 'cDNA library') {
	last;
	}elsif ($type eq 'virtual library'){
	last;
	}else{
	print "Sorry, choose an option (cDNA library or virtual library):\n";
	}
}
if ($type eq 'cDNA library') {
	$type_n=3;
        push @library, $type_n;
} elsif ($type eq 'virtual library') {
	$type_n=7;
        push @library, $type_n;
}

my $library_name = &write_field ('library_name', 'none', 80);
push @library, $library_name;
my $library_shortname =&write_field ('library_shortname', 'yes', 16);
push @library, $library_shortname;
my $authors =&write_field ('authors', 'none', 'none');
push @library, $authors;

  my $organism = &get_library_field('organism', 'none');
  my $query2 = "SELECT organism_id FROM sgn.organism WHERE organism_name =?";
  my $sth=$dbh->prepare($query2);
  $sth->execute($organism);
  my ($organism_id)=$sth->fetchrow_array();
  while (!$organism_id) { 
	print "Sorry, the organism: $organism is not in the sgn.organism table\n";
        print "Do you want abort the process (yes|no, default no)?\n";
        if ( <STDIN> =~ m/^y(es)/i) {
             die "Process aborted.\n";
        } else {
             $organism = &get_library_field('organism', 'none');
             $query2 = "SELECT organism_id FROM sgn.organism WHERE organism_name =?";
             $sth=$dbh->prepare($query2);
             $sth->execute($organism);
             ($organism_id)=$sth->fetchrow_array();
        }
  }
  push @library, $organism;

my $cultivar =&write_field ('cultivar', 'none', 250);
my $tissue =&write_field ('tissue', 'none', 250);
my $dev_stage =&write_field ('development_stage', 'none', 250);
my $cloning_host =&write_field ('cloning_host', 'none', 250);
my $vector =&write_field ('vector', 'none', 250);
my $comments = &write_field ('comments', 'none', 'none');
push @library, ($cultivar, $tissue, $dev_stage, $cloning_host, $vector, $comments);

$dbh->do("CREATE TEMP TABLE $temp_table (clone_name varchar(250), seq text, qscore text, basecall text,
	  library_shortname varchar(16), organism varchar(250), UNIQUE (clone_name));");

if ($opt_f eq 'P') {
     $dbh->do("COPY $temp_table (clone_name, seq, qscore, basecall) FROM '$file'");
} elsif ($opt_f eq 'F') {
     $dbh->do("ALTER TABLE $temp_table DROP COLUMN basecall");
     $dbh->do("COPY $temp_table (clone_name, seq, qscore) FROM '$file'");
}
      
my $update="UPDATE $temp_table SET library_shortname = '$library_shortname', organism = ?";
my $sth=$dbh->prepare($update);
$sth->execute($organism);

$dbh->do("COPY $temp_table TO '$outtab'");
print "\nCopy the sequences datas, library shorname and organism_id to $outtab file.\n";


foreach $lib_field (@library) {
       print $outlib_handle "$lib_field\t";
}
print $outlib_handle "0";
print "\nCopy the library datas into the file: $outlibrary.\n";


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
      sgnpl_create_library_file.pl is a script that write over a tab file (with 3 or 4 columns, like fasta or phred files) another column with the short_name of a library. Also write the datas of the library in another file using STDIN forms. Check in the database some questions, like organism or library name.  
    
    Usage: 
      sgnpl_create_library_file.pl <file> -H <dbhost> -D <dbname> -f <source> 

    Example:
      sgnpl_create_library_file.pl /home/aure/GenBank_sequences/Nicotiana_benthamiana/nb_EST/SAL_AGN/fastaqual.tab -H localhost
      -D sandbox -f F
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -f tab file	      if the tab file comes from Phred (4 columns) or Freefasta (3 columns)
      -h this help

EOF
exit (1);
}

=head2 check_input

  Usage: my ($dbhost, $dbname, $filesource, $filename)=check_input($filename);
  Desc: check the inputs and return them.
  Ret: $dbhost (database host), $dbname (database name), $filesource (the kind of file, so means if the file has 3 o 
       4 columns)
  Args: $filename.
  Side_Effects: if something is wrong, die
  Example: my ($dbhost, $dbname, $filesource, $filename)=check_input($filename);

=cut

sub check_input {
  my $filename= shift;
  my $columns_number=count_columns($filename);
  
  if ($opt_f eq 'P') {
      print "Opt_f ... Ok.\n";
  } elsif ($opt_f eq 'F') {
      print "Opt_f ... Ok.\n";
  } else {
      die "Sorry, -f only can get P or F values.\n";
  }
  
  if ($opt_f =~ m/^P$/ && $columns_number eq 4) {
      print "Opt_f (P) and column file count (4) is right.\n";
  } elsif ($opt_f =~ m/^F$/ && $columns_number eq 3) {
      print "Opt_f (F) and column file count (3) is right.\n";
  } else {
      die "Sorry the Opt_f and the number of columns of the tab file ($columns_number) is wrong!!!.\nRemember, -f P (4 columns) and -f F (3 columns).\n";
  }

  if (!$opt_H || !$opt_D) {
      die "Sorry, you need specify database host (-H) and database name (-D).\n";
  }
  return ($opt_H, $opt_D, $opt_f, $filename);

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

=head2 get_library_field

  Usage: my $field_data=get_library_field($field, $field_constraint_length);
  Desc: Get data for specify field using <STDIN>. Check if the data has more characters than the length constraint.
  Ret: A scalar, the data field.
  Args: $field (field) and $field_constraint_length (maximun length for the data in the specified field)
  Side_Effects: none
  Example: my $library_name=get_library_field('name', 80);

=cut

sub get_library_field {
	my $field = shift;
        my $constraint_length=shift;
	print "Introduce $field ";
	if ($constraint_length ne 'none') {
	    print "with less than $constraint_length characters";
	}
        print ":\n";
	chomp (my $data = <STDIN>);
	my $data_length=length($data);
        if ($constraint_length ne 'none') {
	    while ($data_length > $constraint_length) {
	    	print "The data is not correct.\nIntroduce $field with less than $constraint_length characters:\n";
	        chomp ($data = <STDIN>);
	    }
	}
	return $data; 
}

=head2 check_library_field

  Usage: my $check=check_library_field($field, $data);
  Desc: Check in the database if exists the data in the specified field in the sgn.library table.
  Ret: A scalar, library_id if exists the data.
  Args: $field and $data.
  Side_Effects: none
  Example: my $cultivar_check=check_library_field('cultivar', 'Cherry');

=cut

sub check_library_field {
	my $field = shift;
	my $data = shift;
	my $query = "SELECT library_id FROM sgn.library WHERE $field =?";
	my $sth=$dbh->prepare($query);
  	$sth->execute($data);
  	my ($field_check)=$sth->fetchrow_array();
	if (!$field_check) {
	    $field_check=0;
	}
	return $field_check; 
}

=head2 exists_warning

  Usage: my $warning=exists_warning($field_check);
  Desc: Print warning when exists the data in the sgn.library table and ask to user if want change the data.
  Ret: A scalar, 0 (if want change the STDIN data) or 1 (if wan not change).
  Args: $field_check (if this variable is not equal to 0, the data exists in the database)
  Side_Effects: none
  Example: my $warning=exists_warning($cultivar_check);

=cut

sub exists_warning {
    my $field_check=shift;
    my $warning;
    if ($field_check ne 0) {
	print "\n------------------------------------------------------------\n";
	print "Warning: The data exists with the library_id = $field_check.\n";
	print "Do you want change it (yes|no, default no)?";
	chomp (my $dec = <STDIN>);
	if ($dec =~ m/^y(es)/i ) {
	    $warning = 0;
	}
	print "\n------------------------------------------------------------\n";
    } else {
	$warning = 1;
    }
    return $warning;
}
	
=head2 write_field

  Usage: my $data=write_field($field, $constraint_unique, $constraint_length);
  Desc: Get the data using the get_library_field subroutine and check if exist in the sgn.library table using the 
        check_library_field subroutine. If the warning is 0, reexecute these subroutines.
  Ret: A scalar, the field data.
  Args: $field (field of the sgn.library), $constraint_unique (if the data asociated with the field must be unique ('yes')
        or ('none') ) and $constraint_length (if the data has a length constraint).
  Side_Effects: none
  Example: my $cultivar=write_field('cultivar', 'none', 250); 

=cut

sub write_field {
 	my $field = shift;
        my $constraint_unique = shift;
        my $constraint_length = shift;

	my $data = &get_library_field($field, $constraint_length);
	if ($constraint_unique eq 'yes') {
	    my $check = &check_library_field($field, $data);
	    my $warning=&exists_warning($check);

	    while ($warning eq 0) {
	    	$data = &get_library_field($field);
                $check = &check_library_field($field, $data);
                $warning=&exists_warning($check);
	    }
	}

	return $data;
}
