#!/usr/bin/perl

=head1 NAME

 sgnpl_genbank_sgn_accession_mapping_tool.pl.
 A script to map the GenBank accessions in the database using different fields (version.2.0.).

=head1 SYPNOSIS
  
sgnpl_genbank_sgn_accession_mapping_tool.pl -H <database host> -D <database name> -d <input dataset .tab> -o <organism_name> 
[-p <disable partial match>]
    
=head2 I<Flags:>

=over

=item -H 

database host (mandatory)
      
=item -D 

database name (mandatory)
      
=item -d 

GenBank dataset to compare in tab file  (mandatory) 
      
=item -o 

organism name (with single quotes) (mandatory)
      
=item -p 

disable the partial match option for the clone name
      
=item -h 

print this help

=back

=head1 DESCRIPTION

 This script map the the GenBank accessions in a database. Get the sgn.clone.clone_name and sgn.est.est_id from the database and locus, accession, version and clone from the GenBank input dataset (in tab format). Push in different arrays and compare them. If one element match, print the accession, version and est_id in the output file and the locus, accession, version, clone, sgn.est.est_id and sgn.clone.clone_name in the check output file.

 This script also compare part of the sgn.clone.clone_name with input.clone. To do it, remove the non alphanumeric characters and compare. If they have differents length, remove the begining or the end of the longest (it do both options) and compare.

=head1 AUTHOR

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=head1 METHODS
 
sgnpl_genbank_sgn_accession_mapping_tool.pl


=cut




use strict;

use CXGN::DB::InsertDBH;
use Getopt::Std;
use File::Basename;

our ($opt_D, $opt_H, $opt_d, $opt_o, $opt_P, $opt_h);
getopts("D:H:d:o:Ph");

if (!$opt_D && !$opt_H && !$opt_d && !$opt_o && !$opt_P && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my ($input_tab, $dbname, $dbhost)=validate_input();

our $dbh = CXGN::DB::InsertDBH::connect ({ 	dbname => $dbname, 
						dbhost => $dbhost,
					  });

my $query = "SELECT organism_id FROM sgn.organism WHERE organism_name = ?";
my $sth=$dbh->prepare($query);
$sth->execute($opt_o);
my ($organism_id) = $sth->fetchrow_array();
if (!$organism_id) {
     die ("Sorry, the organism ($opt_o) is not in the $opt_D database.\n");
} else {
     print "\nOrganism check ... Ok ($opt_o => $organism_id)\n";
}

my $temp_table="temp_table_".File::Basename::basename($input_tab);
$temp_table =~ s/\W/_/gi;
my $map_table=$temp_table;
$map_table =~ s/temp_/map_/gi;

print "\nLoading data into temp table ... \n";
$dbh->do("CREATE TEMP TABLE $temp_table 
		(locus varchar(250), accession varchar(250), version varchar(250), primary_id varchar(250), 
                 description text, seq text, qscore_def text, mol_type varchar(250), organism varchar(250), 
                 cell_line varchar(250), cell_type varchar(250), clone varchar(250), clone_lib varchar(250), 
                 cultivar varchar(250), dbxref varchar(250), dev_stage varchar(250), ecotype varchar(250), 
                 environmental_sample varchar(250), lab_host varchar(250),note text, PCR_primers text, 
                 plasmid varchar(250), tissue_lib varchar(250), tissue_type varchar(250), author text, 
                 UNIQUE (accession, version))");

$dbh->do("COPY $temp_table FROM '$input_tab'");
$dbh->do("UPDATE $temp_table SET clone = 'none' WHERE clone IS NULL");


print "\nSearching correspondences between input data and sgn.clone table for organism $opt_o ...\n";
my ($count_match_arrayref, $outfile_arrayref, $checkfile)=map_sgn_clone_name($organism_id, $temp_table);
my @count_match=@$count_match_arrayref;
my @outfiles=@$outfile_arrayref;
my $global_mapping_output=File::Basename::dirname($opt_d)."/OUT_global_dbmap_".File::Basename::basename($opt_d).".tab";
system "cat @outfiles > $global_mapping_output";
my $outputfile;

my $organism_entries_q="SELECT COUNT(sgn.clone.clone_id) FROM sgn.clone 
                        JOIN sgn.library ON sgn.clone.library_id=sgn.library.library_id
                        JOIN sgn.organism ON sgn.library.organism_id=sgn.organism.organism_id 
                        WHERE organism_name = ?";
$sth=$dbh->prepare($organism_entries_q);
$sth->execute($opt_o);
my ($clone_entries) = $sth->fetchrow_array();

print "\n-------------------------------------------------------------\n";
print "MAPPING REPORT:";
print "\n-------------------------------------------------------------\n";
print "Files:";
print "\n-------------------------------------------------------------\n";
print "\tInput_dataset:\t$input_tab\n";
print "\tOutput_check_corresponences: $checkfile\n\n";
print "\tMapping_output_files:\n";
foreach $outputfile (@outfiles) {
    print "\t\t=>\t$outputfile\n";
}
print "\tGlobal_mapping_output_file: $global_mapping_output\n";
print "\n-------------------------------------------------------------\n";
print "Counts:";
print "\n-------------------------------------------------------------\n";
print "\tsgn.clone.clone_name entries: $clone_entries\n";
print "\t|\n\t|\\\n\t| \\\n";
print "\t|   +->\tComplete match with accession+version: $count_match[0]\n";
print "\t|\n\t|\\\n\t| \\\n";
print "\t|   +->\tComplete match with accession: $count_match[1]\n";
print "\t|\n\t|\\\n\t| \\\n";
print "\t|   +->\tComplete match with locus: $count_match[2]\n";
print "\t|\n\t|\\\n\t| \\\n";
print "\t|   +->\tComplete match with library name dash clone name: $count_match[3]\n";
print "\t|\n\t|\\\n\t| \\\n";
print "\t|   +->\tComplete match with clone name: $count_match[4]\n";
if ($opt_P) {
    print "\t|\n\t|\\\n\t| \\\n";
    print "\t|   +->\tPartial match with library name dash clone name: $count_match[5]\n";
    print "\t|\n\t|\\\n\t| \\\n";
    print "\t|   +->\tPartial match with clone: $count_match[6]\n";
}
print "\n-------------------------------------------------------------\n";


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
      This script search correspondences between input dataset and sgn.clone table (using sgn.clone.clone_name field).

    Usage: 
      sgnpl_genbank_sgn_accession_mapping_tool.pl -H <dbhost> -D <dbname> -d <input dataset in tab format> -o <organism_name> [- P]
    
    Example:
      sgnpl_genbank_sgn_accession_mapping_tool.pl -H localhost -D sandbox -d /home/aure/GenBank_sequences/Nicotiana_tabacum/tab_files/nicotiana_tabacum_gbdata.tab -o 'Nicotiana tabacum'

    Flags:
      -H database hostname         for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name             sandbox or cxgn etc (mandatory)
      -d input dataset             genbank file in tab format (with 25 columns) (mandatory) 
      -o organism name             organism name using in the search (mandatory)
      -P disable partial search    disable the option that do a partial match of input.clone and sgn.clone.clone_name
      -h print this help

EOF
exit (1);
}

=head2 validate_input

  Usage: my @var=validate_input();
  Desc:	check the input arguments
  Ret: Three scalars, database name ($opt_D), database host ($opt_H) and input dataset file ($opt_d)
  Args: none
  Side_Effects: die if there are something wrong during the checking (like do not exists mandatory input argument)
  Example: my ($input_tab, $dbname, $dbhost)=validate_input();

=cut

sub validate_input {
    print "\nChecking variables ...\n\n";
    if ($opt_h) {
	help();
    }

    unless ($opt_d) {
	die ("Required argument -d <dataset.tab> was not supplied.\n");
    }

    unless ($opt_o) {
	die ("Required argument -o <organism name> was not supplied.\n");
    }

    my $input_tab_count=count_columns($opt_d);
    if ($input_tab_count != 25) {
	die ("Input tab file has not 25 columns ($input_tab_count).\n");
    } else {
	print "Column number input tab check\t...\tOk.\n\n";
    }
   
    if (!$opt_D || !$opt_H) {
	die ("Required arguments -H <dbhost> or -D <dbname> were not supplied.\n");
    }
    return ($opt_d, $opt_D, $opt_H);

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

=head2 map_sgn_clone_name

  Usage: my ($count_match_arrayref, $outfiles_arrayref, $checkfile) = 
         map_sgn_clone_name($organism_id, $temp_table);
  Desc: this subroutine get the clone_name (from sgn.clone) and the est_id (from sgn.est) and push in paralel in 
        two arrays. Get the locus, accession, version and clone (from input dataset table) and push in paralel in 
	four arrays. Foreach clone_name array element compare with all the elements on the input dataset arrays. If find
	correspondence get the index of array (by count of elements) and print in outfile diferents fields.
  Ret: Three scalars, the first is an array reference of the counts of differents match (clone by clone name...), 
       the second ia an array reference of the output files and the third is the check output file name.
  Args: $organism_id and $temp_table (table with the input dataset)
  Side_Effects: none
  Example: my ($count_match_arrayref, $outfiles_arrayref, $checkfile) = 
           map_sgn_clone_name($organism_id, $temp_table);

=cut

sub map_sgn_clone_name {
    my $organism_id=shift;
    my $temp_table=shift;

    my $map_acc_ver_file=File::Basename::dirname($opt_d)."/dbmap_acc_and_vers_".File::Basename::basename($opt_d).".tab";
    my $map_accession_file=File::Basename::dirname($opt_d)."/dbmap_accession_".File::Basename::basename($opt_d).".tab";
    my $map_locus_file=File::Basename::dirname($opt_d)."/dbmap_locus_".File::Basename::basename($opt_d).".tab";
    my $map_primary_id_file=File::Basename::dirname($opt_d)."/dbmap_primary_id_".File::Basename::basename($opt_d).".tab";
    my $map_libclone_file=File::Basename::dirname($opt_d)."/dbmap_lib_d_clone_".File::Basename::basename($opt_d).".tab";
    my $map_clone_file=File::Basename::dirname($opt_d)."/dbmap_clone_".File::Basename::basename($opt_d).".tab";
   
    my $check_map_file=File::Basename::dirname($opt_d)."/check_correspondences_".File::Basename::basename($opt_d).".tab";
    
    open my $map_acc_ver_filehandle, '>', "$map_acc_ver_file";
    open my $map_accession_filehandle, '>', "$map_accession_file";
    open my $map_locus_filehandle, '>', "$map_locus_file";
    open my $map_primary_id_filehandle, '>', "$map_primary_id_file";
    open my $map_libclone_filehandle, '>', "$map_libclone_file";
    open my $map_clone_filehandle, '>', "$map_clone_file";

    open my $check_map_filehandle, '>', "$check_map_file";

    my ($map_partial_libclone_file, $map_partial_clone_file);

    if ($opt_P) {
        my $map_partial_libclone_file=File::Basename::dirname($opt_d)."/dbmap_partial_lib_d_clone_"
                                     .File::Basename::basename($opt_d).".tab";
        my $map_partial_clone_file=File::Basename::dirname($opt_d)."/dbmap_partial_clone_"
                                  .File::Basename::basename($opt_d).".tab";
        open my $map_partial_libclone_filehandle, '>', "$map_partial_libclone_file";
        open my $map_partial_clone_filehandle, '>', "$map_partial_clone_file";
    }
 


    my $query = "SELECT sgn.est.est_id, sgn.clone.clone_name FROM sgn.est JOIN sgn.seqread ON sgn.est.read_id=sgn.seqread.read_id JOIN sgn.clone ON sgn.seqread.clone_id=sgn.clone.clone_id JOIN sgn.library ON sgn.clone.library_id=sgn.library.library_id WHERE sgn.library.organism_id =?";
    my (@clone_names, @est_ids, @clone_temp, @accession_temp, @version_temp, @accession_vers_temp, @locus_temp, 
        @lib_dash_clone_temp, @primary_id_temp);
    my ($count_pmatch_accver, $count_pmatch_accession, $count_pmatch_locus, $count_pmatch_primary_id, 
        $count_pmatch_c_libclone, $count_pmatch_c_clone, $count_pmatch_p_libclone, $count_pmatch_p_clone)=
       (0, 0, 0, 0, 0, 0, 0, 0);

    print "\tLoading sgn.clone arrays ...";
    my $sth=$dbh->prepare($query);
    $sth->execute($organism_id);
    while (my ($est_id, $clone_name)=$sth->fetchrow_array()) {
	push @clone_names, $clone_name;
	push @est_ids, $est_id;
    }
    my $clonearray_elements=scalar@clone_names;	
 
    print "\t... do ($clonearray_elements elements)\n";
    print "\tLoading input arrays ...";
    my $query2 = "SELECT locus, accession, version, primary_id, clone, clone_lib FROM $temp_table";
    $sth=$dbh->prepare($query2);
    $sth->execute();
    while (my ($locus_temp, $accession_temp, $version_temp, $primary_id_temp, $clone_temp, $library_temp)=
    $sth->fetchrow_array()) {
	push @clone_temp, $clone_temp;
	push @accession_temp, $accession_temp;
	push @version_temp, $version_temp;
	my $acc_ver_temp = $accession_temp.".".$version_temp;
	push @accession_vers_temp, $acc_ver_temp;
	push @locus_temp, $locus_temp;
        my $library_dash_clone_temp= $library_temp."-".$clone_temp;
	push @lib_dash_clone_temp, $library_dash_clone_temp;
        push @primary_id_temp, $primary_id_temp;
    }
    my $accarray_elements=scalar@accession_temp;
    print "\t... do ($accarray_elements elements)\n\n";
    my $clonearray_count=0;
    my $clone_number=scalar @clone_names;
    my $sgn_clone_name;
    foreach $sgn_clone_name (@clone_names) {
	my $count_clone=$clonearray_count+1;
	my $sgn_est_id=$est_ids[$clonearray_count];
        print "\n=> Analyzing SGN-E$sgn_est_id [$sgn_clone_name] ($count_clone of $clone_number)\n";
        my $match_acc_ver=match_accession_version($sgn_clone_name, $sgn_est_id, \@accession_temp, \@version_temp,
						  $check_map_filehandle, \@accession_vers_temp, $map_acc_ver_filehandle);
	if ($match_acc_ver > 0) {
	    $count_pmatch_accver++;
	} else {
	    my $match_accession=match_only_accession($sgn_clone_name, $sgn_est_id, \@accession_temp, \@version_temp,
						     $check_map_filehandle, $map_accession_filehandle);
	    if ($match_accession > 0) {
		$count_pmatch_accession++;
	    } else {
		my $match_locus=match_locus($sgn_clone_name, $sgn_est_id, \@accession_temp, \@version_temp,
					    $check_map_filehandle, \@locus_temp, $map_locus_filehandle);
	        if ($match_locus > 0) {
		    $count_pmatch_locus++;
		} else {
		    my $match_primary_id=match_primary_id($sgn_clone_name, $sgn_est_id, \@accession_temp, \@version_temp,
						          $check_map_filehandle, \@primary_id_temp, 
                                                          $map_primary_id_filehandle);
		    if ($match_primary_id > 0) {
			$count_pmatch_primary_id++;
		    } else {
			my $match_c_libclone=complete_match_library_with_clone($sgn_clone_name, $sgn_est_id, 
                                                                               \@accession_temp, \@version_temp,
					 	                               $check_map_filehandle, 
                                                                               \@lib_dash_clone_temp, 
                                                                               $map_libclone_filehandle);
			if ($match_c_libclone > 0) {
			    $count_pmatch_c_libclone++;
			} else {
			    my $match_c_clone=complete_match_clone_name($sgn_clone_name, $sgn_est_id, \@accession_temp, 
                                                                 \@version_temp, $check_map_filehandle, 
                                                                 \@clone_temp, $map_clone_filehandle);
			    if ($match_c_clone > 0) {
				$count_pmatch_c_clone++;
			    } elsif (!$opt_P) {
				my $match_p_libclone=partial_match_library_with_clone($sgn_clone_name, $sgn_est_id, 
                                                                                      \@accession_temp, \@version_temp,
					 	                                      $check_map_filehandle, 
                                                                                      \@lib_dash_clone_temp, 
                                                                                      $map_libclone_filehandle);
			        if ($match_p_libclone > 0) {
			            $count_pmatch_p_libclone++;
			        } else {
			            my $match_p_clone=partial_match_clone($sgn_clone_name, $sgn_est_id, \@accession_temp, 
                                                                          \@version_temp, $check_map_filehandle, 
                                                                          \@clone_temp, $map_clone_filehandle);
			            if ($match_p_clone > 0) {
				        $count_pmatch_p_clone++;
			            }
				}
			    }
			}
		    }
		}
	    }
	}
	$clonearray_count++;
    }
    my @count_pmatch=($count_pmatch_accver, $count_pmatch_accession, $count_pmatch_locus, $count_pmatch_primary_id, 
	              $count_pmatch_c_libclone, $count_pmatch_c_clone, $count_pmatch_p_libclone, $count_pmatch_p_clone);
    my @output_files=($map_acc_ver_file, $map_accession_file, $map_locus_file, $map_primary_id_file, $map_libclone_file,
                      $map_clone_file);
    if ($opt_P) {
	push @output_files, $map_partial_libclone_file;
        push @output_files, $map_partial_clone_file;
    }
   
    return (\@count_pmatch, \@output_files, $check_map_file);
}

=head2 match_accession_version

  Usage: my $match=match_accession_version($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                           $check_map_filehandle, $accession_version_arrayref, $map_acc_ver_filehandle);
  Desc: Search matches between sgn.clone.clone_name and input.version+accession data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $accession_version_arrayref (array reference for the array with the input.accession"."input.versions),
        $map_acc_ver_filehandle (filehandle of the file where the script is going to write the accession, version and 
        est_id of the match)
  Side_Effects: none
  Example: my $match=match_accession_version($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                             $check_map_filehandle, $accession_version_arrayref, $map_acc_ver_filehandle);

=cut

sub match_accession_version {
    my $sgn_clone_name=shift;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $accession_version_arrayref=shift;
    my $map_acc_ver_filehandle=shift;

    my @acc_ver=@$accession_version_arrayref;
    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $acc_ver;


    foreach $acc_ver (@acc_ver) {
	if ($sgn_clone_name eq $acc_ver) {
	    print "\tComplete match accession+version <=> sgn.clone.clone_name: $acc_ver = $sgn_clone_name\n";
            print $map_acc_ver_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
            print $check_map_filehandle "$acc_ver\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch ACCESSION+VERSION\n";
	    $possitive_match++;
	}
	$array_count++;
    }
    if ($possitive_match eq 0) {
	print "\tNone match with accession+version\n";
    }
    return $possitive_match;
}

=head2 match_only_accession

  Usage: my $match=match_only_accession($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                        $check_map_filehandle, $accession_arrayref, $map_accession_filehandle);
  Desc: Search matches between sgn.clone.clone_name and input.accession data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $accession_arrayref (array reference for the array with the input.accession),
        $map_accession_filehandle (filehandle of the file where the script is going to write the accession, version and 
        est_id of the match)
  Side_Effects: none
  Example: my $match=match_only_accession($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                          $check_map_filehandle, $accession_arrayref, $map_accession_filehandle);

=cut

sub match_only_accession {    
    my $sgn_clone_name=shift;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $map_accession_filehandle=shift;

    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $accession;


    foreach $accession (@accession) {
	if ($sgn_clone_name eq $accession) {
	    print "\tComplete match accession <=> sgn.clone.clone_name: $accession = $sgn_clone_name\n";
            print $map_accession_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
            print $check_map_filehandle "$accession\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch ACCESSION\n";
	    $possitive_match++;
	}
	$array_count++;
    } 
    if ($possitive_match eq 0) {
	print "\tNone match with accession\n";
    }
    return $possitive_match;
}

=head2 match_locus

  Usage: my $match=match_locus($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                               $check_map_filehandle, $locus_arrayref, $map_locus_filehandle);
  Desc: Search matches between sgn.clone.clone_name and input.locus data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $locus_arrayref (array reference for the array with the input.locus),
        $map_locus_filehandle (filehandle of the file where the script is going to write the accession, version and 
        est_id of the match)
  Side_Effects: none
  Example: my $match=match_locus($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                 $check_map_filehandle, $locus_arrayref, $map_locus_filehandle);

=cut

sub match_locus {
    my $sgn_clone_name=shift;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $locus_arrayref=shift;
    my $map_locus_filehandle=shift;

    my @locus=@$locus_arrayref;
    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $locus;


    foreach $locus (@locus) {
	if ($sgn_clone_name eq $locus) {
	    print "\tComplete match locus <=> sgn.clone.clone_name: $locus = $sgn_clone_name\n";
            print $map_locus_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
            print $check_map_filehandle "$locus\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch LOCUS\n";
	    $possitive_match++;
	}
	$array_count++;
    }
    if ($possitive_match eq 0) {
	print "\tNone match with locus\n";
    }
    return $possitive_match;
}

=head2 match_primary_id

  Usage: my $match=match_primary_id($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                    $check_map_filehandle, $primary_id_arrayref, $map_primary_id_filehandle);
  Desc: Search matches between sgn.clone.clone_name and input.primary_id data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $primary_id_arrayref (array reference for the array with the input.primary_id),
        $map_primary_id_filehandle (filehandle of the file where the script is going to write the accession, version and 
        est_id of the match)
  Side_Effects: none
  Example: my $match=match_primary_id($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                      $check_map_filehandle, $primary_id_arrayref, $map_primary_id_filehandle);

=cut

sub match_primary_id {
    my $sgn_clone_name=shift;
    $sgn_clone_name =~ s/gi//gi;
    $sgn_clone_name =~ s/\W//gi;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $primary_id_arrayref=shift;
    my $map_primary_id_filehandle=shift;

    my @primary_id=@$primary_id_arrayref;
    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $primary_id;


    foreach $primary_id (@primary_id) {
	if ($sgn_clone_name eq $primary_id) {
	    print "\tComplete match primary_id <=> sgn.clone.clone_name: $primary_id = $sgn_clone_name\n";
            print $map_primary_id_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
            print $check_map_filehandle "$primary_id\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch PRIMARY_ID\n";
	    $possitive_match++;
	}
	$array_count++;
    }
    if ($possitive_match eq 0) {
	print "\tNone match with primary id\n";
    }
    return $possitive_match;
}

=head2 complete_match_library_with_clone

  Usage: my $match=complete_match_library_with_clone($sgn_clone_name, $sgn_est_id, $accession_arrayref, 
                                                     $version_arrayref, $check_map_filehandle, 
                                                     $lib_dash_clone_arrayref, $map_lib_dash_clone_filehandle);
  Desc: Search matches between sgn.clone.clone_name and input.library-input.clone data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $lib_dash_clone_arrayref (array reference for the array with the input.library."-".input.clone),
        $map_lib_dash_clone_filehandle (filehandle of the file where the script is going to write the accession, 
        version and est_id of the match)
  Side_Effects: none
  Example: my $match=complete_match_library_with_clone($sgn_clone_name, $sgn_est_id, $accession_arrayref, 
                                                       $version_arrayref, $check_map_filehandle, 
                                                       $lib_dash_clone_arrayref, $map_lib_dash_clone_filehandle);

=cut

sub complete_match_library_with_clone {
    my $sgn_clone_name=shift;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $lib_dash_clone_arrayref=shift;
    my $map_lib_dash_clone_filehandle=shift;

    my @lib_dash_clone=@$lib_dash_clone_arrayref;
    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $lib_dash_clone;


    foreach $lib_dash_clone (@lib_dash_clone) {
	if ($sgn_clone_name eq $lib_dash_clone) {
	    print "\tComplete match library-clone <=> sgn.clone.clone_name: $lib_dash_clone = $sgn_clone_name\n";
            print $map_lib_dash_clone_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
            print $check_map_filehandle "$lib_dash_clone\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch LIBRARY+CLONE\n";
	    $possitive_match++;
	}
	$array_count++;
    }
    if ($possitive_match eq 0) {
	print "\tNone match with library-clone\n";
    }
    return $possitive_match;
}

=head2 complete_match_clone_name

  Usage: my $match=complete_match_clone_name($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                             $check_map_filehandle, $clone_arrayref, $map_compl_clone_filehandle);
  Desc: Search matches between sgn.clone.clone_name and input.clone data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $clone_arrayref (array reference for the array with the input.clone),
        $map_compl_clone_filehandle (filehandle of the file where the script is going to write the accession, 
        version and est_id of the match)
  Side_Effects: none
  Example: my $match=complete_match_clone_name($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                               $check_map_filehandle, $clone_arrayref, $map_comp_clone_filehandle);

=cut

sub complete_match_clone_name {
    my $sgn_clone_name=shift;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $clone_arrayref=shift;
    my $map_compl_clone_filehandle=shift;

    my @clone=@$clone_arrayref;
    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $clone;


    foreach $clone (@clone) {
	if ($sgn_clone_name eq $clone) {
	    print "\tComplete match clone_name <=> sgn.clone.clone_name: $clone = $sgn_clone_name\n";
            print $map_compl_clone_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
            print $check_map_filehandle "$clone\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch CLONE_NAME (complete)\n";
	    $possitive_match++;
	}
	$array_count++;
    }
    if ($possitive_match eq 0) {
	print "\tNone match with clone name\n";
    }
    return $possitive_match;
}

=head2 partial_match_library_with_clone

  Usage: my $match=partial_match_library_with_clone($sgn_clone_name, $sgn_est_id, $accession_arrayref, 
                                                    $version_arrayref, $check_map_filehandle, 
                                                    $lib_dash_clone_arrayref, $map_partlib_dash_clone_filehandle);
  Desc: Search PARTIAL matches between sgn.clone.clone_name and input.library-input.clone data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $lib_dash_clone_arrayref (array reference for the array with the input.library."-".input.clone),
        $map_partlib_dash_clone_filehandle (filehandle of the file where the script is going to write the accession, 
        version and est_id of the match)
  Side_Effects: none
  Example: my $match=partial_match_library_with_clone($sgn_clone_name, $sgn_est_id, $accession_arrayref, 
                                                      $version_arrayref, $check_map_filehandle, 
                                                      $lib_dash_clone_arrayref, $map_partlib_dash_clone_filehandle);

=cut


sub partial_match_library_with_clone {
    my $sgn_clone_name=shift;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $lib_dash_clone_arrayref=shift;
    my $map_partlib_clone_filehandle=shift;

    my @lib_dash_clone=@$lib_dash_clone_arrayref;
    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $lib_dash_clone;
    my $best_match = 0;
    my ($best_partial_lib_clone, $best_partial_sgn_clone, $best_array_count);

    foreach $lib_dash_clone (@lib_dash_clone) {
        my ($partial_lib_clone, $partial_sgn_clone, $match)=compare_partial_strings($lib_dash_clone, $sgn_clone_name);
	if ($match == 4) {
	     $possitive_match++;
	     $best_match=3;
	 } elsif ($match == 3 && $best_match < 3) {
	     $possitive_match++;
	     $best_match=2;
	 } elsif ($match == 2 && $best_match < 2) {
             $possitive_match++;
             $best_match=1;
	 }
	 if ($best_match < $match) {
	     print "\tPartial match library-clone <=> sgn.clone.clone_name: $partial_lib_clone = $partial_sgn_clone\n";
             print $map_partlib_clone_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
             print $check_map_filehandle "$lib_dash_clone\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch PARTIAL LIBRARY+CLONE\n";
	 }
         $array_count++;
       
    }
    if ($possitive_match == 0) {
	print "\tNone partial match with library-clone\n";
    } elsif ($possitive_match > 1) {
	print "\tAttention, the sgn_clone: $sgn_clone_name (SGN-E$sgn_est_id) match with more than one partial library+clone.\n";
    }
    return $possitive_match;
}

=head2 partial_match_clone

  Usage: my $match=partial_match_clone($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                       $check_map_filehandle, $clone_arrayref, $map_partial_clone_filehandle);
  Desc: Search PARTIAL matches between sgn.clone.clone_name and input.clone data and write them in a file
  Ret: A scalar (the number of matches for a sgn.clone.clone_name)
  Args: $sgn_clone_name (sgn.clone.clone_name), $sgn_est_id, (sgn.est.est_id associated to clone_name),
        $accession_arrayref (array reference for the array with the input.accessions),
        $version_arrayref (array reference for the array with the input.versions), 
        $check_map_filehandle (filehandle of the file where the script is going to write the matches), 
        $clone_arrayref (array reference for the array with the input.clone),
        $map_partial_clone_filehandle (filehandle of the file where the script is going to write the accession, 
        version and est_id of the match)
  Side_Effects: none
  Example: my $match=partial_match_clone($sgn_clone_name, $sgn_est_id, $accession_arrayref, $version_arrayref, 
                                         $check_map_filehandle, $clone_arrayref, $map_partial_clone_filehandle);

=cut

sub partial_match_clone {
    my $sgn_clone_name=shift;
    my $sgn_est_id=shift;
    my $accession_arrayref=shift;
    my $version_arrayref=shift;
    
    my $check_map_filehandle=shift;
    my $clone_arrayref=shift;
    my $map_partial_clone_filehandle=shift;

    my @clone=@$clone_arrayref;
    my @accession=@$accession_arrayref;
    my @version=@$version_arrayref;

    my ($array_count, $possitive_match)= (0, 0);
    my $clone;
    my $best_match = 0;
    my ($best_partial_clone, $best_partial_sgn_clone, $best_array_count);

    foreach $clone (@clone) {
        my ($partial_clone, $partial_sgn_clone, $match)=compare_partial_strings($clone, $sgn_clone_name);
	if ($match == 4) {
	     $possitive_match++;
	     $best_match=3;
	 } elsif ($match == 3 && $best_match < 3) {
	     $possitive_match++;
	     $best_match=2;
	 } elsif ($match == 2 && $best_match < 2) {
             $possitive_match++;
             $best_match=1;
	 }
	 if ($best_match < $match) {
	     print "\tPartial match clone <=> sgn.clone.clone_name: $partial_clone = $partial_sgn_clone\n";
             print $map_partial_clone_filehandle "$accession[$array_count]\t$version[$array_count]\t$sgn_est_id\n";
             print $check_map_filehandle "$clone\t$sgn_clone_name\tSGN-E$sgn_est_id\tmatch PARTIAL CLONE\n";
	 }
         $array_count++;
       
    }
    if ($possitive_match == 0) {
	print "\tNone partial match with clone\n";
    } elsif ($possitive_match > 1) {
	print "\tAttention, the sgn_clone: $sgn_clone_name (SGN-E$sgn_est_id) match with more than one partial clone.\n";
    }
    return $possitive_match;
}

=head2 compare_partial_strings

  Usage: my ($final_string_A, $final_string_B, $match_type)=partial_process($string_A, $string_B);
  Desc: Search partial matches between two strings (calculate different substrings and compare them)
  Ret: Three scalar, $final_string_A and $final_string_B (substring that match) and $match_type (an integer that 
       describe the type of partial match)
  Args: $string_A and $string_B (strings to compare)
  Side_Effects: none
  Example: my ($sub_sgnclone, $sub_inputclone, $part_match_type)=partial_process($sgnclone, $inputclone);

=cut

sub compare_partial_strings {
    my $string_a=shift;
    my $string_b=shift;
    $string_a =~ s/\W//gi;
    $string_b =~ s/\W//gi;
    $string_a =~ s/_//gi;
    $string_b =~ s/_//gi;

    my $length_a=length($string_a);
    my $length_b=length($string_b);
    my $delta_ab=$length_a - $length_b;
    my ($match, $final_a, $final_b);
    
    if ($length_a > 4 && $length_b > 4) {    
        if ($delta_ab == 0) {
	    my $compare_a=$string_a;
            my $compare_b=$string_b;
	    if ($compare_a eq $compare_b) {
	        $match=4;
                $final_a=$compare_a;
                $final_b=$compare_b;
 	    }
        } elsif ($delta_ab > 0) {
	    my $compare_a_begin=substr $string_a, 0, $length_b;
            my $compare_a_end=substr $string_a, $delta_ab, $length_b;
            my $compare_b_short=$string_b;        

	    if ($compare_a_begin eq $compare_b_short) {
	        $match=3;
                $final_a=$compare_a_begin;
                $final_b=$compare_b_short;
	    } elsif ($compare_a_end eq $compare_b_short) {
	        $match=3;
	        $final_a=$compare_a_end;
                $final_b=$compare_b_short;
	    }
        } elsif ($delta_ab < 0) {
	    my $compare_b_begin=substr $string_b, 0, $length_b;
            my $opossite_delta_ab=$length_b - $length_a;
            my $compare_b_end=substr $string_b, $opossite_delta_ab, $length_b;
            my $compare_a_short=$string_a;
	    if ($compare_b_begin eq $compare_a_short) {
	        $match=2;
                $final_a=$compare_a_short;
                $final_b=$compare_b_begin;
	    } elsif ($compare_b_end eq $compare_a_short) {
	        $match=2;
	        $final_a=$compare_a_short;
                $final_b=$compare_b_begin;
	    }
        } else {
	    $match=0;
        }
    } else {
        $match=0;
    }
    return ($final_a, $final_b, $match);
}

