#!/usr/bin/perl

=head1 NAME

sgnpl_genbank_format_processing.pl.
A script to process information from GenBank files (version.2.0.).

=head1 SYPNOSIS
  
sgnpl_genbank_format_processing.pl <input dir with GenBank files> -H <dbhost> -D <dbname> [-m] <GenBank_mol_type> [-c] <cluster_with_entries_per_lib> [-d] [-r]

=head1 DESCRIPTION

This script get the data of GenBank files (.gb) in a genbank format and give us the same data in a different formats:

=over

=item A-
 
Tab format (with the name db_genbank_output.tab), with the follow fields: Locus, accession, version, description, seq, qscore (added by the script, by default 15 for each nt), mol_type, organism, cell_line, cell_type, clone_name, clone_lib, cultivar, dbxref, dev_stage, ecotype, environmental_sample, lab_host, note, PCR_primers, plasmid, tissue_lib, tissue_type, author.

While the script extract the information of genbank file, replace the mRNA mol type for EST when find EST in the keyword field (the mol_type for EST is mRNA in the GenBank file).

For mRNA sequences without clone_lib (usually the core nucleotides download from the GenBank) update the clone_lib like 'virtual mRNA library for genbank sequences from <organism>'.

=item B-

Fasta format

=item C- 

Libraries data in tab format. Give us two files.

=over

=item > 

Crude libraries_output.tab, this script analyze clone_lib (library_name), author, organism, cultivar, tissue_type, dev_stage, lab_host, plasmid, note and number of accessions associated to each library and give a file with one library in each row.

=item > 

Processed libraries_output.tab, over the crude file, cluster the libraries with less than 100 sequences in the 'virtual mRNA library for genbank sequences from <organism>'. Also add the field type and update it with 3 for common libraries and 17 for virtual library. Also add library_shortname with the follow code gbl_o<organism_id from the database used>.<substring 3 char of tissue type><library id for the load, not in the sgn.library table> (the library_shortname must be UNIQUE when you load in the database) for common libraries and vgbl_o<organism_id>.mRNA.

=back

=back

=head1 ADITIONAL INFORMATION

The mol_type field is a mandatory in the GenBank format fields. The options could be:

=over

=item *

'genomic_DNA', 'genomic_RNA', 'mRNA', 'preRNA', 'tRNA', 'rRNA', 'snoRNA', 'snRNA', 'scRNA', 'tmRNA', 'viralcRNA', 'other_RNA', 'other DNA', 'unassigned_DNA' and 'unassigned_RNA'.

=back

This script change the mol_type mRNA for EST when find EST in the keywords, so EST also is a mol_type for this script.
 
=head1 AUTHOR

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=head1 METHODS
 
sgnpl_genbank_format_processing.pl


=cut

use strict;

use Bio::SeqIO;
use CXGN::DB::InsertDBH;
use Getopt::Std;

					 
my $dir=shift(@ARGV);
$dir =~ s/\/$//;

if ($dir eq '-h') {
    help();
}

our ($opt_m, $opt_c, $opt_D, $opt_H, $opt_d, $opt_r, $opt_h);
getopts("m:c:D:H:drh");

if (!$opt_m && !$opt_c && !$opt_D && !$opt_H && !$opt_d && !$opt_r && !$opt_h && !$dir) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my $check_dir=`ls $dir | grep '.gb' | wc -l`;
chomp($check_dir);
if ($check_dir eq 0) {
    die ("Sorry, the $dir dir have not .gb files");
}

our $dbh = CXGN::DB::InsertDBH::connect ({ 	dbname => $opt_D, 
						dbhost => $opt_H,
					  });

my $mol_type = "(\'".$opt_m."\')";
if ($opt_m) {
     $mol_type =~ s/_/ /gi;
     $mol_type =~ s/,/\',\'/gi;
} else {
     $mol_type = "('EST','mRNA')";
}

my $cluster_lib = $opt_c || 100;

my $presufix1 = lc $dir;
$presufix1 =~ m/(\w+)$/i;
our $sufix = "$1";
 
opendir DH, $dir or die "Cannot open $dir to process: $!";

my $file=shift;
my (@tab_files, @source_files);
my ($tab_file, $source_file);
print_separator('Make output directories:');
my $outdir_fasta= create_dir($dir, 'fasta_files');
my $outdir_tab= create_dir($dir, 'tab_files');
my $join_command = "cat ";
print_separator('Process genbanks files');
foreach $file (readdir DH) {
  if ($file =~ /\.gb$/) {
    my $filecomp=$dir."/".$file;
    push @source_files, $filecomp;
    my $sub_entries_count = `grep 'ACCESSION' $filecomp | wc -l`;
    chomp($sub_entries_count);
    print "The $file file has $sub_entries_count accessions.\n";
    print "Processing the $file file...\n";
    my $sub_outfile_tab=process_genbank_file($dir, $file, $outdir_tab, $mol_type);
    my $sub_tab_entries_count = `cut -f1 $sub_outfile_tab | sort -u | wc -l`;
    chomp($sub_tab_entries_count);
    print "The $file file has been processed.\nThe $sub_outfile_tab tab file has $sub_tab_entries_count entries.\n";
    $join_command .= "$sub_outfile_tab ";
    push @tab_files, $sub_outfile_tab;
  }
}

my $output_global_tab = $outdir_tab."/".$sufix."_db_gbseqdata.tab";
$join_command .= "> $output_global_tab";
system "$join_command";
chmod 0666, "$output_global_tab";
my $gb_entries_count=`cut -f1 $output_global_tab | sort -u | wc -l`;
chomp($gb_entries_count);
print "The global file ($output_global_tab) has $gb_entries_count entries.\n";

print_separator('Load data in temp table');
my $temp_table="temp_gb_".$sufix;
print "Create $temp_table temp table in the database.\n";
$dbh->do("CREATE TEMP TABLE $temp_table 
		(locus varchar(250), accession varchar(250), version varchar(250), primary_id varchar(250), 
                 description text, seq text, qscore_def text, mol_type varchar(250), organism varchar(250), 
                 cell_line varchar(250), cell_type varchar(250), clone varchar(250), clone_lib varchar(250), 
                 cultivar varchar(250), dbxref varchar(250), dev_stage varchar(250), ecotype varchar(250),
                 environmental_sample varchar(250), lab_host varchar(250),note text, PCR_primers text, 
                 plasmid varchar(250), tissue_lib varchar(250), tissue_type varchar(250), author text, 
                 UNIQUE (accession, version))");
print "Loading data from the $output_global_tab in the database temp table $temp_table...\n";
$dbh->do("COPY $temp_table FROM '$output_global_tab'");
print "\t\t\t\t...Load complete.\nChecking...\n";
my $query="SELECT COUNT(*) FROM $temp_table";
my $sth=$dbh->prepare($query);
$sth->execute();
my ($count_entries_temp_table)=$sth->fetchrow_array;
print "There are $count_entries_temp_table entries in the temp_table."; 
if ($count_entries_temp_table != $gb_entries_count) {
    print "\n**********\n*WARNING:*\n**********\n";
    print "The entries in the $temp_table temp table are not the same than in the $gb_entries_count tab file.\n";
    possible_abort_process();
}
print_separator('Make the fasta file');
print "Processing the tab file to fasta file...\n";
my $output_global_fasta=tab_to_fasta($outdir_fasta, $temp_table);
my $count_fasta=`grep '>' $output_global_fasta | wc -l`;
chomp($count_fasta);
print "\t...Do\n";
print "The fasta file ($output_global_fasta) has $count_fasta entries.\n";
if ($count_fasta != $gb_entries_count) {
    print "\n**********\n*WARNING:*\n**********\n";
    print "The entries in the $output_global_fasta fasta file are not the same than in the $gb_entries_count tab file.\n";
    possible_abort_process();
}

print_separator('Make the libraries files');
print "Processing the tab file to get the libraries...\n";
my ($crude_lib_file, $processed_lib_file)=create_library_tab($outdir_tab, $temp_table, $cluster_lib);
my $crude_lib_count= `cut -f1 $crude_lib_file | wc -l`;
$dbh->do("COPY $temp_table TO '$output_global_tab'");
print "Update the $output_global_tab file with the processed library data.\n";
chomp($crude_lib_count);
my $processed_lib_count=`cut -f3 $processed_lib_file | sort -u | wc -l`;
chomp($processed_lib_count);
print "...Do\n";
print "There are two files:\n\t$crude_lib_file\t($crude_lib_count entries, extract from tab file without processed).\n\t$processed_lib_file\t($processed_lib_count, where the libraries with less than $cluster_lib accessions have been clustered in a virtual library).\n"; 

if (!$opt_r) {
print_separator('REPORT');
print "\n\tSource files:\n";
foreach $source_file (@source_files) {
      my $accessions_counts=`grep 'ACCESSION' $source_file| wc -l`;
      chomp($accessions_counts);
      print "\t\t$source_file -----> $accessions_counts accessions\n";      
}
print "\n\tOutput .tab files:\n";
foreach $tab_file (@tab_files) {
    my $entries_counts=`cut -f2 $tab_file | sort -u | wc -l`;
    chomp($entries_counts);
    print "\t\t$tab_file -----> $entries_counts accessions\n";
}
print "\tGlobal .tab file:\n\t\t$output_global_tab -----> $gb_entries_count accessions\n";
print "\n\tOutput .fasta file:\n";
print "\t\t$output_global_fasta -----> $count_fasta entries\n";
print "\n\tOutput .tab libraries files:\n";
print "\t\t$crude_lib_file -----> $crude_lib_count entries\n";
print "\t\t$processed_lib_file -----> $processed_lib_count entries\n";
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
      This script get the data of GenBank files (.gb) in a genbank format and give us the same data in a different formats:
         A - tab format (with the name db_genbank_output.tab), with the follow fields:
               Locus, accession, version, primary_id, description, seq, qscore (added by the script, by default 15 
               for each nt), mol_type, organism, cell_line, cell_type, clone_name, clone_lib, cultivar, dbxref, 
               dev_stage, ecotype, environmental_sample, lab_host, note, PCR_primers, plasmid, tissue_lib, tissue_type,
               author.

             while the script extract the information of genbank file, replace the mRNA mol type for EST when find EST 
               in the keyword field (the mol_type for EST is mRNA in the GenBank file).

             for mRNA sequences without clone_lib (usually the core nucleotides download from the GenBank) update the 
               clone_lib like \'virtual mRNA library for genbank sequences from <organism>\'.

         B - fasta format

         C - libraries data in tab format. Give us two files.
	       > crude libraries_output.tab, this script analyze clone_lib (library_name), author, organism, cultivar, 
                 tissue_type, dev_stage, lab_host, plasmid, note and number of accessions associated to each library 
                 and give a file with one library in each row.
	       > processed libraries_output.tab, over the crude file, cluster the libraries with less than 100 sequences
                 in the \'virtual mRNA library for genbank sequences from <organism>\'. Also add the field type and update
                 it with 3 for common libraries and 17 for virtual library. Also add library_shortname with the follow 
                 code gbl_o<organism_id from the database used>.<substring 3 char of tissue type><library id for the load,
                 not in the sgn.library table> (the library_shortname must be UNIQUE when you load in the database) for
                 common libraries and vgbl_o<organism_id>.mRNA.
    
    Usage: 
      sgnpl_genbank_format_processing.pl <input dir with GenBank files> -H <dbhost> -D <dbname> [-m] <mol_type>  [-c] <cluster_with_entries_per_lib> [-d] [-r]
    
    Example:
      sgnpl_genbank_format_processing.pl /home/aure/GenBank_sequences/Nicotiana_sylvestris -H localhost -D sandbox -m 'EST','mRNA'

    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -m specify the mol_type that you can get of the genbank file (separated by commas ,)(An example (\'EST\',\'mRNA\')).
         By default is mRNA and EST.
      -c number of entries max per library to cluster it into a virtual library. By default is 100.
      -d disable the search of previous accessions in the database. It is enable, put this entries in another file. 
      -r disable the report option
      -h print this help

EOF
exit (1);
}

=head2 print_separator
 
  Usage: print_separator(<message>).
  Desc: this subroutine print a separator format in the screen
  Ret: none
  Args: scalar (message)
  Side_Effects: none
  Example: print_separator('Report:');

=cut

sub print_separator {
    my $var=shift;
    print "\n-------------------------------------------------------------\n";
    print "$var\n";
    print "-------------------------------------------------------------\n";
}

=head2 possible_abort_process

  Usage: possible_abort_process()
  Desc: give the option using STDIN (yes|no) of die the process.
  Ret: none
  Args: none
  Side_Efects: none (die the process)
  Example: if ($a > $b) {
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
 
=head2 create_file

  Usage: my($filehandle, $filepath)=create_file($dir, $filename, $fileextension); 
  Desc: Make and open a file in the specify location
  Ret: $filehandle (file handle) and $filepath (complete path of the file)
  Args: $dir (directory), $filename (name of the file) and $fileextension (extension of the file)
  Side_Efects: change the permissions and now all can write it
  Example: my($test_handle, $test_filepath)=create_file($test_dir, 'test', '.txt');

=cut

sub create_file {
    my $dir = shift;
    my $sufix = shift;
    my $filetype = shift;
    my $filename = $sufix . $filetype;
    my $filepath = $dir ."/". $filename;
    print "Create $filename file.\n";
    open my $filehandle, '>', "$filepath" || die "Cannot create the file $filepath\n";
    chmod 0666, "$filepath";
    return ($filehandle, $filepath);
}

=head2 create_dir

  Usage: my $dir=create_dir($path, $dirname);
  Desc: make a directory in the path specified
  Ret: dir name with the complete path
  Args: $path (path where you want create the file) and $dirname (name of the directory)
  Side_Effects: 
  Example: my $trash_dir=create_dir($path, 'garbage');

=cut

sub create_dir {
    my $path = shift;
    my $sufix = shift;
    my $dir= $path . "/" . $sufix;
    mkdir "$dir", 0755 or warn "Cannot make $dir directory: $!";
    return $dir;
}

=head2 process_genbank_file

  Usage: my $filetab=process_genbank_file($input_dir, $input_filename, $output_dir, $mol_type);
  Desc: script that take the $seq object from the genbank file and load in the database. Before it, checks if the locus
        is in the table (remove duplicates and write them in another file). Also check if exits the accession in the 
        public.dbxref table.
  Ret: the filepath for the genbank file in tab format (the output). A scalar.
  Args: $dir (workdir), $filename (genbank file), $outdir (directory where make the output file, the genbank .tab) and 
        $mol_type (type of molecule type tag to select, the rest is remove from the table).
  Side_Effects: if find that the accession exits in the public.dbxref, remove them of the temp table and put in a file
  Example: my $output_tab = process_genbank_file($genbank, $files, $tab_files, $mol_type);

=cut

sub process_genbank_file {
  my $dir=shift;
  my $filename=shift;
  my $outdir=shift;
  my $mol_type=shift;
  my $input_filename_gb = $dir . "/" . $filename;
  my $temp_table = $filename;
  $temp_table =~ s/\.gb/_temp/i;
  $temp_table =~ s/\W/_/gi;

  my ($check_locus, $check_dbxref, $previous_acc_filehandle, $previous_acc_filepath);

  my $input = Bio::SeqIO->new( 	-file   => 	$input_filename_gb,
                       		-format => 	'genbank'
				);

  $dbh->do("CREATE TEMP TABLE $temp_table 
	  (locus varchar(250), accession varchar(250), version varchar(250), primary_id varchar(250), 
           description text, seq text, qscore_def text, mol_type varchar(250), organism varchar(250), 
           cell_line varchar(250), cell_type varchar(250), clone_name varchar(250), clone_lib varchar(250), 
           cultivar varchar(250), dbxref varchar(250), dev_stage varchar(250), ecotype varchar(250), 
           environmental_sample varchar(250), lab_host varchar(250), note text, PCR_primers text, 
           plasmid varchar(250), tissue_lib varchar(250), tissue_type varchar(250), author text, 
           UNIQUE (accession, version))");
  
  my $duplicatefilename=$filename;
  $duplicatefilename =~ s/\.gb/\.duplicates_list/i;
  my ($duplicatefilehandle, $duplicatefilepath)=create_file($outdir, $duplicatefilename, '.tab'); 
  
  if (!$opt_d) {
      my $previous_acc_filename=$filename;
      $previous_acc_filename =~ s/\.gb/\.previous_accession_in_DB/i;
      ($previous_acc_filehandle, $previous_acc_filepath)=create_file($outdir, $previous_acc_filename, '.tab');
  }
  my $n=0;
  
  while (my $seq = $input->next_seq()) {
      if (!$opt_r) {
	  print "Getting the GenBank fields...\n";
      }

      my ($locus, $accession, $version, $primary_id, $descrip, $sequence, $qscore, $mol_type, $organism, $cell_line, 
          $cell_type, $clone_name, $clone_lib, $cultivar, $dbxref, $dev_stage, $ecotype, $env_sample, $lab_host, $note, 
          $pcr_pr, $plasmid, $tissue_lib, $tissue_type, $authors)=get_genbank_fields($seq);
      if (!$opt_r) {
	  print "... GenBank field got for the accession: $accession.\n";
	  print "Checking for previous sequences by accession in the $temp_table temp table...\n";
      }

      my $check1 = "SELECT locus FROM $temp_table WHERE accession = ? AND version = ?";
      my $sth=$dbh->prepare($check1);
      $sth->execute($accession, $version);
      ($check_locus)=$sth->fetchrow_array();
      if (!$opt_r) {
	  print "...checked... for $accession.$version the locus is:\t";
      }
      if ($check_locus) {
	  print "$check_locus.\n";
      } else {
	  print "None.\n";
      }
      if (!$opt_d) { 
          print "Checking for previous sequences by accession in the public.dbxref table....\n";
          my $check2= "SELECT est_id FROM sgn.est_dbxref JOIN public.dbxref 
                       ON sgn.est_dbxref.dbxref_id=public.dbxref.dbxref_id 
                       WHERE accession = ? AND version = ? AND db_id = 
                       (SELECT db_id FROM public.db WHERE name = 'DB:GenBank_Accession')";
	  $sth=$dbh->prepare($check2);
	  $sth->execute($accession, $version);
	  ($check_dbxref)=$sth->fetchrow_array();
          if (!$opt_r) {
	      print "...checked...for $accession.$version this script has found:\t";
	  }
	  if ($check_dbxref) {
	      print "est_id = $check_dbxref.\n";
	  } else {
	      print "None.\n";
	  }
      }
    
      if (!$check_locus) {
	  if (!$check_dbxref) {
	      my $insert = "INSERT INTO $temp_table VALUES (?,?,?,?,?,?, ?,?,?,?,?,?, ?,?,?,?,?,?, ?,?,?,?,?,?, ?)";
	      $sth=$dbh->prepare($insert);
	      $sth->execute($locus, $accession, $version, $primary_id, $descrip, $sequence, $qscore, $mol_type,
                            $organism, $cell_line, $cell_type, $clone_name, $clone_lib, $cultivar, $dbxref, $dev_stage, 
                            $ecotype, $env_sample, $lab_host, $note, $pcr_pr, $plasmid, $tissue_lib, $tissue_type,
                            $authors);
              $n++;
	      if (!$opt_r) {
		  print "Processed Sequence number $n with the accession:  $accession\n\n";
	      }
	  } else {
	      print $previous_acc_filehandle "$check_dbxref\t$accession\t$version\n";
	  }
      } else {
	  print $duplicatefilehandle "$check_locus\t$accession\t$version\n";
      }    
  }
  
  my $query = "SELECT COUNT(*) FROM $temp_table";
  $sth=$dbh->prepare($query);
  $sth->execute();
  my ($count_temp_table)=$sth->fetchrow_array();
  print "There are $count_temp_table entries in the $temp_table temp table.\n";
   
  if (!$opt_d) {
      my $previous_accessions_count=`cut -f1 $previous_acc_filepath | sort -u | wc -l`;
      chomp($previous_accessions_count);
      if ($previous_accessions_count > 0) {
	  print "There are $previous_accessions_count previous accessions in the database.\n";
	  print "These entries have been moved from the tab file to $previous_acc_filepath file.\n";
          my $count_accession_file=`grep 'ACCESSION' $previous_acc_filepath | wc -l`;
          chomp($count_accession_file);
          if ($previous_accessions_count = $count_accession_file) {
	      die "There aren\'t new accessions in the $filename file.\n";
	  }
      }
  }
  if ($mol_type) {
      print ">>>>>\nDELETE FROM $temp_table WHERE mol_type NOT IN $mol_type.\n>>>>>\n";
      $dbh->do("DELETE FROM $temp_table WHERE mol_type NOT IN $mol_type");
  }

  my $gbtab_filename=$filename;
  $gbtab_filename =~ s/\.gb//;
  my ($gbtab_filehandle, $gbtab_filepath)=create_file($outdir, $gbtab_filename, '.tab');
  $dbh->do("COPY $temp_table TO '$gbtab_filepath'");

  my $accessions_count=`cut -f2 $gbtab_filepath | sort -u | wc -l`;
  chomp($accessions_count);
  print "The $gbtab_filename file has $accessions_count accessions.\n";
  
  my $duplicates_count=`cut -f2 $duplicatefilepath | sort -u | wc -l`;
  chomp($duplicates_count); 
  if ($duplicates_count > 0) {
      print "There are $duplicates_count duplicates accessions in $filename.\n They have been removed in the output file $gbtab_filename\n. The duplicate accession list could be found in the file $duplicatefilename\n";
  
  }

  
  check_organism($temp_table);
  return $gbtab_filepath;
  
}

=head2 get_genbank_fields

  Usage: my @gbfields=get_genbank_fields($seq);
  Desc: this scripts get the different fields associate to a seq object (see tags below) 
  Ret: An array with the seq objects
  Args: $seq (seq object)
  Side_Effects: write EST in mol_type if appears in the keywords field and NULL is there aren\'t any object associate
  Example: my @gbfields=get_genbank_fields($seq);

=cut

sub get_genbank_fields {
    my $seq=shift;
    my ($feat_object, $annot_object, $tags, $feat);

    my $seqlocus = $seq->display_id();
    my $seqid = $seq->accession_number();
    my $seqversion = $seq->seq_version();
    my $primary_id = $seq->primary_id();
    my $seqdesc = $seq->description();
    my $sequence = $seq->seq();
    my $qscore_default = $sequence;
       $qscore_default =~ s/\w/15 /g;
    my @seqkeyword = $seq->get_keywords();
        
    my @tags = qw/organism cell_line cell_type clone clone_lib cultivar dbxref dev_stage ecotype
                      environmental_sample lab_host note PCR_primers plasmid tissue_lib tissue_type/;
    my @field = ($seqlocus, $seqid, $seqversion, $primary_id, $seqdesc, $sequence, $qscore_default); 
	
    for my $feat_object ($seq->get_SeqFeatures) {
        my $primary_tag = $feat_object->primary_tag;
	if ($primary_tag eq 'source') {
	    if ($seqkeyword[0] eq 'EST') {
		push @field, $seqkeyword[0];
	    } else {
		my @feat = $feat_object->get_tag_values('mol_type') if ($feat_object->has_tag('mol_type'));
   	        foreach $feat (@feat) {
		    push @field, $feat;
  		}
            }        
	    
	    foreach $tags (@tags) {
                my $validate = $feat_object->has_tag($tags);
       		if ($validate eq 1) {
		    my @feat = $feat_object->get_tag_values($tags);
		    foreach $feat (@feat) {
                        $feat =~ s/;+/,/gi;
			push @field, $feat;
		    }
		}elsif ($validate eq 0) { 
		    my $default=undef;
		    push @field, $default;
		}
	    }
	    	
	    last;
	}
    }
    
    for my $annot_object ($seq->get_Annotations('reference')) {         		
	my $authors = $annot_object->authors();
	if (!$authors) {
	    $authors=undef;
	}
	push @field, $authors;
	last;	
    }

    return @field;
}

=head2 check_organism

  Usage: check_organism($temp_table);
  Desc: this script checks if there are more than one organism in a genbank temp table
  Ret: none
  Args: $temp_table
  Side_Effects: if find more than one organism give the option of abort the process or remove some of them 
                unless only have one
  Example: check_organism($table_for_genbank);

=cut

sub check_organism {
  my $temp_table=shift;
  my $check = "SELECT COUNT(DISTINCT organism) FROM $temp_table";
  my $sth=$dbh->prepare($check);
  $sth->execute();
  my ($organism_count)=$sth->fetchrow_array();
  
  while ($organism_count > 1) {
      print "There are more than one organism in the file:($organism_count)\n";
      my $query = "SELECT organism, COUNT(locus) FROM $temp_table GROUP BY organism";
      $sth=$dbh->prepare($query);
      $sth->execute();
      while (my (@organism_list)=$sth->fetchrow_array()) {
	  my $organism_check = shift(@organism_list);
	  my $accession_count = shift(@organism_list);
	  print "$organism_check\t\twith\t\t$accession_count accessions\n";
      }
      print "\n**********\n*WARNING:*\n**********\n";
      print "Sorry, you only can have one organism.\nThis script can remove the aditional organism but if you want you can abort the process\n";
      possible_abort_process();
      print "Please write the name or names of the organism that you want remove (with quotes <'> and separated by commas <,> for example, 'Potato virus Y','Potato virus X'\n";
      my $preremove_organism=<STDIN>;
      chomp($preremove_organism);
      my $remove_organism="(".$preremove_organism.")";
      $dbh->do("DELETE FROM $temp_table WHERE organism IN $remove_organism");
      my $check3 = "SELECT COUNT(DISTINCT organism) FROM $temp_table";
      my $sth=$dbh->prepare($check3);
      $sth->execute();
      ($organism_count)=$sth->fetchrow_array();
  }
  my $check4 = "SELECT organism_id FROM sgn.organism JOIN $temp_table ON sgn.organism.organism_name=$temp_table.organism";
  $sth=$dbh->prepare($check4);
  $sth->execute();
  my ($organism_id)=$sth->fetchrow_array();
  if (!$organism_id) {
      print "\n**********\n*WARNING:*\n**********\n";
      print "Sorry, the organism used is not in the sgn.organism table.\n";
      print "Some subfunctions like library_shortname constructor use this data.\n";
      possible_abort_process();
  }
}

=head2 tab_to_fasta

  Usage: my $fasta_file=tab_to_fasta($outdir_fasta, $temp_table);
  Desc: get the secuences from the temp table and write a fasta file with them
  Ret: A scalar. The $filename_fasta
  Args: $outdir_fasta (the directory where make the fasta file) and $temp_table (temp table with all the genbank data)
  Side_Effects: none
  Example: my $ns_fasta=tab_to_fasta($fasta_files_dir, $ns_gb_temp);

=cut

sub tab_to_fasta {
  my $outdir_fasta=shift;
  my $temp_table=shift;
  
  $dbh->do("SELECT accession, seq INTO TEMP TABLE tab_format FROM $temp_table");
  my ($seqtab_filehandle, $seqtab_filepath)=create_file($outdir_fasta, 'seq', '.tab');
  $dbh->do("COPY tab_format TO '$seqtab_filepath'");
                               
  my $output_filename_fasta = $outdir_fasta . "/".$sufix."_db_gbseqdata.fasta";

  my $in  = Bio::SeqIO->new(	-file 	=> 	$seqtab_filepath,
                       		-format => 	'tab'
			);
  my $out = Bio::SeqIO->new(	-file 	=> 	">$output_filename_fasta",
                       		-format => 	'fasta'
			);

  while ( my $seq = $in->next_seq() ) {
	$out->write_seq($seq); 
  }
  return $output_filename_fasta;
}

=head2 create_library_tab

  Usage: my ($crude_lib_filepath, $processed_lib_filepath) = create_library_tab ($outdir, $temp_table, $cluster_n);
  Desc: make to files, crude library file with all te libraries and processed library file where the libraries with 
        less than cluster_n accessions are cluster in the virtual library
  Ret: two scalars, the crude library file path and the processed library file path
  Args: $outdir_tab (output directory), $temp_table (temp table with the genbank datas) and $cluster_lib (the number of 
        accessions cutoff to cluster the libraries in the virtual libraries) 
  Side_Effects: change the clone_lib names in the $temp_table (for example if change to \'virtual library\' two libraries,
                change too the accessions files in the $temp_table that have this clone_lib name
  Example: my ($lib_c, $lib_p) = create_library_tab($output_dir, $temp_table, '100');

=cut

sub create_library_tab {
  my $outdir_tab=shift;
  my $temp_table=shift;
  my $cluster_lib=shift;
  my $temp_table_idx1 = $temp_table . "_idx1";
  my $temp_table_idx2 = $temp_table . "_idx2";
  my $temp_table_idx3 = $temp_table . "_idx3";
  my $temp_table_idx4 = $temp_table . "_idx4";

  $dbh->do("UPDATE $temp_table SET clone_lib = 'virtual mRNA library for genbank sequences from '||organism, 
            author =NULL, cultivar =NULL, tissue_type=NULL, dev_stage=NULL, lab_host=NULL, plasmid=NULL, note=NULL  
            WHERE clone_lib IS NULL AND mol_type LIKE 'mRNA'");

  $dbh->do("SELECT DISTINCT(clone_lib), author, organism, cultivar, tissue_type, dev_stage, lab_host, plasmid, note, COUNT($temp_table.accession) INTO TEMP TABLE crude_libraries FROM $temp_table GROUP BY clone_lib, author, organism, cultivar, tissue_type, tissue_type, dev_stage, lab_host, plasmid, note");
   
  my ($crude_lib_filehandle, $crude_lib_filepath) = create_file($outdir_tab, 'crude_library_output', '.tab');
  
  $dbh->do("COPY crude_libraries TO '$crude_lib_filepath'");

  $dbh->do("CREATE TEMP TABLE processed_libraries
		(gb_lib_id SERIAL, clone_lib varchar(250), author text, organism varchar(250), cultivar varchar(250), 
		 tissue_type varchar(250), dev_stage varchar(250), lab_host varchar(250), plasmid varchar(250), 
		 note text, count_acc int)");
  $dbh->do("COPY processed_libraries (clone_lib, author, organism, cultivar, tissue_type , dev_stage, lab_host,
		 plasmid, note, count_acc) FROM '$crude_lib_filepath'");
  $dbh->do("ALTER TABLE $temp_table ADD COLUMN gb_lib_id int");
  $dbh->do("CREATE INDEX $temp_table_idx1 ON $temp_table (clone_lib)");
  $dbh->do("CREATE INDEX $temp_table_idx2 ON $temp_table (author)");
  $dbh->do("CREATE INDEX $temp_table_idx3 ON $temp_table (cultivar)");
  $dbh->do("CREATE INDEX $temp_table_idx4 ON $temp_table (tissue_type)");
  $dbh->do("CREATE INDEX processed_libraries_idx1 ON processed_libraries (clone_lib)");
  $dbh->do("CREATE INDEX processed_libraries_idx2 ON processed_libraries (author)");
  $dbh->do("CREATE INDEX processed_libraries_idx3 ON processed_libraries (cultivar)");
  $dbh->do("CREATE INDEX processed_libraries_idx4 ON processed_libraries (tissue_type)");

  my @null_fields=('clone_lib', 'author', 'cultivar', 'tissue_type', 'dev_stage', 'lab_host', 'plasmid', 'note');
  my $null_values;
  foreach $null_values (@null_fields) {
      $dbh->do("UPDATE $temp_table SET $null_values = 'unk' WHERE $null_values IS NULL");
      $dbh->do("UPDATE processed_libraries SET $null_values ='unk' WHERE $null_values IS NULL");
  }
  $dbh->do("UPDATE $temp_table SET gb_lib_id = (SELECT gb_lib_id FROM processed_libraries WHERE
            processed_libraries.clone_lib||processed_libraries.author||processed_libraries.organism||
            processed_libraries.cultivar||processed_libraries.tissue_type||processed_libraries.dev_stage||
            processed_libraries.lab_host||processed_libraries.plasmid||processed_libraries.note=
            $temp_table.clone_lib||$temp_table.author||$temp_table.organism||$temp_table.cultivar||
            $temp_table.tissue_type||$temp_table.dev_stage||$temp_table.lab_host||$temp_table.plasmid||
            $temp_table.note)");
  foreach $null_values (@null_fields) {
      $dbh->do("UPDATE $temp_table SET $null_values = NULL WHERE $null_values = 'unk'");
      $dbh->do("UPDATE processed_libraries SET $null_values = NULL WHERE $null_values = 'unk'");
  }
  if ($cluster_lib) {
      $dbh->do("UPDATE processed_libraries SET clone_lib = 'virtual mRNA library for genbank sequences from '||
                organism, author =NULL, cultivar =NULL, tissue_type=NULL, dev_stage=NULL,
                lab_host=NULL, plasmid=NULL, note=NULL  WHERE count_acc < $cluster_lib");
  }

  $dbh->do("UPDATE processed_libraries SET clone_lib = substring(organism from 1 for 1)||'.'||
	       split_part(organism, ' ', 2)||'(renamed library_'||gb_lib_id||')' WHERE clone_lib IS NULL");
  $dbh->do("ALTER TABLE processed_libraries ADD COLUMN type int");
  $dbh->do("UPDATE processed_libraries SET type = 17 WHERE clone_lib LIKE 'virtual mRNA library%' ");
  $dbh->do("ALTER TABLE processed_libraries ADD COLUMN library_shortname varchar(100)");
  $dbh->do("UPDATE processed_libraries SET library_shortname = 'vlGB_o'||
	    (SELECT organism_id FROM sgn.organism WHERE sgn.organism.organism_name=processed_libraries.organism)
            ||'mRNA' WHERE clone_lib LIKE 'virtual mRNA library%' ");
  $dbh->do("UPDATE processed_libraries SET type = 3 WHERE library_shortname IS NULL ");
  $dbh->do("UPDATE processed_libraries SET library_shortname = 'gbl_o'||
            (SELECT organism_id FROM sgn.organism WHERE sgn.organism.organism_name=processed_libraries.organism)
            ||'.'||'tl_'||gb_lib_id WHERE library_shortname IS NULL");
  $dbh->do("UPDATE processed_libraries SET clone_lib = substring(clone_lib from 1 for 80) WHERE char_length(clone_lib) > 80");
  $dbh->do("UPDATE processed_libraries SET library_shortname = substring(library_shortname from 1 for 16) WHERE char_length(library_shortname) > 16");
  $dbh->do("UPDATE $temp_table SET clone_lib = (SELECT clone_lib FROM processed_libraries WHERE processed_libraries.gb_lib_id=$temp_table.gb_lib_id)");
  $dbh->do("ALTER TABLE $temp_table DROP COLUMN gb_lib_id");
  $dbh->do("SELECT type, clone_lib, library_shortname, author, organism, cultivar, tissue_type, dev_stage, lab_host, plasmid, note, sum(count_acc) INTO TEMP TABLE gb_libraries_filtered FROM processed_libraries GROUP BY type, clone_lib, library_shortname, author, organism, cultivar, tissue_type, dev_stage, lab_host, plasmid, note");
  
  my ($processed_lib_filehandle, $processed_lib_filepath)=create_file($outdir_tab, 'processed_libraries_output', '.tab');
  $dbh->do("COPY gb_libraries_filtered TO '$processed_lib_filepath'");
  return ($crude_lib_filepath, $processed_lib_filepath);
}
