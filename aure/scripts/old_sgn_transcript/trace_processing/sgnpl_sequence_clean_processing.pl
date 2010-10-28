#!/usr/bin/perl

=head1 NAME
 
 sgnpl_sequence_clean_processing.pl.
 A script to clean a fasta file using the seqclean program (version.3.1.).
 This version parse the output file to combine them without the use of any database.

=cut

=head1 SYNOPSIS

sgnpl_sequence_clean_processing.pl -w <working_dir> -c <seqclean_path> -i <fasta_file> -v <end-trimming_file> -s <contaminant_screening_file> -l <sequence_length_remove_cutoff>-o <output_file> [-C] [-N] [-M] [-A] [-R]
    
=head2 I<Flags:>

=over
            
=item -w

B<Working dir>          directory where seqclean works and creates the cleaning_file folder. (by default \'pwd\').
      
=item -c

B<Seqclean path>        complete path to access to seqclean executable file. (mandatory)
      
=item -i 

B<Input fasta>          input file for clean. (mandatory)
      
=item -v

B<End-trimming file>    comma delimited list of fasta format sequences files to use for end-trimming of sequences (usually vector sequences). These fasta files must be format by formatdb script.
   
=item -s 

B<Contaminat screen>    comma delimited list of fasta format sequences files to use for screening sequences for contamination (cell host genome, mito\/ribo sequences contamination...). These fasta file must be format by formatdb script.
      
=item -l

B<seq_length_cutoff>    sequence length cutoff value to remove a sequence to the dataset (seqclean parameter)(default 100)

=item -o 

B<Output file>          output file with cleaning coordenates (by default <input_fasta_file>.clean_coord.tab).
      
=item -C 

B<Disable script op>    disable cyclic process while trim_number > 0.
      
=item -N 

B<Disable seqclean op>  disable trimming of ends rich in Ns (undetermined bases).
      
=item -M 

B<Disable seqclean op>  disable trashing of low quality sequences.
      
=item -A 

B<Disable seqclean op>  disable trimming of polyA\/T tails.
           
=item -h 

B<Help>                 print the help

=back

=cut

=head1 DESCRIPTION

This script clean a fasta file. Use seqclean tool (http://compbio.dfci.harvard.edu/tgi/software/) to do it, but this script check the output file. If this output has in its report a trimmed sequence number more than 0, rerun the tool over the clean fasta file. Finally take the all the cleanning coordenates files, calculate the final clean coordenates and put in a tab file. To do this use psql. This script make a cleanning files folder and moves all the files to it. 

The input files are:

=over
	  
=item - 

(tag -i): fasta file to clean.
	  
=item - 

(tag -v): comma delimited list of fasta format sequences files to use for end-trimming of sequences (usually vector sequences). These fasta files must be format by formatdb script.
          
=item - 

(tag -s): comma delimited list of fasta format sequences files to use for screening sequences for contamination (cell host genome, mito\/ribo sequences contamination...). These fasta file must be format by formatdb script.

=back
      
The output file (in the cleaning folder in the work dir) are:

=over

=item - 

cleaning_1 folder: Is the folder where the seqclean program put the cleaning slices files before ensambled all. This folder are removed before run seqclean another time.
          
=item - 

.cdix file: binary file with indexes for each run of the program.
          
=item - 

.clean file: fasta file with the cleanned sequences for each run of the program.
          
=item - 

.cln file: report file with the clean coordenates for each run of the program. 
          
=item - 

seqcl<input_fasta_file>.log: report file of the process for each run of the program.
          
=item - 

err_seqcl<input_fasta_file>.log: report file with possible errors for each run of the program
          
=item - 

outparts_cln_sort: file with the all the search slices for the last run of the program.

=item - 

(tag -o): output file with the clean coordenates for the global process.
          
=item - 

global_seqclean_report.txt: report of all seqclean process.
          
=item - 

equiv_clean_codes: .tab file with the equivalences between files used in the screening and the clean_code

=back

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=cut

=head1 METHODS

sgnpl_sequence_clean_processing.pl


=cut

use strict;
use File::Basename;
use Getopt::Std;

our ($opt_w, $opt_c, $opt_i, $opt_v, $opt_s, $opt_l, $opt_o, $opt_C, $opt_N, $opt_M, $opt_A, $opt_h);
getopts("w:c:i:v:s:l:o:CNMAh");

if (!$opt_w && !$opt_c && !$opt_i && !$opt_v && !$opt_s && !$opt_l && !$opt_o && !$opt_C && !$opt_N 
      && !$opt_M && !$opt_A && !$opt_h) {
             print "There are n\'t any tags. Print help\n\n";
             help();
}

print "\nValidating the input...\n";
my ($dbhost, $dbname, $workdir, $seqclean_path, $input_fasta, $output_file)=validate_the_input();

my $cleandir= $workdir . "/cleaning_files";
mkdir "$cleandir", 0755 or die "Cannot make $cleandir directory: $!";
print "\nCreate cleaning directory ($cleandir)\n";
my @globalreport;

our @seqclean_run_report;
my $n=1;
my $precommand="cd $cleandir ; $seqclean_path $input_fasta ";
my $command=command_seqclean($precommand, $cleandir, $n);
my $count = 1;

print "\n---------------------------------------------------------------------\n";
print "Runing system .......... for ... $count time";
print "\n---------------------------------------------------------------------\n";
system "$command";

my $input_fasta_name=File::Basename::basename($input_fasta);
my $log_file=$cleandir."/seqcl_".$input_fasta_name.".log";
my $report_file=$cleandir."/".$input_fasta_name.".cln";
my $output_seqclean_file=$cleandir."/".$input_fasta_name.".clean";
my $trim_n = get_trim($log_file);
print "The number of sequences trimmed is $trim_n\n";
push @globalreport, $report_file;


if (!$opt_C) {
while ($trim_n > 0) {
        $input_fasta = $output_seqclean_file;
	$output_seqclean_file = $input_fasta . ".clean";
        
        $n++;
        my $runing_c=$count+1;
        print "\n---------------------------------------------------------------------\n";
        print "Runing system .......... for ... $runing_c time";
        print "\n---------------------------------------------------------------------\n";
	$command=command_seqclean("cd $cleandir ; $seqclean_path $input_fasta ", $cleandir, $n);
        system "$command";

        $log_file =~ s/\.log/\.clean\.log/i;
        $report_file =~ s/\.cln$/\.clean.cln/i;

	$trim_n = get_trim($log_file);
	print "\nThe number of sequences trimmed is $trim_n\n";
	$count++;
        push @globalreport, $report_file;
	}
}
print "\n---------------------------------------------------------------------\n";
print "\nThe seqclean program has runned $count over the fasta file\n";
print "\n---------------------------------------------------------------------\n";

combine_data(\@globalreport, $output_file);

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
      This script clean a fasta file. Use seqclean tool (http://compbio.dfci.harvard.edu/tgi/software/) to do it, but this script check the output file. If this output has in its report a trimmed sequence number more than 0, rerun the tool over the clean fasta file. Finally take the all the cleanning coordenates files, calculate the final clean coordenates and put in a tab file. To do this use psql. This script make a cleanning files folder and moves all the files to it. 

      The input files are:
	  - (tag -i): fasta file to clean.
	  - (tag -v): comma delimited list of fasta format sequences files to use for end-trimming of sequences
            (usually vector sequences). These fasta files must be format by formatdb script.
          - (tag -s): comma delimited list of fasta format sequences files to use for screening sequences for 
            contamination (cell host genome, mito\/ribo sequences contamination...). These fasta file must be 
            format by formatdb script.

      The output file (in the cleaning folder in the work dir) are:
	  - cleaning_1 folder: Is the folder where the seqclean program put the cleaning slices files before 
            ensambled all. This folder are removed before run seqclean another time.
          - .cdix file: binary file with indexes for each run of the program.
          - .clean file: fasta file with the cleanned sequences for each run of the program.
          - .cln file: report file with the clean coordenates for each run of the program. 
          - seqcl<input_fasta_file>.log: report file of the process for each run of the program.
          - err_seqcl<input_fasta_file>.log: report file with possible errors for each run of the program
          - outparts_cln_sort: file with the all the search slices for the last run of the program.

          - (tag -o): output file with the clean coordenates for the global process.
          - global_seqclean_report.txt: report of all seqclean process.
          - equiv_clean_codes: .tab file with the equivalences between files used in the screening and the clean_code

    Usage: 
      sgnpl_sequence_clean_processing.pl -w <working_dir> -c <seqclean_path> -i <fasta_file> 
      -v <end-trimming_file> -s <contaminant_screening_file> -o <output_file> [-C] [-N] [-M] [-A] [-R]

    Example:
      sgnpl_sequence_clean_processing.pl -H localhost -D sandbox -w /home/aure/GenBank_sequences/Nicotiana_sylvestris 
      -c /home/aure/seqclean/seqclean/seqclean -i /home/aure/GenBank_sequences/Nicotiana_sylvestris/Ns_seq.fasta
      -v /home/aure/screening/UniVec -s /home/aure/screening/ColiBank95.lib,/home/aure/screening/sol_mit_rib_ncds.fasta
      -o /home/aure/GenBank_sequences/Nicotiana_sylvestris/Ns_seq.clean_coord.tab  
    
    Flags:
      -w working dir          directory where seqclean works and creates the cleaning_file folder. (by default \'pwd\').
      -c seqclean path        complete path to access to seqclean executable file. (mandatory)
      -i input fasta          input file for clean. (mandatory)
      -v end-trimming file    comma delimited list of fasta format sequences files to use for end-trimming 
                              of sequences (usually vector sequences). These fasta files must be format by 
                              formatdb script.
      -s contaminat screen    comma delimited list of fasta format sequences files to use for screening sequences 
                              for contamination (cell host genome, mito\/ribo sequences contamination...). These 
                              fasta file must be format by formatdb script.
      -l min sequence length  minimum sequence length to get a sequence as rigth (default 100)
      -o output file          output file with cleaning coordenates (by default <input_fasta_file>.clean_coord.tab).
      -C disable script op    disable cyclic process while trim_number > 0.
      -N disable seqclean op  disable trimming of ends rich in Ns (undetermined bases).
      -M disable seqclean op  disable trashing of low quality sequences.
      -A disable seqclean op  disable trimming of polyA\/T tails.  
      -h help                 print this help

EOF
exit (1);

   }

=head2 validate_the_input

  Usage: my ($dbhost, $dbname, $workdir, $seqclean_path, $input_fasta, $output_file)=validate_the_input();
  Desc: validate the different input variables, like localhost, database name, seqclean program path...
  Ret: six scalars. $dbhost, $dbname, $workdir, $seqclean, $input_fasta and $output_file
  Args: none (the $opt are our variables)
  Side_Effects: die the process if the input variable is not right. For optional variables give a default value.
  Example: my @input_var=validate_the_input();

=cut

sub validate_the_input {
  
    if ($opt_h) {
	help();
    }
   
    unless ($opt_c) {
	die ("Required argument -c <seqclean_path> was not supplied.\n");
    }
    unless ($opt_i) {
	die ("Required argument -i <input_fasta sequence> was not supplied.\n");
    }

    my $seqclean_path = `ls -F $opt_c`;
    chomp($seqclean_path);
    $seqclean_path =~ m/(seqclean)\*$/;
    if ($1 ne 'seqclean') {
	die ("Sorry, the seqclean path is not right. The $opt_c path not is the seqclean executable path.\n");
    } else {
	print "\tCheck seqclean path ... Ok\n";
    }

    my $input_fasta=$opt_i;
    unless (-e $input_fasta) {
	die ("Sorry, the input fasta file do not exists\n");
    } else {
	print "\tCheck .fasta file extension ... Ok\n";
    }

    my $workdir;
    if ($opt_w) {
	$workdir=$opt_w;
    } else {
	$workdir=`pwd`;
	chomp($workdir);
    }
    print "\tThe working dir will be ... $workdir\n";

    my $output_file;
    if ($opt_o) {
	$output_file=$opt_o;
    } else {
        my $input_fasta_name=File::Basename::basename($input_fasta);
        my $filepath=File::Basename::dirname($input_fasta);
        $output_file=$filepath."/cleaning_files/OUT_CLEANCOORD_".$input_fasta_name.".tab";
    }
    print "\tThe output cleaning coordenate file will be ... $output_file.\n";

    my $count_ifasta=`grep '>' $input_fasta | wc -l`;
    chomp($count_ifasta);
    my ($trim_file, $contam_file);
   
    print "---------------------------------------------------------\n";
    print "Input files:\n";
    print "\tFasta file:\n\t\t$opt_i with $count_ifasta sequences.\n";
    
    if ($opt_v) {
	my @trim_files=split(',', $opt_v);
        print "\tTrim files:\n";
	foreach $trim_file (@trim_files) {
	    my $count_itrim=`grep '>' $trim_file | wc -l`;
	    chomp($count_itrim);
	    print "\t\t$trim_file with $count_itrim sequences.\n";
	}
    }
    if ($opt_s) {
	my @contam_files=split(',', $opt_s);
	print "\tContaminant files:\n";
	foreach $contam_file (@contam_files) {
	    my $count_icontam=`grep '>' $contam_file | wc -l`;
	    chomp ($count_icontam);
	    print "\t\t$contam_file with $count_icontam sequences.\n";
	}
    }
    if ($opt_l) {
	unless ($opt_l =~ m/^\d+$/) {
	    die "Sorry, the -l <minimum_sequence_length> must be an integer.\n";
	}
    }
    print "---------------------------------------------------------\n";
    return ($dbhost, $dbname, $workdir, $opt_c, $input_fasta, $output_file);
}

=head2 command_seqclean

  Usage: my $command = command_seqclean ($precommand, $cleandir, $n);
  Desc: this subroutine is the seqclean command constructor. The command depend of the input variables.
  Ret: A scalar. The command used in the system function.
  Args: $precomand (the begining of the command, without options), $cleandir (directory for the outputs) and 
        $n (number on seqclean run) 
  Side_Effects: none 
  Example: $command=command_seqclean("cd $cleandir ; $seqclean_path $input_fasta ", $cleandir, $n);

=cut

sub command_seqclean {
  my $command=shift;
  my $cleandir=shift;
  my $n=shift;

  if ($opt_v) {
    $command .= "-v $opt_v ";
  }
  if ($opt_s) {
    $command .= "-s $opt_s ";
  }
  if ($opt_l) {
      $command .= "-l $opt_l ";
  }
  if ($opt_N) {
    $command .= "-N ";
  }
  if ($opt_M) {
    $command .= "-M ";
  }
  if ($opt_A) {
    $command .= "-A ";
  }
 return $command;
  
} 

=head2 get_trim

  Usage: my $trim_number=get_trim($report_file);
  Desc: get the trim number of the report file (output of seqclean program).
  Ret: a scalar. The number of sequences trimmed during the seqclean process.
  Args: $report_file (the output file of seqclean program)
  Side_Effects: die if can not read the report file
  Example: my $trim_number=get_trim($report_file);

=cut

sub get_trim {
  my $file=shift;
  my $trim_number ='';
  open (my $input, $file) or die "Can't read $file: $!";
  while ( <$input> ) {
	if ($_ =~ m/(\d+) trimmed/) {
  	$trim_number = "$1";
	}
  }
  return $trim_number;
}

=head2 combine_data

  Usage: combine_data($arrayref_with_the_reports, $output_file );
  Desc: this subroutine calculate the final coordenates based in the reports files.
  Ret: none
  Args: $arrayref_with_the_reports, an array reference of an array with the seqclean reports filenames 
        in the same order than they were made and $output_file, the output file name 
  Side_Effects: Make a output file.
  Example: combine_data($aref_cleanreports, $output_file);

=cut

sub combine_data {
  my $aref_cleanreports=shift;
  my $output_file=shift;
  my @files=@$aref_cleanreports;
  my $firstfile=$files[0];

  my $output_filename_db = $output_file;
  open my $out, '>', $output_file || die "Sorry, I can not open the output file $output_file\n";

  print "\nCombining data...\n";

  my (@id, $file, $idref);
  open my $fh, '<', $firstfile || die "I can not open the file $file\n";

  while (<$fh>) {
    chomp($_);
    my @coord=split(/\t/, $_);
    my $id=$coord[0];
    push @id, $id;
  }
  close $fh || die "I can not close the file $file\n";
  my $n_idref=scalar(@id);
  my @vectors=split(/,/, $opt_v);
  my @contam=split(/,/, $opt_s);
  my ($tagvec, $tagcontam);

  my $i=0;
  foreach $idref (@id) {
     $i++;
     print "Processing the id $idref ($i of $n_idref)\r";
     my ($gNpercent,$gbegin, $gend, $grawlength, $gcleantag)=(0,0,0,0);
     my ($glength, $gtrimtag);
     my $n=0;
     foreach $file (@files) {
         my $parsecoord=`grep '$idref' $file`;
         chomp ($parsecoord);
         if ($parsecoord) {
             my @coord=split(/\t/, $parsecoord);
             my $begin=$coord[2]-1;
             my $length=$coord[3]-$coord[2]+1;
             $gNpercent=$coord[1];
             $gbegin += $begin;
             $glength=$length;
             if ($n == 0) {
		 $grawlength=$coord[4];
	     }
             my $n_coord=scalar(@coord);
             
             my ($v,$s)=(0,0);
             foreach $tagvec (@vectors) {
                 my $tagvecname=File::Basename::basename($tagvec);
                 $v++;
		 if ($tagvecname =~ m/$coord[5]/) {
                     $gcleantag='trim_lib'.$v;
		 }
	     }
	     foreach $tagcontam (@contam) {
                 my $tagcontamname=File::Basename::basename($tagcontam);
		 $s++;
		 if ($tagcontamname =~ m/$coord[5]/) {
		     $gcleantag='contam_lib'.$s;
		 }
	     }
             if ($gcleantag == 0) {
		 $gcleantag=$coord[5];
	     }
             $gtrimtag .= $coord[6];
         }
         $n++;
     }
  
   $gbegin += 1;
   $gend=$glength+$gbegin-1;
   

   print $out "$idref\t$gNpercent\t$gbegin\t$gend\t$grawlength\t$gcleantag\t$gtrimtag\n";
 }
 print "\nProcessing done.\nThe output file is: $output_file\n";

}

=head2 make_cleancode_equiv_table

  Usage: my $clean_codes_equiv=make_cleancode_equiv_table($output_file, $trim_aref, $contam_aref);
  Desc: make a table with the equivalences between the clean codes by the file names ($opt_v and $opt_s) and
        the new names. Also give the possible flags and status associated with this clean_tags.
  Ret: A scalar, the equivalences table name.
  Args: $output_file (the name of the output_file), $trim_aref and $contam_aref (the array references of the 
        old clean codes).
  Side_Effects: none. 
  Example: my $clean_codes_equiv=make_cleancode_equiv_table($output_file, $trim_aref, $contam_aref);
 
=cut

sub make_cleancode_equiv_table {
  my $output_file=shift;
  my $trim_lib_ref=shift;
  my $contam_lib_ref=shift;
  my @trim_lib=@$trim_lib_ref;
  my @contam_lib=@$contam_lib_ref;

  my $clean_codes_equiv = $output_file.".equiv_clean_codes";
  open my $clean_codes_equiv_filehandle, '>', "$clean_codes_equiv";
  chmod 0666, "$clean_codes_equiv";
  
  my ($trim_lib, $contam_lib);
  my $p=1;
  my $q=1;
  foreach $trim_lib (@trim_lib) {
      print $clean_codes_equiv_filehandle "$trim_lib\ttrim_lib$p\t20\t1\n";
      $p++;
  }
  foreach $contam_lib (@contam_lib) {
      print $clean_codes_equiv_filehandle "$contam_lib\tcontam_lib$q\t20\t1\n";
      $q++;
  }
  print $clean_codes_equiv_filehandle "short\tshort\t4\t1\nqshort\tqshort\t4\t1\n";
  print $clean_codes_equiv_filehandle "dust\tdust\t20\t1\nlow_qual\tlow_qual\t8\t1\n";
  
  return $clean_codes_equiv;
}

=head2 print_report

 Usage:
 Desc:
 Ret:
 Args:
 Side_Effects:
 Example:

=cut

sub print_report {
  my $output_file=shift;
  
  my $parse_output=`cut -f6 $output_file | sort -u`;
  chomp($parse_output);
  print "\n=>\n$parse_output\n";

}
