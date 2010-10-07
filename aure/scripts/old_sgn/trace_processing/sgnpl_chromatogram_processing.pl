#!/usr/bin/perl

=head1 NAME
 
 sgnpl_chromatogram_processing.pl.
 A script to process chromatogram files using the phred program (version.2.0.)

=cut

=head1 SYNOPSIS

 sgnpl_chromatogram_processing.pl <dirname> -p <phred_program_path> [-S] [-Q] [-D] [-L] [-X] [-Z] [-t]

 <dirname> : input directory with the chromatograms.

=head2 I<Flags:>

=over

=item -p

B<phred_program_path>   complete file path to access to executable program phred file 
      
=item -S

B<enable option>        enable the option write the sequence in a file. Create a dir (<dirname_seq_dir>) with each sequence in fasta format for each chromatogram (seq_name = chromatogram_file_name.seq).    
      
=item -Q 

B<enable option>        enable the option write the quality values in a file. Create a dir (<dirname_qual_dir>) with each quality values in fasta format for each chromatogram (seq_name = chromatogram_file_name.seq).
      
=item -D

B<enable option>        enable the option write the polimorphism in a file. Create a dir (<dirname_poly_dir>) with each polymorphismfor each chromatogram (seq_name = chromatogram_file_name.seq).
      
=item -L

B<enable option>        enable the option, write a report in a .log file
      
=item -X 

B<enable option>        enable the option: Trim the sequences bases and the quality files with less than cutoff value (see -t argument).
      
=item -Z

B<enable option>        enable the option: Trim the base calling in the .phd files too (see -Z argument).
   
=item -t 

B<trim error-p value>   change the default (0.05) error probability value for the trim process  
      
=item -h 

B<Help>                 print the help

=back

=cut

=head1 DESCRIPTION

 This script process chromatogram files from a dir. Use phred tool (http://www.phrap.org/phredphrapconsed.html) to do it. This program produce different output files for each sequence (depend of the arguments) that this script process into a global files. Also make a global tab file with four columns (id (chromatogram_filename), sequence, quality values and basecalling). 

 Important!!!: You need export the phred parameters in the linux console before run the script. 
               > export PHRED_PARAMETER_FILE=phredpar.dat 

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=cut

=head1 METHODS

sgnpl_chromatogram_processing.pl


=cut


use strict;

use Getopt::Std;
use File::Basename;
					 
our $inputdir=shift(@ARGV);
if ($inputdir =~ m/^-h$/ ) {
     help();
} elsif ($inputdir =~ m/^-S|-Q|-D|-L|-X|-Z$/) {
     die  ("Required argument <chromatogram dir> was not supplied.\n");
}

our ($opt_p, $opt_S, $opt_Q, $opt_D, $opt_L, $opt_X, $opt_Z, $opt_t, $opt_h);
getopts("p:SQDLXZt:h");

if (!$opt_p && !$opt_S && !$opt_Q && !$opt_D && !$opt_L && !$opt_X && !$opt_Z && !$opt_t && !$opt_h && !$inputdir) { 
                print "There are n\'t any tags. Print help\n\n";
                help();
              }

my ($dir, $program) = validate_the_input();

my $phred_command=phred_command_constructor($program, $dir);

#my $phredpar = "/home/aure/trace_processing/phredpar.dat";
#system "export PHRED_PARAMETER_FILE=$phredpar";
#print "export PHRED_PARAMETER_FILE=$phredpar\n";

print "Running ...\n\t$phred_command\n\n"; 
system "$phred_command";

print "Processing the phred output files ...\n";

my $dir_phddir=$dir . "/phd_dir";

opendir DH, $dir_phddir or die "Cannot open $dir_phddir to process: $!";
my $file;
my $count=0;
my $filehandle='fh00';
my $output_file = $dir . "/phred_output.tab";
open my $output, '>', "$output_file";
my $output_file2 = $dir . "/phred_output.fasta";
open my $output_fasta, '>', "$output_file2";
foreach $file (readdir DH) {
        if ($file =~ /\.phd\.1$/) {
  		$filehandle++;
		my $seq_info = &extract_seq_qc($dir_phddir, $file, $filehandle);
		print $output "$seq_info";
		my $seq_fasta = &extract_seq_fasta($dir_phddir, $file, $filehandle);
		print $output_fasta "$seq_fasta";
                $count++;
	}
}
print "\tThere are $count chromatogram files process.\n\n";
print "\tThe global output file in tab format is:\n\t\t$output_file\n\n";
print "\tThe global output file with the sequences in fasta format is:\n\t\t$output_file2\n\n";


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
       This script process chromatogram files from a dir. Use phred tool (http://www.phrap.org/phredphrapconsed.html)
       to do it. This program produce different output files for each sequence (depend of the arguments) that this 
       script process into a global files. Also make a global tab file with four columns (id (chromatogram_filename),
       sequence, quality values and basecalling). 

       Important!!!: You need export the phred parameters in the linux console before run the script. 
                 
                     > export PHRED_PARAMETER_FILE=phredpar.dat 

    Usage: 
       sgnpl_chromatogram_processing.pl <dirname> -p <program path> [-S] [-Q] [-D] [-L] [-X] [-Z] [-t]

       <dirname>: Directory name with the chromatograms.
    
    Example:
       export PHRED_PARAMETER_FILE = phredpar.dat
 
       sgnpl_chromatogram_processing.pl /Nicotiana_tabacum/KF8_100 -p /home/aure/programs/phred/bin/phred -S -Q -D -L -X
          
    Flags:
      -p program path         complete path to access the phred program (mandatory)
      -S write seq in file    enable the option write the sequence in a file. Create a dir (<dirname_seq_dir>) 
                              with each sequence in fasta format for each chromatogram (file_name = 
                              chromatogram_file_name.seq).
      -Q write qual in file   enable the option write the quality values in a file. Create a dir (<dirname_qual_dir>) 
                              with each quality values in fasta format for each chromatogram (file_name = 
                              chromatogram_file_name.qual).
      -D write polym in file  enable the option write the polymorphism in a file. Create a dir (<dirname_poly_dir>) 
                              with each posssible polymorphism for each chromatogram (file_name = 
                              chromatogram_file_name.poly).
      -L write report         enable the option write a run report in a .log file.
      -X enable trim option   trim the sequences bases and the quality files with less than cutoff value 
                              (see -t argument).
      -Z enable trim option   trim the base calling in the .phd files too (see -t argument).
      -t change trim cutoff   change the default (0.05) error probability value for the trim process  
      -h help                 print this help

EOF
exit (1);

   }

=head2 validate_the_input

  Usage: my ($dir, $program)=validate_the_input();
  Desc: validate the input variables, like program path or dir...
  Ret: two scalars. $dir (directory with the cromatograms) and $program (phred program path)
  Args: none (the $opt are our variables)
  Side_Effects: die the process if the input variable is not right. For optional variables give a default value.
  Example: my ($dir, $program) =validate_the_input();

=cut

sub validate_the_input {
    print "Validating the inputs ... \n\n";
    if ($opt_h) {
	help();
    }
   
    unless ($opt_p) {
	die ("Required argument -p <program_path> was not supplied.\n");
    }
    unless (-x $opt_p) {
        die ("Program path is wrong, the file is not a binary file.\n");
    } 
    my $programname=File::Basename::basename($opt_p);
    if ($programname ne 'phred') {
	die "Sorry the -p <program_path> argument is not phred file.\n ";
    }

    my @filextensions=('.SCF', '.ABI', '.ESD');
    my $extensions_count=0;
    my ($valid_file_extensions, $checkdir);
    foreach $valid_file_extensions (@filextensions) {
	$checkdir=`ls $inputdir | grep -i '$valid_file_extensions' | wc -l`;
	chomp($checkdir);
	$extensions_count += $checkdir;
    }
    if ($extensions_count eq 0) {
	die "The input directory ($inputdir) have not .SCF, .ABI or .ESD chromatograms files.\n";
    } else {
	print "The input directory has $extensions_count chromatograms files (.SCF, .ABI or .ESD).\n";
    }
	
    return ($inputdir, $opt_p);
}


=head2 phred_command_constructor

  Usage: my $phred_command=phred_command_constructor($program, $dir);
  Desc: this subroutine construct the phred command for the system perl function.
  Ret: A scalar, the phred command ($phred_command)
  Args: $program (phred program path after check) and $dir (chromatogram directory after check) 
  Side_Effects: If exists some options variables (like $opt_S) make a dir.
  Examples: my $phred_command=phred_command_constructor($program, $dir);

=cut

sub phred_command_constructor {
    my $program=shift;
    my $dir=shift;

    my ($dis_S, $dis_Q, $dis_D, $dis_L, $dis_X, $dis_Z, $dis_t, $dis_ta);
    my $dir_phddir = $dir . "/phd_dir";
    mkdir "$dir_phddir", 0755 or warn "Cannot make $dir_phddir directory: $!";

    if ($opt_S) {
	my $dir_seqdir = $dir . "/seq_dir";
	mkdir "$dir_seqdir", 0755 or warn "Cannot make $dir_seqdir directory: $!";
	$dis_S ="-sd $dir_seqdir";
    }

    if ($opt_Q) {
	my $dir_qualdir = $dir . "/qual_dir";
	mkdir "$dir_qualdir", 0755 or warn "Cannot make $dir_qualdir directory: $!";
	$dis_Q ="-qd $dir_qualdir";
    }

    if ($opt_D) {
	my $dir_polydir = $dir . "/poly_dir";
	mkdir "$dir_polydir", 0755 or warn "Cannot make $dir_polydir directory: $!";
	$dis_D ="-dd $dir_polydir";
    }

    if ($opt_L) {
	$dis_L ='-log';
    }
    if ($opt_X) {
	$dis_X ='-trim_fasta';
    }
    if ($opt_Z) {
	$dis_Z ='-trim_phd';
    }
    if ($opt_X || $opt_Z) {
	$dis_ta= "-trim_alt \"\"";
    }
    if ($opt_t) {
	$dis_t = $opt_t;
    } else {
        $dis_t = 0.05;
    }

    my $phred_command = "$program -id $dir -pd $dir_phddir $dis_S $dis_Q $dis_D $dis_L $dis_X $dis_Z $dis_ta -trim_cutoff $dis_t";

    return $phred_command;
}

=head2 extract_seq_qc 

  Usage: my $seq_info = &extract_seq_qc($dir_phddir, $file, $filehandle);
  Desc: subroutine that open phd file produced by phred, and process the its information to get sequence, 
        quality values and call positions.
  Ret: A scalar, the file info in tab format ($seq_info)
  Args: $dir (phred results directory with the .phd files), $file (name of the file with the sequence) and 
        $filehandle (a filehandle name)
  Side_Effects: die if can not open the file
  Example: my $seq_info = &extract_seq_qc($dir_phddir, $file, $filehandle);

=cut

sub extract_seq_qc {
  my $dir = shift;
  my $file = shift;
  my $filename = $dir . "/" . $file; 
  my $filehandle = shift;
  my $id = $file;
  $id =~ s/\.phd\.1$//g;
  my $seq = "";
  my $qscore = "";
  my $call_positions= "";
  open my $filehandle, $filename or die "can't open $file: $!";
  while (<$filehandle>) {
      my $control;
        $_ =~ /^(BEGIN_DNA)$/;
	$control .= $1;
        $_ =~ /^(END_DNA)$/;
	$control .= $1;
	
	if ($_ =~ /^(\w)\s(\d+)\s(\d+)/) {
	    $seq .= $1;
	    $qscore .= "$2 ";
	    $call_positions .= "$3 ";
	}
	if ($control =~ /BEGIN_DNAEND_DNA/) {
	    $seq = 'N';
            $qscore = '0';
	    $call_positions = '0';
            last;
	}
  }
  my $seq_info = $id . "\t" . $seq . "\t" . $qscore . "\t" . $call_positions . "\n";
  return $seq_info; 
}

=head2 extract_seq_fasta

  Usage: my $seq_info = &extract_seq_fasta($dir_phddir, $file, $filehandle);
  Desc: subroutine that open phd file produced by phred, and process the its information to get sequence.
  Ret: A scalar, the sequence info in fasta format ($seq_fasta). If there are not sequence, return only a single n.
  Args: $dir (phred results directory with the .phd files), $file (name of the file with the sequence) and 
        $filehandle (a filehandle name)
  Side_Effects: die if can not open the file
  Example: my $seq_info = &extract_seq_fasta($dir_phddir, $file, $filehandle);
 
=cut

sub extract_seq_fasta {
  my $dir = shift;
  my $file = shift;
  my $filename = $dir . "/" . $file; 
  my $filehandle = shift;
  my $id = $file;
  $id =~ s/\.phd\.1$//g;
  my $seq = "";
  my $control = "";
  my $n = 1; 
  open my $filehandle, $filename or die "can't open $file: $!";
  while (<$filehandle>) {
	$_ =~ /^(BEGIN_DNA)$/;
	$control .= $1;
	if ($n < 50 && $_ =~ /^(\w)\s\d+\s\d+/) {
	$seq .= $1;
	$n++
	} elsif  ($n = 50 && $_ =~ /^(\w)\s\d+\s\d+/) {
	$seq .= "\n";
	$seq .= $1;
	$n =1;
	}
        $_ =~ /^(END_DNA)$/;
	$control .= $1;
	if ($control =~ /BEGIN_DNAEND_DNA/) {
	$seq .= "N";
	last;
	}	
  }
  my $seq_fasta = ">$id" . "\n" . $seq . "\n";
  return $seq_fasta; 
}
	
