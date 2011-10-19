#!/usr/bin/perl

=head1 NAME

sffstats.pl
A script to get some stats using sffinfo over sff datasets

=cut

=head1 SYPNOSIS

sffstats.pl -d <sffdirs> -o <output_basename>

=head2 I<Flags:>

=over


=item -d

B<sffdir>               List of dirs with sff files separated by ',' (mandatory)
 
=item -o

B<output_basename>      An output basename for the reports (mandatory)

=item -N

B<no_trimmed>           get the no-trimmed sequences from sffinfo.

=item -h

B<help>                 print the help

=back

=cut

=head1 DESCRIPTION

 This script run sffinfo program over different sequence .sff sets returning
 some reports for sequence distribution.
 (It runs only in linux OS)

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 sffstats.pl

=cut

use strict;
use warnings;
use autodie;

use Carp qw( croak );

use File::Spec;
use Cwd;

use Getopt::Std;
use Math::BigFloat;
#use R::YapRI::Base;



our ($opt_d, $opt_o, $opt_N, $opt_h);
getopts("d:o:Nh");
if (!$opt_d && !$opt_o && !$opt_N && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}


## Check arguments

my $dirlist = $opt_d || 
    die("DATA ARGUMENT ERROR: -d <sffdirs> WAS NOT SUPPLIED.\n");

my $outbasename = $opt_o ||
    die("DATA ARGUMENT ERROR: -o <output_basename> WAS NOT SUPPLIED");


## Check if the system has sffinfo, R, cut, sed, and grep:

my %prog_exe = (
    sffinfo => 'sffinfo',
    );

my $path = $ENV{PATH};
if (defined $path) {
    my @paths = split(/:/, $path);
    foreach my $req_exec (keys %prog_exe) {
	my $exec_path;
	foreach my $p (@paths) {	
	    my $ufile = File::Spec->catfile($p, $req_exec);
	    if (-e $ufile) {
		$exec_path = $ufile;
	    }	
	}
	unless (defined $exec_path) {
	    my $message = "ERROR: program $req_exec is not in the exec.";
	    $message .= "PATH. Check that this program is in the system.\n";
	    $message .= "POSSIBLE SOLUTIONS if it is in the system:\n";
	    $message .= "\t1) Move the executable to the PATH\n";
	    $message .= "\t2) Add the executable dir to path\n";
	    $message .= "\t3) Create a symlink of the exec. in the PATH\n";
	    croak("$message\n");
	}
    }
}

## Create the output dir.

my $cwd = getcwd();
my $outdir = File::Spec->catdir($cwd, $outbasename);
mkdir($outdir);
my %length_files = ();

## Get the dirs

my @sffdirs = split(/,/, $dirlist);

## Get each file of each dir

foreach my $dir (@sffdirs) {
    
    opendir my ($dh_dir), $dir;
    my @files = readdir($dh_dir);
    foreach my $file (sort @files) {
       
	if (defined $file && $file =~ m/\.sff$/) {
	    my $infile = File::Spec->catfile($dir, $file);
	    my ($vol, $dirs, $endpath) = File::Spec->splitpath($dir);
	   
	    my $outfile = File::Spec->catfile($outdir, $endpath .'-'. $file);
	    $outfile =~ s/\.sff$/\.tab/;

	    my $cmd = "sffinfo -acc -seq ";
	    if ($opt_N) {
		$cmd .= " -notrim ";
	    }
	    $cmd .= "$infile | grep '>' | cut -d ' ' -f1,2 |";
	    $cmd .= " sed -r 's/>//' | sed -r 's/ length=/\t/' > ";
	    $cmd .= " $outfile";
	    
	    print STDERR "\nRUNNING SYSTEM COMMAND:\n";
	    print STDERR "$cmd\n\n";

	    system($cmd);
	    exit();
	    $length_files{$file} = $outfile;
	}
    }
}

## Load this data in R

##my $n_sff = scalar(keys(%length_files));

##my $rbase = R::YapRI::Base->new();


##foreach my $fname (keys %length_files) {
##    my $rcmd = "$fname <- read.delim(\"$length_files{$fname}\", header=FALSE)";
##    $rbase->add_command($rcmd);
##}




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

      This script run sffinfo program over different sequence .sff sets 
    returning some reports for sequence distribution.
    (It runs only in linux OS)
    
    Usage:
        
      sffstats.pl -d <sffdirs> -o <output_basename>

    Example:
	
      sffstats.pl -d 454sht/-o test
      
    Flags:

      -d <sffdir>       List of dirs with sff files separated by ',' (mandatory)
      -o <output>       An output basename for the reports (mandatory)
      -N <no_trimmed>   Get no-trimmed sequences from sffinfo command
      -h <help>         print the help


     
EOF
exit (1);
}




###
1;#
###
