#!/usr/bin/perl

=head1 NAME

 phygomics.pl
 Script to run PhyGomics pipeline.

=cut

=head1 SYPNOSIS

 phygomics.pl [-h] -c config.file -i input.dir -o output.dir [-C]

=head1 EXAMPLE 

 phygomics.pl -c phygo1.txt -i analysis1 -o results1

 ## To print an ampty configuration file:

 phygomics.pl -C
                                  
=head2 I<Flags:>

=over


=item -c

B<configuration.file>          Configuration file printed using -C.

=item -i

B<input.dir>                   Input dir. with input datasets (mandatory)

=item -o

B<output.dir>                  Output dir. to create output datasets (mandatory)

=item -h

B<help>                        Print the help

=item -C

B<create_config.file>         Create an empty configuration file with the name:
                              phygomics.conf 
=back

=cut

=head1 DESCRIPTION

   phygomics.pl is the scripts that executes the PhyGomics pipeline.

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 phygomics.pl

=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;


our ($opt_c, $opt_i, $opt_o, $opt_h, $opt_C);
getopts("c:i:o:h");
if (!$opt_c && !$opt_i && !$opt_o && !$opt_h && !$opt_C) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Mandatory arguments (die if they are not supplied)

my $err_m = "MANDATORY ARGUMENT ERROR:";

my $control = $opt_c 
    || die("\n\n$err_m -c <control_file> argument was not supplied.\n\n");

my $indir = $opt_i
    || die("\n\n$err_m -i <input_dir> argument was not supplied.\n\n");

my $outdir = $opt_o
    || die("\n\n$err_m -o <out_dir> argument was not supplied.\n\n");



##############################################
#### PARSING FILE  ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);



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
 
      phygomics.pl is the scripts that executes the PhyGomics pipeline.       

    Usage:
  
     phygomics.pl [-h] -c config.file -i input.dir -o output.dir [-C]
                              
    Example:

     phygomics.pl -c phygo1.txt -i analysis1 -o results1

    Flags:
     
      -c <configuration.file>  Configuration file printed using -C.
      -i <input.dir>           Input dir. with input datasets (mandatory)
      -o <output.dir>          Output dir. to create output datasets (mandatory)
      -h <help>                Print the help
      -C <create_config.file>  Create an empty configuration file with the name:
                               phygomics.conf 

EOF
exit (1);
}




####
1; #
####
