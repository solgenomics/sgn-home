#!/usr/bin/perl

=head1 NAME

 run_signalp.pl
 A tool to run signalp using Bio::Tools::Run::Signalp

=cut

=head1 SYPNOSIS

 run_signalp.pl -s <specie> -i <input_sequence_file> [-o <output_results>] [-b <sequence_set_breakpoint>] [-a <signalp_arguments>]
 
=head2 I<Flags:>


=over

=item -s

B<specie>                       specie to pass to signalp program -t euk|gram+|gram- (euk by default)

=item -i

B<input_sequence_file>          input sequence file in fasta format (mandatory)

=item -o

B<output_file>                  output_file (by default: <input_file>.outsignalp)

=item -b

B<seq_break_number>             how many sequences should have the temp files to be used by signalP program (500 by default) 

=item -a

B<signalp_arguments>            transfer signalp arguments, for example '-trunc 50'

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

   This script is a perl wrapper to use signalp program. 

   It need set the location of the signalp program using environment variables
    export SIGNALP='location'

   Please check the signalp manual to know how works.
   (http://www.cbs.dtu.dk/cgi-bin/nph-runsafe?man=signalp)

   Note: 
      SignalP give problems for large datasets, so in that case is better use -b option. 
      By default the program said that it can not use more than 4000 sequence but this value can be less depending
      of the sequence lenght.

      -a <-trunc 50 -f summary> are recommended

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 multilib_expression_tool.pl


=cut

use strict;
use warnings;


use Bio::SeqIO;
use File::Basename;
use File::Temp qw/ tempfile tempdir /;
use Getopt::Std;


our ($opt_s, $opt_i, $opt_o, $opt_a, $opt_b, $opt_h);
getopts("s:i:o:a:b:h");
if (!$opt_s && !$opt_i && !$opt_o && !$opt_a && !$opt_b && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my $sp = $opt_s || 'euk';
my $input_file = $opt_i || die("\nMANDATORY ARGUMENT ERROR: -i <input_file> argument was not supplied.\n");
my $output_file = $opt_o || $input_file . '.signalp_results.txt';
my $break = $opt_b || '500';

unless ( defined $ENV{SIGNALP} ) {
    die("ENVIRONMENT VARIABLE ERROR: The SIGNALP env variable is not set.\n");
}

## Create the temp dir

my $tempdir = File::Temp->newdir();

my $in  = Bio::SeqIO->new( 
                           -file   => $input_file,
                           -format => 'Fasta'
                         );



## First split the file in files with no more than 4000 sequences

my $n = 0;
my $t = 1;
my $tmp_fasta;
my @tmp_fasta_files;


print STDERR "\n\n1) Create the subfasta files.\n\n";

while ( my $seq = $in->next_seq() ) {
    my $id = $seq->id();
    my $sequence = $seq->seq();
	
    if ($n == 0) {

	$tmp_fasta = File::Temp->new( 
	                              TEMPLATE => 'temp_XXXXX',
			       	      DIR      => $tempdir,
				      SUFFIX   => '.fasta',
	                              UNLINK   => 0,
	                            );
	
	print $tmp_fasta ">$id\n$sequence\n";
	$n++;
	my $tmp_filename = $tmp_fasta->filename();
	push @tmp_fasta_files, $tmp_filename;

    }
    elsif ($n < $break) {
	print $tmp_fasta ">$id\n$sequence\n";
	$n++;
    }
    else {
	print $tmp_fasta ">$id\n$sequence\n";
	$n = 0;
    }
}

my $tf_n = scalar(@tmp_fasta_files);
print STDERR "\t...done ($tf_n subsequences files has been created).\n";


## Second run the program over each file and put then in another temp file before join them

print STDERR "\n\n2) Execute signalp over the subfasta files.\n\n";

my @tmp_signalp_files;
my $c = 0;

foreach my $tmp_fasta_file (@tmp_fasta_files) {
    $c++;
    my $tmp_signalp = File::Temp->new( 
	                               TEMPLATE => 'temp_XXXXX',
		  	               DIR      => $tempdir,
			               SUFFIX   => '.signalp.txt',
	                               UNLINK   => 0,
	                             );

    my $run_command = $ENV{SIGNALP} . " -t " . $sp;
    
    if (defined $opt_a) {
	$run_command .= ' ' . $opt_a;
    }

    $run_command .= ' ' . $tmp_fasta_file . ' >> ' . $tmp_signalp->filename();

    my $status = system($run_command);
    if ($status != 0) {
	print STDERR "\t($c file of $tf_n)\n\tPROGRAM EXECUTION ERROR for ($run_command) by system error: ${^CHILD_ERROR_NATIVE} \n";
    }
    else {
	print STDERR "\t($c file of $tf_n) Execute:\n\t$run_command.\n";
    }

    push @tmp_signalp_files, $tmp_signalp->filename();
}
print STDERR "\t...done.\n";


## Third join all the results

print STDERR "\n\n3) Join the signalp results files.\n\n";

my $join_command;
if (scalar(@tmp_signalp_files) > 1) {
    $join_command = "cat " . join(' ', @tmp_signalp_files) . ' > ' . $output_file;
}
else {
    $join_command = "cp " . $tmp_signalp_files[0] . " " . $output_file;
}

my $status = system($join_command);
if ($status != 0) {
    print STDERR "\tPROGRAM EXECUTION ERROR for ($join_command) by system error: $? \n";
}
else {
    print STDERR "\t$join_command executed.\n";
}

print STDERR "\t...done.\n";



## Parse the output if -f summary was used

print STDERR "\n\n4) Creation of tab file (if summary argument was used).\n\n";

my $output_tab = $output_file . ".simply.tab";
if ($opt_a =~ m/summary/) {
    open my $ofh, '<', $output_file ||
	die("Sorry, I can not open the file:$output_file.\n");

    open my $new_out, '>', $output_tab ||
	die("Sorry, I can not open the file:$output_tab.\n");

    my $sep = 0;
    my ($seqid, $prediction, $signal_pep_prob, $signal_anchor_prob);

    while (<$ofh>) {
	chomp($_);
	if ($_ =~ m/^-+/) {
	    $sep = 0;
	    if (defined $seqid) {
		print $new_out "$seqid\t$prediction\t$signal_pep_prob\t$signal_anchor_prob\n";
		($seqid, $prediction, $signal_pep_prob, $signal_anchor_prob) = (undef, undef, undef, undef);
	    }
	}
	else {
	    if ($_ =~ m/^>(.+)$/) {
		$seqid = $1;
	    }
	    elsif ($_ =~ m/^Prediction:\s+(.+)$/) {
		$prediction = $1;
	    }
	    elsif ($_ =~ m/^Signal peptide probability:\s+(.+)$/) {
		$signal_pep_prob = $1;
	    }
	    elsif ($_ =~ m/^Signal anchor probability:\s+(.+)$/) {
		$signal_anchor_prob = $1;
	    }
	}
    }
}
    
print STDERR "\t...done.\n\n";

## Finally delete the temp files

print STDERR "\n\n5) Delete the temp files.\n\n";


print STDERR "\t...done.\n\n";



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

      This script is a perl wrapper to use signalp program. 

      It need set the location of the signalp program using environment variables

      Please check the signalp manual to know how works.
      (http://www.cbs.dtu.dk/cgi-bin/nph-runsafe?man=signalp)

    Note: 
      SignalP give problems for large datasets, so in that case is better use -b option. 
      By default the program said that it can not use more than 4000 sequence but this value can be less depending
      of the sequence lenght.

      -a <-trunc 50 -f summary> are recommended

    Usage:

       run_signalp.pl -s <specie> -i <input_sequence_file> [-o <output_results>] [-b <sequence_set_breakpoint>] [-a <signalp_arguments>]
     
    Flags:
      -s <specie>             specie to pass to signalp program -t euk|gram+|gram- (euk by default)
      -i <input_file>         input file in tab format with -f1 contig_name and -f2 sequence_id (mandatory)
      -o <output_file>        output_file (by default: <input_file>.analyzed.tab)
      -b <seq_break_number>   how many sequences should have the temp files to be used by signalP program (500 by default) 
      -h <help>               print the help

EOF
exit (1);
}

