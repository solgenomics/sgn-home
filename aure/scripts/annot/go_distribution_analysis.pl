#!/usr/bin/perl

=head1 NAME

 go_distribution_analysis.pl
 Script to compare different GO Terms distributions

=cut

=head1 SYPNOSIS

 go_distribution_analysis.pl [-h] -r <reference_go_distribution_file> -t <test_go_distribution_files> [-o <output>] [-l <level_filter>] [-n <sequence_number_filter>]

=head1 EXAMPLE 

 go_distribution_analysis.pl [-h] -r global_go_distrib.txt -t group1_go.txt,group2_go.txt -o go_comparisson -l 2 -n 5
                                  
=head2 I<Flags:>

=over


=item -r

B<reference_go_file>           Reference GO Distribution file with 7 columns (mandatory)
                               (Level;GO_ID;Term;Type;#Seqs;Graph_Score;Sequences) 

=item -t

B<test_go_files>               Test GO Distribution files (separated by commas) with 7 columns (mandatory)
                               (Level;GO_ID;Term;Type;#Seqs;Graph_Score;Sequences) 

=item -o

B<output>                      Output basename (go_comparisons.output by default)

=item -l

B<level_filter>                Level used in the comparisson (all by default)

=item -n

B<seq_number_filter>           Minimun number of sequences used in the comparison (1 by default)

=item -h

B<help>                         Print the help

=back

=cut

=head1 DESCRIPTION

  This script parse the GO distribution output produced by Blast2GO (Gotz. et al. 2008), 
  (Analysis > Make Combine Graph > Export Graph Information as text), calculate the percentages
  after filtering and produce a table file (.txt), that can be used by R to represent as a 
  graph

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 go_distribution_analysis.pl

=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;


our ($opt_r, $opt_t, $opt_o, $opt_l, $opt_n, $opt_h);
getopts("r:t:o:l:n:h");
if (!$opt_r && !$opt_t && !$opt_o && !$opt_l && !$opt_n && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Mandatory arguments (die if they are not supplied)

my $ref = $opt_r 
    || die("\n\nMANDATORY ARGUMENT ERROR: -r <reference_go_file> argument was not supplied.\n\n");
my $samples = $opt_t
    || die("\n\nMANDATORY ARGUMENT ERROR: -t <test_go_files> argument was not supplied.\n\n");

## Default arguments (use default values if they are not supplied)

my $output = $opt_o 
    || 'go_comparisons.output';

my $seq_n = $opt_n 
    || 1;

my $level = $opt_l;


##############################################
#### PARSING FILE  ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);

print STDERR "\n\nSCRIPT START POINT ($date)\n\n";

$date = `date`;
chomp($date);
print STDERR "\t1) Parsing reference file: $ref ($date)\n\n";

my %go_ref = ();
my %ordered_ref = ();
my $total_go_ref = 0;
    
my $t = 0;
my $T = `wc -l $ref`;
chomp($T);
$T =~ s/$ref//;
    
open my $ref_fh, '<', $ref;
    
while(<$ref_fh>) {
    chomp($_);
    $t++;
	
    print STDERR "\t\tParsing line $t of $T lines in the file: $ref\r";

    unless ($_ =~ m/^Level/i) {
	
	my @t_data = split(/\t/, $_);
	my @seq_ids = split(/,\s/, $t_data[6]);
	my %data = ( 'level'       => $t_data[0],
		     'go_id'       => $t_data[1],
		     'term'        => $t_data[2],
		     'type'        => $t_data[3],
		     'seq_n'       => $t_data[4],
		     'graph_score' => $t_data[5],
		     'seq_ids'     => \@seq_ids,
	           );
	$total_go_ref += $data{'seq_n'};

	my $selected = 1;
	if (defined $level) {
	    unless ($level == $data{'level'}) {
		$selected = 0;
	    }
	}
	if (defined $seq_n) {
	    if ($data{'seq_n'} < $seq_n) {
		$selected = 0;
	    }
	}   
	if ($selected == 1) {
	    $go_ref{$data{'go_id'}} = \%data;
	    $ordered_ref{$data{'go_id'}} = $data{'level'};
	} 
    }
}
print STDERR "\n\n";


$date = `date`;
chomp($date);
print STDERR "\t2) Parsing files to compare ($date): $samples\n\n";

my @infiles = sort(split(/,/, $samples));

my %go_files = ();

foreach my $infile (@infiles) {

    my %go_data = ();
    my $total_go_file = 0;

    print STDERR "\t\t+ Parsing file to compare: $infile\n\n";

    open my $in_fh, '<', $infile;
    my $l = 0;
    my $L = `wc -l $infile`;
    chomp($L);
    $L =~ s/$infile//;

    while (<$in_fh>) {
	$l++;
	chomp($_);
	print STDERR "\t\tParsing line $l of $L lines in the file: $infile\r";
	
	unless ($_ =~ m/^Level/i) {

	    my @i_data = split(/\t/, $_);
	    my @seq_ids = split(/,\s/, $i_data[6]);
	    my %data = ( 'level'       => $i_data[0],
			 'go_id'       => $i_data[1],
			 'term'        => $i_data[2],
			 'type'        => $i_data[3],
			 'seq_n'       => $i_data[4],
			 'graph_score' => $i_data[5],
			 'seq_ids'     => \@seq_ids
		       );
	    $total_go_file += $data{'seq_n'};

	    my $selected = 1;
	    if (defined $level) {
		unless ($level == $data{'level'}) {
		    $selected = 0;
		}
	    }
	    if (defined $seq_n) {
		if ($data{'seq_n'} < $seq_n) {
		    $selected = 0;
		}
	    }   
	    if ($selected == 1) {
		$go_data{$data{'go_id'}} = \%data;
	    }
	}
    }
    close $in_fh;

    $go_files{$infile} = { 'data_parsed' => \%go_data, 'seq_total_n' => $total_go_file };

    print STDERR "\n\n";

}
print STDERR "\n\n";

$date = `date`;
chomp($date);
print STDERR "\t3) Analyzing Data.\n\n";

open my $o_fh, '>', $output . '.txt';

## Print headers.

print $o_fh "\t\"level\"\t\"term\"\t\"ref_perc\"";

foreach my $ifile (@infiles) {
    my $new_name = $ifile;
    $new_name =~ s/\./_/g;
    print $o_fh "\t\"$new_name\"";
}
print $o_fh "\n";


foreach my $go_id (sort {$ordered_ref{$a} <=> $ordered_ref{$b}} keys %ordered_ref) {
 
    my $lv = $go_ref{$go_id}->{'level'};
    my $term =  $go_ref{$go_id}->{'term'};
    my $seq_n = $go_ref{$go_id}->{'seq_n'};
    my $perc_seq = $seq_n * 100 / $total_go_ref;
    my $perc_seq_obj = Math::BigFloat->new($perc_seq);
    my $rd_perc_seq = $perc_seq_obj->bfround(-3);

    print $o_fh "\"$go_id\"\t$lv\t\"$term\"\t$rd_perc_seq";

    foreach my $comp_file (@infiles) {
	
	my $total_seqn = $go_files{$comp_file}->{'seq_total_n'};
	
	if (defined $go_files{$comp_file}->{'data_parsed'}) {
	    my %comp_godata = %{$go_files{$comp_file}->{'data_parsed'}};

	    if (defined $comp_godata{$go_id}) {
		my $comp_seq_n = $comp_godata{$go_id}->{'seq_n'};
		my $comp_perc_seq = $comp_seq_n * 100 / $total_seqn;
		my $comp_perc_seq_obj = Math::BigFloat->new($comp_perc_seq);
		my $comp_rd_perc_seq = $comp_perc_seq_obj->bfround(-3);
		
		print $o_fh "\t$comp_rd_perc_seq";
	    }
	    else {
		print $o_fh "\t0.000";
	    }
	}
	else {
	    print $o_fh "\t0";
	}
    }
    print $o_fh "\n";
}




$date = `date`;
chomp($date);
print STDERR "\n\nEND OF THE SCRIPT ($date).\n\n";


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
 
      This script parse the GO distribution output produced by Blast2GO (Gotz. et al. 2008), 
    (Analysis > Make Combine Graph > Export Graph Information as text), calculate the percentages
    after filtering and produce a table file (.txt), that can be used by R to represent as a 
    graph

    Usage:
  
     go_distribution_analysis.pl [-h] -r <reference_go_distribution_file> 
                                      -t <test_go_distribution_files> 
                                      [-o <output>] 
                                      [-l <level_filter>] 
                                      [-n <sequence_number_filter>]
                              
    Example:

      go_distribution_analysis.pl -r global_go_distrib.txt 
                                  -t group1_go.txt,group2_go.txt 
                                  -o go_comparisson 
                                  -l 2 
                                  -n 5

    Flags:
    
     -r <reference_go_file>     Reference GO Distribution file with 7 columns (mandatory)
                                 (Level;GO_ID;Term;Type;\#Seqs;Graph_Score;Sequences) 
     -t <test_go_files>         Test GO Distribution files (separated by commas) with 7 columns (mandatory)
                                 (Level;GO_ID;Term;Type;\#Seqs;Graph_Score;Sequences) 
     -o <output>                Output basename (go_comparisons.output by default)
     -l <level_filter>          Level used in the comparisson (all by default)
     -n <seq_number_filter>     Minimun number of sequences used in the comparison (1 by default)
     -h <help>                  Print the help


EOF
exit (1);
}




####
1; #
####
