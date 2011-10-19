#!/usr/bin/perl

=head1 NAME

454metricsparser.pl
A script to parse 454 metrics and produce 454 compatible files.

=cut

=head1 SYPNOSIS

454metricparser.pl -i <inputs>

=head2 I<Flags:>

=over


=item -i

B<input>                Input 454Metrics.txt files, separated by ',' (mandatory)

=item -h

B<help>                 print the help

=back

=cut

=head1 DESCRIPTION

 This script parse 454NewblerMetrics.txt files and produce some outputs 
 compatible with R. 

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 454metricsparser.pl

=cut

use strict;
use warnings;
use autodie;

use Carp qw( croak );

use File::Spec;
use Cwd;

use Getopt::Std;
use Math::BigFloat;
use File::Basename;

#use R::YapRI::Base;



our ($opt_i, $opt_h);
getopts("i:h");
if (!$opt_i && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}


## Check arguments

my $inlist = $opt_i || 
    die("DATA ARGUMENT ERROR: -i <input> WAS NOT SUPPLIED.\n");


## Get the files

my @metr_files = split(/,/, $inlist);

## Get each file of each dir

my %parsfiles = ();
my %sets = ();

## PARSE FILE

print STDERR "\n\nPARSING FILES:\n";

foreach my $in (@metr_files) {

    print STDERR "\n\tParsing file: $in\n";
    
    open my $infh, '<', $in;

    ## Declare my switchs
    my ($readdata_switch, $readdatafile_switch) = (0, 0);
    my ($alndata_switch, $alndatafile_switch) = (0, 0);

    my $fileset = '';

    my %readdata = ( 
	untrimmed_reads => 0,
	untrimmed_bases => 0,
	trimmed_reads   => 0,
	trimmed_bases   => 0,
	paired_reads    => 0,
	);

    my %alndata = (
	aligned_reads            => 0,
	aligned_reads_perc       => 0,
	aligned_bases            => 0,
	aligned_bases_perc       => 0,
	inferred_read_error_perc => 0,
	inferred_read_error      => 0,
	both_mapped              => 0,
	one_mapped               => 0,
	multiple_mapped          => 0,
	both_unmapped            => 0,
	);


    my $l = 0;
    while(<$infh>) {
	chomp($_);
	my $line = $_;
	$l++;

	if ($line =~ m/(runData|pairedReadData)/) {
	    $readdata_switch = 1;
	}
	if ($readdata_switch == 1 && $line =~ m/\s+file/) {
	    $readdatafile_switch = 1;
	}
	if ($readdatafile_switch == 1) {
	    if ($line =~ m/\s+path = "(.+)";/) {
		$fileset = File::Basename::basename($1);
	    }
	    elsif ($line =~ m/\s+numberOfReads\s+=\s+(\d+),\s+(\d+);/) {
		$readdata{untrimmed_reads} = $1;
		$readdata{trimmed_reads} = $2;
	    }
	    elsif ($line =~ m/\s+numberOfBases\s+=\s+(\d+),\s+(\d+);/) {
		$readdata{untrimmed_bases} = $1;
		$readdata{trimmed_bases} = $2;
	    }
	    elsif ($line =~ m/\s+numWithPairedRead\s+=\s+(\d+);/) {
		$readdata{paired_reads} = $1;
	    }
	}
	if ($line =~ m/(readAlignmentResults|pairedReadResults)/) {
	    $alndata_switch = 1;
	}
	if ($alndata_switch == 1 && $line =~ m/\s+file/) {
	    $alndatafile_switch = 1;
	}
	if ($alndatafile_switch == 1) {

	    if ($line =~ m/\s+path = "(.+)";/) {
		$fileset = File::Basename::basename($1);
	    }
	    elsif ($line =~ m/\s+numAlignedReads\s+=\s+(\d+),\s+(.+)%;/) {
		$alndata{aligned_reads} = $1;
		$alndata{aligned_reads_perc} = $2;
	    }
	    elsif ($line =~ m/\s+numAlignedBases\s+=\s+(\d+),\s+(.+)%;/) {
		$alndata{aligned_bases} = $1;
		$alndata{aligned_bases_perc} = $2;
	    }
	    elsif ($line =~ m/\s+inferredReadError\s+=\s+(.+)%,\s+(\d+);/) {
		$alndata{inferred_read_error_perc} = $1;
		$alndata{inferred_read_error} = $2;
	    }
	    elsif ($line =~ m/\s+numberWithBothMapped\s+=\s+(\d+);/) {
		$alndata{both_mapped} = $1;
	    }
	    elsif ($line =~ m/\s+numWithOneUnmapped\s+=\s+(\d+);/) {
		$alndata{one_mapped} = $1;
	    }
	    elsif ($line =~ m/\s+numWithMultiplyMapped\s+=\s+(\d+);/) {
		$alndata{multiple_mapped} = $1;
	    }
	    elsif ($line =~ m/\s+numWithBothUnmapped\s+=\s+(\d+);/) {
		$alndata{both_unmapped} = $1;
	    }
	}



	if ($line =~ m/^}\s*$/) {
	    $readdata_switch = 0;
	    $alndata_switch = 0;
	}
	if (defined $fileset && $line =~ m/^\s+}\s*$/) {
	    
	    if ($readdatafile_switch == 1) {
		$readdatafile_switch = 0;

		## Add to the hash
		unless ($sets{$fileset}) {
		    $sets{$fileset} = { $in => {} };		
		}
		foreach my $key (keys %readdata) {
		    $sets{$fileset}->{$in}->{$key} = $readdata{$key};
		    $readdata{$key} = 0;
		}

		$fileset = undef;
	    }
	    elsif ($alndatafile_switch == 1) {

		$alndatafile_switch = 0;
		
		## The hash should be created in the previous steps
		
		foreach my $alnkey (keys %alndata) {
		    $sets{$fileset}->{$in}->{$alnkey} = $alndata{$alnkey};
		    $alndata{$alnkey} = 0;
		}
		$fileset = undef;
	    }
	}
    }
    close($infh);
    print STDERR "\n\tDone.\n";
}


## PRINT RESULTS

print STDERR "\n\nPRINTING RESULTS:\n";

## Define the order for print results

my @fields = ( 'untrimmed_reads', 'untrimmed_bases', 'trimmed_reads', 
	       'trimmed_bases', 'paired_reads', 'aligned_reads', 
	       'aligned_reads_perc', 'aligned_bases', 'aligned_bases_perc',
	       'inferred_read_error_perc', 'inferred_read_error', 
	       'both_mapped', 'one_mapped', 'multiple_mapped', 'both_unmapped');

my $outfile = '454NewblerMetrics_StatsForR.tab';

open my $ofh, '>', $outfile;

print STDERR "\n\tResult file: $outfile created\n";

## Print the headers

my @headline = ('UTcRD', 'UTcBS', 'TRcRD', 'TRcBS', 'P2cRD', 'ALcRD', 'ALpRD',
		'ALcBS', 'ALpBS', 'IpERR', 'IcERR', 'P2MP2', 'P2MP1', 'P2MPm', 
		'P2UMP');

my @headerline = ();
my $fn = 0;
foreach my $ifle (@metr_files) {
    $fn++;
    foreach my $head (@headline) {
	push @headerline, $head . '_' . $fn;
    }
}
print $ofh join("\t", @headerline);
print $ofh "\n";

foreach my $fileset (sort keys %sets) {
     
    print $ofh "$fileset\t";
    my @datafile;
    foreach my $infile (@metr_files) {	
	
	unless (exists $sets{$fileset}->{$infile}) {
	    foreach my $field (@fields) {
		push @datafile, 'NA';
	    }
	}
	else {
	    foreach my $field (@fields) {

		my $data = $sets{$fileset}->{$infile}->{$field} || 0;
		push @datafile, $data;
	    }
	}	
    }
    my $dataline = join("\t", @datafile);
    print $ofh "$dataline\n";
}

print STDERR "\n\tDone.\n\n";


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

      This script parse 454NewblerMetrics.txt files and produce some outputs 
      compatible with R. 
    
    Usage:
        
      454metricparser.pl [-h] -i <inputs>

    Example:
	
      454metricparser.pl -i 'assembly1metrics.txt,assembly2metrics.txt'
      
    Flags:

      -i <input>      Input 454Metrics.txt files, separated by ',' (mandatory)
      -h <help>       print the help


     
EOF
exit (1);
}




###
1;#
###
