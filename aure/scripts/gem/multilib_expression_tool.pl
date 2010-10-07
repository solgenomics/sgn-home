#!/usr/bin/perl

=head1 NAME

 multilib_expression_tool.pl
 Get the contig composition and calculate expression values

=cut

=head1 SYPNOSIS

 multilib_expression_tool.pl -i <input_file> -l <library_files> -o <output_file>

=head1 EXAMPLE 

 multilib_expression_tool.pl -i contig_read.tab -l lib_H01.tab,lib_H02.tab,lib_P01.tab,lib_P02.tab

=head2 I<Flags:>

=over


=item -i

B<input_file>                   input file, two columns -f1 contig name and -f2 sequence_id (mandatory)

=item -l

B<library_file>                 library files with -f1 sequence_id and -f2 library id in the second (mandatory)

=item -o

B<output_file>                  output_file (by default: <input_file>.analyzed.tab)

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

    This script get the composition of each contig, by libraries and calculate the abundance of them. 

    Produce an output file with:
      -f1: contig_name
      -f2: n_total_reads
      -f(3,4,5)xL: library name, n_reads_in_this_contig_for_this_library and normalized value in ppm 
                   (N_contigreads_lib / N_read_lib * 1000000)

      -f_last: R value calculated according: 

       Stekel DJ, Git Y and Faciani F.
       The comparison of gene expression from multiple cDNA libraries.
       Genome Res. 2000 Dec;10(12):2055-61

    According data randomization:
    
        +====+=================+
	| R  |	Believability  |  				     
	+====+=================+
	| 13 |           99.9% |
	| 12 |           99.9% |
	| 11 |           99.8% |
	| 10 |           99.7% |
	|  9 |           99.0% |
	|  8 |           98.2% |
	|  7 |           97.0% |
	|  6 |           91.5% |
	|  5 |           86.3% |
	|  4 |           82.2% |
	|  3 |           57.8% |
	|  2 |           26.8% |
        |  1 |           46.8% |
	+====+=================+


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

use File::Basename;
use Getopt::Std;
use Math::BigFloat;

our ($opt_i, $opt_l, $opt_o, $opt_h);
getopts("i:l:o:h");
if (!$opt_i && !$opt_l && !$opt_o && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my $input_file = $opt_i || die("\nMANDATORY ARGUMENT ERROR: -i <input_file> argument was not supplied.\n");
my $library_files = $opt_l || die("\nMANDATORY ARGUMENT ERROR: -l <library_files> argument was not supplied.\n");

my $output_file = $opt_o || $input_file . '.expression_analysis.tab';

## Parse the library files

my %lib_reads = ();
my %lib_counts = ();

my @library_files = split(/,/, $library_files);

print STDERR "\nPARSING LIBRARY FILES:\n\n";

foreach my $library_file (@library_files) {

    open my $lfh, '<', $library_file || die("\nSorry, I can not open the library file: $library_file.\n");

    my $l = 0;
    while (<$lfh>) {
	chomp($_);
	my @data = split(/\t/, $_);
	$l++;

	print STDERR "\tParsing line:$l from file:$library_file                                   \r";
	## Skip the row if it have more than two columns and print a warning message 

	if (scalar(@data) != 2) {
	    print STDERR "\nPARSING WARNING: line=$l in the file=$library_file have more than two columns. Skip parse.\n";
	}
	else {

	    ## Add the data to the lib_read hash

	    unless (exists $lib_reads{$data[1]}) {
		$lib_reads{$data[1]} = $data[0];

		## Add the data to the lib_count reads

		unless (exists $lib_counts{$data[0]}) {
		    $lib_counts{$data[0]} = 1;
		}
		else {
		    $lib_counts{$data[0]}++;
		}
 
	    }
	    else {
		print STDERR "\nPARSING WARNING: Seq_id=$data[1] was parsed before. Skip parse for line:$l file:$library_file.\n";
	    }
	}
    }
    print STDERR "\n\n";
    close $lfh;
}


my %contig_lib_count = ();
foreach my $lib (keys %lib_counts) {
    $contig_lib_count{$lib} = 0;
}


print STDERR "\n\nPARSING CONTIG COMPOSITION FILE:\n\n";

open my $ifh, '<', $input_file || die("\nSorry, I can not open the input file: $input_file.\n");

my %contig_read_count = ();
my %contig_comp = ();

my $m = 0;
while (<$ifh>) {
    chomp($_);
    my @idata = split(/\t/, $_);
    $m++;
    
    print STDERR "\tParsing line:$m from file:$input_file                                   \r";
    my $ilib = $lib_reads{$idata[1]}; 

    if (defined $ilib) {
	unless (exists $contig_comp{$idata[0]}) {      
	    $contig_comp{$idata[0]} = { $ilib => 1 };
	    $contig_read_count{$idata[0]} = 1;
	}
	else {
	    my $ilib_href = $contig_comp{$idata[0]};
	    $contig_read_count{$idata[0]}++;

	    if (defined $ilib_href->{$ilib}) {
		$ilib_href->{$ilib}++;
	    }
	    else {
		$ilib_href->{$ilib} = 1;
	    }
	}
    }
    else {
	print STDERR "\nPARSING ERROR: read:$idata[1] in line:$m have not parsed in the library files.\n";
    }
}

print STDERR "\n\nCHECKING THE LIBRARY PARAMETERS TO APPLY THE Stekel Formula:\n\n";

my $non_valid = 0;
foreach my $t_lib (sort keys %lib_counts) {
    my $t_lib_n = $lib_counts{$t_lib};
    
    if ($t_lib_n < 51) {
	$non_valid = 1;
	print STDERR "\nThe library: $t_lib has less than 50 reads. Please remove from the analysis.\n";
    }

    my $lib20_perc = $t_lib_n/5;

    foreach my $c_name (sort keys %contig_comp) {
	my $t_lib_href = $contig_comp{$c_name};
	my $c_read_count = $t_lib_href->{$t_lib} || 0;
	
	if ($c_read_count > $lib20_perc) {
	    $non_valid = 1;
	    print STDERR "\nThe contig:$c_name represent more than 20 percentage of the library:$t_lib. Please remove from analysis.\n";
	}
    }    
}

if ($non_valid == 1) {
    print STDERR "\n\tWARNING: One or some of the libraries do not satisfy the conditions to apply Sketel Formula for abundance.\n";
}
else {
     print STDERR "\n\tCHEKING OK: All the libraries satisfy the conditions to apply Sketel Formula for abundance.\n";
}


print STDERR "\n\nCALCULATING AND PRINTING OUTPUT FILE:\n\n";



my @lib_analyzed = sort keys %lib_counts;

my $total_reads = 0;
foreach my $lib (@lib_analyzed) {
    my $lib_read_n = $lib_counts{$lib};
    $total_reads += $lib_read_n;
}
my $N = Math::BigFloat->new($total_reads);

open my $out_fh, '>', $output_file ||
    die("Sorry, I can not open the file: $output_file.\n");

my $C = scalar(keys %contig_comp);
my $c = 0;

foreach my $contig_comp_name (sort keys %contig_comp) {
    my $c_lib_href = $contig_comp{$contig_comp_name};

    my ($R, $r_value);
    $c++;

    my $total_read_count_by_contig = $contig_read_count{$contig_comp_name};
    my $Xi = Math::BigFloat->new($total_read_count_by_contig);

    print $out_fh "$contig_comp_name\t$total_read_count_by_contig\t";
    print STDERR "\tprocessing contig: $contig_comp_name ($c of $C)                            \r";

    foreach my $lib_analyzed_name (@lib_analyzed) {
	my $contig_read_count = $c_lib_href->{$lib_analyzed_name} || 0;
	my $Xij = Math::BigFloat->new($contig_read_count);

	my $total_read_count_by_lib = $lib_counts{$lib_analyzed_name};
	my $Ni = Math::BigFloat->new($total_read_count_by_lib);

	my $ppm_coef = Math::BigFloat->new('1000000');
	my $Xij_c = $Xij->copy();
	$Xij_c->bdiv($Ni);
	$Xij_c->bmul($ppm_coef);
	my $norm_value = $Xij_c->bfround(-2);

	print $out_fh "$lib_analyzed_name\t$contig_read_count\t$norm_value\t";

	## Now it will calculate the values to get R (Stekel DJ et al. Genome Res. 2000)

	my $fi = $Xi->copy();
	$fi->bdiv($N);

	if ($non_valid == 0) {
	    my $log_Cij = $Xij->copy();
	    unless ($log_Cij->is_zero) {
		$log_Cij->bdiv($Ni);
		$log_Cij->bdiv($fi);
		$log_Cij->blog('10');

		my $Xij_log_Cij = $Xij->copy();
		$Xij_log_Cij->bmul($log_Cij);	
       	
		if (defined $R) {
		    $R = $Xij_log_Cij->badd($R);
		}
		else {
		    $R = $Xij_log_Cij->copy();
		}
	    }
	}
    }
    if ($non_valid == 0) {
	$r_value = $R->bfround(-2);
    }
    else {
	$r_value = 'NA';
    }
    print $out_fh "$r_value\n";
}
print STDERR "\n\nDONE\n\n";


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
      This script get the composition of each contig, by libraries and calculate the abundance of them. 

      Produce an output file with:
      -f1: contig_name
      -f2: n_total_reads
      -f(3,4,5)xL: library name, n_reads_in_this_contig_for_this_library and normalized value in ppm 
                   (N_contigreads_lib / N_read_lib * 1000000)

      -f_last: R value calculated according: 

       Stekel DJ, Git Y and Faciani F.
       The comparison of gene expression from multiple cDNA libraries.
       Genome Res. 2000 Dec;10(12):2055-61

      According data randomization:
        +====+=================+
	| R  |	Believability  |  				     
	+====+=================+
	| 13 |           99.9% |
	| 12 |           99.9% |
	| 11 |           99.8% |
	| 10 |           99.7% |
	|  9 |           99.0% |
	|  8 |           98.2% |
	|  7 |           97.0% |
	|  6 |           91.5% |
	|  5 |           86.3% |
	|  4 |           82.2% |
	|  3 |           57.8% |
	|  2 |           26.8% |
        |  1 |           46.8% |
	+====+=================+

    Usage:
      multilib_expression_tool.pl [-h] -i <input_file> -l <library_file> [-o <output_file>]

    Example:
      multilib_expression_tool.pl -i contig_read.tab -l lib_H01.tab,lib_H02.tab,lib_P01.tab,lib_P02.tab

    Flags:
      -i <input_file>         input file in tab format with -f1 contig_name and -f2 sequence_id (mandatory)
      -l <library_file>       library files with -f1 sequence_id and -f2 library id in the second (mandatory)
      -o <output_file>        output_file (by default: <input_file>.analyzed.tab)
      -h <help>               print the help

EOF
exit (1);
}

