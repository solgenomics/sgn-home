#!/usr/bin/perl

=head1 NAME

 interpro2go.pl
 Script to extract go terms from interpro result raw files.

=cut

=head1 SYPNOSIS

 interpro2go.pl [-h] -i <input_file> -o <output_file> -D
                     [-m <filter_by_method>]
                     [-l <filter_by_length_domain_match>]
                     [-e <filter_by_evalue>]

=head1 EXAMPLE 

 go_distribution_analysis.pl -i myfile.txt -o go.txt -m 'HMMPIR,HMMPanther'
                                  
=head2 I<Flags:>

=over


=item -i

B<input_file>                  input file in raw format (mandatory)

=item -o

B<output_file>                 output filename (mandatory) 

=item -m

B<filter_by_method>            method name (or names) to filter the file

=item -l

B<filter_by_length_match>      length of the domain match to use as filter

=item -e

B<filter_by_evalue>            evalue used as filtering cutoff

=item -D

B<print_go_descriptions>       add go description to the go terms

=item -h

B<help>                        print the help

=back

=cut

=head1 DESCRIPTION

  This script parse the interpro result files in rawn format and produces a 
  output file with each of the sequences and a list (separated by semicolon)
  of GO terms (<ID><tab><GOID1><;><GOID2><;>...).

  If -D option is used it will print: <ID><tab><GODescrip><space><GOID1><;>...

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 interpro2go.pl

=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;


our ($opt_i, $opt_o, $opt_m, $opt_l, $opt_e, $opt_D, $opt_h);
getopts("i:o:m:l:e:Dh");
if (!$opt_i && !$opt_o && !$opt_m && !$opt_l && !$opt_e && !$opt_D && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my $input = $opt_i ||
    die("\nARGUMENT ERROR: No input file was supplied to this script.\n");

my $output = $opt_o ||
    die("\nARGUMENT ERROR: No output filename was supplied to this script.\n");

my $f_meth = $opt_m;
my $f_leng = $opt_l;
my $f_eval = $opt_e;


##############################################
#### PARSING FILE  ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);

print STDERR "\n\n**************************************************";
print STDERR "\nSCRIPT START POINT ($date)\n";
print STDERR "**************************************************\n";

my %golist = ();
my %id_idx = ();
my $idx = 1;

print STDERR "\nParsing file ($date)\n\n";

    
my $t = 0;
my $T = `cut -f1 $input | wc -l`;
chomp($T);
    
open my $ifh, '<', $input;
    
while(<$ifh>) {
    chomp($_);
    $t++;
	
    print STDERR "\t\tParsing line $t of $T lines in the file: $input\r";

    my @data = split(/\t/, $_);
    my $id = $data[0];
    my $goline = $data[13];

    if (defined $id) {
	unless (exists $id_idx{$id}) {
	    $id_idx{$id} = $idx;
	    $idx++;
	}
    }

    if (defined $id && defined $goline) {

	## Go terms are separated by comma, but sometimes there are some commas
	## inside the go description. It will replace the commas close to 
	## the definition with semicolons

	$goline =~ s/, Molecular Function:/; Molecular Function:/g;
	$goline =~ s/, Biological Process:/; Biological Process:/g;
	$goline =~ s/, Cellular Component:/; Cellular Component:/g;
    
	my @go = split(/; /, $goline);

	my $selected = 1;
    
	## Use filters
    
	if (defined $f_meth) {
	    my $match = 0;
	    my @meths = split(/,/, $f_meth);
	    foreach my $meth (@meths) {
		if ($meth =~ m/$data[3]/i) {
		    $match = 1;
		}
	    }
	    if ($match == 0) {
		$selected = 0;
	    }
	}
	if (defined $f_leng) {
	    if ($data[7] =~ m/^\d+$/ && $data[6] =~ m/^\d+$/) { ## Verify line
		my $domlen = $data[7] - $data[6];
		if ($domlen < $f_leng) {
		    $selected = 0;
		}
	    }
	}
	if (defined $f_eval) {
	    my $e_val = Math::BigFloat->new($data[8]);
	    my $cutoff = Math::BigFloat->new($f_eval);
	    my $test = $e_val->bcmp($cutoff);

	    if ($e_val->is_pos()) {
		if ($e_val->bcmp($cutoff) > 0 ) {
		    $selected = 0;
		}
	    }
	}
	
	## If selected still have a value of 1

	if ($selected == 1) {
	    foreach my $go (@go) {

		if ($go =~ m/^(\w+\s+\w+):(.+)\s+\((GO:\d+)\)/) {
		    my $gotype = $1;
		    my $descrip = $2;
		    my $goid = $3;
		    
		    if (exists $golist{$id}) {
			if (exists $golist{$id}->{$gotype}) {
			    $golist{$id}->{$gotype}->{$goid} = $descrip;
			}
			else {
			    $golist{$id}->{$gotype} = {$goid => $descrip};
			}
		    }
		    else {
			$golist{$id} = {$gotype => { $goid => $descrip } };
			
		    }
		}
	    }
	}
    }
}

$date = `date`;
chomp($date);

print STDERR "\n\nPrinting output file ($date)\n\n";

open my $ofh, '>', $output;

foreach my $id (sort {$id_idx{$a} <=> $id_idx{$b}} keys %id_idx) {
    if (defined $golist{$id}) {
	
	my @ord_types = ( 'Biological Process', 
			  'Molecular Function', 
			  'Cellular Component' );

	my @list = ();

	foreach my $type (@ord_types) {
	
	    if (exists $golist{$id}->{$type}) {
		
		my %goclass = %{$golist{$id}->{$type}};
		foreach my $go (keys %goclass) {
		    if ($opt_D) {
			push @list, "$goclass{$go} $go";
		    }
		    else {
			push @list, "$go";
		    }
		}
	    }
	}
	my $list_line = join('; ', @list);
	print $ofh "$id\t$list_line\n";
    }
}


$date = `date`;
chomp($date);
print STDERR "\n**************************************************";
print STDERR "\nEND OF THE SCRIPT ($date)\n";
print STDERR "**************************************************\n\n";


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
 
       This script parse the interpro result files in rawn format and produces 
     an output file with each of the sequences and a list (separated by 
     semicolon) of GO terms (<ID><tab><GOID1><;><GOID2><;>...).

        If -D option is used it will print: 
        <ID><tab><GODescrip><space><GOID1><;>...

    Usage:
  
      interpro2go.pl [-h] -i <input_file> -o <output_file> -D
                     [-m <filter_by_method>]
                     [-l <filter_by_length_domain_match>]
                     [-e <filter_by_evalue>]
    
    Example:

      go_distribution_analysis.pl -i myfile.txt -o go.txt -m 'HMMPIR,HMMPanther'
    
    Flags:
    
     -i <input_file>                 input file in raw format (mandatory)
     -o <output_file>                output filename (mandatory) 
     -m <filter_by_method>           method name (or names) to filter the file
     -l <filter_by_length_match>     length of the domain match to use as filter
     -e <filter_by_evalue>           evalue used as filtering cutoff
     -D <print_go_descriptions>      add go description to the go terms
     -h <help>                       print the help



EOF
exit (1);
}




####
1; #
####
