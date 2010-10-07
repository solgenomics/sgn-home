#!/usr/bin/perl

=head1 NAME

 frg_format_tool.pl
 Multiple functions to work with frg files from wgs-assembler (version.0.1)

=cut

=head1 SYPNOSIS

 frg_format_tool.pl [-h] -i <input_files> -o <output_basename> [-f <filter_file>] [-t <tag_file>] [-S] [-Q]

=head2 I<Flags:>

=over


=item -i

B<input_files>            frg input files inside quotes and separated by commas (mandatory)

=item -o

B<output_basename>        basename for the output files (mandatory)

=item -f

B<filter_file>            filter file with a list of sequences id to remove from the input files

=item -t

B<tag_file>               tag file in tab format (-f1 id, -f2 tag, -f3 value)

=item -S

B<get_sequence>           get sequence from the frg files to create a fasta file

=item -Q

B<get_qscores>            get qscores from the frg files to create a qual file

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script can be used to modify the information in frg files. This is the format used by wgs-assembler 
 (http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=FRG_Files). Check this page to get
 more informatio about frg format and the fields/tags used.

 To convert fasta, qual, sff or trace files to frg files use wgs scripts (convert-fasta-to-v2.pl, 
 tracedb-to-frg.pl, sffToCA, or fastqToCA) 

 Examples: 

    ## To remove a list of sequences

    frg_format_tool.pl -i test.frg -o outputname -f filter_test.txt

    ## To change tags in the frg file (example of tag_test.txt line: TEST_ID   sta    V)

    frg_format_tool.pl -i test.frg -o outputname -t tag_test.txt

    ## To get sequences

    frg_format_tool.pl -i test.frg -o outputname -S

    ## To get quality

    frg_format_tool.pl -i test.frg -o outputname -Q

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 frg_format_tool.pl


=cut

use strict;
use warnings;

use Getopt::Std;
use Bio::Seq::Quality;
use Bio::SeqIO;

our ($opt_i, $opt_o, $opt_f, $opt_t, $opt_S, $opt_Q, $opt_h);
getopts("i:o:f:t:SQh");
if (!$opt_i && !$opt_o && !$opt_f && !$opt_t && !$opt_S && !$opt_Q && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
elsif (!$opt_f && !$opt_t && !$opt_S && !$opt_Q) {
    die("DATA ARGUMENT ERROR: one of this argument should be supplied to this script: -f, -t, -S, or -Q.\n");
}

if ($opt_h) {
    help();
}

my $input_files = $opt_i || die("DATA ARGUMENT -i input_file WAS NOT SUPPLIED.\n");
my $output_name = $opt_o || die("DATA ARGUMENT -o output_name WAS NOT SUPPLIED.\n");


## First get the names of the input files as an array

my @infiles = split(/,/, $input_files);

## Get the data nad put in the memory as hashes spend a lot of memory, so it will
## process the data on the fly, storing only the filter data to be used in the parsing
## of frg files

## Filter option:

print STDERR "\n\nFILTER BY ID:";

## Define the excluided ids hash

my %excl_ids;
my $outf_fh;

if ($opt_f) {

    print STDERR " [Enabled]\n\n";

    open my $ffh, '<', $opt_f ||
	die("OPEN FILE ERROR: $opt_f file can not be opened (system error: $!)\n");

    my $F = `wc -l $opt_f`;
    chomp($F);

    my $f = 0;

    while(<$ffh>) {
	chomp($_);
	$f++;

	## Use a filter, do not parse the lines that begin with #

	print STDERR "\tProcessing the filter file: $opt_f, line:$f of $F                \r";

	unless ($_ =~ m/^#/) {

	    my @data = split(/\t/, $_);

	    ## Take the first column and forget the rest

	    $excl_ids{$data[0]} = 1;

	}
    }

    ## It also will create the output file

    my $outf_name = $output_name . 'filtered.frg';

    open $outf_fh, '>', $outf_name ||
	die("OPEN FILE ERROR: The output file:$outf_name (system error: $!)\n");
    
    print STDERR "\n";
}
else {
    print STDERR " [Disabled]\n\n";
}


## Tag option:

my %edit_tag_data;
my $outt_fh;

print STDERR "\n\nEDIT DATA BY TAGS:";

if ($opt_t) {

    print STDERR " [Enabled]\n\n";

    open my $tfh, '<', $opt_t ||
	die("OPEN FILE ERROR: $opt_t file can not be opened (system error: $!)\n");

    my $T = `wc -l $opt_t`;
    chomp($T);

    my $t = 0;

    while(<$tfh>) {
	chomp($_);
	$t++;

	## Get the first three columnar data and forget the rest

	print STDERR "\tProcessing the tag file: $opt_t, line:$t of $T                    \r";

	my ($id, $tag, $value) = split(/\t/, $_);

	## Store tag/value in a hash reference of a hash with key=id

	if (exists $edit_tag_data{$id}) {
	    
	    $edit_tag_data{$id}->{$tag} = $value;
	}
	else {
	    $edit_tag_data{$id} = { $tag => $value };
	}
	
    }
    print STDERR "\n";

    my $outt_name = $output_name . 'tag_edited.frg';

    open $outt_fh, '>', $outt_name ||
	die("OPEN FILE ERROR: The output file:$outt_name (system error: $!)\n");


}
else {
    print STDERR " [Disabled]\n\n";
}


## Seq and Qual options

my ($seqio, $qualio);

print STDERR "\n\nEXTRACT SEQUENCES AND QSCORES:";

if ($opt_S || $opt_Q) {

    if ($opt_S) {
	
	print STDERR "\n\tSequences [Enabled]";

	$seqio = Bio::SeqIO->new( -file   => ">$output_name".'.seq', 
	                          -format => 'fasta' );
    }
    
    if ($opt_Q) {

	print STDERR "\n\tQscores [Enabled]";

	$qualio = Bio::SeqIO->new( -file   => ">$output_name".'.qual', 
	                           -format => 'qual' );
    }  
}
else {
    print STDERR " [Disabled]\n\n";
}


print STDERR "\n\nPROCESSING FRG INPUT FILES:\n\n";

foreach my $infile (@infiles) {

    my $tag_data_href;

    ## Open the file

    open my $ifh, '<', $infile ||
	die("OPEN FILE ERROR: $infile file can not be openned (error system: $!).\n");

    my ($class);    
    my ($switch_comment, $switch_source, $switch_feature, $switch_seq, $switch_qlt, $switch_hps) = (0, 0, 0, 0, 0, 0);
    my ($comment, $source, $feature, $seq, $qlt, $hps);

    my $L = `wc -l $infile | cut -d ' ' -f1`;
    chomp($L);

    my $l;
    while(<$ifh>) {
	chomp($_);
	$l++;
	print STDERR "\tProcessing line:$l of $L lines for file:$infile                             \r";

	if ( $_ =~ m/^\{(\w+)/) {
	    $class = $1;	    
	}
	elsif ( $_ =~ m/^\}/) {
	   	    	   
	    ## Here it will print into the output files the data with the previous filtering
	    
	    if ($opt_f) {
	       	if ($class ne 'FRG') {
		    print_frg_entry($outf_fh, $class, $tag_data_href);
		}
		else {

		    if (exists $tag_data_href->{'acc'}) {
			unless (defined $excl_ids{$tag_data_href->{'acc'}}) {
			    
			    print_frg_entry($outf_fh, $class, $tag_data_href);
			}
		    }
		}
	    }

	    if ($opt_t) {
		if ($class ne 'FRG') {
		    print_frg_entry($outt_fh, $class, $tag_data_href);
		}
		else {

		    if (exists $tag_data_href->{'acc'}) {

			if (defined $edit_tag_data{$tag_data_href->{'acc'}}) {
			    my %edittags = %{$edit_tag_data{$tag_data_href->{'acc'}}}; 
			    
			    foreach my $tag (keys %edittags) {
				$tag_data_href->{$tag} = $edittags{$tag};
			    }
			    print_frg_entry($outt_fh, $class, $tag_data_href);
			}
		    }
		}
	    }

	    if ($class eq 'FRG') {
		    my $id = $tag_data_href->{'acc'};
		    my $seq = $tag_data_href->{'seq'};
		    $seq =~ s/\n//g;
		    my $qlt = $tag_data_href->{'qlt'};
		    $qlt =~ s/\n//g;
		    my $qual = frg_qlt_to_qual($qlt);

		    my $seqobj = Bio::Seq::Quality->new( -id => $id, -seq => $seq, -qual => $qual );

		    if ($opt_S) {
			$seqio->write_seq($seqobj);
		    }
		    if ($opt_Q) {
			$qualio->write_seq($seqobj);
		    }
	    }

	    
	    ## Also it will overwrite with undef values these variables

	    $class = '';

	    ## And delete all the entries for the array (better delete to spend less memory)

	    foreach my $k (keys %{$tag_data_href}) {
		delete $tag_data_href->{$k};
	    }
	}

	if ($class eq 'BAT') {
	    if ($_ =~ m/^bna:(.+)$/) {
		$tag_data_href->{'bna'} = $1;
	    }
	    elsif ($_ =~ m/^acc:(.+)$/) {
		$tag_data_href->{'acc'} = $1;
	    }
	    elsif ($_ =~ m/^com:/) {
		$switch_comment = 1;
	    }
	    elsif ($_ =~ m/^\.$/) {
		$switch_comment = 0;
		$tag_data_href->{'com'} = $comment;
		$comment = '';
	    }
	    
	    if ($switch_comment == 1) {
		unless ($_ =~ m/com:/) {
		    $comment .= $_ . "\n";
		}
	    }
	}
	elsif ($class eq 'VER') {
	    if ($_ =~ m/^ver:(.+)$/) {
		$tag_data_href->{'ver'} = $1;
	    }
	}
	elsif ($class eq 'LIB') {
	    if ($_ =~ m/^act:(.+)$/) {
		$tag_data_href->{'act'} = $1;
	    }
	    elsif ($_ =~ m/^acc:(.+)$/) {
		$tag_data_href->{'acc'} = $1;
	    }
	    elsif ($_ =~ m/^ori:(.+)$/) {
		$tag_data_href->{'ori'} = $1;
	    }
	    elsif ($_ =~ m/^mea:(.+)$/) {
		$tag_data_href->{'mea'} = $1;
	    }
	    elsif ($_ =~ m/^std:(.+)$/) {
		$tag_data_href->{'std'} = $1;
	    }
	    elsif ($_ =~ m/^src:(.+)$/) {
		$tag_data_href->{'src'} = $1;
	    }
	    elsif ($_ =~ m/^nft:(.+)$/) {
		$tag_data_href->{'nft'} = $1;
	    }
	    elsif ($_ =~ m/^src:/) {
		$switch_source = 1;
	    }
	    elsif ($_ =~ m/^fea:/) {
		$switch_feature = 1;
	    }
	    elsif ($_ =~ m/^\.$/) {
		if ($switch_source == 1) {
		    $switch_source = 0;
		    $tag_data_href->{'src'} = $source;
		    $source = '';
		}
		elsif ($switch_feature == 1) {
		    $switch_feature = 0;
		    $tag_data_href->{'fea'} = $feature;
		    $feature = '';
		}
	    }
	    
	    if ($switch_source == 1) {
		unless ($_ =~ m/src:/) {
		    $source .= $_ . "\n";
		}
	    }
	    elsif ($switch_feature == 1) {
		unless ($_ =~ m/fea:/) {
		    $feature .= $_ . "\n";
		}
	    }
	}
	elsif ($class eq 'FRG') {
	    if ($_ =~ m/^act:(.+)$/) {
		$tag_data_href->{'act'} = $1;
	    }
	    elsif ($_ =~ m/^acc:(.+)$/) {
		$tag_data_href->{'acc'} = $1;
	    }
	    elsif ($_ =~ m/^rnd:(.+)$/) {
		$tag_data_href->{'rnd'} = $1;
	    }
	    elsif ($_ =~ m/^sta:(.+)$/) {
		$tag_data_href->{'sta'} = $1;
	    }
	    elsif ($_ =~ m/^lib:(.+)$/) {
		$tag_data_href->{'lib'} = $1;
	    }
	    elsif ($_ =~ m/^pla:(.+)$/) {
		$tag_data_href->{'pla'} = $1;
	    }
	    elsif ($_ =~ m/^loc:(.+)$/) {
		$tag_data_href->{'loc'} = $1;
	    }
	    elsif ($_ =~ m/^clv:(.+)$/) {
		$tag_data_href->{'clv'} = $1;
	    }
	    elsif ($_ =~ m/^clq:(.+)$/) {
		$tag_data_href->{'clq'} = $1;
	    }
	    elsif ($_ =~ m/^clr:(.+)$/) {
		$tag_data_href->{'clr'} = $1;
	    }
	    elsif ($_ =~ m/^src:/) {
		$switch_source = 1;
	    }
	    elsif ($_ =~ m/^seq:/) {
		$switch_seq = 1;
	    }
	    elsif ($_ =~ m/^qlt:/) {
		$switch_qlt = 1;
	    }
	    elsif ($_ =~ m/^hps:/) {
		$switch_hps = 1;
	    }
	    elsif ($_ =~ m/^\.$/) {
		if ($switch_source == 1) {
		    $switch_source = 0;
		    $tag_data_href->{'src'} = $source;
		    $source = '';
		}
		elsif ($switch_seq == 1) {
		    $switch_seq = 0;
		    $tag_data_href->{'seq'} = $seq;
		    $seq = '';
		}
		elsif ($switch_qlt == 1) {
		    $switch_qlt = 0;
		    $tag_data_href->{'qlt'} = $qlt;
		    $qlt = '';
		}
		elsif ($switch_hps == 1) {
		    $switch_hps = 0;
		    $tag_data_href->{'hps'} = $hps;
		    $hps = '';
		}
	    }
	    
	    if ($switch_source == 1) {
		unless ($_ =~ m/src:/) {
		    $source .= $_ . "\n";
		}
	    }
	    elsif ($switch_seq == 1) {
		unless ($_ =~ m/seq:/) {
		    $seq .= $_ . "\n";
		}
	    }
	    elsif ($switch_qlt == 1) {
		unless ($_ =~ m/qlt:/) {
		    $qlt .= $_ . "\n";
		}
	    }
	    elsif ($switch_hps == 1) {
		unless ($_ =~ m/hps:/) {
		    $hps .= $_ . "\n";
		}
	    }
	}
	elsif ($class eq 'LKG') {
	    if ($_ =~ m/^act:(.+)$/) {
		$tag_data_href->{'act'} = $1;
	    }
	    elsif ($_ =~ m/^frg:(.+)$/) {
		if (exists $tag_data_href->{'frg'}) {
		    my @frg = @{$tag_data_href->{'frg'}};
		    push @frg, $1;
		    $tag_data_href->{'frg'} = \@frg;
		}
		else {
		    $tag_data_href->{'frg'} = [$1];
		}
	    }
	}
	elsif ($class eq 'PLC') {
	    if ($_ =~ m/^act:(.+)$/) {
		$tag_data_href->{'act'} = $1;
	    }
	    elsif ($_ =~ m/^frg:(.+)$/) {
		if (exists $tag_data_href->{'frg'}) {
		    my @frg = @{$tag_data_href->{'frg'}};
		    push @frg, $1;
		    $tag_data_href->{'frg'} = \@frg;
		}
		else {
		    $tag_data_href->{'frg'} = [$1];
		}
	    }
	}
    }
    print STDERR "\n\n";
    
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
  
     This script can be used to modify the information in frg files. This is the format used by wgs-assembler 
     (http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=FRG_Files). Check this page to get
     more informatio about frg format and the fields/tags used.

     To convert fasta, qual, sff or trace files to frg files use wgs scripts (convert-fasta-to-v2.pl, 
     tracedb-to-frg.pl, sffToCA, or fastqToCA) 

    Examples: 

      ## To remove a list of sequences

        frg_format_tool.pl -i test.frg -o outputname -f filter_test.txt

      ## To change tags in the frg file (example of tag_test.txt line: TEST_ID   sta    V)

        frg_format_tool.pl -i test.frg -o outputname -t tag_test.txt

      ## To get sequences

        frg_format_tool.pl -i test.frg -o outputname -S

      ## To get quality

        frg_format_tool.pl -i test.frg -o outputname -Q


    Usage:
     
     frg_format_tool.pl [-h] -i <input_files> -o <output_basename> [-f <filter_file>] [-t <tag_file>] [-S] [-Q]
  
    Flags:

     -i <input_files>            frg input files inside quotes and separated by commas (mandatory)
     -o <output_basename>        basename for the output files (mandatory)
     -f <filter_file>            filter file with a list of sequences id to remove from the input files
     -t <tag_file>               tag file in tab format (-f1 id, -f2 tag, -f3 value)
     -S <get_sequence>           get sequence from the frg files to create a fasta file
     -Q <get_qscores>            get qscores from the frg files to create a qual file
     -h <help>                   print the help

     

EOF
exit (1);
}


=head2 print_frg_entry

  Usage: print_frg_entry($fh, $class, $frg_data_href);

  Desc: Print into a file the data in frg format

  Ret: None

  Args: $fh, a file handle writable
        $class, frg format class
        $frg_data_href, a hash reference with
          key=tag_name
          value=data

  Side_Effects: die if is supplied none argument

  Example: print_frg_entry($fh, $class, $frg_data_href);

=cut


sub print_frg_entry {
    
    my $fh = shift ||
	die("ARGUMENT ERROR: None argument was supplied to print_frg_file function.\n");
    my $class = shift ||
	die("ARGUMENT ERROR: FRG class data was not supplied to print_frg_file function.\n");
    my $frg_data_href = shift ||
	die("ARGUMENT ERROR: FRG tag hash ref was not supplied to print_frg_file function.\n");

    ## We define the order for each data type

    my %order = (
	          'BAT' => ['bna', 'acc', 'com'],
	          'VER' => ['ver'],
	          'LIB' => ['act', 'acc', 'ori', 'mea', 'std', 'src', 'nft', 'fea'],
	          'FRG' => ['act', 'acc', 'rnd', 'sta', 'lib', 'pla', 'loc', 'src', 'seq', 'qlt', 'hps', 'clv', 'clq', 'clr'],
	          'LKG' => ['act', 'frg'],
	          'PLC' => ['act', 'frg']
	        );
    
    ## Also we defined the line and the multitags

    my %line_tags = (
	              'com' => 1,
	              'src' => 1, 
	              'fea' => 1, 
	              'seq' => 1, 
	              'qlt' => 1, 
	              'hps' => 1
	            );

    my %multi_tags = (
	               'frg' => 1
	             );


    print $fh "\n{$class\n";
	    
    my %tag_data = %{ $frg_data_href };

    my @tag_ordered = @{ $order{$class} };

    foreach my $tag (@tag_ordered) {
	my $data = $tag_data{$tag} || '';

	if (exists $line_tags{$tag}) {
	    print $fh "$tag:\n$data.\n";
	}
	elsif (exists $multi_tags{$tag}) {
	    my @multidata = @{$data};
	    foreach my $multidata (@multidata) {
		print $fh "$tag:$multidata\n";
	    }
	}
	else {
	    print $fh "$tag:$data\n";
	}
    }
    print $fh "}\n";
}

=head2 frg_qlt_to_qual

  Usage: my $qual = frg_qlt_to_qual($qlt);

  Desc: change a qlt format to qual

  Ret: $qual, a string

  Args: $qlt, a scalar

  Side_Effects: none

  Example: my $qual = frg_qlt_to_qual($qlt);

=cut


sub frg_qlt_to_qual {
    my $qlt = shift;

    my @qual;
    
    if (defined $qlt) {

	my @ascii = split(//, $qlt);
	
	foreach my $ascii (@ascii) {
	    my $q = ord($ascii)-48;
	    push @qual, $q;
	}
    }

    my $qual = join(' ', @qual);

    return $qual;    
}
