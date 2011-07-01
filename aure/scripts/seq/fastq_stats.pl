#!/usr/bin/perl

=head1 NAME

 fastq_stats.pl
 A script to get some stats from a fastq file.

=cut

=head1 SYPNOSIS

 fastq_stats.pl [-h] -p fastq_file1 [-q fastq_file2] - t <seq_technology> [-V]
                [-c <cutoff_value>]               

=head2 I<Flags:>

=over


=item -p

B<fastq_file1>            fastq file (mandatory)

=item -q

B<fastq_file2>            if a pair end file is analyzed, second fastq file.

=item -t 

B<seq_technology>         sequencing technology used (see descrip., mandatory)

=item -c

B<tagging_cutoff_value>   q15 index value for tag a seq as poor quality (90)

=item -V

B<verbose>                print running messages

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script parse the fastq file and return the following files:

 1) input_filename.stats.txt, a file with 5 columns:
     - ID, 
     - raw_length, 
     - Number_of_Ns, 
     - Q15_trim_length, 
     - Q20_trim_length, 
     - Q25_trim_length,
     - Q15_index (% of bases with Qscore >= 15).
     - Q20_index (% of bases with Qscore >= 20).
     - Q20_index (% of bases with Qscore >= 25).

 2) If pair ends are supplied, a pairends_cmp.txt file with the double columns
    for the same values than the first file.
    
 Note: This scripts supports the following sequencing technologies for 
       qscore calculation.

       S or sanger or 454, phred+33  (range 0,40)
       X or solexa,        solexa+64 (range -5,40)
       I or illumina1.3,   phred+64, (range 0,40)
       J or illumina1.5,   phred+64, (range 3,40, 0 and 1 unused, 2 qual.ctrl)


=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 fastq_stats.pl


=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::SeqIO;
use Math::BigFloat;

our ($opt_p, $opt_q, $opt_t, $opt_c, $opt_V, $opt_h);
getopts("p:q:t:c:Vh");
if (!$opt_p && !$opt_q && !$opt_t && !$opt_c && !$opt_V && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check them

my $fastq1 = $opt_p || 
    die("INPUT ARGUMENTS ERROR: -p <fastq_file1> has not been supplied.\n");

my $fastq2 = $opt_q;

my $seqtech = $opt_t ||
    die("INPUT ARGUMENTS ERROR: -t <seq_technology has not been supplied.\n");

my $q15_cutoff = $opt_c || 90;

my %permtech = (
		'S'           => 'sanger',
		'454'         => 'sanger',
		'X'           => 'solexa',
		'I'           => 'illumina',
		'illumina'    => 'illumina',
		'illumina1.3' => 'illumina',
		'J'           => 'illumina',
		'illumina1.5' => 'illumina',
		);

unless (exists $permtech{$seqtech}) {
    die("-t $seqtech is not a permited technology. See -help.\n");
}

## Define the order for the printed columns

my @col = ('raw_length', 'n_count', 'q15_length', 'q20_length', 
	   'q25_length', 'q15_idx', 'q20_idx', 'q25_idx');


## Parse the file and get the stats

my $dt = `date`;
chomp($dt);
print STDERR "\n\n===> Script Starts [$dt].<===\n";

$dt = `date`;
chomp($dt);
print STDERR "\n1) Parsing first fastq file: $fastq1 [$dt]...\n";

my ($seqstats1_href, $seqidx1_href) = seqstats($fastq1, $permtech{$seqtech});

$dt = `date`;
chomp($dt);
print STDERR "\t\tDone [$dt].\n";


## Print the first output

my $outfile1 = $fastq1 . ".stats.txt";
open my $out_fh1, '>', $outfile1;

$dt = `date`;
chomp($dt);
print STDERR "\n2) Printing first status file: $outfile1 [$dt].\n";


foreach my $idx1 (sort keys %{$seqidx1_href}) {

    my $id1 = $seqidx1_href->{$idx1};

    print $out_fh1 "$id1/1\t";
    
    my @data = ();
    foreach my $c (@col) {
	push @data, $seqstats1_href->{$id1}->{$c};
    }
    
    my $dataline = join("\t", @data);
    print $out_fh1 "$dataline\n";
}

$dt = `date`;
chomp($dt);
print STDERR "\t\tDone [$dt].\n";

## Repeat the same action if exists the second file

if (defined $fastq2) {

    $dt = `date`;
    chomp($dt);
    print STDERR "\n===> Pair file was supplied. <===\n";
    print STDERR "\n3) Parsing second fastq file: $fastq2 [$dt]...\n";

    my ($seqstats2href, $seqidx2href) = seqstats($fastq2, $permtech{$seqtech});

    $dt = `date`;
    chomp($dt);
    print STDERR "\t\tDone [$dt].\n";

    ## Print the second output

    my $outfile2 = $fastq2 . ".stats.txt";
    open my $out_fh2, '>', $outfile2;

    $dt = `date`;
    chomp($dt);
    print STDERR "\n4) Printing second status file: $outfile2 [$dt].\n";

    foreach my $idx2 (sort keys %{$seqidx2href}) {

	my $id2 = $seqidx2href->{$idx2};

	print $out_fh2 "$id2/2\t";
    
	my @data2 = ();
	foreach my $c (@col) {
	    push @data2, $seqstats2href->{$id2}->{$c};
	}
    
	my $dataline2 = join("\t", @data2);
	print $out_fh2 "$dataline2\n";
    }

    $dt = `date`;
    chomp($dt);
    print STDERR "\t\tDone [$dt].\n";

    ## Now it will compare both files and print the comparison

    my $outfile3 = "pairends_cmp.txt";
    open my $out_fh3, '>', $outfile3;

    ## Print headers

    print $out_fh3 "##        SEQUENCE_PAIR_ID        ##\t#RawL#\t#NCnt#\t";
    print $out_fh3 "#Q15L#\t#Q20L#\t#Q25L#\t";
    print $out_fh3 "#(Q<15)%#\t#(Q<20)%#\t#(Q<25)%#\t";
    print $out_fh3 "#QUALITY_TAG#\n";

    $dt = `date`;
    chomp($dt);
    print STDERR "\n5) Printing comparison file: $outfile3 [$dt].\n";

    foreach my $idx1 (sort keys %{$seqidx1_href}) {

	my $id1 = $seqidx1_href->{$idx1};
	
	my @data3 = ($id1);
	if (exists $seqstats2href->{$id1}) {
	    
	    my ($bad1, $bad2) = (0, 0);

	    foreach my $c (@col) {
		my $coldata1 = $seqstats1_href->{$id1}->{$c};
		my $coldata2 = $seqstats2href->{$id1}->{$c};
		push @data3, $coldata1 . "|" . $coldata2;

		if ($c eq 'q15_idx') {
		    if ($coldata1 >= $q15_cutoff) {
			$bad1 = 1;
		    }
		    if ($coldata2 >= $q15_cutoff) {
			$bad2 = 1;
		    }
		}
	    }
	    my $line = join("\t", @data3);
	    print $out_fh3 "$line\t";

	    if ($bad1 == 1 && $bad2 == 0) {
		print $out_fh3 "PAIR1 POOR QUALITY\n";
	    }
	    elsif ($bad1 == 0 && $bad2 == 1) {
		print $out_fh3 "PAIR2 POOR QUALITY\n";
	    }
	    elsif ($bad1 == 1 && $bad2 == 1) {
		print $out_fh3 "BOTH PAIRS POOR QUALITY\n";
	    }
	    else {
		print $out_fh3 "GOOD QUALITY\n";
	    }
	}
	else {
	    foreach my $c (@col) {
		my $coldata1 = $seqstats1_href->{$id1}->{$c};
		push @data3, $coldata1 . "|NA";
	    }
	    my $line = join("\t", @data3);
	    print $out_fh3 "$line\tABSENT PAIR2\n";
	}
    }

    foreach my $idx2 (sort keys %{$seqidx2href}) {

	my $id2 = $seqidx2href->{$idx2};
	
	my @data4 = ($id2);
	
	unless (exists $seqstats1_href->{$id2}) {
	    
	    foreach my $c (@col) {
		
		my $coldata2 = $seqstats2href->{$id2}->{$c};
		push @data4, "NA|" . $coldata2;
	    }

	    my $line = join("\t", @data4);
	    print $out_fh3 "$line\tABSENT PAIR1\n";	    
	}
    }

    $dt = `date`;
    chomp($dt);
    print STDERR "\t\tDone [$dt].\n";
}

$dt = `date`;
chomp($dt);
print STDERR "\n\n===> Script Ends [$dt].<===\n\n";



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
      
       This script parse the fastq file and return the following files:

         1) input_filename.stats.txt, a file with 5 columns:
            - ID, 
            - raw_length, 
            - Number_of_Ns, 
            - Q15_trim_length, 
            - Q20_trim_length, 
            - Q25_trim_length,
            - Q15_index (% of bases with Qscore >= 15).
            - Q20_index (% of bases with Qscore >= 20).
            - Q20_index (% of bases with Qscore >= 25).

         2) If pair ends are supplied, a pairends_cmp.txt file with the double 
            columns for the same values than the first file.
    
    Note: 
  
       This scripts supports the following sequencing technologies for qscore 
       calculation.

         S or sanger or 454, phred+33  (range 0,40)
         X or solexa,        solexa+64 (range -5,40)
         I or illumina1.3,   phred+64, (range 0,40)
         J or illumina1.5,   phred+64, (range 3,40, 0 and 1 unused,2 qual.ctrl)

    Usage:
     
      fastq_stats.pl [-h] -p fastq_file1 [-q fastq_file2] - t <seq_technology> 
                     [-V] [-c <cutoff_value>]     

    Flags:

      -p <fastq_file1>            first fastq file (mandatory)
      -q <fastq_file2>            second fastq file (pair end)(optional) 
      -t <seq_technology>         sequence technology used (mandatory)
      -c <tagging_cutoff_value>   q15 index value for tag a seq as poor quality
                                  (90 by default)
      -V <verbose>                print parsing status messages 
      -h <help>                   print the help
     

EOF
exit (1);
}

=head2 seqstats

  Usage: my ($seqstats_href, $seqidx_href) = seqstats($filename, $filevariant);
  Desc: Parse the file and calculate the stats
  Ret: A hashref with key=ID and value=hashref. with key=(raw_length, n_count, 
       q15_length, q20_length, q25_length, q15_idx, q20_idx, q25_idx)and 
       value=value
       Another hashref. with key=sequence index and value=sequence id
  Args: $filename, filename
        $filevariant, s scalar (sanger, solexa or illumina)
  Side_Effects: none
  Example: my ($seqstats, $seqidx) = seqstats($filename, $filevariant);

=cut

sub seqstats {

    my $seqfile = shift;
    my $variant = shift;

    ## Define variable

    my %seq = ();
    my %seqidx = ();  ## To preserve the sequence order in the file

    ## If verbose is enabled it will count how many sequence are and
    ## how many sequences are processing.

    my $tseq;
    if ($opt_V) {
	
	$tseq = `grep -c '^\@' $seqfile`;
	chomp($tseq);
    }

    ## Open the file and parse

    my $seqio1 = Bio::SeqIO->new( -format  => 'fastq',
				  -variant => $variant,
				  -file    => $seqfile,
				  );


    my $idx = 0;
    while (my $seqobj = $seqio1->next_seq()) {
	
	$idx++;

	if ($opt_V) {
	    print STDERR "\t\tParsing sequence $idx of $tseq     \r";
	}

	my $id = $seqobj->id();

	## It needs to remove the last two digits (pair end index) from the id

	$id =~ s/..$//;

	$seqidx{$idx} = $id;

	## Basic stats

	my $seq = $seqobj->seq();    
	my $length = length($seq);
	my $seq_woutn = $seq;
	$seq_woutn =~ s/N//gi;
	my $n_count = $length - length($seq_woutn);
	
	## Qscores

	my @qual = @{$seqobj->qual()};
    
	my ($q15, $q20, $q25) = (0, 0, 0);                          ## Counts
	my ($i_q15, $i_q20, $i_q25) = (1, 1, 1);                    ## Initial
	my ($e_q15, $e_q20, $e_q25) = ($length, $length, $length);  ## End

	my $p = 0;
    
	foreach my $q (@qual) {
	    $p++;                     ## First value is 1

	    if ($q < 25) {
		$q25++;
		if ($i_q25 == $p) {
		    $i_q25++;
		}
	    
		if ($q < 20) {
		    $q20++;
		    if ($i_q20 == $p) {
			$i_q20++;
		    }
	    
		    if ($q < 15) { 
			$q15++;
			if ($i_q15 == $p) {	    
			    $i_q15++;
			}
		    }
		}
	    }
	}
	
	my $r = $length;
	foreach my $qr (reverse @qual) {
	    
	    if ($qr < 25) {
		if ($e_q25 == $r) {
		    $e_q25--;
		}
	    
		if ($qr < 20) {
		    if ($e_q20 == $r) {
			$e_q20--;
		    }
		
		    if ($qr < 15) {
			if ($e_q15 == $r) {
			    $e_q15--;
			}
		    }
		}
	    }
	    $r--;  ## Reduce 1 the reverse
	}

	## Store the data into a hash

	my $q15_length = $e_q15 - $i_q15 + 1;
	if ($i_q15 > $e_q15) {
	    $q15_length = 0;
	}
	my $q20_length = $e_q20 - $i_q20 + 1;
	if ($i_q20 > $e_q20) {
	    $q20_length = 0;
	}
	my $q25_length = $e_q25 - $i_q25 + 1;
	if ($i_q25 > $e_q25) {
	    $q25_length = 0;
	}
	
	my $q15_perc = Math::BigFloat->new(( $q15 / $length ) * 100)
	                             ->bfround(-2);
	my $q20_perc = Math::BigFloat->new(( $q20 / $length ) * 100)
	                             ->bfround(-2);
	my $q25_perc = Math::BigFloat->new(( $q25 / $length ) * 100)
	                             ->bfround(-2);


	my %stats = (
		     raw_length => $length,
		     n_count    => $n_count,
		     q15_length => $q15_length,
		     q20_length => $q20_length,
		     q25_length => $q25_length,
		     q15_idx    => $q15_perc,
		     q20_idx    => $q20_perc,
		     q25_idx    => $q25_perc,
		     );
    
	$seq{$id} = \%stats;
    }
    
    if ($opt_V) {
	print STDERR "\n";
    }

    return (\%seq, \%seqidx);
}
