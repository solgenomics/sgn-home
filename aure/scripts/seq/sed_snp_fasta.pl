#!/usr/bin/perl

=head1 NAME

 sed_snp_fasta.pl
 
 Regenerate a new fasta file replacing snps in a reference 
 sequence (version.0.1).

=cut

=head1 SYPNOSIS

 sed_snp_fasta.pl [-h] -s <snp_file> -f <fasta_input_file> 
                       [-r <region_extract_file>] -o <output_file> [-S]

=head2 I<Flags:>

=over


=item -f

B<fasta_file>           fasta file (mandatory)

=item -s

B<snp_file>             snp file with 4 columns (-f1 seq_id,
                        -f2 position, -f3 ref_nt, -f4 snp )

=item -r

B<region_extract_file>  coordinates file to extract sequence
                        regions from the fasta file

=item -o

B<output_file>          output filename (outfile by default)

=item -S

B<check_snp_ref>        check snp reference. If it is not 
                        the same skip change

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script replace the nucleotide in the specified region 
 with the nucleotide from the snp file (id, position, 
 ref_nt and snp_nt)

 Also have the option of get a sequence between a start and
 end position detailed in a tab file (id, start, end)

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 sed_snp_fasta.pl

=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::SeqIO;
use Math::BigFloat;


our ($opt_f, $opt_s, $opt_r, $opt_o, $opt_S, $opt_h);
getopts("f:s:r:o:Sh");
if (!$opt_f && !$opt_s && !$opt_r && !$opt_o && !$opt_S && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}


## Check arguments

my $fastafile = $opt_f 
    || die("DATA ARGUMENT ERROR: -f fasta file WAS NOT SUPPLIED.\n");

my $snpfile = $opt_s;
my $regionfile = $opt_r;

my $outfile = $opt_o
    || 'outfile';

## Check the goal of the script

if (!$snpfile && !$regionfile) {
    die("DATA ARGUMENT ERROR: None -s <snp_file> or -r <region_file> were supplied.\n");
}


## Parse SNPs file

my %snps;
 
if (defined $snpfile) {

    ## Use of autodie
    open(my $snp_fh, '<', $snpfile); 

    print STDERR "\nPARSING snp file:$snpfile\n\n";
    my $l = 0;

    while(<$snp_fh>) {
	chomp($_);
	$l++;

	print STDERR "\tParsing line $l from snp_file.\r";

	## Ignore comments
	## First column => seq_id
	## Second column => snp_position
	## Third column => reference sequence nt
	## Forth column => snp nt

	unless ($_  =~ m/^#/) {
	    my @data = split(/\t/, $_);
	
	    my $seq_id = $data[0];
	    
	    if (scalar(@data) >= 4) {
		if (exists $snps{$seq_id}) {
		    $snps{$seq_id}->{$data[1]} = [ $data[2], $data[3] ];
		}
		else {
		    $snps{$seq_id} = { $data[1] => [ $data[2], $data[3] ] };
		}
	    }
	}
    }
    print STDERR "\n\tDone.\n\n";
}

## Parse region file

my %regions;
 
if (defined $regionfile) {
    
    ## Use of autodie
    open(my $reg_fh, '<', $regionfile); 

    print STDERR "\nPARSING region file:$regionfile\n\n";
    my $n = 0;

    while(<$reg_fh>) {
	chomp($_);
	$n++;

	print STDERR "\tParsing line $n from region_file.\r";

	## Ignore comments
	## First column => seq_id
	## Second column => start position
	## Third column => end position

	unless ($_  =~ m/^#/) {
	    my @data = split(/\t/, $_);
	
	    my $seq_id = $data[0];

	    if (exists $regions{$seq_id}) {
		push @{$regions{$seq_id}}, [ $data[1], $data[2]];
	    }
	    else {
		$regions{$seq_id} = [ [$data[1], $data[2] ] ];
	    }
	}	
    }
    print STDERR "\n\tDone.\n\n";
}

## Now it will open the fasta file and it will do the nt changes

my $seqio = Bio::SeqIO->new( -file => $fastafile, -format => 'fasta' );
my $outseqio = Bio::SeqIO->new( -file => ">$outfile", -format => 'fasta' );

my %changes = (); ## To count the changes per sequence
my %extract = (); ## To count the extractions per sequence

print STDERR "SNPs changes and sequences extractions.\n\n";

## If opt_S is used it will store the SNPs differences between ref. and
## the fasta file in a file

my $diff_fh;
if (defined $opt_S) {
    my $snp_file_diff = $outfile . '.snp_dif.tab';
    open ($diff_fh, '>', $snp_file_diff);
}

while( my $seqobj = $seqio->next_seq() ) {
    
    ## It will print in the end of the script all the sequences from this parsing
    my @seqprint;

    my $id = $seqobj->id();
    my $seq = $seqobj->seq();
    my $length = $seqobj->length();
    
    if (defined $snps{$id}) {

	my @nts = split(//, $seq);
	my %snp_seq = %{$snps{$id}};
	
	foreach my $snp_p (sort {$a <=> $b} keys %snp_seq) {
	
	    if ($snp_p <= $length) {
		my $ref = $snp_seq{$snp_p}->[0];
		my $new = $snp_seq{$snp_p}->[1];
		my $seqnt = $nts[$snp_p-1];
	    
		if ($seqnt ne $ref && $opt_S) { ## Means that it is the nt to change
		    print STDERR "\nWARNING: The nt ($seqnt) for the position $snp_p in the sequence $id is not $ref. Skipping change.\n";
		    print $diff_fh "$id\t$snp_p\t$seqnt\t!=\t$ref\n";
		}
		else {
		    $nts[$snp_p-1] = $new;
		    print STDERR "\tChanging $seqnt for $new in the sequence $id ($snp_p)\r";
		    $changes{$id}++;		    
		}
	    }
	    else {
		print STDERR "\nWARNING: The position $snp_p is out the sequence (length=$length). Skipping position.\n";
	    }
	}
	$seqobj->seq(join('', @nts));
    }

    if (defined $regions{$id}) {
	
	my @regions_aref = @{$regions{$id}};
	
	foreach my $reg_aref (@regions_aref) {
	    my $st = $reg_aref->[0];
	    my $en = $reg_aref->[1];
	    print STDERR "\tExtracting region [$st-$en] from the sequence $id\r";
	    my $rev = 0;
	    if ($st > $en) {
		$rev = 1;
		$st = $reg_aref->[1];
		$en = $reg_aref->[0];
	    }

	    my $new_seqobj = $seqobj->trunc($st,$en);
	    if ($rev == 1) {
		$new_seqobj = $new_seqobj->revcom();
	    }

	    ## Also it will need rewrite the id with the region
	    $new_seqobj->id($id . '_[' . $st . '-' . $en. ']' );
	    $extract{$id}++;
	    push @seqprint, $new_seqobj;
	}  
    }
    

    ## After the changes it will create the new seq and store into the new file    
    
    if (defined $regionfile) {
	foreach my $seqobj_reg (@seqprint) {
	    $outseqio->write_seq($seqobj_reg);
	}
    }
    else {
	$outseqio->write_seq($seqobj);
    }
}

## Finally it will print the results

print STDERR "\n\n** Script Results:\n\n";
my ($total_ch, $total_ex) = (0, 0);

foreach my $ch (keys %changes) {
    $total_ch += $changes{$ch};
    print STDOUT "seq_id=$ch\tchanges=$changes{$ch}\n";
}
print STDOUT "Total changes in the fasta file: $total_ch\n\n";

foreach my $ex (keys %extract) {
    $total_ex += $extract{$ex};
    print STDOUT "seq_id=$ex\tregion extractions=$extract{$ex}\n";
}
print STDOUT "Total region extractions in the fasta file: $total_ex\n\n";


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

        This script replace the nucleotide in the specified region 
      with the nucleotide from the snp file (id, position, 
      ref_nt and snp_nt)

      Also have the option of get a sequence between a start and
      end position detailed in a tab file (id, start, end)
    
    Usage:
        sed_snp_fasta.pl [-h] -s <snp_file> -f <fasta_input_file> 
                              [-r <region_extract_file>] -o <output_file> [-S]
      
    Example:
	sed_snp_fasta.pl -s tomato_cv_ta496.snp.tab -f tomato_cv_m86.fasta 
      
    Flags:

	-f <fasta_file>           fasta file (mandatory)
        -s <snp_file>             snp file with 4 columns (-f1 seq_id,
                                  -f2 position, -f3 ref_nt, -f4 snp )
	-r <region_extract_file>  coordinates file to extract sequence
                                  regions from the fasta file
        -o <output_file>          output filename (outfile by default)
	-S <check_snp_ref>        check snp reference. If it is not 
                                  the same skip change
	-h <help>                 print the help

     
EOF
exit (1);
}




###
1;#
###
