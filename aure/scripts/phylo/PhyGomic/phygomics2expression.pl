#!/usr/bin/perl

=head1 NAME

 phygomics2expression.pl
 Script to calculate RPKM expression for PhygOmics pipeline

=cut

=head1 SYPNOSIS

 phygomics2expression.pl [-h] -c <read_count_files> -s <strain_file> 
                              -m <map_file> -r <reference_sequence_file>
                              

=head1 EXAMPLE 

 phygomics2expression.pl -c a_count.tab,b_count.tab -s str_ab.tab
                         -m map_ab2ref,tab -r ref.fasta

=head2 I<Flags:>

=over


=item -c

B<read_count_files>             read count per library files with two columns:
                                -f1, contig_id, -f2, read count. (mandatory)

=item -s

B<strain_file>                  strain file with two columns:
                                -f1, contig_id, -f2, library/strain (mandatory)

=item -m

B<map_file>                     map file with at least two columns: 
                                -f1, contig_id, -f2, reference_id (mandatory)

=item -r

B<ref_sequence>                 reference sequence fasta file (mandatory)

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

    This script calculate the expression of reads in an assembly and compare
    them along different species/strains based in the map with a reference.
    
    It produces a file with:
    -f1: reference_id,
    -f2: contig_id,
    -f3: strain/specie/library
    { 
       -f4: n_reads
       -f5: RPKM  
    } x times as samples are.
    -f last: R value calculated according:

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

 phygomics2expression.pl


=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;
use Bio::SeqIO;

our ($opt_c, $opt_s, $opt_m, $opt_r, $opt_h);
getopts("c:s:m:r:h");
if (!$opt_c && !$opt_s && !$opt_m && !$opt_r && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my $date = `date`;
chomp($date);
print STDERR "\n************************************************************\n";
print STDERR "** phygomics2exp.pl starts ($date) **\n";
print STDERR "************************************************************\n";


## CHECK MANDATORY ARGUMENTS:

my $countfiles = $opt_c || 
    die("MANDATORY ARGUMENT ERROR: -c <read_count_files> was not supplied.\n");

my $strfile = $opt_s ||
    die("MANDATORY ARGUMENT ERROR: -s <strain_file> was not supplied.\n");

my $mapfile = $opt_m ||
    die("MANDATORY ARGUMENT ERROR: -m <map_file> was not supplied.\n");

my $reffile = $opt_r ||
    die("MANDATORY ARGUMENT ERROR: -r <referencefile> was not supplied.\n");


## PARSE FILES:

print STDERR "\n\t1) PARSING READ COUNT FILES:\n";

my @rcfiles = split(/,/, $countfiles);
my $rc_files = scalar(@rcfiles);
print STDERR "\n\t\t $rc_files have been detected.\n";

my %libs = ();

foreach my $rcfile (@rcfiles) {

    print STDERR "\n\t\tParsing file: $rcfile\n\n";

    open my $rcfh, '<', $rcfile;

    my $total = 0;
    my $line = 0;
    my $L = `cut -f1 $rcfile | wc -l`;
    chomp($L);
    
    my %readsc = ();

    while(<$rcfh>) {
    
	$line++;
	print STDERR "\t\tParsing line $line of $L                \r";

	chomp($_);
	my @data = split(/\t/, $_);

	$readsc{$data[0]} = $data[1];

	if ($data[1] !~ m/^\d+$/) {
	    die("-f2 read count is not a numeric data ($rcfile, line $line)\n");
	}
	$total += $data[1];
    }

    ## Add data to the hash
    $libs{$rcfile} = { total => $total, counts => \%readsc };
    print STDERR "\n\t\tDone.\n";
}


## Parse the reference file and store the length

print STDERR "\n\t2) PARSING REFERENCE FILE:\n";
my $refio = Bio::SeqIO->new( -file => $reffile, -format => 'fasta' );
print STDERR "\n\t\tDone\n";


my %refs = ();

while (my $refseq = $refio->next_seq()) {

    my $len = $refseq->length();
    my $id =  $refseq->id();

    print STDERR "\t\tProcessing ID:$id ($len)              \r";
    $refs{$id} = $len; 
}
print STDERR "\n\n\t\tDone.\n";


## Parse the map file

print STDERR "\n\t3) PARSING MAP FILE:\n";

open my $mapfh, '<', $mapfile;

my $map_l = 0;
my $map_L = `cut -f1 $mapfile | wc -l`;
chomp($map_L);
    
my %map = ();

while(<$mapfh>) {
    
    $map_l++;
    print STDERR "\t\tParsing line $map_l of $map_L              \r";
    
    chomp($_);
    my @data = split(/\t/, $_);
    
    $map{$data[0]} = $data[1];
}
print STDERR "\n\t\tDone\n";

## Finally parse the strain/library file

print STDERR "\n\t4) PARSING STRAIN FILE:\n";

open my $strfh, '<', $strfile;

my $str_l = 0;
my $str_L = `cut -f1 $strfile | wc -l`;
chomp($str_L);
    
my %str = ();

while(<$strfh>) {
    
    $str_l++;
    print STDERR "\t\tParsing line $str_l of $str_L              \r";
    
    chomp($_);
    my @data = split(/\t/, $_);
    
    $str{$data[0]} = $data[1];
}
print STDERR "\n\t\tDone\n";


## Define a na hash to store contigs were Sketel formula can not be applied

my %na = ();

## Now it will calculate the RPKM for each contig

my %exp = ();
my %strcount = ();

print STDERR "\n\t5) CALCULATING RPKM\n";

my @ord_contigs = sort { $map{$a} cmp $map{$b} } keys %map;

foreach my $contig_id (@ord_contigs) {

    my $length = $refs{$map{$contig_id}};
    my $strain = $str{$contig_id};
    
    print STDERR "\t\tProcessing contig_id $contig_id ($strain)($length)    \r";

    foreach my $src (sort keys %libs) {

	my $count = $libs{$src}->{counts}->{$contig_id};
	my $total = $libs{$src}->{total};

	## Calculate the 20%, it it is higher, add to the na hash

	if (defined $count && defined $total && defined $length) {
	   
	    my $perc = $count * 100 / $total;
	    if ($perc > 20 ) {
		if (exists $na{$contig_id}) {
		    push @{$na{$contig_id}}, $src;
		}
		else {
		    $na{$contig_id} = [$src];
		}
	    }
	    
	    my $rpkm = ( $count * 1000000 * 1000 ) / ( $total * $length ); 
	    my $rpkmobj = Math::BigFloat->new($rpkm);
	    my $form_rpkm = $rpkmobj->bfround(-2);

	    my %data = ( rpkm          => $form_rpkm, 
			 count         => $count, 
			 totallibrary  => $total ); 

	    if (exists $exp{$map{$contig_id}}) {
		if (exists  $exp{$map{$contig_id}}->{$strain}) {
		    $exp{$map{$contig_id}}->{$strain}->{$contig_id} = \%data;
		}
		else {
		    $exp{$map{$contig_id}}->{$strain} = {$contig_id => \%data};
		}
	    }
	    else {
		$exp{$map{$contig_id}} = { $strain => {$contig_id => \%data }}; 
	    }

	    if (exists $strcount{$strain}) {
		$strcount{$strain} += $count;
	    }
	    else {
		$strcount{$strain} = $count;
	    }
	}
    }
}
print STDERR "\n\n\t\tDone.\n";

## First check if the Stekel formula can be applied.

print STDERR "\n\t6) CHECKING THE LIBRARY PARS. TO APPLY THE Stekel Formula:\n";

## 1) The Library should have more than 50 reads.

foreach my $libfile (keys %libs) {

    my $total = $libs{$libfile}->{total};
    if ($total < 51) {
	print STDERR "\nThe file: $libfile has less than 50 total reads.";
	die("Please remove from the analysis.\n");
    }
}

## 2) No contig should represent more than 20% of the library
##    They have been stored at na hash

if (scalar(keys %na) > 0) {
    print STDERR "\n\nThere are some contigs were the library representation";
    print STDERR "is higher than 20%:\n\t";
    
    my $na_contigs = join('\n\t', keys %na);
    print STDERR "$na_contigs\n";
    die("Please, revise them before run the pipeline.\n\n");
}

print STDERR "\n\t\tCHEKING OK: All the libraries satisfy the conditions";
print STDERR " to apply Sketel Formula for abundance.\n";



## Finally it will calculate R and print the file

print STDERR "\n\n\t7) CALCULATING R DATA AND PRINTING OUTPUT FILE:\n\n";

my $output = "phygomics2expression.output.tab";
open my $outfh, '>', $output;

my $total_reads = 0;
foreach my $libfile (keys %libs) {
    $total_reads += $libs{$libfile}->{total};
}

my $N = Math::BigFloat->new($total_reads);

my @strains = sort keys %strcount;

foreach my $mapctg_id (sort keys %exp) {
    
    my %mapstr = %{$exp{$mapctg_id}};
    my @data = ();

    ## Calculate the total number of reads per cluster

    my $totalcls = 0;
    foreach my $str (@strains) {
	if (exists $mapstr{$str}) {
	    my %ctgdata = %{$mapstr{$str}};
	    foreach my $ctg_id (sort keys %ctgdata) {
		$totalcls += $ctgdata{$ctg_id}->{count};
	    }
	}
    }

    my $Xi = Math::BigFloat->new($totalcls);

    ## Now it will calculate R for each cluster

    my $r_value = 'NA';
    my $R;
    
    
    foreach my $str (@strains) {
	my $ctg_count = 0;
	if (exists $mapstr{$str}) {
	    my %ctgdata = %{$mapstr{$str}};
	    foreach my $ctg_id (sort keys %ctgdata) {
		$ctg_count += $ctgdata{$ctg_id}->{count};
	    }
	}
	my $Xij = Math::BigFloat->new($ctg_count);
	my $Ni = Math::BigFloat->new($strcount{$str});
	
	my $fi = $Xi->copy();
	$fi->bdiv($N);

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
    if (defined $R) {
	$r_value = $R->bfround(-2);
    }

    ## Now it will print the output

    foreach my $str (@strains) {
	if (exists $mapstr{$str}) {
	    my %ctgdata = %{$mapstr{$str}};
	    foreach my $ctg_id (sort keys %ctgdata) {
		my %spdata = %{$ctgdata{$ctg_id}};
		print $outfh "$mapctg_id\t$str\t$ctg_id\t$spdata{count}\t";
		print $outfh "$spdata{totallibrary}\t$spdata{rpkm}\t";
		print $outfh "$r_value\n";
	    }
	}
	else {
	    print $outfh "$mapctg_id\t$str\tNA\tNA\tNA\tNA\tNA\n";
	}
    }
}

print STDERR "\n\t\tDONE\n\n";


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
     
      This script calculate the expression of reads in an assembly and compare
    them along different species/strains based in the map with a reference.
    
    It produces a file with:
    -f1: reference_id,
    -f2: contig_id,
    -f3: strain/specie/library
    -f4: n_reads, 
    -f5: total reads
    -f7: RPKM
    -f6: R value calculated according:

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
 
       phygomics2expression.pl [-h] -c <read_count_files> -s <strain_file> 
                                    -m <map_file> -r <reference_sequence_file>

    Example:
    
        phygomics2expression.pl -c a_count.tab,b_count.tab -s str_ab.tab
                                -m map_ab2ref,tab -r ref.fasta


    Flags:
      
      -c <read_count_files>     read count per library files with two columns:
                                -f1, contig_id, -f2, read count. (mandatory)
      -s <strain_file>          strain file with two columns:
                                -f1, contig_id, -f2, library/strain (mandatory)
      -m <map_file>             map file with at least two columns: 
                                -f1, contig_id, -f2, reference_id (mandatory)
      -r <ref_sequence>         reference sequence fasta file (mandatory)
      -h <help>                 print the help

EOF
exit (1);
}

