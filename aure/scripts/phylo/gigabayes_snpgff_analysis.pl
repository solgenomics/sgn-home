#!/usr/bin/perl

=head1 NAME

 gigabayes_snpgff_analysis.pl
 Script to parse Gigabayes GFF files

=cut

=head1 SYPNOSIS

 gigabayes_snpgff_analysis.pl [-h] -g <gff_file> -s <sequence_file> [-p <ploidy>] [-f <parental_filter>] [-o <output>]

=head1 EXAMPLE 

 gigabayes_snpgff_analysis.pl -g SolSp.gff -s sol.fasta -p 2 -f 'L1=AA,L2=BB'
                                  
=head2 I<Flags:>

=over


=item -g

B<gff_file>                     GFF file produced by GigaBayes (mandatory)

=item -s

B<sequence_file>                Sequence file in fasta format (mandatory)

=item -p

B<ploidy>                       Ploidy used in Gigabayes (1 by default)

=item -f

B<strain_filter>                Strain filter detailed as Strain_Name=Allele1[Allele2] separated by commas

=item -o

B<output_filename>              Output filename (gigabayes_gff_analysis by default)

=item -h

B<help>                         Print the help

=back

=cut

=head1 DESCRIPTION

   Scripts with some functions to analyze the GFF file produced by GigaBayes.

   1) Number of SNP and INDEL per sequence.

   2) Average of SNP/INDEL per sequence and strain.

  A SNP filter can be used to ignore some SNPs. Here some examples:

   + Ignoring SNP when one of the strains is heterozygous -f Str1=AB

   + Ignoring SNP when two of the strains are heterozygous without relation
     between them. -f Str1=AB,Str2=CD

   + Ignoring SNP when two of the strains are heterozygous with relation
     between them. -f Str1=AB,Str2=AB

   + Ignoring SNP when two strains are different -f Str1=A,Str2=B

   + Ignoring SNP when two strains are the same -f Str1=A,Str2=A


=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 gigabayes_snpgff_analysis.pl


=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;
use Bio::SeqIO;
use Bio::Tools::SeqStats;

our ($opt_g, $opt_p, $opt_s, $opt_f, $opt_o, $opt_h);
getopts("g:p:s:f:o:h");
if (!$opt_g && !$opt_p && !$opt_s && !$opt_f && !$opt_o && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Mandatory arguments (die if they are not supplied)

my $gff = $opt_g 
    || die("\nMANDATORY ARGUMENT ERROR: -g <gff_file> argument was not supplied.\n");
my $fasta = $opt_s 
    || die("\nMANDATORY ARGUMENT ERROR: -g <sequence_file> argument was not supplied.\n");

## Default arguments (use default values if they are not supplied)

my $ploidy = $opt_p 
    || 1;

unless ($ploidy =~ m/\d+/i) {
    die("OPTIONAL ARGUMENT ERROR: -p <ploidy> only can take numeric values.\n\n")
}

my $output = $opt_o 
    || 'gigabayes_gff_analysis';


## Now it will parse the filter if it exists

my %str_filter = ();
my %grp_filter = ();

if (defined $opt_f) {
    my @ftr_elements = split(/,/, $opt_f);

    foreach my $element (@ftr_elements) {
	if ($element =~ m/(.+?)=(.+)/) {
	    $str_filter{$1} = $2;
	    if (exists $grp_filter{$2}) {
		push @{$grp_filter{$2}}, $1;
	    }
	    else {
		$grp_filter{$2} = [$1];
	    }
	}    
    }
}



##############################################
#### PARSING FILE  ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);

print STDERR "\n\nSCRIPT START POINT ($date)\n\n";

print STDERR "\t1) Parsing GFF file: $gff\n\n";

open my $ifh, '<', $gff;
my $l = 0;
my $L = `cut -f1 $gff | wc -l`;
chomp($L);

my %seq = ();
my %seqfilter = ();

my $disc_fh;
if ($opt_f) {
    my $discarted_file = $output . '.removed_by_filter.txt';
    open $disc_fh, '>', $discarted_file;
}

while (<$ifh>) {
    $l++;
    chomp($_);

    print STDERR "\t\tParsing line $l of $L lines in the file:$gff\r";

    ## Ignore the lines that start with #

    unless ($_ =~ m/^\s*#/) {

	my ($seqid, $source, $type, $start, $end, $score, $strand, $phase, $attline) = split(/\t/, $_); 
    
	## Parsing attributes ##

	my @attributes = split(/;/, $attline);

	## Now it is specific from Gigabayes format

	my %attr;

	foreach my $single_attr (@attributes) {
	    if ($single_attr =~ m/^(.+?)=(.+)$/) {
		$attr{$1} = $2;
	    }
	}

	my $data_href = { 
	                   'seqid'      => $seqid,
	                   'source'     => $source,
	                   'type'       => $type, 
	                   'start'      => $end, 
	                   'end'        => $end, 
	                   'score'      => $score, 
	                   'strand'     => $strand, 
	                   'phase'      => $phase,
			   'attr_line'  => $attline,
	                   'attributes' => \%attr 
	                };

	## Now it will filter if exists $opt_f

	my $add_switch = 1;
	if (defined $opt_f) {
	    
	    ## Get the individual genotypes to apply the filter

	    if (defined $attr{'individualGenotypes'}) {
		my @genotype_list = split(/,/,  $attr{'individualGenotypes'});
		my %genotypes;

		foreach my $gtype (@genotype_list) {
		    if ($ploidy == 2) {
			if ($gtype =~ m/^(.+?):(.)(.)\|(.+)$/) {
			    my $strain = $1;
			    my $allele1 = $2;
			    my $allele2 = $3;
			    my $a_score = $4;

			    $genotypes{$strain} = $allele1 . $allele2;

			    if (defined $str_filter{$strain}) {
				my @f_alleles = split(//, $str_filter{$strain}); 
				if ($f_alleles[0] eq $f_alleles[1]) {
				    if ($allele1 eq $allele2) {
					$add_switch = 0;
				    }
				}
				else {
				    if ($allele1 ne $allele2) {
					$add_switch = 0;
				    }
				}
			    }
			}
			elsif ($gtype =~ m/^(.+?):(.)\|(.+)$/) {
			    $genotypes{$1} = $2;
			}
		    }
		}

		foreach my $grp (keys %grp_filter) {
		    my @grp_strains = @{$grp_filter{$grp}};

		    my $rep_allele;
		    foreach my $gprstr (@grp_strains) {
			my $allele = $genotypes{$gprstr};
			
			unless (defined $rep_allele) {
			    $rep_allele = $allele;
			}
			else {
			    unless ($rep_allele ne $allele) {
				$add_switch = 0;
			    }
			}
		    }
		}
	    }
	}
       
	if ($add_switch == 1) {
	    if (exists $seq{$seqid}) {
		$seq{$seqid}->{$l} = $data_href;
	    }
	    else {
		$seq{$seqid} = { $l => $data_href};
	    }
	}
	else {
	    print $disc_fh "$seqid\t$attline\n";
	    if (exists $seqfilter{$seqid}) {
		$seqfilter{$seqid}->{$l} = $data_href;
	    }
	    else {
		$seqfilter{$seqid} = { $l => $data_href};
	    }
	}
    }
}
print STDERR "\n\n";

$date = `date`;
chomp($date);
print STDERR "\t2) Parsing sequence file: $fasta ($date)\n\n";

my %seqlength;
my %seqorder;
my %seqcomp;

my $seqio = Bio::SeqIO->new( -file => "<$fasta", -format => 'fasta');
my $s = 0;
my $S = `grep -c '>' $fasta`;
chomp($S);


while (my $seqobj = $seqio->next_seq()) {
    $s++;
    my $id = $seqobj->id();
    my $length = $seqobj->length();
    
    my $nt_stats_href = Bio::Tools::SeqStats->count_monomers($seqobj);

    print STDERR "\t\tExtracting length from sequence $id ($s of $S)\r";
    $seqlength{$id} = $length;
    $seqorder{$id} = $s;
    $seqcomp{$id} = $nt_stats_href;
}
print STDERR "\n\n";


##############################################
#### DATA ANALYSIS ###########################
##############################################


## First create the outputfiles.


my $sequence_file = $output . '.by_sequence.txt';
my $report_file = $output . '.report.txt';

open my $o1_fh, '>', $sequence_file;

## Print header

print $o1_fh "#ID#\t#Length#\t#%GC#\t#changesN#\t#changesNp100bp#\t#snpN#\t#snpNp100bp#\t#indelN#\t#indelNp100bp#\n";

$date = `date`;
chomp($date);
print STDERR "\t3) Data analysis:\n\n";

## Total variables

my ($T_bp, $T_changes, $T_snp, $T_indel) = (0, 0, 0, 0);
my %T_nt = ( 'A' => 0, 'T' => 0, 'C' => 0, 'G' => 0, 'N' => 0);
my %T_alleles = ( 'AT' => 0, 'AC' => 0, 'AG' => 0, 'CT' => 0, 'GT' => 0, 'CG' => 0, 'ACG' => 0, 'ACT' => 0, 'AGT' => 0, 'CGT' => 0);
my %T_indels = ( '*A' => 0, '*C' => 0, '*G' => 0, '*T' => 0, '*AT' => 0, '*AC' => 0, '*AG' => 0, '*CT' => 0, '*GT' => 0, '*CG' => 0);

my $t = 0;
my $T = scalar(keys %seqorder);

foreach my $seq_id (sort {$seqorder{$a} <=> $seqorder{$b}} keys %seqorder) {
    $t++;
    print STDERR "\t\tAnalyzing data for seq_id=$seq_id ($t of $T)\r";

    my $seqlength = $seqlength{$seq_id};
    $T_bp += $seqlength;

    my %seqnt = %{$seqcomp{$seq_id}};
    foreach my $nt (keys %seqnt) {
	$T_nt{$nt} += $seqnt{$nt};
    }
    
    my $gc_perc = ( $seqnt{'G'} + $seqnt{'C'} ) * 100 / $seqlength;
    my $gc_numobj = Math::BigFloat->new($gc_perc);
    my $r_gc_perc = $gc_numobj->bround(3);
    
    if (exists $seq{$seq_id}) {
	my %snps = %{$seq{$seq_id}};
	my $changes_n = scalar(keys %snps);
	$T_changes += $changes_n;

	my ($seq_snp, $seq_indel) = (0, 0);

	foreach my $ln (sort keys %snps) {
	    my %snpdata = %{$snps{$ln}};
	    my $alleles = $snpdata{'attributes'}->{'alleles'};

	    if ($snpdata{'type'} eq 'SNP') {
		$seq_snp++;
		$T_snp++;		
		my $ord_alleles = join('', sort(split(/,/, $alleles)));
		$T_alleles{$ord_alleles}++;				
	    }
	    elsif ($snpdata{'type'} eq 'INDEL') {
		$seq_indel++;
		$T_indel++;		
		my $ord_indels = join('', sort(split(/,/, $alleles)));
		$T_indels{$ord_indels}++;
	    } 
	}

	## Finally it will calculate snp/100bp
	my $changes_p100bp = $changes_n * 100 / $seqlength;
	my $ch_numobj = Math::BigFloat->new($changes_p100bp);
	my $rp_changes = $ch_numobj->bround(3);
	my $snps_p100bp = $seq_snp * 100 / $seqlength;
	my $snp_numobj = Math::BigFloat->new($snps_p100bp);
	my $rp_snps = $snp_numobj->bround(3);
	my $indels_p100bp = $seq_indel * 100 / $seqlength;
	my $ind_numobj = Math::BigFloat->new($indels_p100bp);
	my $rp_indels = $ind_numobj->bround(3);
	
	print $o1_fh "$seq_id\t$seqlength\t$r_gc_perc\t$changes_n\t$rp_changes\t$seq_snp\t$rp_snps\t$seq_indel\t$rp_indels\n";
    }
    else {
	print $o1_fh "$seq_id\t$seqlength\t$r_gc_perc\t0\t0\t0\t0\t0\t0\n";
    }
}
print STDERR "\n\n";

## Print report file

open my $o2_fh, '>', $report_file;

print $o2_fh "Sequences Analyzed:\t$T\n";
print $o2_fh "Total length:\t$T_bp\n";

## GC content

my $ntG = $T_nt{'G'};
my $ntC = $T_nt{'C'};
my $GC = ($ntG + $ntC) * 100 / $T_bp;
my $nobj_gc = Math::BigFloat->new($GC);
my $r_GC = $nobj_gc->bround(3);

print $o2_fh "GC content:\t$r_GC\n";

## Composition

foreach my $comp (sort {$T_nt{$b} <=> $T_nt{$a} } keys %T_nt) {
    print $o2_fh "$comp count:\t$T_nt{$comp}\n";
}
print $o2_fh "\n\n";

## SNPs

print $o2_fh "Total changes:\t$T_changes\n";
print $o2_fh "Total SNPs:\t$T_snp\n";
print $o2_fh "Total INDELs:\t$T_indel\n\n\n";

## Finally SNPs composition

foreach my $snp_type (sort {$T_alleles{$b} <=> $T_alleles{$a}} keys %T_alleles) {

    print $o2_fh "Allele_type $snp_type: $T_alleles{$snp_type}\n";
}
print $o2_fh "\n\n";

foreach my $ind_type (sort {$T_indels{$b} <=> $T_indels{$a}} keys %T_indels) {

    print $o2_fh "Indel_type $ind_type: $T_indels{$ind_type}\n";
}
print $o2_fh "\n\n";


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
 
     Scripts with some functions to analyze the GFF file produced by GigaBayes.

       1) Number of SNP and INDEL per sequence.

       2) Average of SNP/INDEL per sequence and strain.

     A SNP filter can be used to ignore some SNPs. Here some examples:

       + Ignoring SNP when one of the strains is heterozygous -f Str1=AB

       + Ignoring SNP when two of the strains are heterozygous without relation
         between them. -f Str1=AB,Str2=CD

       + Ignoring SNP when two of the strains are heterozygous with relation
         between them. -f Str1=AB,Str2=AB

       + Ignoring SNP when two strains are different -f Str1=A,Str2=B

       + Ignoring SNP when two strains are the same -f Str1=A,Str2=A

    Usage:
 
      gigabayes_snpgff_analysis.pl [-h] -g <gff_file> -s <sequence_file> [-p <ploidy>] [-f <parental_filter>] [-o <output>]
                              
    Example:

       gigabayes_snpgff_analysis.pl -g SolSp.gff -s sol.fasta -p 2 -f 'L1=AA,L2=BB'

    Flags:
     
       -g <gff_file>                GFF file produced by GigaBayes (mandatory)
       -s <sequence_file>           Sequence file in fasta format (mandatory)
       -p <ploidy>                  Ploidy used in Gigabayes (1 by default)
       -f <strain_filter>           Strain filter detailed as Strain_Name=Allele1[Allele2] separated by commas
       -o <output_filename>         Output filename (gigabayes_gff_analysis by default)
       -h <help>                    Print the help


EOF
exit (1);
}




####
1; #
####
