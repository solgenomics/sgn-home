#!/usr/bin/perl

=head1 NAME

 get_snp_from_mira_acefile.pl
 Script to parse and get snp data from ace assembly file from MIRA (version.0.1).

=cut

=head1 SYPNOSIS

 get_snp_from_mira_acefile.pl [-h] -i <ace_file> [-o <output_basename>] -s <strain_file> [-c <strain_compare>]

=head2 I<Flags:>

=over


=item -i

B<ace_file>                     assembly file in ace format format (mandatory)

=item -o

B<output_basename>              output basename (input_file.output by default)

=item -s

B<strain_file>                  file with two columns: -f1 read_id -f2 strain (mandatory) 

=item -c

B<strain_compare>               compare the strains as ref={strain_names},target={single_strain} (optional)

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

 This script parse the assembly file in ace format that mira can produce.
 It get the snp data based in the tags: 

  SAOr => SNPs intra-organism
  SROr => SNPs inter-organism
  SIOr => SNPs inter- and intra-organism 

 Get the type of the SNP and link with the strains.

 Finally should print:
    -f1 Contig_id
    -f2 SNP tag
    -f3 SNP position
    -f4..X SNP tag asspciated to the strain.
     (STRAIN:{type_a}xNa+{type_b}xNb...)
     Where type_a and type_b are different snps types
     for that possition and Na and Nb the count of them

 When the option -c is used, it will compare the reference strains with the target.

 For example, ref={strain1,strain2},target={strain3}, will compare each snp
 for each sequence, giving a value.

 contig_1   15   strain1={A}   strain2={C}   strain3={C}  ## strain3=strain2 ## 0 1 0
 contig_1   18   strain1={T}   strain2={T}   strain3={*}  ## !=              ## 0 0 1
 contig_1   38   strain1={T}   strain2={C}   strain3={C}  ## strain3=strain2 ## 0 1 0
 contig_1   75   strain1={*}   strain2={A}   strain3={A}  ## strain3=strain2 ## 0 1 0
 contig_1   81   strain1={C}   strain2={A}   strain3={A}  ## strain3=strain2 ## 0 1 0
 contig_1   90   strain1={C}   strain2={G}   strain3={C}  ## strain3=strain1 ## 1 0 0

 Finally it will give the values
   contig_1 strain3 strain1=rt{1/6} strain2=rt{4/6} diff=rt{1/6}
   contig_1 strain3 strain1=%{16.6} strain2=%{66.6} diff=%{16.6} 

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 get_snp_from_mira_acefile.pl

=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Math::BigFloat;

our ($opt_i, $opt_o, $opt_s, $opt_c, $opt_h);
getopts("i:o:s:c:h");
if (!$opt_i && !$opt_o && !$opt_s && !$opt_c && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check the mandatory ones

my $ace = $opt_i
    || die("MANDATORY ARGUMENT ERROR: -i <input_acefile> was not supplied.\n");
my $baseout = $opt_o 
    || $opt_i . '_output';
my $strain = $opt_s
    || die("MANDATORY ARGUMENT ERROR: -s <strain_file> was not supplied.\n");

## Now it will check if the opt_c argument is right

my @str_refs;
my $str_target;

if (defined $opt_c) {
    if ($opt_c =~ m/^ref=\{(.+?)\}\s*,\s*target=\{(.+?)\}\s*$/ ) {
	@str_refs = split(/,/, $1);
	$str_target = $2;  

	print STDERR "\nCalculate Strain SNP Ratio option: Enabled.\n";
	print STDERR "\tTarget=$str_target\n";
	my $str_references = join(',', @str_refs);
	print STDERR "\tReferences=$str_references";
    }
    else {
	die("ARGUMENT ERROR: -c <compare_strain> argument have not the structure: ref={strain1,strain2},target={strain3}.\n");
    }
}


## First, parse the file and get the data

print STDERR "\n\n1) PARSING ACE FILE.\n\n";

open my $ace_fh, '<', $ace 
    or die("OPEN FILE ERROR: file=$ace could not be openned (system error=$!).\n");

my $blast_results = {};
my $l = 0;

## To reduce the amount of memory used it will store only the data that it needs
## It will not use Bioperl to parse the ace because it parse and 
## everything (and for big files have memory problems)

## Define the interline variables

my $tag = '';
my $seq = '';
my $contig_id = '';
my $snp_id = '';
my $read_id = '';
my %snp;

my ($n_contigs, $n_reads) = (0, 0);

my (%co_seq, %re_seq);

while (<$ace_fh>) {
    $l++;

    print STDERR "\tParsing line=$l for file=$ace. ($n_contigs contigs and $n_reads reads have been parsed).      \r";

    chomp($_);

    ## Ace file works with tags. The tags used by this scripts are:
    ##
    ## CO contig_id

    ## Catch the contig_id

    if ($_ =~ m/^CO\s+(.+?)\s+/) { 
	$contig_id = $1;
	$tag = 'CO';
    }

    ## When the tag is CO it will get the sequence
    ## With an empty line, it will load the data into a hash
    ## and undef seq and tag variables

    if ($tag eq 'CO') {
	if ($_ =~ m/^$/) {
	    my $seqobj = Bio::Seq->new( -id => $contig_id, -seq => $seq );
	    $co_seq{$contig_id} = { seq => $seqobj };
	    
	    $seq = '';
	    $tag = '';
	    $n_contigs++;
	}
	elsif ($_ !~ m/^CO/) {
	    $seq .= $_;
	}
    }
    
    ## Get the alignments from AF 

    if ($_ =~ m/^AF\s+(\w+)\s+(\w)\s+(\d+)/) {
	
	$co_seq{$contig_id}->{reads}->{$1} = {c_align_start => $3, direction => $2};
    }    

    ## Next, catch the tags for the SNPs
    ##
    ## CT{
    ## contig_id tag coordinate
    ## COMMENT{
    ## confidence
    ## nucleotide}
    ## }
    ##

    if ($_ =~ m/^CT\{/) { 
	$tag = 'CT';
    }

    if ($tag eq 'CT') {
	if ($_ =~ m/^$/) {	    	    	    
	    $tag = '';
	    $snp_id = '';
	}
	elsif ($_ =~ /$contig_id\s+(S\w+)\s+MIRA\s+(\d+)\s+(\d+)/) {
	 
	    ## This only will take tags that start with S as SAOr, SIOr or SROr
	    
	    $snp_id = $2 . '-' . $1 . '-' .$contig_id;

	    if (exists $co_seq{$contig_id}->{'SNP'} ) {
		push @{$co_seq{$contig_id}->{'SNP'}}, $snp_id;
	    }
	    else {
		$co_seq{$contig_id}->{'SNP'} = [$snp_id]; 
	    }

	    $snp{$snp_id} = { type => $1, start => $2, end => $3 };
	}
	elsif ($_ =~ m/COMMENT{/) {
	    $tag = 'CT-COM';
	}
    }

    if ($tag eq 'CT-COM') {
	unless ($_ =~ m/^\}$/) {
	    $snp{$snp_id}->{comment} .= $_;
	}
	else {
	    $tag = '';
	    if ($snp{$snp_id}->{comment} =~ m/COMMENT\{(highC|mediumC|weakC)\}/) {
		$snp{$snp_id}->{comment} = $1;
	    }
	    my $c = $snp{$snp_id}->{comment};
	}
    }

    ## RD read_id
    ## sequence   

    if ($_ =~ m/^RD\s+(.+?)\s+/) { 
	$read_id = $1;
	$tag = 'RD';
	$n_reads++;
    }

    if ($tag eq 'RD') {
	if ($_ =~ m/^$/) {
	    my $seqobj = Bio::Seq->new( -id => $read_id, -seq => $seq );
	    
	    ## The read_id exists from the AF tag, so it only will add the sequence

	    $co_seq{$contig_id}->{reads}->{$read_id}->{seq} = $seqobj;

	    $seq = '';
	    $tag = '';
	}
	elsif ($_ !~ m/^RD/) {
	    $seq .= $_;
	}
    }

    ## QA read1 qual_start qual_end align_start align_end
    ##

    if ($_ =~ m/^QA\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)/) {
	
	$co_seq{$contig_id}->{reads}->{$read_id}->{r_align_start} = $3;
	$co_seq{$contig_id}->{reads}->{$read_id}->{r_align_end} = $4;

	$read_id = '';
    }    
}

## After parse the file with the strains

print STDERR "\n\n2) PARSING STRAIN FILE.\n\n";

open my $st_fh, '<', $strain
    or die("OPEN FILE ERROR: file=$strain could not be openned (system error=$!).\n");

my %strain = ();

my $ls =  0;
while (<$st_fh>) {
    
    $ls++;
    chomp($_);
    
    if ($_ =~ m/^\w+/) {
	my @data = split(/\t/, $_);

	print STDERR "\tParsing line=$ls ($data[0] => $data[1]).       \r";
	
	if (defined $data[0] && defined $data[1]) {
	    $strain{$data[0]} = $data[1];
	}
    }
}
print "\n";

## After parse the file it will get every SNP

print STDERR "\n\n3) GETTING SNPS.\n\n";

## Creation of different output_files.

my $snplist_file = $baseout . ".snplist.tab";
open my $out1_fh, '>', $snplist_file
    or die("OPEN FILE ERROR: file=$snplist_file could not be openned (system error=$!).\n");

my $snpstrain_file = $baseout . ".snp_strain.tab";
open my $out2_fh, '>', $snpstrain_file
    or die("OPEN FILE ERROR: file=$snpstrain_file could not be openned (system error=$!).\n");

my $snpratio_file = $baseout . ".snp_ratio_analysis.tab";
my $out3_fh;

if ($opt_c) {
    open $out3_fh, '>', $snpratio_file
	or die("OPEN FILE ERROR: file=$snpratio_file could not be openned (system error=$!).\n");
}

print STDERR "Output_files:\n\n\t* $snplist_file (snp list by read)\n\t* $snpstrain_file (snp composition per contig position)\n";
if ($opt_c) {
    print STDERR "\t* $snpratio_file (snp ration per contig)\n";
}

my $p = 0;
print STDERR "\n";

my %co_strain_ratio = ();

foreach my $contigid (sort keys %co_seq) {

    my %contigdata = %{$co_seq{$contigid}};

    ## It will create the count hash for the -c option

    my %strain_ratio = ( diff => 0, total => 0, subtotal => 0 );
    foreach my $refs (@str_refs) {
	$strain_ratio{$refs} = 0;
    }
    
    if (defined $contigdata{'SNP'}) {
    
	my @snps = @{$contigdata{'SNP'}};

	foreach my $snp_id (@snps) {


	    my @check = split(/-/, $snp_id);

	    if ($check[2] ne $contigid) {
		print STDERR "WARNING: contig_id=$contigid is not contained in snp_id=$snp_id.\n";
	    }

	    print STDERR "\tProcessing contig_id=$contigid, snp=$snp_id     \r";

	    if (defined $snp{$snp_id}) {

		my $snptype = $snp{$snp_id}->{type};
		my $position = $snp{$snp_id}->{start};
		my $comment = $snp{$snp_id}->{comment} || '';
		
		my $snp = $contigdata{'seq'}->subseq($position, $position); 

		my %reads = %{$contigdata{'reads'}};
		my $c_reads = scalar(keys %reads);

		## define the strain count as two dimension STRAIN+CHANGE

		my %str_snp = ();
		
		foreach my $read_id (keys %reads) {
		    my $contig_st = $reads{$read_id}->{c_align_start};
		    my $align_dir = $reads{$read_id}->{direction};
		    my $read_st = $reads{$read_id}->{r_align_start};
		    my $read_en = $reads{$read_id}->{r_align_end};

		    my $read_strain = $strain{$read_id} || 'without_strain';

		    if ( $contig_st + $read_st - 1 <= $position && $position <= $contig_st + $read_en - 1) {
			my $rel_position = $position - $contig_st + $read_st;
			my $snp_read = $reads{$read_id}->{seq}->subseq($rel_position, $rel_position);
			print $out1_fh "$contigid\t$snptype\t$position\t$comment\t$snp\t$read_id\t$snp_read\t$align_dir\n";
		    
			my $str_snp_code = $read_strain . ':{' . $snp_read . '}';
			
			if (exists $str_snp{$str_snp_code}) {
			    $str_snp{$str_snp_code}++;
			}
			else {
			    $str_snp{$str_snp_code} = 1;
			}
		    }
		}

		print $out2_fh "$contigid\t$snptype\t$position\t$comment\t$snp";
	
		my %fsort_str = ();
		foreach my $str_code (sort keys %str_snp) {
		
		    my ($str, $snpc) = split(/:/, $str_code);

		    if (exists $fsort_str{$str}) {
			$fsort_str{$str} .= '+'. $snpc . 'x' . $str_snp{$str_code};
		    }
		    else {
			$fsort_str{$str} = $snpc . 'x' . $str_snp{$str_code};
		    }
		}
		foreach my $fstr_code (sort keys %fsort_str) {
		    print $out2_fh "\t$fstr_code" . '=' . "$fsort_str{$fstr_code}";
		}

		## Now if option_c is enabled, it will calculate the ratios of the reference
		
		if ($opt_c) {
		
		    my $snp_sim_true = 0;
		    my $target_snp = $fsort_str{$str_target};
		    		    
		    ## Sometimes the snps is not defined for the target and the references

		    if (defined $target_snp) {
			my @tsnps = split(/\+/, $target_snp);
			
			## check that the snps exists for all the references and the targets

			my $check_refs = 1;
			foreach my $check_str_ref (@str_refs) {
				    
			    unless (defined $fsort_str{$check_str_ref}) { 
				$check_refs = 0;
			    }
			}
			
			if ($check_refs == 1) { ## Only if the snps is defined for all the strains it will calculate the perc.

			    foreach my $tsnp (@tsnps) {
			    
				if ($tsnp =~ m/\{(.)\}x\d+/) {
				    my $ts = $1; 

				    foreach my $str_ref (@str_refs) {
					
					if (defined $fsort_str{$str_ref}) {
					    my $ref_snp = $fsort_str{$str_ref};
					    my @rsnps = split(/\+/, $ref_snp);
				   
					    foreach my $rsnp (@rsnps) {
					    
						if ($rsnp =~ m/\{(.)\}x\d+/) {
						    my $rs = $1;
						    if ($rs eq $ts) {
							$strain_ratio{$str_ref}++;
							$strain_ratio{'total'}++;
							$strain_ratio{'subtotal'}++;
							$snp_sim_true = 1;
						    }
						}
					    }
					}
				    }
				}
			    }			
			
			    if ($snp_sim_true == 0) {
				$strain_ratio{'subtotal'}++;
			    }
			}
		    }
		    if ($snp_sim_true == 0) {
			$strain_ratio{'diff'}++;
			$strain_ratio{'total'}++;
		    }
		}
  					
		print $out2_fh "\n";
	    }
	    else {
		print STDERR "\nWARNING ($contig_id): snp_id=$snp_id was not parsed.\n";
	    }
	}
    
    
	## Finally it will add the hash to the contig
	## And print the percentages
	
	if ($opt_c ) {
    
	    ## Print Headers

	    if ($p == 0) {
		print $out3_fh "#contig_id\t#Tot_C\t#SubtC\t#DiffC";
		foreach my $strref (@str_refs) {
		    print $out3_fh "\t#$strref" . "C";
		}
		print $out3_fh "\t#DiffP";
		foreach my $strref (@str_refs) {		
		    print $out3_fh "\t#$strref" . "P";
		}
		$p++;
		print $out3_fh "\n";
	    }
	    

	    print $out3_fh "$contigid\t$strain_ratio{'total'}\t$strain_ratio{'subtotal'}\t$strain_ratio{'diff'}";
	
	    foreach my $strref (@str_refs) {
		print $out3_fh "\t$strain_ratio{$strref}";
	    }	
	
	    my $diff_perc = ($strain_ratio{'diff'}/$strain_ratio{'total'})*100;
	    my $diff_perc_obj = Math::BigFloat->new($diff_perc);
	    my $f_diff_perc = $diff_perc_obj->bfround(-2);

	    print $out3_fh "\t$f_diff_perc%";
	
	    foreach my $strref (@str_refs) {
		my $str_perc = ($strain_ratio{$strref}/$strain_ratio{'total'})*100;
		my $str_perc_obj = Math::BigFloat->new($str_perc);
		my $f_str_perc = $str_perc_obj->bfround(-2);
		print $out3_fh "\t$f_str_perc%";
	    }
	    print $out3_fh "\n";
	}
    }    
}

print STDERR "\n\n";



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

     This script parse the assembly file in ace format that mira can produce.
     It get the snp data based in the tags: 

      SAOr => SNPs intra-organism
      SROr => SNPs inter-organism
      SIOr => SNPs inter- and intra-organism 

     Get the type of the SNP and link with the strains.

     Finally it should print:
      -f1 Contig_id
      -f2 SNP tag
      -f3 SNP position
      -f4..X SNP tag asspciated to the strain.
             (STRAIN:{type_a}xNa+{type_b}xNb...)
             Where type_a and type_b are different snps types
             for that possition and Na and Nb the count of them    
      When the option -c is used, it will compare the reference strains with the target.

      For example, ref={strain1,strain2},target={strain3}, will compare each snp
      for each sequence, giving a value.

      contig_1   15   strain1={A}   strain2={C}   strain3={C}  ## strain3=strain2 ## 0 1 0
      contig_1   18   strain1={T}   strain2={T}   strain3={*}  ## !=              ## 0 0 1
      contig_1   38   strain1={T}   strain2={C}   strain3={C}  ## strain3=strain2 ## 0 1 0
      contig_1   75   strain1={*}   strain2={A}   strain3={A}  ## strain3=strain2 ## 0 1 0
      contig_1   81   strain1={C}   strain2={A}   strain3={A}  ## strain3=strain2 ## 0 1 0
      contig_1   90   strain1={C}   strain2={G}   strain3={C}  ## strain3=strain1 ## 1 0 0

      Finally it will give the values
      contig_1 strain3 strain1=rt{1/6} strain2=rt{4/6} diff=rt{1/6}
      contig_1 strain3 strain1=%{16.6} strain2=%{66.6} diff=%{16.6} 
	           
    Usage:
  
     get_snp_from_mira_acefile.pl [-h] -i <ace_file> -s <strain_file> [-o <output_basename>]
       
    Flags:
      
     -i <ace_file>              assembly file in ace format format (mandatory)
     -o <output_basename>       output basename (input_file.output by default)
     -s <strain_file>           file with two columns (-f1 read_id, -f2 strain) (mandatory)
     -h <help>                  print the help
      

EOF
exit (1);
}
