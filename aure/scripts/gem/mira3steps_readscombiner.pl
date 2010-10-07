#!/usr/bin/perl

=head1 NAME

 mira3steps_readscombiner.pl
 Map the reads between different assemblies

=cut

=head1 SYPNOSIS

 mira3steps_readscombiner.pl [-h] -a assembly_file1 -b assembly_files2 -c assembly_file3 -l library_file -o output_name [-m allele_file] [-Q] 

=head1 EXAMPLE 

 mira3steps_readscombiner.pl [-h] -a step1_info_contigreadlist.txt 
                                  -b step2_Nsyl_info_contigreadlist.txt,step2_Ntom_info_contigreadlist.txt,step2_Ntab_info_contigreadlist.txt,step2_remain_info_contigreadlist.txt 
                                  -c step3_info_contigreadlist.txt 
                                  -l Nspecies_straindata_in.txt
                                  -o assembly_reads_mapping.txt
                                  -m Ntabacum_put_alleles.txt 
                                  

=head2 I<Flags:>

=over


=item -a

B<step1_contigreadslist>        step1_contigreadslist file, contig reads file produced by mira in the first step (mandatory)

=item -b

B<step2_contigreadslist>        step2_contigreadslist files, contig reads file produced by mira in the second step separated by commas (mandatory)

=item -c

B<step3_contigreadslist>        step3_contigsreadslist file, contig read file produced by mira in the third step (mandatory)

=item -l 

B<library_file>                 library file with reads in the first column and library in the second (mandatory)

=item -o

B<output_filename>              output filename, (output_mira3stepscombiner.txt by default)

=item -m

B<allele_file>                  file with contig_id (step2 with strain names) first column, putative allele in second and note in third (optional)

=item -Q

B<run_quiet>                    run the script in quiet mode, without produce any progress message (optional)

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

    This script is a simple mapping tool between differnt columnar files oriented for mira output.

    Basically combine all these files in a columnar file.

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 mira3steps_readscombiner.pl


=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Math::BigFloat;

our ($opt_a, $opt_b, $opt_c, $opt_l, $opt_o, $opt_m, $opt_Q, $opt_h);
getopts("a:b:c:l:o:m:Qh");
if (!$opt_a && !$opt_b && !$opt_c && !$opt_l && !$opt_o && !$opt_m && !$opt_Q && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Mandatory arguments (die if they are not supplied)

my $step1_file = $opt_a 
    || die("\nMANDATORY ARGUMENT ERROR: -a <step1_contigreadslist> argument was not supplied.\n");
my $step2_file = $opt_b 
    || die("\nMANDATORY ARGUMENT ERROR: -b <step2_contigreadslist> argument was not supplied.\n");
my $step3_file = $opt_c 
    || die("\nMANDATORY ARGUMENT ERROR: -c <step3_contigreadslist> argument was not supplied.\n");
my $library_file = $opt_l 
    || die("\nMANDATORY ARGUMENT ERROR: -l <library_file> argument was not supplied.\n");

## Default arguments (use default values if they are not supplied)

my $output_file = $opt_o 
    || 'output_mira3stepscombiner.txt';

## Optional arguments

my $allele_file = $opt_m;



##############################################
#### PARSING FILES ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);

unless ($opt_Q) { print STDERR "\n\nSCRIPT START POINT ($date)\n\nSTEP 1: PARSING FILES.\n\n"; }

## The parse order will be: library, step1, step2, step3 and allele.

## Parsing library file and storing data in %read_lib and %libs_count

my %read_lib = ();
my %libs_count = ();

$date = `date`;
chomp($date);
unless ($opt_Q) { print STDERR "\tParsing library_file=$library_file ($date)\n\n"; }

open my $lib_fh, '<', $library_file;   ## It is using autodie

my $lib_l = 0;

while(<$lib_fh>) {
    chomp($_);
    $lib_l++;

    unless ($opt_Q) { print STDERR "\t\tParsing line=$lib_l for library file=$library_file\r"; }

    ## Ignore lines that start with pounds

    unless ($_ =~ m/^#/) {
    
	my @data = split('\t', $_);

	## It will take always the first and second column and it will ignore the rest

	my $read = $data[0]
	    || die("\n\nLIBRARY FILE PARSING ERROR: line=$lib_l has not read data.\n\n");
	my $library = $data[1]
	    || die("\n\nLIBRARY FILE PARSING ERROR: line=$lib_l has not library data.\n\n");

	## store data in hashes

	$read_lib{$read} = $library;
	
	if (exists $libs_count{$library}) {
	    $libs_count{$library}++;
	}
	else {
	    $libs_count{$library} = 1;
	}
    }
}
unless ($opt_Q) { print STDERR "\n\n"; }


## Parsing first contigreadslist and store data into 
## it will store the data in both ways 
## %step1_readcontig where keys=read_id and value=step1_contig_id
## %step1_contigread where keys=step1_contig_id and value=array ref. with reads_id

my %step1_readcontig = ();
my %step1_contigread = ();

$date = `date`;
chomp($date);
unless ($opt_Q) { print STDERR "\tParsing step1_contigreadlist_file=$step1_file ($date)\n\n"; }

open my $s1_fh, '<', $step1_file;   ## It is using autodie

my $s1_l = 0;

while(<$s1_fh>) {
    chomp($_);
    $s1_l++;

    unless ($opt_Q) { print STDERR "\t\tParsing line=$s1_l for step1 file=$step1_file\r"; }

    ## Ignore lines that start with pounds

    unless ($_ =~ m/^#/) {

	my @data = split('\t', $_);

	## It will take always the first and second column and it will ignore the rest

	my $step1_contig_id = $data[0]
	    || die("\n\nSTEP1_CONTIGREADSLIST FILE PARSING ERROR: line=$s1_l has not contig_id data.\n\n");
	my $step1_read_id = $data[1]
	    || die("\n\nSTEP1_CONTIGREADSLIST FILE PARSING ERROR: line=$s1_l has not read_id data.\n\n");

	$step1_readcontig{$step1_read_id} = $step1_contig_id;

	if (exists $step1_contigread{$step1_contig_id}) {
	    push @{$step1_contigread{$step1_contig_id}}, $step1_read_id
	}
	else {
	    $step1_contigread{$step1_contig_id} = [$step1_read_id]
	}
    }
}
unless ($opt_Q) { print STDERR "\n\n"; }


## Parsing second contigreadslist and store data into 
## it will store the data in both ways 
## %step2_readcontig where keys=read_id and value=step2_contig_id
## %step2_contigread where keys=step2_contig_id and value=array ref. with reads_id

my %step2_readcontig = ();
my %step2_contigread = ();

$date = `date`;
chomp($date);
unless ($opt_Q) { print STDERR "\tParsing step2_contigreadlist_file=$step2_file ($date)\n\n"; }

my @step2_filelist = split(',', $step2_file);

## Check the number of files

unless (scalar(@step2_filelist) == scalar(keys %libs_count) + 1) {
    unless (scalar(@step2_filelist) == scalar(keys %libs_count)) { ## Ignoring remain strain
	warn("STEP2_CONTIGREADSLIST have a different number of files than libraries or libraries + remain library.\n\n")
    }
}

foreach my $step2_f (@step2_filelist) {
 
    open my $s2_fh, '<', $step2_f;   ## It is using autodie

    my $s2_l = 0;

    while(<$s2_fh>) {
	chomp($_);
	$s2_l++;

	unless ($opt_Q) { print STDERR "\t\tParsing line=$s2_l for step1 file=$step2_f\r"; }

	## Ignore lines that start with pounds

	unless ($_ =~ m/^#/) {

	    my @data = split('\t', $_);

	    ## It will take always the first and second column and it will ignore the rest

	    my $step2_contig_id = $data[0]
		|| die("\n\nSTEP1_CONTIGREADSLIST FILE PARSING ERROR: line=$s2_l for file $step2_f has not contig_id data.\n\n");
	    my $step2_read_id = $data[1]
		|| die("\n\nSTEP1_CONTIGREADSLIST FILE PARSING ERROR: line=$s2_l for file $step2_f has not read_id data.\n\n");

	    $step2_readcontig{$step2_read_id} = $step2_contig_id;

	    if (exists $step2_contigread{$step2_contig_id}) {
		push @{$step2_contigread{$step2_contig_id}}, $step2_read_id
	    }
	    else {
		$step2_contigread{$step2_contig_id} = [$step2_read_id]
	    }
	}
    }
    close($s2_fh);
    unless ($opt_Q) { print STDERR "\n\n"; }
}

## Parsing third contigreadslist and store data into 
## it will store the data in both ways 
## %step3_readcontig where keys=read_id and value=step3_contig_id
## %step3_contigread where keys=step3_contig_id and value=array ref. with reads_id

my %step3_readcontig = ();
my %step3_contigread = ();

$date = `date`;
chomp($date);
unless ($opt_Q) { print STDERR "\tParsing step3_contigreadlist_file=$step3_file ($date)\n\n"; }

open my $s3_fh, '<', $step3_file;   ## It is using autodie

my $s3_l = 0;

while(<$s3_fh>) {
    chomp($_);
    $s3_l++;

    unless ($opt_Q) { print STDERR "\t\tParsing line=$s3_l for step1 file=$step3_file\r"; }

    ## Ignore lines that start with pounds

    unless ($_ =~ m/^#/) {

	my @data = split('\t', $_);

	## It will take always the first and second column and it will ignore the rest

	my $step3_contig_id = $data[0]
	    || die("\n\nSTEP3_CONTIGREADSLIST FILE PARSING ERROR: line=$s3_l has not contig_id data.\n\n");
	my $step3_read_id = $data[1]
	    || die("\n\nSTEP3_CONTIGREADSLIST FILE PARSING ERROR: line=$s3_l has not read_id data.\n\n");

	$step3_readcontig{$step3_read_id} = $step3_contig_id;

	if (exists $step3_contigread{$step3_contig_id}) {
	    push @{$step3_contigread{$step3_contig_id}}, $step3_read_id
	}
	else {
	    $step1_contigread{$step3_contig_id} = [$step3_read_id]
	}
    }
}
unless ($opt_Q) { print STDERR "\n\n"; }

## Finally it will parse the allele file defining two variables
## %allele_names with keys=step2_contig_id and value=allele_name
## %allele_notes with keys=step2_contig_id and value=allele_note

my %allele_name = ();
my %allele_note = ();

if (defined $opt_m) {
    $date = `date`;
    chomp($date);
    unless ($opt_Q) { print STDERR "\tParsing allele_file=$allele_file ($date)\n\n"; }
    
    open my $a_fh, '<', $allele_file;   ## It is using autodie
    
    my $a_l = 0;

    while(<$a_fh>) {
	chomp($_);
	$a_l++;

	unless ($opt_Q) { print STDERR "\t\tParsing line=$a_l for allele file=$allele_file\r"; }

	## Ignore lines that start with pounds

	unless ($_ =~ m/^#/) {

	    my @data = split('\t', $_);

	    ## It will take always the first, second and third column and it will ignore the rest

	    my $step2_contig_id = $data[0]
		|| die("\n\nALLELE FILE PARSING ERROR: line=$s3_l has not contig_id data.\n\n");
	    my $allele_name = $data[1]
		|| die("\n\nALLELE FILE PARSING ERROR: line=$s3_l has not allele_name data.\n\n");
	    my $allele_note = $data[2];
	    
	    $allele_name{$step2_contig_id} = $allele_name;
	    
	    if (defined $allele_note) {
		$allele_note{$step2_contig_id} = $allele_note;
	    }
	}
    }
    unless ($opt_Q) { print STDERR "\n\n"; }
}

##############################################
#### COMBINING DATA ##########################
##############################################


## 2) Combine data using library as base

$date = `date`;
chomp($date);

unless ($opt_Q) { print STDERR "\n\nSTEP 2: COMBINING DATA ($date).\n\n"; }

## It will take all the reads_id from library file and get from each 
## hash the data, replacing by 'NA' (NoAvailable) when none data exists
## for that hash

my $reads_n = scalar(keys %read_lib);
my $read_c = 0;

open my $out_fh, '>', $output_file;

foreach my $readid (sort keys %read_lib) {
    $read_c++;

    unless ($opt_Q) { print STDERR "\t\tCombining data for read_id=$readid ($read_c of $reads_n)\r"; }

    my $libname = $read_lib{$readid};

    ## First get step1_contig and step2_contig

    my $step1_ctgid = $step1_readcontig{$readid} || 'NA';
    my $step2_ctgid = $step2_readcontig{$readid};

    ## If exists step2_contig, it will take the rest of the data

    my ($step3_ctgid, $allelename, $allelenote) = ('NA', 'NA', 'NA');

    if (defined $step2_ctgid) {
	$step2_ctgid =~ s/step2_//;
	$step3_ctgid = $step3_readcontig{$step2_ctgid} || 'NA';
	$allelename = $allele_name{$step2_ctgid} || 'NA';
	$allelenote = $allele_note{$step2_ctgid} || 'NA';	
    }
    else {
	$step2_ctgid = 'NA';
    }

    print $out_fh "$readid\t$libname\t$step1_ctgid\t$step2_ctgid\t$step3_ctgid\t$allelename\t$allelenote\n";
}
unless ($opt_Q) { print STDERR "\n\n"; }

$date = `date`;
chomp($date);
unless ($opt_Q) { print STDERR "\n\nEND OF THE SCRIPT ($date).\n\n"; }


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
  
      This script is a simple mapping tool between 
    different columnar files oriented for mira output.

    Basically combine all these files in a columnar file.
 
    Usage:

      mira3steps_readscombiner.pl [-h] -a assembly_file1 -b assembly_files2 -c assembly_file3 -l library_file -o output_name [-m allele_file] [-Q]

    Example:
       
       mira3steps_readscombiner.pl [-h] -a step1_info_contigreadlist.txt
                                         -b step2_Nsyl_info_contigreadlist.txt,step2_Ntom_info_contigreadlist.txt,step2_Ntab_info_contigreadlist.txt,step2_remain_info_contigreadlist.txt
                                         -c step3_info_contigreadlist.txt
                                         -l Nspecies_straindata_in.txt
                                         -o assembly_reads_mapping.txt
                                         -m Ntabacum_put_alleles.txt

    Flags:
       -a  step1_contigreadslist        step1_contigreadslist file, contig reads file produced by mira in the first step (mandatory)
       -b  step2_contigreadslist        step2_contigreadslist files, contig reads file produced by mira in the second step separated by commas (mandatory)
       -c  step3_contigreadslist        step3_contigsreadslist file, contig read file produced by mira in the third step (mandatory)
       -l  library_file                 library file with reads in the first column and library in the second (mandatory)
       -o  output_filename              output filename, (output_mira3stepscombiner.txt by default)
       -m  allele_file                  file with contig_id (step2 with strain names) first column, putative allele in second and note in third (optional)
       -Q  run_quiet                    run the script in quiet mode, without produce any progress message (optional)
       -h  help                         print the help

      

EOF
exit (1);
}




####
1; #
####
