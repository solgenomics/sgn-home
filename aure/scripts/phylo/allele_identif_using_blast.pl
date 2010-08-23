#!/usr/bin/perl

=head1 NAME

 allele_identif_using_blast.pl
 Identification of allele based in blast result with parent species

=cut

=head1 SYPNOSIS

 allele_identif_using_blast.pl [-h] -b <blast_result_files> -s <dir_with_sequence_files> [-v <variable_for_selection>] [-o <output_filename>]

=head1 EXAMPLE 

 allele_identif_using_blast.pl -b 'Sp1=mysample.blastn.sp1.m8;Sp2=mysample.blastn.sp2.m8' -v score -o palleles.tab
                                  
=head2 I<Flags:>

=over


=item -b

B<blast_result_files>           Parent species and blast result files in m8 noted as SpecieName=BlastResultFilename separated by semicolon ';' (mandatory)

=item -s

B<dir_sequence_files>           Dir with the sequence files (mandatory)

=item -v

B<variable_for_selection>       Variable choose for allele asignment using blast result ('identity', 'evalue' or 'score') (score by default)

=item -o

B<output_filename>              Output filename (putative_alleles.tab by default)

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

    This script identify the best match between some blast results and assign that as putative allele.

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

use Bio::SeqIO;
use Bio::DB::Fasta;

our ($opt_b, $opt_v, $opt_o, $opt_s, $opt_h);
getopts("b:v:s:o:h");
if (!$opt_b && !$opt_v && !$opt_s && !$opt_o && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Mandatory arguments (die if they are not supplied)

my $blast_list = $opt_b 
    || die("\nMANDATORY ARGUMENT ERROR: -b <blast_result_list> argument was not supplied.\n");
my $seqdir = $opt_s 
    || die("\nMANDATORY ARGUMENT ERROR: -s <dir_sequence_list> argument was not supplied.\n");

## Default arguments (use default values if they are not supplied)

my $select_var = $opt_v 
    || 'score';

unless ($select_var =~ m/[identity|evalue|score]/i) {
    die("OPTIONAL ARGUMENT ERROR: -v <variable_for_selection> only can take the following values: identity, evalue or score.\n\n")
}

my $output_file = $opt_o 
    || 'putative_alleles_by_blast.tab';



##############################################
#### PARSING FILES ###########################
##############################################


## 1) Print start message and start to parse

my $date = `date`;
chomp($date);

print STDERR "\n\nSCRIPT START POINT ($date)\n\n";

## Divide the files and put into a hashes

my %blastfiles;
my @parents;

my @blastR = split(/;/, $blast_list);

foreach my $blast_r (@blastR) {
    
    my @blast_pair = split(/=/, $blast_r);

    $blastfiles{$blast_pair[0]} = $blast_pair[1];
    
    ## Also it will take the parents names
    
    push @parents, $blast_pair[0];
}

$date = `date`;
chomp($date);

print STDERR "STEP 0: CREATING BIOPERL SEQUENCE DB FOR SEQUENCE FILES ($date).\n\n";

my $seqdb = Bio::DB::Fasta->new($seqdir);

## Now it will parse each blast file taking always the first match

$date = `date`;
chomp($date);

print STDERR "STEP 1: PARSING BLAST FILES ($date).\n\n";

my %blast_parse = ();

foreach my $parent (sort keys %blastfiles) {

    my $blast_file = $blastfiles{$parent};

    print STDERR "\t\t+ Parsing blast file:$blast_file for parent $parent.\n\n";

    my $b_l = 0;
    my $B_l = `cut -f1 $blast_file | wc -l`;
    chomp($B_l);
    open my $b_fh, '<', $blast_file;
    
    while(<$b_fh>) {
	chomp($_);
	$b_l++;

	my @data = split(/\t/, $_);
	
	print STDERR "\t\tParsing line=$b_l of $B_l (obj_id=$data[0], sub_id=$data[1] and parent=$parent) for blast file=$blast_file.\r";
	
	## Take only the first hit (generally with best score)

	if (defined $blast_parse{$data[0]}) {
	    my %blast_parents = %{$blast_parse{$data[0]}};

	    unless (defined $blast_parents{$parent}) {
		$blast_parse{$data[0]}->{$parent} = { 'subject_id' => $data[1], 
						      'identity'   => $data[2], 
						      'alignment'  => $data[3], 
						      'mismatches' => $data[4], 
						      'gaps'       => $data[5], 
						      'qstart'     => $data[6], 
						      'qend'       => $data[7], 
						      'sstart'     => $data[8], 
						      'send'       => $data[9], 
						      'evalue'     => $data[10], 
						      'score'      => $data[11] 
		                                    };
	    }
	}
	else {
	    $blast_parse{$data[0]} = { $parent =>  { 'subject_id' => $data[1], 
						      'identity'   => $data[2], 
						      'alignment'  => $data[3], 
						      'mismatches' => $data[4], 
						      'gaps'       => $data[5], 
						      'qstart'     => $data[6], 
						      'qend'       => $data[7], 
						      'sstart'     => $data[8], 
						      'send'       => $data[9], 
						      'evalue'     => $data[10], 
						      'score'      => $data[11] 
		                                    }
	    };
	}
    }	
    print STDERR "\n\n";
    close($b_fh);
}

## Now it will compare different blast results to identify possible alleles

$date = `date`;
chomp($date);

print STDERR "STEP 2: CHOOSING BEST BLAST HIT ($date).\n\n";

my $selected_seqs = 'allele_selected_seqs';
mkdir($selected_seqs);

my %files = ( 'multifasta' => {} ) ;

open my $out_fh, '>', $output_file;

my $s_l = 0;
my $s_L = scalar(keys %blast_parse);

foreach my $seqid (sort keys %blast_parse) {

    print STDERR "\t\tProcessing seq_id=$seqid ($s_l of $s_L sequences).\r";

    my %parent_blast = %{$blast_parse{$seqid}};

    ## Now it will take all the parents from the file, if the blast match haven't all of them...
    ## it will write MissingBlastHit() in note

    my $filename = $selected_seqs . '/' . $seqid . '.fasta';
    $files{'multifasta'}->{$seqid} = $filename;
    my $fasta_io = Bio::SeqIO->new( -file => ">$filename", -format => 'fasta' );

    my $target_seqobj = $seqdb->get_Seq_by_id($seqid);
    
    if (defined $target_seqobj) {
	$fasta_io->write_seq($target_seqobj);
    }
    else {
	print STDERR "\n\tWARNING: Sequence_id=$seqid does not exist in the sequence dir.\n";
    }

    my %hits = ();
    my %par = ();

    foreach my $parent (@parents) {
    
	if (defined $parent_blast{$parent}) {
	    my %blast_data = %{$parent_blast{$parent}};

	    my $selected = $blast_data{$select_var};
	    $hits{$blast_data{'subject_id'}} = $selected;
	    $par{$blast_data{'subject_id'}} = $parent;
	  
	    my $parent_obj = $seqdb->get_Seq_by_id($blast_data{'subject_id'});
	    if (defined $parent_obj) {
		$fasta_io->write_seq($parent_obj);
	    }
	    else {
		print STDERR "\n\tWARNING: Sequence_id=$blast_data{'subject_id'} does not exist in the sequence dir.\n";
	    }
	}
	else {
	    print STDERR "\n\tWARNING: Sequence_id=$seqid has not blast hit for parent=$parent.\n";
	}	
    }

    ## Now it will choose the best hit

    my @best_hits = ();
    
    if (scalar(keys %hits) > 1) {
	if ($select_var =~ m/evalue/) {
	    @best_hits = sort {$hits{$b} <=> $hits{$a}} keys %hits;
	}
	else {
	    @best_hits = sort {$hits{$a} <=> $hits{$b}} keys %hits;
	}
	
	my $best = shift(@best_hits);
	my $better = shift(@best_hits);
	my @best_hits_list = ($best);
	
	if ($hits{$best} == $hits{$better}) {
	    push @best_hits_list, $better;
	    while ($hits{$best} == $hits{$better} && scalar(@best_hits) > 0) {
		$best = $better;
		$better = shift(@best_hits);
		push @best_hits_list, $better;
	    } 
	    my $best_list = join(',', @best_hits_list);
	    print $out_fh "$seqid\t$best_list\tBlast hit by $select_var ($best => $hits{$best}) with two or more parents\n";
	}
	else {
	    print $out_fh "$seqid\t$par{$best}\tBest blast hit by $select_var ($best => $hits{$best})\n";
	}
    }
    elsif (scalar(keys %hits) == 1) {
	my ($unique) = keys %hits;
	print $out_fh "$seqid\t$par{$unique}\tUnique blast hit ($unique => $hits{$unique})\n";
    }
    else {
	print $out_fh "$seqid\tNo hits\tNo hits\n";
    }
}

print STDERR "\n\n";

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
  
       This script identify the best match between some blast results and assign that as putative allele.
 
    Usage:
 
     allele_identif_using_blast.pl [-h] -b <blast_result_files> -s <sequence_files> [-v <variable_for_selection>] [-o <output_filename>]
                              
    Example:
      
     allele_identif_using_blast.pl -b 'Sp1=mysample.blastn.sp1.m8;Sp2=mysample.blastn.sp2.m8' -v score -o palleles.tab

    Flags:
     
     -b <blast_result_files>           Parent species and blast result files in m8 noted as SpecieName=BlastResultFilename separated by semicolon ';' (mandatory)
     -s <dir_sequence_files>           Dir with sequence files (mandatory)
     -v <variable_for_selection>       Variable choose for allele asignment using blast result ('identity', 'evalue' or 'score') (score by default)
     -o <output_filename>              Output filename (putative_alleles.tab by default)
     -h <help>                         print the help


EOF
exit (1);
}




####
1; #
####
