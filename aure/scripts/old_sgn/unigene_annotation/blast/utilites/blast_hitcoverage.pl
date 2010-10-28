#!/usr/bin/perl

=head1 NAME

 blast_hitcoverage.pl
 Script to filter sequences based in the coverage and identity of a match (version.0.1).

=cut

=head1 SYPNOSIS

 blast_hitcoverage.pl [-h] -r <blast_result_file> -s <sequence_file> 
                      [-c <coverage>] [-i <percentage_identity>] [-g <gaps>] [-H] > output.tab

=head2 I<Flags:>

=over


=item -r

B<blast_result_file>            blast result file in m8 format (mandatory)

=item -s

B<sequence_file>                sequence file in fasta format (mandatory)

=item -c

B<coverage>                     percentage of the coverage (75 by default)

=item -i

B<identity>                     percentage of the identity (90 by default)

=item -g

B<gaps_number>                  number of gaps permited in the coverage (10 by default) 

=item -H

B<print header>                 print header as stdout

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

 This script parse the blast result file, get the length of the sequence file and calculate the coverage for each
 sequence.

 Filter the sequences and the blast file based in a coverage and identity cutoff value.

 It will produce as STDOUT: $query_id\t$subject_id\t$coverage\t$identity

 Note: This script uses a lot of memory. For big blast results is better if you reduce the size of the file by
       filterind using other scripts as blast_handle_tool.pl

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 blast_result_handle_tool.pl


=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use Bio::SeqIO;
use Math::BigFloat;

our ($opt_r, $opt_s, $opt_c, $opt_i, $opt_g, $opt_H, $opt_h);
getopts("r:s:c:i:g:Hh");
if (!$opt_r && !$opt_s && !$opt_c && !$opt_i && !$opt_g && !$opt_H && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check the mandatory ones

my $blast = $opt_r 
    || die("MANDATORY ARGUMENT ERROR: -r <blast_result_file> was not supplied.\n");
my $seqfile = $opt_s 
    || die("MANDATORY ARGUMENT ERROR: -s <sequences_file> was not supplied.\n");
my $cutoff_cov = $opt_c 
    || 75;
my $cutoff_ident = $opt_i 
    || 90;
my $cutoff_gaps = $opt_g
    || 10;

## Parse the blast file clustering by query_id and subject_id with the
## following structure:
##    $blast_result = { query_id => { subject_id => [{ blast_field => value }] } }
##

print STDERR "\n\n1) PARSING BLAST FILE.\n\n";

open my $bl_fh, '<', $blast 
    || die("OPEN FILE ERROR: file=$blast could not be openned (system error=$!).\n");

my $blast_results = {};
my $l = 0;

while (<$bl_fh>) {
    $l++;
    print STDERR "\tParsing line=$l for blast file=$blast.                  \r";

    chomp($_);
    my @data = split(/\t/, $_);

    ## To reduce the amount of memory used it will store only the data that it needs

    my $blast_data = { 
		       identity     => $data[2],
		       align_length => $data[3],
		       q_start      => $data[6],
		       q_end        => $data[7],
                     };

    if (exists $blast_results->{$data[0]}) {
	my $subj_results = $blast_results->{$data[0]};
	
	if (exists $subj_results->{$data[1]}) {
	    push @{$subj_results->{$data[1]}}, $blast_data;
	}
	else {
	    $subj_results->{$data[1]} = [$blast_data];
	}
    }
    else {
	$blast_results->{$data[0]} = { $data[1] => [$blast_data] };
    }
}


## Now it will get the length for each sequence from the sequence file
## and store in:
##   $seq_length = { $query_id => $length };

print STDERR "\n\n2) EXTRACTING SEQUENCE LENGTHS.\n\n";

my $seq_length = {};
my $s = 0;

my $seqio = Bio::SeqIO->new( -file => "$seqfile", -format => 'fasta' ); 

while( my $seqobj = $seqio->next_seq() ) {

    my $id = $seqobj->primary_id();
    my $length = $seqobj->length();

    $seq_length->{$id} = $length;

    $s++;
    print STDERR "\tExtracting length ($length) for id=$id (seq $s).                  \r";
}


## It will calculate the percentage of the coverage calculating the global alignment length
## for each query_id/subject_id pair.

print STDERR "\n\n3) CALCULATING COVERAGE AND FILTERING RESULTS.\n\n";

my $t = -1;

foreach my $query_id (sort keys %{$blast_results}) {
    
    my $subj_results = $blast_results->{$query_id};

    my $query_length = $seq_length->{$query_id} ||
	print STDERR "WARNING: sequence_id=$query_id have not sequence length from $seqfile. Skipping sequence.\n\n";

    if (defined $query_length) { 
	
	foreach my $subject_id (sort keys %{$subj_results}) {
	
	    my $global_align = 0;
	    my $global_identity = 0;

	    my %match_regions = ();

	    foreach my $hit_results ( @{$subj_results->{$subject_id}} ) {
		
		my $new_match_region;
		 
		my $qst = $hit_results->{'q_start'};
		my $qen = $hit_results->{'q_end'};
		my $qid = $hit_results->{'identity'};
		my $qname = $query_id . '{m}' . $subject_id . '{region}[' . $qst . '-' . $qen . ']';

		print STDERR "\tCalculating overlappings for $qname                \r";

		## Switch the coordinate for - strands

		if ($qst > $qen) {
		    $qst = $hit_results->{'q_end'};
		    $qen = $hit_results->{'q_start'};
		}

		## Compare with previous match regions and overwrite the overlapping

		my $old = 0;

		if (scalar(keys %match_regions) > 0) {
		    foreach my $mname (keys %match_regions) {

			my $match_region = $match_regions{$mname};

			my $mst = $match_region->{'start'};
			my $men = $match_region->{'end'};
			my $mid = $match_region->{'identity'};

			## If exists the overlapping it will replace the old match region for the new one

			if ($men+1 >= $qst && $mst <= $qen+1) { ## Exists overlapping

			    if ($mst <= $qst && $men >= $qen) {
				## Do nothing, just ignore the match region is in the previous one
				$old += 1;
			    }
			    elsif ($mst < $qst  && $men < $qen) { ## overlap in the end
				$match_region->{'end'} = $qen;
				$qst = $mst;

				if ($mid < $qid) {
				    $qid = $mid;
				}

				delete $match_regions{$qname};
				delete $match_regions{$mname};
				$qname = $query_id . '{m}' . $subject_id . '{region}[' . $qst . '-' . $qen . ']';
				$match_regions{$qname} = { start => $qst, end => $qen, identity => $qid };

				$old += 1;
			    } 
			    elsif ($mst > $qst  && $men > $qen) { ## overlap in the start
				$match_region->{'start'} = $qst;
				$qen = $men;
				
				if ($mid < $qid) {
				    $qid = $mid;
				}

				delete $match_regions{$qname};
				delete $match_regions{$mname};
				$qname = $query_id . '{m}' . $subject_id . '{region}[' . $qst . '-' . $qen . ']';
				$match_regions{$qname} = { start => $qst, end => $qen, identity => $qid };

				$old += 1;
			    }
			    elsif ($mst > $qst && $men < $qen) { ## overlap both of them
				$match_region->{'end'} = $qen;
				$match_region->{'start'} = $qst;

				if ($mid < $qid) {
				    $qid = $mid;
				}
				
				delete $match_regions{$qname};
				delete $match_regions{$mname};
				$qname = $query_id . '{m}' . $subject_id . '{region}[' . $qst . '-' . $qen . ']';
				$match_regions{$qname} = { start => $qst, end => $qen, identity => $qid };

				$old += 1;
			    }			
			}
		    }
		}
		
		## If has been defined a new match region, it will added to the array
		
		if ($old == 0) {
		    $match_regions{$qname} = { start => $qst, end => $qen, identity => $qid };
		}
	    }

	    ## After parse all the hits, it will sum all the match_regions

	    foreach my $m_name (keys %match_regions) {
		my $align_l = $match_regions{$m_name}->{'end'} - $match_regions{$m_name}->{'start'} + 1;
		
		$global_align += $align_l;
		$global_identity += $match_regions{$m_name}->{'identity'};
	    }

	    ## Coverage will be:

	    my $cov = ($global_align / $query_length) * 100;
	    my $covobj = Math::BigFloat->new($cov);

	    ## Global identity will be:

	    my $ident = $global_identity / scalar(keys %match_regions);
	    my $identobj = Math::BigFloat->new($ident);

	    ## Gaps will be:

	    my $gaps = scalar(keys %match_regions) - 1;

	    if ($opt_H && $t == -1) {
		$t++;
		print STDOUT "##query_id##\t##subject_id##\t##query_length##\t##form_cov##\t##form_ident##\t##gaps##\n";
	    }
	    if ($cov > $cutoff_cov && $ident > $cutoff_ident && $gaps <= $cutoff_gaps) {
		$t++;
		my $form_cov = $covobj->bfround(-2);
		my $form_ident = $identobj->bfround(-2);
		print STDOUT "$query_id\t$subject_id\t$query_length\t$form_cov\t$form_ident\t$gaps\n";
	    }	    
	}
    }
}


print STDERR "\n\n\tDone... $t matches has been printed as STDOUT.\n\n";


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
    
      This script parse the blast result file, get the length of the sequence file and 
      calculate the coverage for each sequence.

      Also can filter the sequences and the blast file based in a coverage and identity cutoff value.

      It will produce as STDOUT: query_id subject_id query_length coverage identity gaps

    Usage:
  
       blast_hitcoverage.pl [-h] -r <blast_result_file> -s <sequence_file> 
                            [-c <coverage>] [-i <percentage_identity>] [-g <gaps_permited>] [-H] > output.tab
      
    Flags:
      
       -r <blast_result_file>            blast result file in m8 format (mandatory)
       -s <sequence_file>                sequence file in fasta format (mandatory)
       -c <coverage>                     percentage of the coverage (75 by default)
       -i <identity>                     percentage of the identity (90 by default)
       -g <gaps>                         number of gaps permited in the coverage analysis (10 by default)
       -H <print_header>                 print header
       -h <help>                         print this help


EOF
exit (1);
}
