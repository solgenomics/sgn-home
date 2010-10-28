#!/usr/bin/perl

=head1 NAME

 full_length_screen.pl.
 A script to estimate how many full length sequence there are into a fasta file using a protein file (version.1.0.).

=cut

=head1 SYPNOSIS

full_length_screen.pl -d <full_length_dataset> -b <blastx_results_in_m8> [-o <output_filename>] [-i <identity_cutoff] [-l
<length_percentage_cutoff>] 
    
=head2 I<Flags:>

=over
      
=item -o

B<output_file>              output filename (by default: possible_fulllength.tab)

=item -d

B<full_length_dataset>      full length dataset in fasta format (mandatory)

=item -b

B<blastx_results_in_m8>     blastx results file in m8 format (mandatory)

=item -i

B<identity_cutoff>          miminum identity percentage to count the match as positive

=item -l 

B<length_percentage_cutoff> minimum length percentage to tag a sequence as full length 

=item -h 

B<help>                   show the help

=back

=cut

=head1 DESCRIPTION

 This script analyze the possible full length in a sequence based in the results of a blast against swissprot
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

full_length_screen.pl


=cut

use strict;
use warnings;

use Getopt::Std;
use File::Basename;
use Bio::SeqIO;

our ($opt_q, $opt_d, $opt_b, $opt_o, $opt_i, $opt_l, $opt_S, $opt_s, $opt_h);
getopts("q:d:b:o:i:l:Ss:h");
if (!$opt_q && !$opt_d && !$opt_b && !$opt_o && !$opt_i && !$opt_l && !$opt_S && !$opt_s && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
my $input_query=$opt_q;
my $fulllength_dataset=$opt_d || die "Requested argument -d <full_length_dataset> was not supplied\n";
my $blast_result=$opt_b || die "Requested argument -b <blastx_result> was not supplied\n";
my $minimum_identity_cutoff=$opt_i || 80;
my $minimum_aligm_perc=$opt_l || 99;
my $output_file=$opt_o || 'possible_full_length.txt';
my $p_fl_species_cutoff=$opt_s || 5;

if ($opt_s) {
    $opt_S=1;
}
    
open my $out, '>', $output_file || die "I can not open the output file: $output_file\n";
print "\nPARAMETERS:\n\tfull_length_dataset_file:\t\t\t$fulllength_dataset\n\tblast_result_file:\t\t\t\t$blast_result\n";
print "\tminimum_identity_percentage_cutoff:\t\t$minimum_identity_cutoff\n\tminimum_percentage_of_full_length_aligment:\t$minimum_aligm_perc\n\toutput file:\t\t\t\t\t$output_file\n\n";


print "LOADING the full length dataset file...\n";
my ($length_ds_href,$desc_ds_href)=get_length_from_input($fulllength_dataset);

my %swp_organism;
if ($opt_S) {
    print "PARSING SWISSPROT ANNOTATIONS...\n";
    my ($swp_name_href,$swp_os_href,$swp_gn_href,$swp_ec_href,$swp_sv_href)=parse_swissprot_description($desc_ds_href);
    %swp_organism=%$swp_os_href;
}

print "PARSING AND LOADING the blast result file...\n";
my $best_match_href=parse_blast_result($blast_result);
my @query_ids=sort keys %$best_match_href;
my (%full_lenght_dist,%organism_FL_count,%organism_FL_ident);
foreach my $ids (@query_ids) {
    my $best_match_aref=$$best_match_href{$ids};
    my $subject_id=$$best_match_aref[1];
    my $identity_percentage=$$best_match_aref[2];
    if ($identity_percentage >= $minimum_identity_cutoff) {
	my $subject_fulllength=$$length_ds_href{$subject_id};
	my $percent_aligment_fulllength=get_full_length_percentage($subject_fulllength,$best_match_aref);
        
        if ($percent_aligment_fulllength >= $minimum_aligm_perc) {
            if ($opt_S) {
		my $organism=$swp_organism{$subject_id};
		if (exists $organism_FL_count{$organism}) {
		    my $oc=$organism_FL_count{$organism};
		    $organism_FL_count{$organism}=$oc+1;
		    my @org_ident=@{$organism_FL_ident{$organism}};
		    push @org_ident, $identity_percentage;
		    $organism_FL_ident{$organism}=\@org_ident;
		} else {
		    $organism_FL_count{$organism}=1;
		    my @org_ident=($identity_percentage);
		    $organism_FL_ident{$organism}=\@org_ident;
		}
		print $out "$ids\t$subject_id\t$organism\t";
		print $out "possible_full_length ($percent_aligment_fulllength \% subject_length)\n";
	    } else {
		print $out "$ids\t$subject_id\tpossible_full_length ($percent_aligment_fulllength \% subject_length)\n";
	    }
        }       
	$full_lenght_dist{$ids}=$percent_aligment_fulllength;
    }
}
my ($total,$t,$c10,$c20,$c30,$c40,$c50,$c60,$c70,$c80,$c90,$c100,$c_full_length)=(0,0,0,0,0,0,0,0,0,0,0,0,0);
my @fulllengths=values %full_lenght_dist;
foreach my $value (@fulllengths) {
    $t++;
    $total += $value;
    if ($value < 10) {
	$c10++;
    } elsif ($value >= 10 && $value < 20) {
	$c20++;
    } elsif ($value >= 20 && $value < 30) {
	$c30++;
    } elsif ($value >= 30 && $value < 40) {
	$c40++;
    } elsif ($value >= 40 && $value < 50) {
	$c50++;
    } elsif ($value >= 50 && $value < 60) {
	$c60++;
    } elsif ($value >= 60 && $value < 70) {
	$c70++;
    } elsif ($value >= 70 && $value < 80) {
	$c80++;
    } elsif ($value >= 80 && $value < 90) {
	$c90++;
    } else {
	$c100++;
    }
    if ($value >= $minimum_aligm_perc) {
	$c_full_length++
    }
}

my ($first_met_c, $stop_codon_c, $pred_full_length_c, $without_predict_begin_or_end_c)=(0,0,0,0); 
if ($opt_q) {
    print "PROCESSING THE QUERY SEQUENCES ACOORDING WALKING METHODS...\n";
    my ($tag_begin_href, $tag_end_href)=predict_fulllength_walking_method($input_query, $best_match_href);
    my @id_list=sort keys %$tag_begin_href;
    my %tag_begin=%$tag_begin_href;
    my %tag_end=%$tag_end_href;
    my $output_tag_file='walking_predict_method_tag_file.txt';
    open my $tag_fh, '>', $output_tag_file || die "I can not open the file $output_tag_file\n";
    my $full_length_file='possible_full_length.tab';
    open my $full_length_fh, '>', $full_length_file || die "I can not open the file $full_length_file\n";

    foreach my $id_of_list (@id_list) {
	my $tag_begin_id=$tag_begin{$id_of_list};
	my $tag_end_id=$tag_end{$id_of_list};
	print $tag_fh "$id_of_list\t$tag_begin_id\t$tag_end_id\n";
	if ($tag_begin_id =~ m/first_methionine_position/ && $tag_end_id =~ m/stop_codon_position/) {
	    $pred_full_length_c++;
	    my $begin=$tag_begin_id;
	    $begin =~ s/first_methionine_position_//;
	    my $end=$tag_end_id;
	    $end =~ s/stop_codon_position_//;
	    print $full_length_fh "$id_of_list\t$begin\t$end\n";
	} elsif ($tag_begin_id =~ m/first_methionine_position/ && $tag_end_id =~ m/unknown_stop_codon/ ) {
	    $first_met_c++;
	} elsif ($tag_begin_id =~ m/unknown_first_methionine/ && $tag_end_id =~ m/stop_codon_position/ ) {
	    $stop_codon_c++;
	} else {
	    $without_predict_begin_or_end_c++;
	}
    }
}    


my $average=$total/$t;
my $query_n=`cut -f1 $blast_result | sort -u | wc -l`;
my $s=$c_full_length*$p_fl_species_cutoff/100;
chomp($query_n);
print "REPORT:\n";
print "\tquery_ids number in blast report file = $query_n\n";
print "\tquery_ids number with identity percentage >= $minimum_identity_cutoff = $t\n\n";
print "\tFull Length Cover Distribution (percentage_range = n_ids)\n";
print "\t [ 0 - 10] = $c10\n\t [10 - 20] = $c20\n\t [20 - 30] = $c30\n\t [30 - 40] = $c40\n\t [40 - 50] = $c50\n\t [50 - 60] = $c60\n\t [60 - 70] = $c70\n\t [70 - 80] = $c80\n\t [80 - 90] = $c90\n\t [90 -100] = $c100\n\t(possible full lengths: $c_full_length)\n\n";

if ($opt_S) {
    print "\tPossible Full Length for Species Distribution:\n";
    
    my @species=sort keys %organism_FL_count;
    my $minor_count;
    my %print_by_ident;
    foreach my $sp (@species) {
	my $count=$organism_FL_count{$sp};
	my @identity=@{$organism_FL_ident{$sp}};
	my $a=scalar(@identity);
	my ($g,$decimal);
	foreach my $ident_value (@identity) {
	    $g+=$ident_value;
	}
	my $i=$g/$a;
	my $identity_average=sprintf("%.3f",$i);
	my $specie=sprintf("%50s",$sp);
	if ($count >= $s) {
	    my $print = "$specie\t$count\t$identity_average % of identity (average)\n";
	    $print_by_ident{$identity_average}=$print;
	} else {
	    $minor_count += $count;
	}
    }
    my @list_asc=sort keys %print_by_ident;
    my @list_desc=reverse @list_asc;
    foreach my $item (@list_desc) {
	my $item_to_print=$print_by_ident{$item};
	print "$item_to_print";
    }
    print "\n\tother species (with less than $p_fl_species_cutoff % of the total of possible full length matches)\t$minor_count\n";
}
print "\n";

if ($opt_q) {
    print "Predict full length based in codon walking:\n";
    print "Sequences with first Met: $first_met_c\n";
    print "Sequences with stop codon: $stop_codon_c\n";
    print "Full length sequences (First Met and Stop codon): $pred_full_length_c\n";
    print "Sequences without first met or stop codon: $without_predict_begin_or_end_c\n";
}
print "\n";

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
         This script analyze the possible full length in a sequence based in the results of a blast against swissprot. 
         
         If the sequence match have not the complete full length, this script walk in 5 and 3 direction until find the first
      metionine or the stop codon.

    Usage:
      full_length_screen.pl [-h] -d <full_length_dataset> -b <blastx_results_in_m8> [-q <query_sequence_file>]
      [-o <output_filename>] [-i <identity_cutoff>] [-l <length_percentage_cutoff>]

    Flags:
       -q <query_sequence_file>          query sequence file in fasta format

       -d <full_length_dataset>          full length dataset (protein) (mandatory)

       -b <blastx_result_with_dataset>   blastx result in m8 format  (mandatory)
       
       -o <output_filename>              output filename (by default: possible_full_length.txt) 

       -i <identity_perc_cutoff>         minimum identity percentage to match a sequence (default: 80) 

       -l <length_percentage_cutoff>     minimum percentage of the length to tag as full length (default: 99) 

       -S <species_distribution>         print the distribution of the full length based in the swissprot anotation

       -h <help>                         print the help

EOF
exit (1);
}


=head2 get_length_from_input

  Usage: my ($length_href, $description_href)=get_length_from_input($input_fasta_seq);
  Desc: get the length and the descriptions for sequences in a fasta file
  Ret: Two hash references ($length_href and $description_href) with keys=query_id and values=sequence_length and description 
       respectively
  Args: $input_fasta_seq, a sequence file in fasta format
  Side_Effects: print the progress in the process
  Example: my ($length_href, $description_href)=get_length_from_input($input_fasta_seq);

=cut

sub get_length_from_input {
  my $inputfastaseq=shift;

  my (%length,%descrip);
  my $seq_n=`grep '>' $inputfastaseq | cut -d ' ' -f1 | wc -l`;
  chomp($seq_n);
  my $n=0;
  my $input_seq = Bio::SeqIO->new(  -file   =>      $inputfastaseq,
                                    -format =>      'fasta'    );
  while (my $seq_obj=$input_seq->next_seq() )  {
      $n++;
      my $id=$seq_obj->display_id();
      my $seqlength=$seq_obj->length();
      my $description=$seq_obj->desc();
      print "\tloading $id with length: $seqlength Aa (seq:$n of $seq_n)                      \r";
      $length{$id}=$seqlength;
      $descrip{$id}=$description;
  }
  print "\n";
  return (\%length,\%descrip);
}

=head2 parse_blast_result

  Usage: my $best_match_href=parse_blast_results($input);
  Desc: parse the blast results getting the best match
  Ret: $best_match_href, an hash reference with keys=query_id and value=array_reference with the best match
  Args: $input, the input file
  Side_Effects: print the progress of the parsing
  Example:  my $best_match_href=parse_blast_results($input);

=cut 

sub parse_blast_result {
  my $input=shift;

  my (%best_hit_score, %best_match);
  my $lines_n=`cut -f1 $input | wc -l`;
  chomp($lines_n);
  my $n=1;  
  open my $in_fh, '<', $input || die "I can not open the blast result input file: $input\n";
  while (<$in_fh>) {
      chomp($_);
      print "\tprocessing the line ($n) of $lines_n lines from the file:$input\r";
      my @line=split(/\t/, $_);
      my ($query_id, $subject_id, $identity_perc, $aligment_length, $mismatches, $gaps, 
	  $q_start, $q_end, $s_start, $s_end, $e_value, $hit_score)=@line;
      if (exists $best_hit_score{$query_id}) {
	  if ($best_hit_score{$query_id} < $hit_score) {
	      $best_hit_score{$query_id}=$hit_score;
	      $best_match{$query_id}=\@line;
	  }
      } else {
	  $best_hit_score{$query_id}=$hit_score;
	  $best_match{$query_id}=\@line;
      }
      $n++;
  }
  print "\n";
  return \%best_match;
}

=head2 get_full_length_percentage

  Usage: my $aligment_full_length_percentage=get_full_length_percentage($subject_length, $blast_best_match_aref);
  Desc: Calculate the percentage of the full length cover by an aligment
  Ret:  A scalar, $aligment_fulllenght_perc, with the percentage of full length in a aligment
  Args: $subject_lentgh, a scalar, the subject length and an array reference with the best match for the query with the subject
  Side_Effects: none
  Example:  my $aligment_full_length_percentage=get_full_length_percentage($subject_length, $blast_best_match_aref);

=cut

sub get_full_length_percentage {
  my $subject_length=shift;
  my $blast_best_match_aref=shift;

  my ($query_id, $subject_id, $identity_perc, $aligment_length, $mismatches, $gaps, 
      $q_start, $q_end, $s_start, $s_end, $e_value, $hit_score)=@$blast_best_match_aref;

  my $aligment_fulllength_perc=($aligment_length/$subject_length)*100;
  
  return $aligment_fulllength_perc;
}

=head2 parse_swissprot_description

  Usage: my ($name_href, $organism_href, $gene_name_href, $prot_code_href, $version_href)=parse_swissprot_description($description_href)
  Desc: get an array of hash references with part of the descriptions of the swissprot defline
  Ret: Five hash references with keys=swissprot_accessions and values=name, organism, gene_name, protein_existence_code and 
       sequence_version respectively
  Args: An hash reference, $description_href with keys=swissprot_accession and value=defline
  Side_Effects: none
  Example:  my ($name_href,$organism_href,$gene_name_href,$protcode_href,$version_href)=parse_swissprot_description($description_href)

=cut

sub parse_swissprot_description {
  my $description_href=shift;
  
  my %description=%$description_href;
  my (%name, %organism, %gene_name, %protein_existence_code, %sequence_version);
  my @ids= keys %description;
  my $n_ids=scalar(@ids);
  my $n=1;
  foreach my $id (@ids) {
      print "\tparsing description for id:$id ($n of $n_ids)          \r";
      $n++;
      my $ind_description=$description{$id};
      if ($ind_description =~ m/(.+)\s+OS=(.+)\s+GN=(.+)\s+PE=(\d+)\s+SV=(\d+)/) {
	  $name{$id}=$1;
	  $organism{$id}=$2;
	  $gene_name{$id}=$3;
	  $protein_existence_code{$id}=$4;
	  $sequence_version{$id}=$5;
      } elsif ($ind_description =~ m/(.+)\s+OS=(.+)\s+PE=(\d+)\s+SV=(\d+)/) {
	  $name{$id}=$1;
	  $organism{$id}=$2;
	  $protein_existence_code{$id}=$3;
	  $sequence_version{$id}=$4;
      } else {
	  print "\nI can not parse the description for the id:$id\n";
      }
  }
  print "\n";
  return (\%name, \%organism, \%gene_name, \%protein_existence_code, \%sequence_version);
}


=head2 methionine_walking_5coord

  Usage: my ($firstaa, $first_triplet_codon)=methionine_walking_5coord($seq_obj, $best_match_blast_href);
  Desc: Get a sequence and walk in 5 direction until find the first aa
  Ret: Two scalars, $first_aa, the first aa and $first_triplet_codon and first triplet condon
  Args: $seq_obj, a Bioperl sequence object and $best_match_blast_href, the best match blast href
  Side_Effect: none
  Example: my ($first_aa, $first_codon)=methionine_walking_5coord($seq_obj, $best_match_blast_href);

=cut

sub methionine_walking_5coord {
    my $seq_obj=shift;
    my $best_match_blast_href=shift;
    my %best_match_blast=%$best_match_blast_href;

    my $id=$seq_obj->display_id();
    if ($best_match_blast{$id}) {
        my $best_match_for_id_aref=$best_match_blast{$id};
	my ($query_id, $subject_id, $identity_perc, $aligment_length, $mismatches, $gaps, $q_start, $q_end, 
	    $s_start, $s_end, $e_value, $hit_score)=@$best_match_for_id_aref;
	
	my $first_triplet_start=$q_start;
	my $first_triplet_end=$first_triplet_start+2;
        if ($q_start > $q_end) {
	    $first_triplet_start=$q_end;
	    $first_triplet_end=$first_triplet_start+2;
	}
	#print "\tfirst_check_met_parameters: $first_triplet_start, $first_triplet_end\n";
	my $first_aa_obj=$seq_obj->trunc($first_triplet_start,$first_triplet_end)->translate();
	my $first_aa=uc $first_aa_obj->seq();

	while ($first_aa ne 'M' && $first_triplet_start > 3) {
	    $first_triplet_start -= 3;
	    $first_triplet_end = $first_triplet_start+2;
	    $first_aa_obj= $seq_obj->trunc($first_triplet_start,$first_triplet_end)->translate();
	    $first_aa=uc $first_aa_obj->seq();   
	}
        #print "\tcheck_met_parameters: $first_triplet_start, $first_triplet_end, $first_aa\n";
        return ($first_aa, $first_triplet_start);	
    }
}

=head2 stop_codon_3walking

  Usage: my ($last_aa, $last_stop_codon)=stop_codon_3walking($seq_obj, $best_match_blast_href);
  Desc: find, if is possible the last aa and the last stop codon in a sequence object
  Ret: $last_aa, last aminoacid and $last_stop_codon, last stop codon
  Args: $seq_obj, a Bioperl sequence object and $best_match_blast_href, the best match blast href
  Side_Effects: none
  Example: my ($last_aa, $last_stop_codon)=stop_codon_3walking($seq_obj, $best_match_blast_href);


=cut 

sub stop_codon_3walking {
    my $seq_obj=shift;
    my $best_match_blast_href=shift;
    my %best_match_blast=%$best_match_blast_href;

    my $id=$seq_obj->display_id();
    my $length=$seq_obj->length();
    my $length_limit=$length-2;
    if (exists $best_match_blast{$id}) {
	my $best_match_for_id_aref=$best_match_blast{$id};
	my ($query_id, $subject_id, $identity_perc, $aligment_length, $mismatches, $gaps, $q_start, $q_end, 
	    $s_start, $s_end, $e_value, $hit_score)=@$best_match_for_id_aref;

	my $last_triplet_start=$q_end-2;
	my $last_triplet_end=$q_end;
	if ($q_start > $q_end) {
	    $last_triplet_start=$q_start-2;
	    $last_triplet_end=$q_start;
	}
	#print "\tfirst_check_stop_parameters: $length_limit, $last_triplet_start, $last_triplet_end\n";
	my $last_aa_obj=$seq_obj->trunc($last_triplet_start,$last_triplet_end)->translate();
        my $last_aa=uc $last_aa_obj->seq();        
	
	while ($last_aa ne '*' && $last_triplet_end < $length_limit) {
	    $last_triplet_start += 3;
	    $last_triplet_end = $last_triplet_start+2;
	    $last_aa_obj=$seq_obj->trunc($last_triplet_start,$last_triplet_end)->translate();
	    $last_aa=uc $last_aa_obj->seq();
	}
	#print "\tcheck_stop_parameters: $length_limit, $last_triplet_start, $last_triplet_end, $last_aa\n";
        return ($last_aa, $last_triplet_end);
    }
}

=head2 predict_fulllength_walking_method

  Usage: my ($five_extreme_tags_href, $three_extreme_tags_href)=predict_fulllength_walking_method($query_file, $best_match_href);
  Desc: Get the sequence from a file and walk over the sequence in 5 and 3 extreme until find Met and Stop
  Ret: two hash references with keys=sequence_id and values=tags
  Args: $query_file, input file with the sequences and $best_match_href, best match hash reference
  Side_Effects: none
  Example: my ($5tags_href, $3tags_href)=predict_fulllength_walking_method($query_file, $best_match_href);

=cut

sub predict_fulllength_walking_method {
    my $query_file=shift;
    my $best_match_href=shift;
    my %best_match=%$best_match_href;

    my ($five_coord_tag, $three_coord_tag);
    my (%fl5_tag, %fl3_tag);
    my $n=0;
    my $total=`grep '>' $query_file | cut -d ' ' -f1 | wc -l`;
    chomp($total);
    my $input_seq = Bio::SeqIO->new(  -file   =>      $query_file,
                                      -format =>      'fasta'    );
    while (my $seq_obj=$input_seq->next_seq() )  {
      $n++;
      print "\tprocessing sequence $n of $total         \r";
      my $id=$seq_obj->display_id();
      my $length=$seq_obj->length();
      if (exists $best_match{$id}) {
	  my ($first_aa, $first_triplet_begin)=methionine_walking_5coord($seq_obj, $best_match_href);
	  my ($last_aa, $last_triplet_end)=stop_codon_3walking($seq_obj, $best_match_href);

     
	  if ($first_aa eq 'M') {
	      $five_coord_tag='first_methionine_position_'.$first_triplet_begin;
	  } else {
	      $five_coord_tag='unknown_first_methionine';
	  }
	  if ($last_aa eq '*') {
	      $three_coord_tag='stop_codon_position_'.$last_triplet_end;
	  } else {
	      $three_coord_tag='unknown_stop_codon';
	  }
	  $fl5_tag{$id}=$five_coord_tag;
	  $fl3_tag{$id}=$three_coord_tag;
      } else {
	  $fl5_tag{$id}="$id"."_is_not_in_the_blast_result_file";
	  $fl3_tag{$id}="$id"."_is_not_in_the_blast_result_file";
      }
  }
  print "\n";  
  return (\%fl5_tag, \%fl3_tag);
}
