#!/usr/bin/perl

=head1 NAME

 get_unigeneseq_from_sgndb.pl
 This script get unigene data from the SGN database (version.1.0.).

=cut

=head1 SYPNOSIS

 get_unigeneseq_from_sgndb.pl [-h] -H <database_host> -D <database_name> [-I | -G | -u <unigene_build> | -o <organism_name>]
 [-b <basename>] [-F | -Q | -M | -A] [-C | -S]  [-l] [-a]

=head2 I<Flags:>

=over


=item -H

B<database_host>                database host (mandatory)

=item -D

B<database_name>                database name (mandatory)

=item -G

B<get_all_uniegene>             this option get all the current unigene build... so be patient.

=item -I

B<get_info_about_unigene_build> get information about all the unigene builds from the database

=item -u

B<unigene_build>                unigene_build_id of the unigene_build to get

=item -o

B<organism_name>                organism_name of the CURRENT unigene_build to get

=item -b

B<basename>                     basename to the output files

=item -l

B<min_seq_length>               minimum sequence length to get a sequence

=item -F

B<get_sequence_in_fasta>        get the DNA sequence in fasta format (option by default).

=item -Q

B<get_quality_in_fasta>         get the qscores in fasta format

=item -M

B<get_member_mapping>           get the member mapping in tab format (-f1 unigene_id, -f2 est_id, -f3 consensus_seq_start, -f4 consensus_seq_stop, -f5 est_start, -f6 est_stop, -f7 direction)

=item -A

B<get_defline>                  get a defline for unigene as SGN-Udddddd\tannotation\n

=item -C

B<get_only_contigs>             get only the contigs (nr_member > 1)

=item -S

B<get_only_singlets>            get only the singlets (nr_member = 1)

=item -a

B<annotation_code>              see annotation note ('A' or combination of 'UB', 'MB', 'DB:dataset' or 'DFL:int')

=item -h

B<print_help>                   print this help

=back

=cut

=head1 DESCRIPTION

 This script get unigene data in many ways:

 * 1- Info option (-I): Get information about Unigene_build_id, Organism, Member_Nr, Contig_Nr, Singlet_Nr and Date.
      (for more information you can use REPORT_about_unigene.pl)

 * 2- Get sequence (-S) and/or qscore (-Q) in fasta format without annotation, with complete annotations (-A) or with annotations filtered (-a <filter>). There are two types of filter... an integer that cut the length of the annotation and a database name (as genbank).

 * 3- Get membership unigene mapping (-M) as a file with seven columns (-f1 unigene_id, -f2 est_id, -f3 consensus_seq_start, -f4 consensus_seq_stop, -f5 est_start, -f6 est_stop, -f7 direction)

 Annotation Note:

   The annotation parameter (-a) is a predefined argument that you can use to get specific annotations over the unigenes. There some options, that you can combine (except 'all' that you can not combine with others).

   + 'A': Get all the annotations, that means:
	        
          - Unigene_ID as SGN-Udddddd
	  - Unigene_build as Name_Unigene_build\#d
	  - Number of members as [d ESTs aligned]
	  - One or more dataset annotations from blast results, separated by ;. Appears the best hit as:
	    (*DBCode) DBAccession (e_value=dd) Defline_Annotation;
            
    + 'UB', get only the unigene_build as Name_Unigene_build\#d

    + 'MB', get only the unigene_member_nr as [d ESTs sligned]

    + 'DB:database_name', get the blast result annotation for a specific dataset as 'genbank'.
      The current possibilities are: 'genbank', 'arabidopsis' or 'swissprot'.

    + 'DFL:integer', limits the number of characters of a defline (it could be really long in genbank)(only for DB:tag).  
         
      So, it means that if you want an annotations as unigene_build, and annotation for genbank in no more than 60 characters, you should use: -A 'UB,DB:genbank,DFL:60'.
    
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 get_unigeneseq_from_sgndb.pl


=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use CXGN::DB::InsertDBH;
use Math::BigFloat;


our ($opt_H, $opt_D, $opt_I, $opt_G, $opt_u, $opt_o, $opt_b, $opt_F, $opt_Q, $opt_M, $opt_A, $opt_C, $opt_S, $opt_l,
     $opt_a, $opt_h);
getopts("H:D:IGu:o:b:FQMACSl:a:h");
if (!$opt_H && !$opt_D && !$opt_I && !$opt_G && !$opt_u && !$opt_o && !$opt_b && !$opt_F && !$opt_Q && !$opt_M && !$opt_A 
    && !$opt_C && !$opt_S && !$opt_l && !$opt_a && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
} elsif ($opt_h) {
    help();
}

print "\n\nValidating the input data...\n";
my ($dbh, $min_length, $unigene_build_id_aref)=validate_input();
my @unigene_build_ids=@$unigene_build_id_aref;
print "\t\t... input data is okay\n\n";

my $unigene_build_description=join(',', @unigene_build_ids);
print "Parameters:\n";
print "\tUnigene_builds: $unigene_build_description\n";
print "\tMinimum_sequence_length: $min_length\n";
print "\tOutput_files:\n";
my ($seqFH, $qualFH, $membermapFH, $deflineFH);
my $seqfile=$opt_b.'.seq';
my $qualfile=$opt_b.'.qual';
my $membermapfile=$opt_b.'.membermapping.tab';
my $deflinefile=$opt_b.'.annotations.tab';

if (!$opt_F && !$opt_Q && !$opt_M && !$opt_A) {
    open $seqFH, '>', $seqfile || die "I can not open the file:$seqfile\n\n";
    print "\t\t+ Sequence_file: $seqfile\n";
}
if ($opt_F) {
    open $seqFH, '>', $seqfile || die "I can not open the file:$seqfile\n\n";
    print "\t\t+ Sequence_file: $seqfile\n";
}
if ($opt_Q) {
    open $qualFH, '>', $qualfile || die "I can not open the file:$qualfile\n\n";
    print "\t\t+ Quality_file: $qualfile\n";
}
if ($opt_M) {
    open $membermapFH, '>', $membermapfile || die "I can not open the file:$membermapfile\n\n";
    print "\t\t+ Membermap_file: $membermapfile\n";
}
if ($opt_A) {
    open $deflineFH, '>', $deflinefile || die "I can not open the file:$deflinefile\n";
    print "\t\t+ Annotation_file: $deflinefile\n";
}
print "\n\tAnnotations_options:";
if (!$opt_a) {
    print "\tnone\n\n";
} else {
    print "\t$opt_a\n\n";
}

foreach my $unigene_build_id (@unigene_build_ids) {
    my $check="SELECT COUNT(unigene_id), build_nr, sgn.groups.comment FROM sgn.unigene 
               JOIN sgn.unigene_build USING(unigene_build_id)
               JOIN sgn.groups ON sgn.unigene_build.organism_group_id=sgn.groups.group_id 
               WHERE sgn.unigene_build.unigene_build_id=? GROUP BY build_nr, sgn.groups.comment";
    my $sth=$dbh->prepare($check);
    $sth->execute($unigene_build_id);
    my ($unigene_n, $build_nr, $comment)=$sth->fetchrow_array;
    my $u=0;
    
    print "\nGetting data for unigene_build_id = $unigene_build_id ($comment#$build_nr)...\n";
    my $query="SELECT unigene_id, consensi_id, nr_members FROM sgn.unigene WHERE unigene_build_id=? ORDER BY unigene_id";
    $sth=$dbh->prepare($query);
    $sth->execute($unigene_build_id);
    while (my ($unigene_id, $consensi_id, $nr_members)=$sth->fetchrow_array() ) {
	$u++;
	print "\tExtracting unigene: SGN-U$unigene_id ($u of $unigene_n)\t\t\t\r";
	if ($nr_members > 1) {
	    my $query2="SELECT seq, qscores FROM sgn.unigene_consensi WHERE consensi_id=?";
	    my $sth2=$dbh->prepare($query2);
	    $sth2->execute($consensi_id);
	    my ($seq, $qscores)=$sth2->fetchrow_array();
	    $qscores =~ s/\s+$//;
	    my $seqlength=length($seq);
	    my @qscores=split(/ /, $qscores);
	    my $qscorelength=scalar(@qscores);
	    if ($seqlength != $qscorelength ) {
		print STDERR "\nError!!!\tThe unigene:SGN-U$unigene_id (consensi_id=$consensi_id) has seq length and ";
		print STDERR "qscore length different ($seqlength | $qscorelength)\n";
	    }
	    if ($seqlength >= $min_length) {
		my $annotation;
		if (!$opt_S) {
		    if ($opt_a) {
			$annotation=get_annotations($unigene_id, $dbh, $build_nr, $comment);
		    }
		    if ($seqFH) {
			if ($annotation) {
			    print $seqFH ">SGN-U$unigene_id $annotation\n$seq\n";
			} else {
			    print $seqFH ">SGN-U$unigene_id\n$seq\n";
			}
		    }
		    if ($qualFH) {
			if ($annotation) {
			    print $qualFH ">SGN-U$unigene_id $annotation\n$qscores\n";
			} else {
			    print $qualFH ">SGN-U$unigene_id\n$qscores\n";
			}
		    }
		    if ($membermapFH) {
			my $query3="SELECT est_id, start, stop, qstart, qend, dir FROM sgn.unigene_member 
                                    WHERE unigene_id=?";
			my $sth3=$dbh->prepare($query3);
			$sth3->execute($unigene_id);
			while (my ($est_id, $start, $stop, $qstart, $qend, $dir)=$sth3->fetchrow_array) {
			    print $membermapFH "SGN-U$unigene_id\tSGN-E$est_id\t$start\t$stop\t$qstart\t$qend\t$dir\n";
			}
		    }
		    if ($deflineFH) {
			my $print_annotation=$annotation || ' ';
			print $deflineFH "SGN-U$unigene_id\t$print_annotation\n";
		    }
		}
	    }
	} elsif ($nr_members == 1) {
	    my $query2="SELECT seq, qscore, hqi_start, hqi_length, est_id FROM sgn.est 
                        JOIN sgn.qc_report USING(est_id) JOIN sgn.unigene_member USING(est_id) WHERE unigene_id=?";
	    my $sth2=$dbh->prepare($query2);
	    $sth2->execute($unigene_id);
	    my ($seq, $qscore, $start, $length, $est_id)=$sth2->fetchrow_array();
            $qscore =~ s/\s+$//;
            my $seqlength=length($seq);
	    my @qscore=split(/ /, $qscore);
	    my $qscorelength=scalar(@qscore);
	    if ($seqlength != $qscorelength ) {
		print STDERR "\nError!!!\tThe unigene:SGN-U$unigene_id before trimming (est_id=$est_id) has seq length ";
		print STDERR "and qscore length different ($seqlength | $qscorelength)\n";
	    }

	    my $trimseq=substr $seq, $start, $length;
	    my $trimseqlength=length($trimseq);
	    my @trimqscore;
	    if ($length > 1) {
		my $end=$length+$start-1;
		@trimqscore=@qscore[$start..$end];
	    } else {
		@trimqscore=$qscore[$start];
	    }
	    my $trimqscorelength=scalar(@trimqscore);
	    my $trimqscore="@trimqscore";
	    my $annotation;
	    if ($trimseqlength != $trimqscorelength ) {
		print STDERR "\nError!!!\tThe unigene:SGN-U$unigene_id after trimming (est_id=$est_id) has seq length";
		print STDERR " and qscore length different ($trimseqlength | $trimqscorelength)\n";
	    }
	    if ($trimseqlength >= $min_length) {
		if (!$opt_C) {
		    if ($opt_a) {
			$annotation=get_annotations($unigene_id, $dbh, $build_nr, $comment);
		    }
		    if ($seqFH) {
			if ($annotation) {
			    print $seqFH ">SGN-U$unigene_id $annotation\n$trimseq\n";
			} else {
			    print $seqFH ">SGN-U$unigene_id\n$trimseq\n";
			}
		    }
		    if ($qualFH) {
			if ($annotation) {
			    print $qualFH ">SGN-U$unigene_id $annotation\n$trimqscore\n";
			} else {
			    print $qualFH ">SGN-U$unigene_id\n$trimqscore\n";
			}
		    }
		    if ($membermapFH) {
			my $query3="SELECT est_id, start, stop, qstart, qend, dir FROM sgn.unigene_member 
                                    WHERE unigene_id=?";
			my $sth3=$dbh->prepare($query3);
			$sth3->execute($unigene_id);
			while (my ($est_id, $start, $stop, $qstart, $qend, $dir)=$sth3->fetchrow_array) {
			    print $membermapFH "SGN-U$unigene_id\tSGN-E$est_id\t$start\t$stop\t$qstart\t$qend\t$dir\n";
			}
		    }
		    if ($deflineFH) {
			my $print_annotation=$annotation || ' ';
			print $deflineFH "SGN-U$unigene_id\t$print_annotation\n";
		    }
			
		}
	    }
	}
    }
}
print "\n\ndone!!!\n\n";

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
       This script get unigene data in many ways:

          * 1- Info option (-I): Get information about Unigene_build_id, Organism, Member_Nr, Contig_Nr, 
               Singlet_Nr and Date. (for more information you can use REPORT_about_unigene.pl)

          * 2- Get sequence (-S) and/or qscore (-Q) in fasta format without annotation, with complete 
               annotations (-A) or with annotations filtered (-a <filter>). There are two types of filter... 
               an integer that cut the length of the annotation and a database name (as genbank). Other option is
               -G, that means get ALL the CURRENT unigenes... so be patient.

          * 3- Get membership unigene mapping (-M) as a file with seven columns (-f1 unigene_id, -f2 est_id, 
               -f3 consensus_seq_start, -f4 consensus_seq_stop, -f5 est_start, -f6 est_stop, -f7 direction).

	Annotation Note:

	       The annotation parameter (-a) is a predefined argument that you can use to get specific annotations
	  over the unigenes. There some options, that you can combine (except 'all' that you can not combine with others).

	    + 'A': Get all the annotations, that means:
	        
	         - Unigene_ID as SGN-Udddddd
		 - Unigene_build as Name_Unigene_build\#d
		 - Number of members as [d ESTs aligned]
		 - One or more dataset annotations from blast results, separated by ;. Appears the best hit as:
		   (*DBCode) DBAccession (e_value=dd) Defline_Annotation;
            
            + 'UB', get only the unigene_build as Name_Unigene_build\#d

	    + 'MB', get only the unigene_member_nr as [d ESTs sligned]

	    + 'DB:database_name', get the blast result annotation for a specific dataset as 'genbank'.
	      The current possibilities are: 'genbank', 'arabidopsis' or 'swissprot'.

	    + 'DFL:integer', limits the number of characters of a defline (it could be really long in genbank).
              (only for DB:tag)  
         
	      So, it means that if you want an annotations as unigene_build, and annotation for genbank in no more than
	  60 characters, you should use: -A 'UB,DB:genbank,DFL:60'.
    
    Usage: 

       get_unigeneseq_from_sgndb.pl -H <database_host> -D <database_name> [-I | -B | -u <unigene_build> | -o <organism_name>]
       [-b <basename>] [-F | -Q | -M | -A] [-C | -S] [-l] [-a <annotation_code>] [-h]
           
    Flags:

       -H <database_host>                database host (mandatory)

       -D <database_name>                database name (mandatory)

       -I <info_about_unigene>           get information about all the unigene builds from the database

       -G <get_all_unigenes>             get all the CURRENT unigene builds... so be patient  

       -u <unigene_build>                unigene_build_id of the unigene_build to get

       -o <organism_name>                organism_name of the CURRENT unigene_build to get

       -b <basename>                     basename to the output files (mandatory)

       -F <get_sequence_in_fasta>        get the sequence in fasta format (option by default).

       -l <min_sequence_length>          minimum length value to get a sequence

       -Q <get_quality_in_fasta>         get the qscores in fasta format

       -M <get_member_mapping>           get the member mapping in tab format

       -A <get_defline>                  get defline for unigene (as SGN-Udddddd\tannotation)

       -C <get_only_contigs>             get only the contigs (nr_member > 1)

       -S <get_only_singlets>            get only the singlets (nr_member = 1)

       -a <annotation_code>              see annotation note ('A' or combination of 'UB', 'MB', 'DB:dataset' or 'DFL:int')

       -h <print_help>                   print this help    

EOF
exit (1);
}

=head2 validate_input

  Usage: my ($dbh, $min_length, $unigene_build_aref)=validate_input()
  Desc: check if exists some mandatory inputs. Also check the compatibility of them.
  Ret: $dbh, database conection object and $min_length, minimum length to get a sequence and $unigene_build_aref
       a list of unigene builds
  Args: none
  Side_Effects: die if there are something wrong. It do a new connection with the database.
  Example: my ($dbh, $min_length, $unigene_build_aref)=validate_input();

=cut

sub validate_input {

  if (!$opt_H || !$opt_D) {
      die "Sorry, database paramaters -H <hostname> or/and -D <databasename> were not supplied.\n\n";
  }
  my $i=0;
  if ($opt_I) {
      $i++;
  }
  if ($opt_G) {
      $i++;
  }
  if ($opt_u) {
      $i++;
  }
  if ($opt_o) {
      $i++;
  }
  if ($i == 0) {
      print "Sorry, input parameters were not supplied. You must use one of these input parameters\n";
      print "\t-I <get_info_about_all_unigene_builds>\n";
      print "\t-G <get_all_the_CURRENT_unigene_builds>\n";
      print "\t-u <get_the_unigene_build_id=X>\n";
      print "\t-o <get_the_current_unigene_buil_for_organism=Y>\n";
      die "For more information use -help or perldoc\n\n";
  } elsif ($i > 1) {
      die "Sorry, you only can use one of these options -I, -G, -u or -o.\n\n";
  }

  if ($opt_A && !$opt_a) {
      die "Sorry, -A <get_deflines> only can be used with -a <annotation> option\n\n";
  }

  if (!$opt_I && !$opt_b) {
      die "Sorry, -b <output_basename> parameter was not supplied.\n\n";
  }
  my $min_length=$opt_l || 50 ;
  
  if ($opt_C && $opt_S) {
      die "Sorry, you can not use -C <get_only_contig> and -S <get_only_singlet> at the same time.\n\n";
  }

  my $dbh = CXGN::DB::InsertDBH::connect({       dbname=>$opt_D ,
                                                 dbhost=>$opt_H ,
                                                 dbargs=>{ RaiseError => 1, PrintError => 1 }
                                        })
        or die "Failed to connect to database ($DBI::errstr)\n\n";
  print "\n\n";

  my @unigene_build_list;
  if ($opt_I) {
      print_unigene_report($dbh);
      die "\n\n";
  } elsif ($opt_G) {
      my $query="SELECT unigene_build_id FROM sgn.unigene_build WHERE status=? ORDER BY unigene_build_id";
      my $sth=$dbh->prepare($query);
      $sth->execute('C');
      while (my $unigene_build_id=$sth->fetchrow_array) {
	  push @unigene_build_list, $unigene_build_id;
      }
  } elsif ($opt_u) {
      my @unigenes_ids = split(/,/, $opt_u);
      foreach my $unigene_id (@unigenes_ids) {
	  if ($unigene_id =~ m/^\d+$/) {
	      my $check="SELECT unigene_build_id FROM sgn.unigene_build WHERE unigene_build_id=?";
	      my $sth=$dbh->prepare($check);
	      $sth->execute($unigene_id);
	      my $right_check=$sth->fetchrow_array;
	      if ($right_check) {
		  push @unigene_build_list, $right_check;
	      } else {
		  die "Sorry, the unigene_build_id=$opt_u doesn't exists in the database\n\n";
	      }
	  } else {
	      die "Sorry, -u <unigene_build_id> must be an integer\n\n";
	  }
      }
  } elsif ($opt_o) {
      my $check="SELECT unigene_build_id FROM sgn.unigene_build 
                 JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id
                 JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id
                 WHERE organism_name ILIKE ? AND status=? ORDER BY build_date DESC LIMIT 1";
      my $organism="%".$opt_o."%";
      my $sth=$dbh->prepare($check);
      $sth->execute($organism, 'C');
      my $unigene_build_id=$sth->fetchrow_array();
      if (!$unigene_build_id) {
	  die "Sorry, the organism=$opt_o haven't any unigene_build_id with status 'current'.\n\n";
      } else {
	  push @unigene_build_list, $unigene_build_id;
      }
  }
  return ($dbh, $min_length, \@unigene_build_list);
}


=head2 print_unigene_report

  Usage: print_unigene_report($dbh);
  Desc: get unigene data and print in a friedly format (as table)
  Ret: none
  Args: $dbh, a database conection object
  Side_Effects: none
  Example: print_unigene_report($dbh);

=cut 

sub print_unigene_report {
  my $dbh=shift;

  ## Option I means get all the information, print them in a user friendly format and die, so first get data
  my (%organisms, %build_nr, %build_date, %status, %unigene_n, %est_n); 
  my $report1="SELECT unigene_build_id, build_nr, build_date, status, COUNT(DISTINCT sgn.unigene.unigene_id), 
               COUNT(est_id) FROM sgn.unigene_build JOIN sgn.unigene USING(unigene_build_id) JOIN sgn.unigene_member 
               USING(unigene_id) GROUP BY unigene_build_id, build_nr, build_date, status";
  my $sth1=$dbh->prepare($report1);
  $sth1->execute();
       
  while (my ($unigene_build_id, $build_nr, $build_date, $status, $unigene_n, $est_n)=$sth1->fetchrow_array() ) {
      $build_nr{$unigene_build_id}=$build_nr;
      $build_date{$unigene_build_id}=$build_date;
      $status{$unigene_build_id}=$status;
      $unigene_n{$unigene_build_id}=sprintf('%9s', $unigene_n);
      $est_n{$unigene_build_id}=sprintf('%9s', $est_n);
	
      my @organism_list;
      my $report2="SELECT organism_name FROM sgn.organism 
                   JOIN sgn.group_linkage ON sgn.organism.organism_id=sgn.group_linkage.member_id 
                   JOIN sgn.unigene_build ON sgn.group_linkage.group_id=sgn.unigene_build.organism_group_id
                   WHERE unigene_build_id=? ORDER BY organism_name";
      my $sth2=$dbh->prepare($report2);
      $sth2->execute($unigene_build_id);
      while (my $organism_name=$sth2->fetchrow_array() ) {
	  $organism_name =~ s/\(.+\)$//g;
	  my $organism_name_f=sprintf('%50s', $organism_name);
	  push @organism_list, $organism_name_f;
      }
      $organisms{$unigene_build_id}=\@organism_list;
  }
      
  ## Now print them in a friendly format 
  my @unigene_build_ids_list=sort {$a <=> $b} keys %organisms;
     
  ## Print the header
  print "+----------+----------------------------------------------------+----------+------------+";
  print "--------+-----------+-----------+\n";
  print "| BUILD_ID | ORGANISM_NAME                                      | BUILD_NR | BUILD_DATE |";
  print " STATUS | UNIGENE_N | EST_N     |\n";
  print "+==========+====================================================+==========+============+";
  print "========+===========+===========+\n";

  ## Print the table
  foreach my $unigene_build_id (@unigene_build_ids_list) {
      my @organisms_list=@{$organisms{$unigene_build_id}};
      my $first_organism=shift(@organisms_list);
      my $rest_count=scalar(@organisms_list);
      my $print_build_id=sprintf('%8s', $unigene_build_id);
      my $print_build_nr=sprintf('%8s', $build_nr{$unigene_build_id});
      my $print_build_date=sprintf('%10s', $build_date{$unigene_build_id});
      my $print_status=sprintf('%6s', $status{$unigene_build_id});
      print "| $print_build_id | $first_organism | $print_build_nr | $print_build_date |";
      print " $print_status | $unigene_n{$unigene_build_id} | $est_n{$unigene_build_id} |\n";
      if ($rest_count > 0) {
          foreach my $rest_organism (@organisms_list) { 
     	  print "|          | $rest_organism |          |            |        |           |           |\n";
          }
      }
      print "+----------+----------------------------------------------------+----------+------------+";
      print "--------+-----------+-----------+\n";
  }
}

=head2 get_annotations

  Usage: my $annotation=get_annotations($unigene_id, $dbh, $build_nr, $comment_group)
  Desc: get and parse annotations from the sgn.blast tables (get the annotation with highest score)
  Ret: $annotation, a scalar with the annotations
  Args: $dbh, a database conection object, $unigene_id, an integer, $build_nr, and $comment_group used in UB
  Side_Effects: none
  Example: my $annotation=get_annotation($unigene_id, $dbh, $build_nr, $comment_group)

=cut

sub get_annotations {
  my $unigene_id=shift;
  my $dbh=shift;
  my $build_nr=shift;
  my $comment_group=shift;

  my ($switch_ub, $switch_um)=(0,0);
  my $defline_charlimit;
  my @selected_datasets;
  my @parsed_a=split(/,/, $opt_a);
  foreach my $parsed_a (@parsed_a) {
      if ($parsed_a eq 'A') {
	  $switch_ub=1;
	  $switch_um=1;
	  @selected_datasets=("'genbank\/nr'","'arabidopsis\/peptide'","'swissprot'");
      } elsif ($parsed_a eq 'UB') {
	  $switch_ub=1;
      } elsif ($parsed_a eq 'MB') {
	  $switch_um=1;
      } elsif ($parsed_a =~ m/DFL:/i) {
	  if ($parsed_a =~ m/^DFL:(\d+)$/) {
	      $defline_charlimit=$1;
	  } else {
	      print "/nWarning!!! The -a DFL value is not an integer, it will not use DFL option/n";
	  }
      } elsif ($parsed_a =~ m/^DB:(\w+)/) {
	  if ($1 =~ m/genbank/i) {
	      my $dataset="'genbank\/nr'";
	      push @selected_datasets, $dataset;
	  } elsif ($1 =~ m/arabidopsis/i) {
	      my $dataset="'arabidopsis\/peptide'";
	      push @selected_datasets, $dataset;
	  } elsif ($1 =~ m/swissprot/i) {
	      my $dataset="'swissprot'";
	      push @selected_datasets, $dataset;
	  } else {
	      print "/nWarning!!! the DB:$1 is not in the database. DB only can be: genbank, arabidopsis and swissprot/n";
	  }
      } else {
	  die "\nError!!! You are using -a option with wrong tags ($parsed_a) instead 'A', 'UB', 'MB', 'DB' or 'DFL'\n\n";
      }
  }
  my $desclength=$defline_charlimit || '250';
  my $selected_datasets_nr=scalar(@selected_datasets);  
  my $selected_datasets='('.join(',', @selected_datasets).')';
  my $annotations;
  if ($switch_ub == 1) {
      $annotations .= $comment_group.'#'.$build_nr.' ';
  }
  
  if ($switch_um == 1) {
      my $query0="SELECT nr_members FROM sgn.unigene WHERE unigene_id=?";
      my $sth0=$dbh->prepare($query0);
      $sth0->execute($unigene_id);
      my $nr_member=$sth0->fetchrow_array();
      $annotations .= '['.$nr_member.' ESTs aligned] ';
  }
  if ($selected_datasets_nr > 0) {
      my $query="SELECT blast_target_id, db_name FROM sgn.blast_targets
                 JOIN sgn.blast_annotations USING(blast_target_id)
                 JOIN sgn.unigene ON sgn.blast_annotations.apply_id=sgn.unigene.unigene_id
                 WHERE unigene_id=? AND db_name IN $selected_datasets
                 GROUP BY db_name, sgn.blast_targets.blast_target_id
                 ORDER BY sgn.blast_targets.blast_target_id;";
      my $sth=$dbh->prepare($query);
      $sth->execute($unigene_id); 
      while(my ($blast_target_id, $defline_db)=$sth->fetchrow_array()) {
         my ($dbtag, $presingle_description, $single_description, $trimmed_desc);
	 my $query2="SELECT sgn.blast_defline.target_db_id, sgn.blast_defline.defline, sgn.blast_hits.evalue
                     FROM sgn.blast_defline
                     JOIN sgn.blast_hits USING(defline_id)
                     JOIN sgn.blast_annotations USING(blast_annotation_id)
                     WHERE sgn.blast_annotations.apply_id=? AND sgn.blast_annotations.blast_target_id=?
                     ORDER BY sgn.blast_hits.score DESC LIMIT 1";
	 my $sth2=$dbh->prepare($query2);
	 $sth2->execute($unigene_id, $blast_target_id);
	 my ($accession, $defline, $e_value)=$sth2->fetchrow_array();
	 if ($accession) {
	     if ($defline) {
		 if ($defline_db =~ m/genbank/i) {
		     $dbtag='GB';
		     $defline =~ s/>gi/gi/g;
		     my @gi=split(/gi\|/, $defline);
		     my $gi_count=scalar(@gi);
		     if ($gi_count > 1) {
			 $presingle_description='gi|'.$gi[1];
		     } else {
			 $presingle_description=$defline;
		     }
		     my @parse=split(/ /, $presingle_description);
		     if ($parse[0] =~ m/^gi\|/i) {
			 my $id=shift(@parse);
		     }
		     $single_description="@parse";
		 } elsif ($defline_db =~ m/arabidopsis/i) {
		     $dbtag='ATH';
		     $presingle_description=$defline;
		     my @parse=split(/ /, $presingle_description);
		     if ($parse[0] =~ m/^At/i) {
			 my $id=shift(@parse);
		     }
		     $single_description="@parse";
		 } elsif ($defline_db =~ m/swissprot/i) {
		     $dbtag='SWP';
		     $presingle_description=$defline;
		     my @parse=split(/ /, $presingle_description);
		     if ($parse[0] =~ m/^sp\|/i) {
			 my $id=shift(@parse);
		     }
		     $single_description="@parse";
		 }
		 if ($single_description) {
		     my $length_desc=length($single_description);
		     if ($length_desc > $desclength) {
			 $trimmed_desc=substr($single_description, 1, $desclength).'...';
			 $trimmed_desc .= ' ; ';
		     } else {
			 $trimmed_desc=$single_description.'; ';
		     }
		 }
	     } else {
		 $trimmed_desc='; ';
	     }
	     my $complete_desc='(*'.$dbtag.') '.$accession.' (e_value='.$e_value.') '.$trimmed_desc;
	     $annotations .= $complete_desc;
	 }
     }
  }
  return $annotations;
}


