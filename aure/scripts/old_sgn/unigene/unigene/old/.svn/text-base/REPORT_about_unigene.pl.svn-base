=head1 NAME

 REPORT_organism_for_unigene.pl.
 A script to get information about unigene (build and annotations) tables for each organism from SGN database  (version.1.0.).

=head1 SYPNOSIS
  
 REPORT_organism_for_unigene.pl [-h] -H <dbhost> -D <dbname> [-o <organism_name>]

=head1 DESCRIPTION

This script get a table from the database with unigene information.

=head1 AUTHOR

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=head1 METHODS
 
REPORT_organism_for_trace_processing.pl


=cut

use strict;

use CXGN::DB::InsertDBH;
use Getopt::Std;

our ($opt_D, $opt_H, $opt_o, $opt_C, $opt_h);
getopts("D:H:o:Ch");

if (!$opt_D && !$opt_H && !$opt_o && !$opt_C && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
   help();
}

if (!$opt_D || !$opt_H) { 
      die "\nSorry, this script needs -D (database name) and -H (host) arguments.\n";
  }

our $dbh = CXGN::DB::InsertDBH::connect ({ 	dbname => $opt_D, 
						dbhost => $opt_H,
					  });

my $var;
my (@unigene_build_global, @unigene_cds_annot_global, @unigene_blast_annot_global);
if ($opt_C) {
    $var="WHERE status='C'";
}
my $query;
if ($opt_o) {
    $var =~ s/WHERE/AND/gi;
    my $organism_q="'%".$opt_o."%'";
    $query="SELECT unigene_build_id, build_date FROM sgn.unigene_build 
               JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id 
               JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id 
               WHERE sgn.organism.organism_name ILIKE $organism_q $var ORDER BY unigene_build_id";
} else {
    $query="SELECT unigene_build_id, build_date FROM sgn.unigene_build $var ORDER BY unigene_build_id";
}

my $sth=$dbh->prepare($query);
$sth->execute();

while (my ($unigene_build_id, $build_date)=$sth->fetchrow_array() ) {
    if (!$unigene_build_id) {
        die "There are unigene_build_id in the database for the arguments specified\n";
    }
    my (@unigene_build_report, @unigene_cds_annot_report, @unigene_blast_annot_report);
    push @unigene_build_report, $unigene_build_id;
    push @unigene_cds_annot_report, $unigene_build_id;
    push @unigene_blast_annot_report, $unigene_build_id;

    push @unigene_build_report, $build_date;
    my $organism_est_count=0;
    my @organismname=get_organism_names($dbh, $unigene_build_id);
    my $organismname_aref=\@organismname;
    
    my $est_count=count_ests($dbh, $unigene_build_id);
    push @unigene_build_report, $est_count;

    my $unigene_members_count=count_members($dbh, $unigene_build_id);
    push @unigene_build_report, $unigene_members_count;

    my $unigenes_count=count_unigenes($dbh, $unigene_build_id);
    push @unigene_build_report, $unigenes_count;

    my $singletons_count=count_singletons($dbh, $unigene_build_id);
    push @unigene_build_report, $singletons_count;

    my $contigs_count=count_contigs($dbh, $unigene_build_id);
    push @unigene_build_report, $contigs_count;
    push @unigene_build_report, $organismname_aref;   
    
    #my @cds_reports=count_cds($dbh, $unigene_build_id);
    #if (!@cds_reports) {
    #    @cds_reports=('None prediction cds');
    #}
    #push @unigene_cds_annot_report, \@cds_report;
    
    my @domain_matches_report=count_cds_with_domain_match($dbh, $unigene_build_id);
    if (!@domain_matches_report) {
        @domain_matches_report=('None', 0, 0, 'None', 0);
    }
    push @unigene_cds_annot_report, \@domain_matches_report;

    my @blast_reports=count_unigene_with_blast_annotations($dbh, $unigene_build_id);
    if (!@blast_reports) {
        @blast_reports=('None', 0, 0, 0);
    }
    push @unigene_blast_annot_report, \@blast_reports;
   
    push @unigene_build_global, \@unigene_build_report;
    push @unigene_cds_annot_global, \@unigene_cds_annot_report;
    push @unigene_blast_annot_global, \@unigene_blast_annot_report;
}

my @ub_header= ('UB_ID', 'BUILD_DATE', 'EST_N', 'MEMBER_N', 'UNIGENE_N', 'SINGLETS_N', 'CONTIGS_N', 'ORGANISM');  
unshift @unigene_build_global, \@ub_header;
print_simple_aoa_table(\@unigene_build_global);
print "\n\n";

my @uca_header= ('UB_ID', 'METHOD', 'CDS_RUN_ID', 'CDS_N', 'DOMAIN_RUN_ID', 'CDS_WITH_DOMAIN');
unshift @unigene_cds_annot_global, \@uca_header;
print_complex_aoa_table(\@unigene_cds_annot_global);
print "\n\n";

my @uba_header= ('UB_ID', 'DB_TARGET', 'LAST_UPDATED_DB-FIELD', 'LAST_UPDATED_TIME', 'UNIGENE_BLASTMATCH_N');
unshift @unigene_blast_annot_global, \@uba_header;
print_complex_aoa_table(\@unigene_blast_annot_global);            
print "\n\n";

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
      This script check different data in the trace tables related with the unigene build.

    Usage: 
      REPORT_about_unigene.pl [-h] -H <dbhost> -D <dbname> [-o <organism_name>] [-C]
    
    Example:
       REPORT_about_unigene.pl -H localhost -D sandbox -C

    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -o organism name        the organism name inside single quotes ('xxxx')
      -C enalble option       get only the current unigene builds
      -h print this help

EOF
exit (1);
}

=head2 get_organism_names

  Usage: my @organism_name=get_organism_names($dbh, $unigene_build_id);
  Desc: get the organism names associated to an unigene build.
  Ret: An array, @organism_name, with the organisms name associated to an unigene build.
  Args: $dbh, database connection object and $unigene_build_id, unigene_build_id, an integer.
  Side_Effects: none
  Example: my @organism_names=get_organism_names($dbh, $unigene_build_id);

=cut

sub get_organism_names {
  my $dbh=shift;
  my $unigene_build_id=shift;
  my @organism_names;

  my $query="SELECT sgn.organism.organism_name FROM sgn.unigene_build 
             JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id
             JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id
             WHERE sgn.unigene_build.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  while (my ($organism_name) = $sth->fetchrow_array() ) {  
      push @organism_names, $organism_name;
  }
  return @organism_names;
}

=head2 count_ests

  Usage: my $est_count=count_est($dbh, $unigene_build_id);
  Desc: count the ALL the est (also the est with status and flags != 0) for the organisms used in a unigene build. 
  Ret: A scalar, a est count.
  Args: $dbh, database connection object and $unigene_build_id, unigene_build_id, an integer.
  Side_Effects: none
  Example: my $est_count=count_est($dbh, $unigene_build_id);

=cut

sub count_ests {
  my $dbh=shift;
  my $unigene_build_id=shift;

  my $query="SELECT COUNT(sgn.est.est_id) FROM sgn.unigene_build 
             JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id
             JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id
             JOIN sgn.library ON sgn.organism.organism_id=sgn.library.organism_id
             JOIN sgn.clone ON sgn.library.library_id=sgn.clone.library_id
             JOIN sgn.seqread ON sgn.clone.clone_id=sgn.seqread.clone_id
             JOIN sgn.est ON sgn.seqread.read_id=sgn.est.read_id
             WHERE sgn.unigene_build.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  my ($est_count) = $sth->fetchrow_array();
  return $est_count;
}  

=head2 count_members

  Usage: my $members_count=count_members($dbh, $unigene_build_id);
  Desc: count unigene members to an unigene build
  Ret: A scalar, $member_count, an integer
  Args: $dbh, database connection object and $unigene_build_id, unigene_build_id, an integer.
  Side_Effects: none
  Examples: my $member_count=count_members($dbh, $unigene_build_id);

=cut

sub count_members {
  my $dbh=shift;
  my $unigene_build_id=shift;

  my $query="SELECT COUNT(sgn.unigene_member.unigene_member_id) FROM sgn.unigene
             JOIN sgn.unigene_member ON sgn.unigene.unigene_id=sgn.unigene_member.unigene_id
             WHERE sgn.unigene.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  my ($members_count) = $sth->fetchrow_array();
  return $members_count;
}

=head2 count_clusters

  Usage: my $clusters_count=count_clusters($dbh, $unigene_build_id);
  Desc: count the number of cluster for an unigene build (the unigenes with cluster_no = -1 
      are sequences without real cluster)
  Ret: A scalar, $clusters_count (An integer)
  Args: $dbh (database connection object) and $unigene_build_id;
  Side_Effects: none
  Examples: my $clusters_count=count_clusters($dbh, $unigene_build_id);

=cut

sub count_clusters {
  my $dbh=shift;
  my $unigene_build_id=shift;

  my $query="SELECT COUNT(DISTINCT(sgn.unigene.cluster_no)) FROM sgn.unigene 
             WHERE sgn.unigene.cluster_no != -1 AND sgn.unigene.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  my ($clusters_count) = $sth->fetchrow_array();
  return $clusters_count;
}  

=head2 count_seq_outside_cluster

  Usage: my $seq_outside_cluster_count=count_seq_outside_cluster($dbh, $unigene_build_id);
  Desc: count the number of unigene with cluster_no = -1 (sequences that are not in the clusters)
  Ret: A scalar, $seq_outside_cluster_count (an integer)
  Args: $dbh (database connection object) and $unigene_build_id
  Side_Effects: none
  Examples: my $seq_outside_cluster_count=count_seq_outside_cluster($dbh, $unigene_build_id);

=cut

sub count_seq_outside_cluster {
  my $dbh=shift;
  my $unigene_build_id=shift;

  my $query="SELECT COUNT(sgn.unigene.unigene_id) FROM sgn.unigene 
             WHERE sgn.unigene.cluster_no = 1 AND sgn.unigene.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  my ($seq_outside_cluster_count) = $sth->fetchrow_array();
  return $seq_outside_cluster_count;   
}  

=head2 count_unigenes

  Usage: my $unigenes_count=count_unigene($dbh, $unigene_build_id);
  Desc: count the number of unigene for the specified unigene_build_id
  Ret: A scalar, $unigene_count (an integer)
  Args: $dbh (database connection object) and $unigene_build_id 
  Side_Effects: none
  Examples: my $unigenes_count=count_unigene($dbh, $unigene_build_id);

=cut

sub count_unigenes {
  my $dbh=shift;
  my $unigene_build_id=shift;

  my $query="SELECT COUNT(sgn.unigene.unigene_id) FROM sgn.unigene WHERE sgn.unigene.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  my ($unigenes_count) = $sth->fetchrow_array();
  return $unigenes_count;
}

=head2 count_singletons

  Usage: my $singletons_count=count_singletons($dbh, $unigene_build_id);
  Desc: count the number the singletons (unigenes composed by only one sequence, these could be in a cluster).
  Ret: A scalar, $singletons_count (an integer)
  Args: $dbh (database connection object) and $unigene_build_id
  Side_Effects: none 
  Example: my $singletons_count=count_singletons($dbh, $unigene_build_id);

=cut

sub count_singletons {
  my $dbh=shift;
  my $unigene_build_id=shift;

  my $query="SELECT COUNT(sgn.unigene.unigene_id) FROM sgn.unigene 
             WHERE sgn.unigene.nr_members=1 AND sgn.unigene.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  my ($singletons_count) = $sth->fetchrow_array();
  return $singletons_count;
}

=head2 count_contigs

  Usage: my $contigs_count=count_contigs($dbh, $unigene_build_id);
  Desc: count the number of contigs (unigenes composed by more than one sequences)
  Ret: A scalar, $contigs_count (an integer)
  Args: $dbh (database connection object) and $unigene_build_id
  Side_Effects: none 
  Example: my $contigs_count=count_contigs($dbh, $unigene_build_id);

=cut

sub count_contigs {
  my $dbh=shift;
  my $unigene_build_id=shift;

  my $query="SELECT COUNT(DISTINCT(sgn.unigene.consensi_id)) 
             FROM sgn.unigene WHERE sgn.unigene.unigene_build_id=?";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  my ($contigs_count) = $sth->fetchrow_array();
  return $contigs_count;
}

=head2 count_cds

  Usage: my @cds_reports=count_cds($dbh, $unigene_build_id);
  Des: count the cds_id for the diferent protein prediction methods asociated to unigene build
  Ret: An array of array_references with three scalars per array of these references.
        => @cds_reports=($cds_run1_array_ref, $cds_run2_array_ref...)
        => @cds_run1=($method1, $run_id, $cds_count)
  Args: $dbh (database connection object) and $unigene_build_id
  Side_Effects: none
  Example: my @cds_report=count_cds($dbh, $unigene_build_id);

=cut

sub count_cds {
  my $dbh=shift;
  my $unigene_build_id=shift;
  my @cds_reports; 

  my $query="SELECT sgn.cds.method, sgn.cds.run_id, COUNT(sgn.cds.cds_id) FROM sgn.cds
             RIGHT JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id 
             WHERE sgn.unigene.unigene_build_id=? GROUP BY sgn.cds.method, sgn.cds.run_id";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  while (my @cds_run=$sth->fetchrow_array() ) {
      push @cds_reports, \@cds_run;
  }
}

=head2 count_cds_with_domain_match

  Usage: my @domain_matches_report=count_domain_matches($dbh, $unigene_build_id);
  Desc: count the domain ids associated to a concrete protein prediction of an unigene build
  Ret: An array of array_references with four scalars per array of these references.
        => @domain_matches_reports=($domain_matches_run1_array_ref, $domain_matches_run2_array_ref...)
        => @domain_matches_run1=($method1, $cds_run_id, $domain_match_run_id, $domain_matches_count)
  Args: $dbh (database connection object) and $unigene_build_id
  Side_Effects: none
  Example:my @domain_matches_report=count_domain_matches($dbh, $unigene_build_id);

=cut

sub count_cds_with_domain_match {
  my $dbh=shift;
  my $unigene_build_id=shift;
  my @domain_match_report;

  my $query="SELECT sgn.cds.method, sgn.cds.run_id AS cds_run_id, COUNT(sgn.cds.cds_id) AS cds_n, 
             sgn.domain_match.run_id AS domain_match_run_id, COUNT(DISTINCT(sgn.domain_match.cds_id)) 
             FROM sgn.domain_match 
             RIGHT JOIN sgn.cds ON sgn.domain_match.cds_id=sgn.cds.cds_id
             RIGHT JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id
             WHERE sgn.unigene.unigene_build_id=? GROUP BY method, cds_run_id, domain_match_run_id";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  while (my @domain_match_run=$sth->fetchrow_array() ) {
      if ($domain_match_run[2]==0) {
          $domain_match_run[0]='None';
      }
      push @domain_match_report, \@domain_match_run;
  }
  return @domain_match_report;
}

=head2 count_unigene_with_blast_annotations

  Usage: my @blast_reports=count_unigene_with_blast_annotations($dbh, $unigene_build_id);
  Desc: count the unigene that have blast match for a concrete blast target database
  Ret: An array of array_references with four scalars per array of these references.
        => @blast_reports=($blast_run1_array_ref, $blast_run2_array_ref...)
        => @blast_run1=($blast_target_db, $last_updated, $last_updated_time, $unigene_count);
        $last_updated is the time in linux format when the blast results were loaded into the database, 
        this could be usefull to get/manipulate a specific unigene annotation but not to know when the blast was runned.
        To know it, you can use $last_upated_time, but remember that this time do not appears in the database 
        as a field.  
  Args: $dbh (database connection object) and $unigene_build_id
  Side_Effects: none
  Example: my @blast_reports=count_unigene_with_blast_annotations($dbh, $unigene_build_id);

=cut

sub count_unigene_with_blast_annotations {
  my $dbh=shift;
  my $unigene_build_id=shift;
  my @blast_reports;
  my $query="SELECT sgn.blast_annotations.blast_target_id, sgn.blast_annotations.last_updated, 
             COUNT(DISTINCT(sgn.blast_annotations.apply_id)) 
             FROM sgn.blast_annotations LEFT JOIN sgn.unigene ON sgn.unigene.unigene_id=sgn.blast_annotations.apply_id
             WHERE sgn.unigene.unigene_build_id = ? AND sgn.blast_annotations.n_hits > 0
             GROUP BY sgn.blast_annotations.blast_target_id, sgn.blast_annotations.last_updated";
  my $sth=$dbh->prepare($query);
  $sth->execute($unigene_build_id);
  while (my ($blast_target_db, $last_updated, $unigene_count) = $sth->fetchrow_array() ) {
      $blast_target_db =~ s/^1$/GenBank/;
      $blast_target_db =~ s/^2$/Arabidopsis/;
      $blast_target_db =~ s/^3$/Swissprot/;
      my $last_updated_time=convert_linux_time($last_updated);
      my @blast_run=($blast_target_db, $last_updated, $last_updated_time, $unigene_count);
      push @blast_reports, \@blast_run;
  }
  return @blast_reports;
} 
      

=head2 simple_columns_sizes

  Usage: my @column_sizes=simple_columns_sizes($aoa_simple_table_aref);
  Desc: count the maximum length for the elements in an aoa_table. 
  Ret: An array with the fields sizes for a aoa simple table (see below)
  Args: $aoa_table_aref, an array of arrays where the first array contains the array references of the rows.
  Side_Effects: none
  Example: my @column_sizes=simple_columns_sizes($aoa_table_aref);

=cut

sub simple_columns_sizes {
  my $table_aref=shift;
  my @table=@$table_aref;

  my $row_aref;
  my $field;
  my @field_size;

  foreach $row_aref (@table) {
     my @row=@$row_aref;
     my $field_posit=0;
     foreach $field (@row) {
         my $field_length; 
         if ($field =~ m/ARRAY/i ) {
             my $elements;
             my $element_size=0;
             my @field_array=@$field;
             foreach $elements (@field_array) {
		 my $element_length=length($elements);
		 if ($element_size < $element_length) {
		     $element_size = $element_length;
		     $field_length=$element_size;
		 }
	     }
	 } else {
	     $field_length=length($field);
	 }
         $field_length += 2;
         if ($field_size[$field_posit] < $field_length) {
             $field_size[$field_posit]= $field_length;
         }
         
         $field_posit++;
     }
 }
  return @field_size;
}

=head2 print_simple_aoa_table

  Usage: print_simple_aoa_table($aoa_table_aref);
  Desc: print a table with an array of arrays like argument.
  Ret: none
  Args: $aoa_table_aref, an array reference of array references. 
        $aoa_table_aref = \@aoa_table;
        @aoa_table= (\@header, \@row1, \@row2 ...);
        @row1= ($field1, $field2, $field3, \@field4 ...);
        @field4=($entry4_1, $entry4_2, $entry_4_3 ...);

        +----------+----------+----------+-----------+-----+
        | $header1 | $header2 | $header3 | $header4  | ... |==> @header
        +----------+----------+----------+-----------+-----+
        | $field1  | $field2  | $field3  | $entry4_1 | ... |==> @row1
        |          |          |          | $entry4_2 |     |
        |          |          |          | $entry4_3=|=====|==> @field4
        |          |          |          |    ...    |     |
        +----------+----------+----------+-----------+-----+
        | $field1  | $field2  | $field3  | $entry4_1 | ... |==> @row2
        +----------+----------+----------+-----------+-----+
        |   ...    |   ...    |   ...    |    ...    | ... |
        +----------+----------+----------+-----------+-----+ 

  Side_Effects: none
  Example: print_simple_aoa_table(\@unigene_build_global);

=cut

sub print_simple_aoa_table {
  my $aoa_table_aref=shift;
  my @aoa_table=@$aoa_table_aref;

  my @columns_sizes=simple_columns_sizes($aoa_table_aref);
  my (@row, $row_aref, $line);

  foreach $row_aref (@aoa_table) {
     @row=@$row_aref;
     my ($field, $line);
     my $n=0;
     foreach $line (@columns_sizes) {
        print "+-";
        my $m=0;
        while ($m < $line) {
           print "-";
           $m++;
        }
     } 
     my @empty_lines;
     print "+\n";
     my $t=scalar(@row);
     foreach $field (@row) {
        if ($field =~ m/array/i) {
            my @array_field=@$field;
            my $array_element;
            my $first_element=shift(@array_field);
            printf "| %-*s", $columns_sizes[$n], $first_element;
            
            foreach $array_element (@array_field) {
               my $a=0;
               my @empty_row; 
               while ($a < $n) {
                  my $empty_field=" ";
                  push @empty_row, $empty_field;
                  $a++;
               }
               push @empty_row, $array_element;
               $a++;
               while ($a < $t) {
                  my $empty_field=" ";
                  push @empty_row, $empty_field;
                  $a++;
               }
               push @empty_lines, \@empty_row;
            }
        } else {
          printf "| %-*s", $columns_sizes[$n], $field;
        }
        $n++;       
     }
     printf "|\n";
     if (@empty_lines) {
         my $emptyrow_aref;
         foreach $emptyrow_aref(@empty_lines) {
            my @empty_row=@$emptyrow_aref;
            my $o=0;
            my $empty_field;
            foreach $empty_field (@empty_row) { 
               printf "| %-*s", $columns_sizes[$o], $empty_field;
               $o++;
            }
            printf "|\n";
         }
     }
  }
  foreach $line (@columns_sizes) {
     print "+-";
     my $m=0;
     while ($m < $line) {
         print "-";
         $m++;
     }
  }
  print "+\n"; 
}

=head2 complex_columns_sizes
 
  Usage: my @column_sizes=simple_columns_sizes($aoa_complex_table_aref);
  Desc: count the maximum length for the elements in an aoa_complex_table. 
  Ret: An array with the fields sizes for a aoa simple table (see below)
  Args: $aoa_table_aref, an array of arrays where the first array contains the array references of the rows.
  Side_Effects: none
  Example: my @column_sizes=simple_columns_sizes($aoa_table_aref);

=cut

sub complex_columns_sizes {
  my $table_aref=shift;
  my @table=@$table_aref;
  my $row_aref;
  my $field;
  my @field_size;

  foreach $row_aref (@table) {
     my @row=@$row_aref;
     my $field_posit=0;
     my $opt_pos;
     foreach $field (@row) {
         my $field_length;
         if ($field =~ m/array/i ) {
             my $elements;
             my $element_size=0;
             my @field_array=@$field; 
             foreach $elements (@field_array) {
                 if ($elements =~ m/array/i ) {
		     my @elements_array=@$elements;
                     my $unit;
                     my ($unit_size, $opt_pos)=(0, 0);
                     my $opt_pos += $field_posit;
                     foreach $unit (@elements_array) {
			 my $unit_length=length($unit);
                         if ($unit_size < $unit_length) {
                             $unit_size = $unit_length;
                             $field_length=$unit_size;
			 }
                         if ($field_size[$opt_pos] < $field_length) {
                             $field_size[$opt_pos]= $field_length;
			 }
                         $opt_pos++; 
		     }
		 } else {
                     my $element_length=length($elements);
		     if ($element_size < $element_length) {
			 $element_size = $element_length;
			 $field_length=$element_size;
		     }
		 }
	     }
	 } else {
	     $field_length=length($field);
	 }
         $field_length += 2;
         if ($field_size[$field_posit] < $field_length) {
             $field_size[$field_posit]= $field_length;
         }
         $field_posit++;
	 
     }
 }
  return @field_size;
}

=head2 print_complex_aoa_table

  Usage: print_complex_aoa_table($aoa_table_aref);
  Desc: print a table with an array of arrays like argument.
  Ret: none
  Args: $aoa_table_aref, an array reference of array references. 
        $aoa_table_aref = \@aoa_table;
        @aoa_table= (\@header, \@row1, \@row2 ...);
        @row1= ($field1, \@field2 ...);
        @field2=(\@array_entry1, \@array_entry2 ...);
        @array_entry1=($entry1_1, $entry1_2, $entry1_3 ...);

        +----------+-----------+-----------+-----------+-----+
        | $header1 | $header2  | $header3  | $header4  | ... |==> @header
        +----------+-----------+-----------+-----------+-----+
        | $field1  | $entry1_1 | $entry1_2 | $entry1_3 | ... |==> @row1
        |          | $entry2_1 | $entry2_2 | $entry2_2 |     |
        |          |     +=====|=====+=====|=====+=====|=====|==> @array_entry1 ...
        +----------+-----------+-----------+-----------+-----+
        | $field1  | $entry1_1 | $entry1_2 | $entry1_3 | ... |==> @row2
        +----------+-----------+-----------+-----------+-----+
        |   ...    |   ...     |   ...     |    ...    | ... |
        +----------+-----------+-----------+-----------+-----+ 

  Side_Effects: none
  Example: print_complex_aoa_table(\@unigene_cds_annot_global);

=cut

sub print_complex_aoa_table {
  my $aoa_table_aref=shift;
  my @aoa_table=@$aoa_table_aref;

  my @columns_sizes=complex_columns_sizes($aoa_table_aref);
  my (@row, $field_array, $row_aref, $line);

  foreach $row_aref (@aoa_table) {
     @row=@$row_aref;
     my ($field, $line);
     my $n=0;
     foreach $line (@columns_sizes) {
        print "+-";
        my $m=0;
        while ($m < $line) {
           print "-";
           $m++;
        }
     } 
     my @empty_lines;
     print "+\n";
     my $t=scalar(@row);
     foreach $field_array (@row) {
        if ($field_array =~ m/array/i) {
            my @array_field=@$field_array;
            my $field;
            foreach $field (@array_field) {
		if ($field =~ m/array/i) {
		    my $array_element;
		    my $first_element_aref=shift(@array_field);
                    my @first_elements=@$first_element_aref;
                    my $first_element;
                    my $s=$n;
                    foreach $first_element (@first_elements) {
           	       printf "| %-*s", $columns_sizes[$s], $first_element;
                       $s++;
		    }
		    foreach $array_element (@array_field) {
			my $a=0;
			my @empty_row; 
			while ($a < $n) {
			    my $empty_field=" ";
			    push @empty_row, $empty_field;
			    $a++;
			}
			if ($array_element =~ m/array/i) {
			    my @array_element=@$array_element;
                            my $unit;
                            foreach $unit (@array_element) {
				push @empty_row, $unit;
            			$a++;
			    }
			} else {
			    push @empty_row, $array_element;
            		    $a++;
			}
			while ($a < $t) {
			    my $empty_field=" ";
			    push @empty_row, $empty_field;
			    $a++;
			}
			push @empty_lines, \@empty_row;
		    }
		} else {
		    printf "| %-*s", $columns_sizes[$n], $field;
                    $n++;
		}
	    }
	} else {
	    printf "| %-*s", $columns_sizes[$n], $field_array;
            $n++
	}    
     }
     printf "|\n";
     if (@empty_lines) {
         my $emptyrow_aref;
         foreach $emptyrow_aref(@empty_lines) {
            my @empty_row=@$emptyrow_aref;
            my $o=0;
            my $empty_field;
            foreach $empty_field (@empty_row) { 
               printf "| %-*s", $columns_sizes[$o], $empty_field;
               $o++;
            }
            printf "|\n";
         }
     }
  }
  foreach $line (@columns_sizes) {
     print "+-";
     my $m=0;
     while ($m < $line) {
         print "-";
         $m++;
     }
  }
  print "+\n"; 
}

=head2 convert_linux_time

  Usage: my $time=convert_linux_time($linux_time);
  Desc: convert the linux time to the format XXXX-YY-ZZ where XXXX=> year, YY=> month and ZZ=> day. 
  Ret: A scalar, $time (in human)
  Args: $linux_time
  Side_Effects: none
  Example: my $time=convert_linux_time($linux_time);

=cut

sub convert_linux_time {
  my $linux_time=shift;
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($linux_time);
 
  my $c_year=$year + 1900;
  my @months=('01', '02', '03', '04', '05', '06', '07', '08', '09', '10', '11', '12');
  my $c_month=$months[$mon];
  my $c_day;
  if ($mday < 10) {
      $c_day='0'.$mday;
  } else {
      $c_day=$mday;
  }
  my $time=$c_year."-".$c_month."-".$c_day;
  return $time;
}
