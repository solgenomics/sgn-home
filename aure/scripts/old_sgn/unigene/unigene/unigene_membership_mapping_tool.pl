#!/usr/bin/perl

=head1 NAME

 unigene_membership_mapping_tool.pl
 A script to map sequences ID's in the SGN databases.(version.1.1.).

=cut

=head1 SYPNOSIS

unigene_membership_mapping_tool.pl -H <dbhost> -D <dbname> -b <basename> [-o <organism_name> || -u <unigene_build>]
 -L -C -T -E -G -U -P [-m <method>] [-j <join_kind>] [-f <flags>] [-s <status>] [-h]
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>     for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>         sandbox or cxgn etc (mandatory)

=item -b

B<output_basename>       basename for the output (mandatory)
      
=item -o

B<organism_name>         organism name (co-mandatory with -u <unigene_build>)
      
=item -u 

B<unigene_build>         unigene build (co-mandatory with -o <organism_name>)
      
=item -L

B<get_library_shortname> get the library_shortname in the mapping

=item -C 

B<get clone_id>          get the clone_id (as SGN-C <clone_id>) in the mapping
      
=item -T 

B<get trace_id>          get the trace_id (as SGN-T <read_id>) in the mapping

=item -E 

B<get est_id>            get the est_id (as SGN-E <est_id>) in the mapping
      
=item -U 

B<get unigene_id>        get the unigene_id (as SGN-U <unigene_id>) in the mapping

=item -P

B<get protein_id>        get the protein_id (or cds, as SGN-P <cds_id>) in the mapping     

=item -m

B<protein pred method>   specify the protein prediction method for the SGN-P mapping (always take the last run_id). 

=item -j

B<specify join>          specify the join ('R' for RIGHT JOIN, 'L' for LEFT JOIN) (JOIN by default)

=item -f

B<specify flags>         specify flags, like 0 (dependient of -E)

=item -s

B<specify status>        specify status, like 0 (dependient of -E) 

=item -h 

B<help>                  show the help

=back

=cut

=head1 DESCRIPTION

This script get SGN-ID's and map with other SGN-ID's. The file output has a tab format. You can choose more than one ID to map, for exmaple you can choose GenBank accessions (-G), SGN-E <est_id> and SGN-U <unigene_id> for an organism. The use of -U with -o <organism_name> produce one column per build (from oldest to newest). The use of -P produce more than one column order by run_id and method.

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

unigene_membership_mapping_tool.pl


=cut

use strict;
use Getopt::Std;
use CXGN::DB::InsertDBH;
use File::Basename;

our ($opt_H, $opt_D, $opt_b, $opt_o, $opt_u, $opt_L, $opt_C, $opt_T, $opt_E, $opt_G, $opt_U, $opt_P, $opt_m, $opt_j, $opt_f, $opt_s, $opt_h);
getopts("H:D:b:o:u:LCTEGUPm:j:f:s:h");
if (!$opt_H && !$opt_D && !$opt_b && !$opt_o && !$opt_u && !$opt_L && !$opt_C && !$opt_T && !$opt_E && !$opt_G && !$opt_U && !$opt_P && !$opt_m && !$opt_j && !$opt_f && !$opt_s && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

if ($opt_h) {
    help();
}

unless ($opt_H && $opt_D) {
    die ("Required argument -H <dbhost> AND -D <dbname> were not supplied.\n");
}
unless ($opt_b) {
    die ("Required argument -b <basename> was not supplied.\n");
}
unless ($opt_o || $opt_u) {
    die ("Required argument -o <organism_name> OR -u <unigene_build> (not both) was not supplied.\n");
}
if ($opt_o && $opt_u) {
    die ("Incompatible arguments -o <organism_name> and -u <unigene_build> were supplied. Please, use only one of them.\n");
}
if (!$opt_L && !$opt_C && !$opt_T && !$opt_E && !$opt_G && !$opt_U && !$opt_P) {
    print "There are not any mapping option. Please rerun the script with at least one of them.\n\n";
    help();
}
if ($opt_P && ($opt_L || $opt_C || $opt_T || $opt_E || $opt_G)) {
    die "\nIncompatible arguments were chosen.\nPlease, use help (-h) for more information\n\n";
} elsif ($opt_U && $opt_G && ($opt_L || $opt_C || $opt_T)) {
    die "\nIncompatible arguments were chosen.\nPlease, use help (-h) for more information\n\n";
}

if (!$opt_P && $opt_m) {
    die "Argument -m <protein_prediction_method> only can be used with -P <get_protein_id> argument.\n";
}

if ($opt_u) {
    unless ($opt_u =~ m/\d+/) {
        die ("Argument -u must be an integer\n");
    }
}

if ($opt_j) {
    unless ($opt_j =~ m/^[R|L]/i) {
        die "Argument -j only can have 'R' or 'L' values\n";
    }
}

if ($opt_E) {
    if ($opt_f) {
	unless ($opt_f =~ m/^[\d+]$/i) {
            die "Argument -f <flags> only can have integer values (numbers)\n";
	}
    }
    if ($opt_s) {  
        unless ($opt_s =~ m/^[\d+]$/i) {
            die "Argument -s <status> only can have integer values (numbers)\n";
	}
    }
}

our $dbh = CXGN::DB::InsertDBH::connect({	dbname=>$opt_D ,
						dbhost=>$opt_H ,
					  	dbargs=>{ RaiseError => 1, PrintError => 1 }
					})
	or die "Failed to connect to database ($DBI::errstr)";

my $global_results_aref;
if ($opt_o) {
    my $organism="%".$opt_o."%";
    my $check_organism="SELECT organism_id FROM sgn.organism WHERE organism_name ILIKE ?";
    my $sth=$dbh->prepare($check_organism);
    $sth->execute($organism);
    my ($organism_id) = $sth->fetchrow_array();

    if (!$organism_id) {
        die ("The organism -o $opt_o is not in the $opt_D database.\n");
    } else {
        print "Processing search for $opt_o (organism_id: $organism_id) ...\n";
        $global_results_aref=organism_query_results_temp_table($dbh);
    }
}
if ($opt_u) {
    my $check_unigene="SELECT build_date FROM sgn.unigene_build WHERE unigene_build_id = ?";
    my $sth=$dbh->prepare($check_unigene);
    $sth->execute($opt_u);
    my ($build_date) = $sth->fetchrow_array();

    if (!$build_date) {
        die ("The build id -u $opt_u is not in the $opt_D database.\n");
    }
        print "Processing search for $opt_u (unigene_build_date: $build_date)\n";
        $global_results_aref=unigene_build_query_results_temp_table($dbh);
}

my $output_filename=$opt_b;
if ($opt_L) {
    $output_filename .= "_L";
}
if ($opt_C) {
    $output_filename .= "_C";
}
if ($opt_T) {
    $output_filename .= "_T";
}
if ($opt_E) {
    $output_filename .= "_E";
}
if ($opt_G) {
    $output_filename .= "_G";
}
if ($opt_U) {
    $output_filename .= "_U";
}
if ($opt_P) {
    $output_filename .= "_P";
}
$output_filename .= ".tab";

my @global_results=@$global_results_aref;
my $n_entries=scalar(@global_results);
print "Elements processed: $n_entries\n";

if ($n_entries > 0) {
    my $complete_filename;
    open my $filehandle, '>', "$output_filename" || die "Cannot create the file $output_filename\n";
    chmod 0666, "$output_filename";

    my $row_result;
    foreach $row_result (@global_results) {
       my @results=@$row_result;
       my $single_result;
       my $m=scalar(@results);
       my $n=0;
       foreach $single_result (@results) {
          print $filehandle "$single_result";
          $n++;
          if ($n < $m) {
              print $filehandle "\t";
          }
       }
       print $filehandle "\n";
    }
} else {
    print "The query results has $n_entries results.\n";
}
 
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
      This script map different SGN identifiers for a concrete organism or unigene build. The available options are:
      
        * Using organism_name in the search:
             -L (sgn.library_shortname) (incompatible with -P, used with -U give the libraries for unigene members)
             -C (SGN-C<clone_id>) (incompatible with -P, used with -U give the libraries for unigene members)
             -T (SGN-T<read_id>) (incompatible with -P, used with -U give the libraries for unigene members)
             -E (SGN-E<est_id>) (used with -U give the unigene members)
             -G (GenBank accessions) (used with -U give the GenBank accession for the unigene members, 
                 in that case it is incompatible with -L, -C or -T)
             -U (SGN-U<unigene_id>)
             -P (SGN-P<cds_id>) (incompatible with -L, -C, -T, -E and -G)

        * Using the unigene build in the search:            
             -U (SGN-U<unigene_id>) (mandatory with this option)
             -E (SGN-E<est_id>) (used with -U give the unigene members)
             -G (GenBank accessions) (used with -U give the GenBank accession for the unigene members)
             -P (SGN-P<cds_id>) (incompatible with -E and -G)  

    Usage: 
      unigene_membership_mapping_tool.pl -H <dbhost> -D <dbname> -b <basename> [-o <organism_name> OR -u <unigene_build>] 
      [-L] [-C] [-T] [-E] [-G] [-U] [-P] [-m <method>] [-j <join_kind>] [-f <flags>] [-s <status>] [-h]
    
    Flags:
      -H database hostname      for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name          sandbox or cxgn etc (mandatory)
      -b basename               basename for the output_file (it will be <basename>+map+<options>+.tab) (mandatory)
      -o organism_name          organism_name for the search (co-mandatory with -u)
      -u unigene_build_id       unigene_build_id for the search (co-mandatory with -o)
      -L get library            get library_shortname for the mapping
      -C get clone id           get clone id as SGN-C<clone_id> for the mapping
      -T get trace id           get trace id as SGN-T<read_id> for the mapping
      -E get sequence id        get sequence id as SGN-E<est_id> for the mapping
      -G get genbank accession  get GenBank accession for the mapping
      -U get unigene id         get unigene id as SGN-U<unigene_id> for the mapping
      -P get protein id         get protein id as SGN-P<cds_id>, method and run_id if do not specify the method 
                                with -m, and only SGN-P<cds_id> if it is specified (get the id for the last run_id)
      -m prediction method      specify a protein prediction method for the search (always it get the last run_id)
      -j specify join           specify the join type ('R'->right join or 'L'-> left join) (simple 'join' by default)
      -f specify flag           specify flags for sgn.est table 
      -s specify status         specify status for sgn.est table
      -h this help

EOF
exit (1);
}

=head2 organism_query_results_temp_table

  Usage: my $global_results_aref=organism_query_results_temp_table($dbh);
  Desc: select into a temp table the diferent arguments for organism name argument
  Ret: A array reference of arrays references (array with row and each row as another array), $global_results_aref;
  Args: $dbh, the database conection object and $basename, a basename for the temp table
  Side_Effects: none
  Example: my $global_results_aref=organism_query_results_temp_table($dbh);

=cut

sub organism_query_results_temp_table {
    my $dbh=shift;
    my @args;

    my $query="SELECT ";
    if ($opt_U) {
	$query .= "sgn.unigene.unigene_id AS SGN_U, ";
        push @args, 'SGN_U';
    }
    if ($opt_L) {
        $query .= "sgn.library.library_shortname AS lib, ";
        push @args, 'lib';
    }
    if ($opt_C) {
	$query .= "sgn.clone.clone_id AS SGN_C, ";
        push @args, 'SGN_C';
    }
    if ($opt_T) {
	$query .= "sgn.seqread.read_id AS SGN_T, ";
        push @args, 'SGN_T';
    }
    if ($opt_E) {
        $query .= "sgn.est.est_id AS SGN_E, ";
        push @args, 'SGN_E';
    }
    if ($opt_G) {
	$query .= "public.dbxref.accession AS gb, ";
        push @args, 'gb';
    }
    if ($opt_P && !$opt_m) {
	$query .= "sgn.cds.cds_id as SGN_P, sgn.cds.method as ppm, sgn.cds.run_id as ppr, ";
        push @args, 'SGN_P', 'ppm', 'ppr';
    } elsif ($opt_P && $opt_m) {
        $query .= "sgn.cds.cds_id as SGN_P, ";
        push @args, 'SGN_P';
    }

    $query =~ s/,\s$/ /gi;     
    if (!$opt_U) {
        $query .= "FROM sgn.organism 
                   JOIN sgn.library ON sgn.organism.organism_id=sgn.library.organism_id 
                   JOIN sgn.clone ON sgn.library.library_id=sgn.clone.library_id 
                   JOIN sgn.seqread ON sgn.clone.clone_id=sgn.seqread.clone_id
                   JOIN sgn.est ON sgn.seqread.read_id=sgn.est.read_id
                   JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id ";
	if ($opt_G) {
	    $query .= "JOIN sgn.est_dbxref ON sgn.est.est_id=sgn.est_dbxref.est_id
                       JOIN public.dbxref ON sgn.est_dbxref.dbxref_id=public.dbxref.dbxref_id ";
	}
    } else {
        if ($opt_P) {
            $query .= "FROM sgn.organism
                   JOIN sgn.group_linkage ON sgn.organism.organism_id=sgn.group_linkage.member_id
                   JOIN sgn.unigene_build ON sgn.group_linkage.group_id=sgn.unigene_build.organism_group_id
                   JOIN sgn.unigene ON sgn.unigene_build.unigene_build_id=sgn.unigene.unigene_build_id
                   JOIN sgn.cds ON sgn.unigene.unigene_id=sgn.cds.unigene_id ";
        } elsif ($opt_G) {
            $query .= "FROM sgn.organism
                       JOIN sgn.group_linkage ON sgn.organism.organism_id=sgn.group_linkage.member_id
                       JOIN sgn.unigene_build ON sgn.group_linkage.group_id=sgn.unigene_build.organism_group_id
                       JOIN sgn.unigene ON sgn.unigene_build.unigene_build_id=sgn.unigene.unigene_build_id
                       JOIN sgn.unigene_member ON sgn.unigene.unigene_id=sgn.unigene_member.unigene_id
                       JOIN sgn.est ON sgn.unigene_member.est_id=sgn.est.est_id
                       JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id
                       JOIN sgn.est_dbxref ON sgn.est.est_id=sgn.est_dbxref.est_id
                       JOIN public.dbxref ON sgn.est_dbxref.dbxref_id=public.dbxref.dbxref_id ";
        } else {
            $query .= "FROM sgn.organism
                       JOIN sgn.group_linkage ON sgn.organism.organism_id=sgn.group_linkage.member_id
                       JOIN sgn.unigene_build ON sgn.group_linkage.group_id=sgn.unigene_build.organism_group_id
                       JOIN sgn.unigene ON sgn.unigene_build.unigene_build_id=sgn.unigene.unigene_build_id
                       JOIN sgn.unigene_member ON sgn.unigene.unigene_id=sgn.unigene_member.unigene_id
                       JOIN sgn.est ON sgn.unigene_member.est_id=sgn.est.est_id
                       JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id
                       JOIN sgn.seqread ON sgn.est.read_id=sgn.seqread.read_id
                       JOIN sgn.clone ON sgn.seqread.clone_id=sgn.clone.clone_id
                       JOIN sgn.library ON sgn.clone.library_id=sgn.library.library_id ";
        }
    }
    my $organism="'%".$opt_o."%'";
    $query .= "WHERE organism_name ILIKE $organism ";
    if ($opt_f) {
        $query .= "AND sgn.est.flags=$opt_f ";
    } elsif ($opt_f ==0) {
        $query .= "AND sgn.est.flags=0 ";
    }
    if ($opt_s) {
        $query .= "AND sgn.est.status=$opt_s ";
    } elsif ($opt_s ==0) {
	$query .= "AND sgn.est.status=0 ";
    }
    if ($opt_m) {
        my $method="'%".$opt_m."%'";
        my @run_ids;
        my $check1="SELECT DISTINCT(sgn.cds.run_id) FROM sgn.cds 
                   JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id
                   JOIN sgn.unigene_build ON sgn.unigene.unigene_build_id=sgn.unigene_build.unigene_build_id
                   JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id
                   JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id
                   WHERE organism_name ILIKE $organism AND sgn.cds.method ILIKE $method ORDER BY sgn.cds.run_id ASC";
        my $sth=$dbh->prepare($check1);
        $sth->execute();
        while (my ($run_id) = $sth->fetchrow_array() ) {
	    push @run_ids, $run_id;
	}
        my $last_run_id=pop(@run_ids);
        $query .="AND sgn.cds.method ILIKE $method ";
        if ($last_run_id) {
            $query .= "AND sgn.run_id = $last_run_id ";
        } 
    }
    if ($opt_L && $opt_U && $opt_E) {
	$query .= "ORDER BY sgn.library.library_shorname, sgn.unigene.unigene_id, sgn.est.est_id";
    } elsif (!$opt_L && $opt_U && $opt_E && !$opt_P) {
        $query .= "ORDER BY sgn.unigene.unigene_id, sgn.est.est_id ";
    } elsif (!$opt_L && $opt_U && !$opt_E && !$opt_P) { 
	$query .= "ORDER BY sgn.unigene.unigene_id ";
    } elsif (!$opt_L && !$opt_U && $opt_E) {
        $query .= "ORDER BY sgn.est.est_id ";
    } elsif ($opt_U && $opt_P && !$opt_m) {
        $query .= "ORDER BY sgn.cds.method, sgn.unigene.unigene_id, sgn.cds.cds_id ";
    }

    if ($opt_j =~ m/^L/i) {
        $query =~ s/JOIN/LEFT JOIN/gi;
    } elsif ($opt_j =~ m/^R/i) {
        $query =~ s/JOIN/RIGHT JOIN/gi;
    }

    print "QUERY TEST:\n$query\n";
    my @global;
    my $sth1=$dbh->prepare($query);
    $sth1->execute();
    while (my @results=$sth1->fetchrow_array() ) {
       my ($args_f, $result_f);
       my $args_c=scalar(@args);
       print "----------------\nGet data: @results\nArgument array: @args ($args_c)\n";
       my $ap=0;
       foreach $args_f (@args) {
	  print "arguments list inside foreach: $args_f\t";
          if ($args_f =~ m/SGN/i) {
              my $subs=$args_f;
              $subs =~ s/_/-/gi;
              my $var=$results[$ap];
              $result_f=$subs.$var;
              $results[$ap]=$result_f;
              print "check ap: $ap, check result_f:$result_f, check_subs: $subs\n";
          }
          $ap++;
       }
       print "Data after change: @results\nArgument array: @args\n--------------------\n";
       push @global, \@results;
    }   
    return \@global;
}

=head2 unigene_build_query_results_temp_table

  Usage: my $temp_table=unigene_build_query_results_temp_table($dbh, $basename);
  Desc: select into a temp table the diferent arguments for unigene build argument
  Ret: A scalar, $temp_table the name of a temp table that stores the results of the query
  Args: $dbh, the database conection object and $basename, a basename for the temp table
  Side_Effects: none
  Example: my $temp_table=unigene_build_query_results_temp_table($dbh, $basename);

=cut

sub unigene_build_query_results_temp_table {
    my $dbh=shift;
    my @args;

    my $query="SELECT ";
    if ($opt_L) {
        $query .= "sgn.library.library_shortname AS lib, ";
        push @args, 'lib';
    }
    if ($opt_C) {
	$query .= "sgn.clone.clone_id AS SGN_C, ";
        push @args, 'SGN_C';
    }
    if ($opt_T) {
	$query .= "sgn.seqread.read_id AS SGN_T, ";
        push @args, 'SGN_T';
    }
    if ($opt_E) {
        $query .= "sgn.est.est_id AS SGN_E, ";
        push @args, 'SGN_E';
    }
    if ($opt_G) {
	$query .= "public.dbxref.accession AS gb, ";
        push @args, 'gb';
    }
    if ($opt_U) {
	$query .= "sgn.unigene.unigene_id AS SGN_U, ";
        push @args, 'SGN_U';
    }
    if ($opt_P && !$opt_m) {
	$query .= "sgn.cds.cds_id as SGN_P, sgn.cds.method as ppm, sgn.cds.run_id as ppr, ";
        push @args, 'SGN_P', 'ppm', 'ppr';
    } elsif ($opt_P && $opt_m) {
        $query .= "sgn.cds.cds_id as SGN_P, ";
        push @args, 'SGN_P';
    }
    $query =~ s/,\s$/ /gi;
    if ($opt_P) {
        $query .= "FROM sgn.unigene_build
                   JOIN sgn.unigene ON sgn.unigene_build.unigene_build_id=sgn.unigene.unigene_build_id
                   JOIN sgn.cds ON sgn.unigene.unigene_id=sgn.cds.unigene_id ";
    } elsif ($opt_G) {
        $query .= "FROM sgn.unigene_build
                  JOIN sgn.unigene ON sgn.unigene_build.unigene_build_id=sgn.unigene.unigene_build_id
                  JOIN sgn.unigene_member ON sgn.unigene.unigene_id=sgn.unigene_member.unigene_id
                  JOIN sgn.est ON sgn.unigene_member.est_id=sgn.est.est_id
                  JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id
                  JOIN sgn.est_dbxref ON sgn.est.est_id=sgn.est_dbxref.est_id
                  JOIN public.dbxref ON sgn.est_dbxref.dbxref_id=public.dbxref.dbxref_id ";
    } else {
        $query .= "FROM sgn.unigene_build
                   JOIN sgn.unigene ON sgn.unigene_build.unigene_build_id=sgn.unigene.unigene_build_id
                   JOIN sgn.unigene_member ON sgn.unigene.unigene_id=sgn.unigene_member.unigene_id
                   JOIN sgn.est ON sgn.unigene_member.est_id=sgn.est.est_id
                   JOIN sgn.qc_report ON sgn.est.est_id=sgn.qc_report.est_id
                   JOIN sgn.seqread ON sgn.est.read_id=sgn.seqread.read_id
                   JOIN sgn.clone ON sgn.seqread.clone_id=sgn.clone.clone_id
                   JOIN sgn.library ON sgn.clone.library_id=sgn.library.library_id ";
    }
    $query .= "WHERE sgn.unigene_build.unigene_build_id=$opt_u ";
    if ($opt_m) {
        my $method="'%".$opt_m."%'";
        my @run_ids;
        my $check1="SELECT DISTINCT(sgn.cds.run_id) FROM sgn.cds 
                   JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id
                   JOIN sgn.unigene_build ON sgn.unigene.unigene_build_id=sgn.unigene_build.unigene_build_id
                   WHERE sgn.unigene_build.unigene_build_id = $opt_u AND sgn.cds.method ILIKE $method 
                   ORDER BY sgn.cds.run_id ASC";
        my $sth=$dbh->prepare($check1);
        $sth->execute();
        while (my ($run_id) = $sth->fetchrow_array() ) {
	    push @run_ids, $run_id;
	}
        my $last_run_id=pop(@run_ids);
        $query .="AND sgn.cds.method ILIKE $method ";
        if ($last_run_id) {
            $query .= "AND sgn.run_id = $last_run_id ";
        } 
    }
    if ($opt_L && $opt_U && $opt_E) {
	$query .= "ORDER BY sgn.library.library_shorname, sgn.unigene.unigene_id, sgn.est.est_id";
    } elsif (!$opt_L && $opt_U && $opt_E && !$opt_P) {
        $query .= "ORDER BY sgn.unigene.unigene_id, sgn.est.est_id ";
    } elsif (!$opt_L && $opt_U && !$opt_E && !$opt_P) { 
	$query .= "ORDER BY sgn.unigene.unigene_id ";
    } elsif (!$opt_L && !$opt_U && $opt_E) {
        $query .= "ORDER BY sgn.est.est_id ";
    } elsif ($opt_U && $opt_P && !$opt_m) {
        $query .= "ORDER BY sgn.cds.method, sgn.unigene.unigene_id, sgn.cds.cds_id ";
    }
    
    if ($opt_j =~ m/^L/i) {
        $query =~ s/JOIN/LEFT JOIN/gi;
    } elsif ($opt_j =~ m/^R/i) {
        $query =~ s/JOIN/RIGHT JOIN/gi;
    }
 
    print "QUERY TEST:\n$query\n";
    my @global;
    my $sth1=$dbh->prepare($query);
    $sth1->execute();
    while (my @results=$sth1->fetchrow_array() ) {
       my ($args_f, $result_f);
       my $ap=0;
       foreach $args_f (@args) {
          if ($args_f =~ m/SGN/i) {
              my $subs=$args_f;
              $subs =~ s/_/-/gi;
              my $var=$results[$ap];
              $result_f=$subs.$var;
              $results[$ap]=$result_f;
          }
          $ap++;
       }
       push @global, \@results;
    }   
    return \@global;
}
