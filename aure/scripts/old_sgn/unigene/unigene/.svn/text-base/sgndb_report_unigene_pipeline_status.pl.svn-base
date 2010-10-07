#!/usr/bin/perl
=head1 NAME

 sgndb_report_unigene_pipeline_status.pl
 A script to get different data for unigene build (version.1.0.).

=cut

=head1 SYPNOSIS

 sgndb_report_unigene_pipeline_status.pl [-h] -H <dbhost> -D <dbname> [-A | -o <organism_name> | -b <build_id>] [-s <status>]
  
    
=head2 I<Flags:>

=over

=item -H

B<database_host>                database host (mandatory if you want check the relations in the database)

=item -D

B<database_name>                database name (mandatory if you want check the relations in the database)

=item -A 

B<get ALL>                      get all the unigene builds (co-mandatory with -o or -b)

=item -o 

B<organism name>                get the organism name or (organism names, separated by commas)(co-mandatory with -A or -b)
      
=item -b 

B<build id>                     get the unigene build id (or build ids, separated by commas)(co-mandatory with -A or -o)
      
=item -s 

B<status>                       filter, get only the unigenes with this status (permited values: current, previous and deprecated)

=item -h

B<help>                         print the help  

=back

=cut

=head1 DESCRIPTION

     This script get unigene data and print this information in the screen (you can store in a file if you want using > ). This script
      give four tables:

	  UNIGENE_ASSEMBLY_STATS, with nine columns:
	     -f1: build_id, id for the unigene build in the SGN database.
	     -f2: organism, organism name if the sequences used in the unigene assembly.
	     -f3: build_nr, number of unigene build for these dataset (2 means that its the second unigene build done for the dataset)
	     -f4: build_date, date when was done the assembly and load in the database
	     -f5: status of the unigene build, could be current, previous and deprecated
	     -f6: unigene_nr, number of unigene that have the unigene build
	     -f7: contig_nr, number of consensus sequences that have the unigene build (unigenes composed by more than one sequence)
	     -f8: member_nr, number of sequences used in the assembly
	     -f9: redundancy, percentage of the unigene sequence that are contigs [(contig_nr/unigene_nr)*100]

	  UNIGENE_IDENTIFIER_RANGE, it is a table with three columns, the first for the unigene build id, and the second and the third
	     for the first and last SGN unigene accession respectively.

	  UNIGENE_BLAST_ANNOTATION_STATS, is a table with the number of unigene annotated using blast homology searches.

	  UNIGENE_PROTEIN_PREDICTION_AND_ANNOTATION_STATS, is a table with number of cds (and protein) predicted for the unigene
	     build dataset and the number of different domains found in these sequences.
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

sgndb_report_unigene_pipeline_status.pl


=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;
use Math::BigInt;
use CXGN::DB::InsertDBH;

our ($opt_H, $opt_D, $opt_A, $opt_o, $opt_b, $opt_s, $opt_h);
getopts("H:D:Ao:b:s:h");
if (!$opt_H && !$opt_D && !$opt_A && !$opt_o && !$opt_b && !$opt_s && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
} elsif ($opt_h) {
    help();
}

my ($dbh, $unigene_build_id_list)=validate_input();
my (@data_assembly1, @data_assembly2, @data_blast_annotation, @data_prot_annotation);
foreach my $unigene_build_id (@$unigene_build_id_list) {
    print "Getting unigene assembly stats for unigene_build_id: $unigene_build_id ...\n";
    my ($data_assembly_href, $data_assembly_ext_href)=get_unigene_assembly_stats($dbh, $unigene_build_id);
    print "\t\t\t\t\t\t\t\t... done\n\n";
    push @data_assembly1, $data_assembly_href;
    push @data_assembly2, $data_assembly_ext_href;
    print "Getting unigene annotations stats for unigene_build_id: $unigene_build_id ...\n";
    my ($data_annotation_blast_href, $data_annotation_prot_href)=get_unigene_annotation_stats($dbh, $unigene_build_id);
    print "\t\t\t\t\t\t\t\t... done\n\n";
    push @data_blast_annotation, $data_annotation_blast_href;
    push @data_prot_annotation, $data_annotation_prot_href;
}
print "\n";
if (scalar(@data_assembly1) > 0) {
    print_table(\@data_assembly1, 'unigene_assembly_stats');
    print "\n\n";
}
if (scalar(@data_assembly2) > 0) {
    print_table(\@data_assembly2, 'unigene_identifiers_range');
    print "\n\n";
}
if (scalar(@data_blast_annotation) > 0) {
    print_table(\@data_blast_annotation, 'unigene_blast_annotations_stats');
    print "\n\n";
}
if (scalar(@data_prot_annotation) > 0) {
    print_table(\@data_prot_annotation, 'unigene_protein_prediction_and_annotation_stats');
    print "\n\n";
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
      This script get unigene data and print this information in the screen (you can store in a file if you want using > ). This script
      give four tables:

	  UNIGENE_ASSEMBLY_STATS, with nine columns:
	     -f1: build_id, id for the unigene build in the SGN database.
	     -f2: organism, organism name if the sequences used in the unigene assembly.
	     -f3: build_nr, number of unigene build for these dataset (2 means that its the second unigene build done for the dataset)
	     -f4: build_date, date when was done the assembly and load in the database
	     -f5: status of the unigene build, could be current, previous and deprecated
	     -f6: unigene_nr, number of unigene that have the unigene build
	     -f7: contig_nr, number of consensus sequences that have the unigene build (unigenes composed by more than one sequence)
	     -f8: member_nr, number of sequences used in the assembly
	     -f9: redundancy, percentage of the unigene sequence that are contigs [(contig_nr/unigene_nr)*100]

	  UNIGENE_IDENTIFIER_RANGE, it is a table with three columns, the first for the unigene build id, and the second and the third
	     for the first and last SGN unigene accession respectively.

	  UNIGENE_BLAST_ANNOTATION_STATS, is a table with the number of unigene annotated using blast homology searches.

	  UNIGENE_PROTEIN_PREDICTION_AND_ANNOTATION_STATS, is a table with number of cds (and protein) predicted for the unigene
	     build dataset and the number of different domains found in these sequences.
    
    Usage: 
      sgndb_report_unigene_pipeline_status.pl [-h] -D <dbname> -H <dbhost> [-A | -o <organism_name> | -b <build_id>] [-s <status>]

    Example:
      sgndb_report_unigene_pipeline_status.pl -D sandbox -H localhost -o 'Nicotiana' -s 'current' > Nicotiana_species_unigene_report.txt
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory for check domains in the database)(mandatory)
      -D database name        sandbox or cxgn etc (mandatory for check domains in the database)(mandatory)
      -A get ALL              get all the unigene builds (co-mandatory with -o or -b)
      -o organism name        get the organism name or (organism names, separated by commas)(co-mandatory with -A or -b)
      -b build id             get the unigene build id (or build ids, separated by commas)(co-mandatory with -A or -o)
      -s status               filter, get only the unigenes with this status (permited values: current, previous and deprecated)
      -h this help

EOF
exit (1);
}

=head2 validate_input

  Usage: my ($dbh, $unigene_build_aref)=validate_input();
  Desc: Validate the input parameters supplied by getopt
  Ret: $dbh, a database conection object and $unigene_build_aref, an array reference with the unigene_build_ids
  Args: none
  Side_Effects: die if there are something wrong
  Example: my ($dbh, $unigene_build_aref)=valiadate_input();

=cut

sub validate_input {
    my $dbh;
    my @unigene_build_ids;

    if (!$opt_H || !$opt_D) {
	die "Sorry, the database conection paramaters (-H <dbhost> and -D <dbname>) were not supplied.\n\n";
    } else {
	$dbh = CXGN::DB::InsertDBH::connect({      dbname=>$opt_D,
						   dbhost=>$opt_H,
						   dbargs=>{ RaiseError => 1, PrintError => 1 }
                                        })
	or die "Failed to connect to database ($DBI::errstr)\n\n";
    }
    print STDERR "\n\n";
    if (!$opt_A && !$opt_o && !$opt_b) {
	print STDERR "Sorry, any dataset information was supplied.";
	die " Please use: -A (get_ALL_unigene_builds), -o <organism_name> or -b <unigene_build_id> options\n\n";
    } elsif ($opt_A && $opt_o) {
	print STDERR "Sorry, the parameter -A (get_ALL_unigene_builds) is incompatible with -o <organism_name>.";
	die " Please use: -A (get_ALL_unigene_builds) OR -o <organism_name>, don't use both.\n\n";
    } elsif ($opt_A && $opt_b) {
	print STDERR "Sorry, the parameter -A (get_ALL_unigene_builds) is incompatible with -b <unigene_build_ids>.";
	die " Please use: -A (get_ALL_unigene_builds) OR -b <unigene_build_ids>, don't use both.\n\n";
    } elsif ($opt_o && $opt_b) {
	print STDERR "Sorry, the parameter -o <organism_name> is incompatible with -b <unigene_build_ids>.)";
	die " Please use: -o <organism_name> OR -b <unigene_build_id>, don't use both.\n\n";
    } else {
        my $st="'C', 'P', 'D'";
	if ($opt_s) {
	    my @st;
	    my @status=split(/,/, $opt_s);
	    foreach my $status (@status) {
		if ($status =~ m/^c$|^current$/i) {
		    push @st, "'C'";
		} elsif ($status =~ m/^p$|^previous$/i) {
		    push @st, "'P'";
		} elsif ($status =~ m/^d$|^deprecated$/i) {
		    push @st, "'D'";
		} else {
		    die "Sorry, the status used ($opt_s) are not valid. Permited values are: Current, Previous or/and Deprecated\n\n";
		}
	    }
	    $st=join(', ', @st);
	}
	my ($query, $sth);
	print STDERR "\nDATASET_SELECTED:\n";
	if ($opt_A) {
	    print STDERR "\tAll the Unigene Builds:";
	    $query="SELECT unigene_build_id FROM sgn.unigene_build WHERE status IN ($st) ORDER BY organism_group_id, build_date";
	    $sth=$dbh->prepare($query);
	    $sth->execute();
	    while (my $unigene_build_id=$sth->fetchrow_array() ) {
		push @unigene_build_ids, $unigene_build_id;
	    }
	} elsif ($opt_o) {
	    print STDERR "\t Unigene Builds for organism ($opt_o):";
	    my @organism=split(/,/, $opt_o);
	    foreach my $organism (@organism) {
		my $c=0;
		my $org="'%".$organism."%'";
		$query="SELECT unigene_build_id FROM sgn.unigene_build 
                            JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id 
                            JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id 
                            WHERE organism_name ILIKE $org AND status IN ($st) ORDER BY organism_group_id, build_date";
		$sth=$dbh->prepare($query);
		$sth->execute();
		
		while (my ($unigene_build_id_o)=$sth->fetchrow_array() ) {
		    $c++;
		    push @unigene_build_ids, $unigene_build_id_o;
		}
		if ($c < 1) {
		    print STDERR "\nWarning!!!,\n";
		    print STDERR "The organism:$org have not any unigene_build_id associated to it with status($st).\n\n";
		}
	    }
	} elsif ($opt_b) {
	    print STDERR "\tUnigene Builds for unigene build ids:";
	    my @unigene_build=split(/,/, $opt_b);
	    foreach my $unigene_build (@unigene_build) {
		$query="SELECT unigene_build_id FROM sgn.unigene_build WHERE unigene_build_id=? AND status IN ($st)";
		$sth=$dbh->prepare($query);
		$sth->execute($unigene_build);
		my ($unigene_build_id_b)=$sth->fetchrow_array();
		if (!$unigene_build_id_b) {
		    print STDERR "\nWarning!!!, ";
		    print STDERR "Do not exists any unigene_build_id ($unigene_build) associated to it with status($opt_s).\n\n";
		} else {
		    push @unigene_build_ids, $unigene_build_id_b;
		}
	    }
	}		
	my $unig_build_nr=scalar(@unigene_build_ids);
	if ($unig_build_nr < 1) {
	    print STDERR "Sorry, there are enought unigene_build_ids to do the report.\n";
	    die " Please check that the -o <organism_name> or -b <unigene_build> detailed exists in the database used.\n\n ";
	}
	my $n=0;
	while ($n < $unig_build_nr) {
	   $unig_build_nr -= 1; 
	   my $m;
	   if ($unig_build_nr < 50) {
	       $m=$unig_build_nr;
	   } else {
	       $m=49;
	   }
	   my $list=join(',', @unigene_build_ids[$n..$m]);
	   print STDERR "\t\t$list\n";
	   $n += 50;
	   $unig_build_nr -= 49;
	}
    }
    print STDERR "\n\n";
    return ($dbh, \@unigene_build_ids);
}

=head2 get_unigene_assembly_stats

  Usage: my ($unigene_stats_href, $unigebe_stats_ext_href)=get_unigene_assembly_stats($dbh, $unigene_build_id);
  Desc: Get then unigene stats for a concrete unigene id
  Ret: Two hash reference with key=field and value=value. (The second contains additional data as unigene_ids)
  Args: A database conection object, $dbh and an unigene build id, $unigene_build_id
  Side_Effects: the data is formated
  Example: my ($unigene_stats_href, $unigene_stats_href_ext)=get_unigene_assembly_stats($dbh, $unigene_build_id);

=cut

sub get_unigene_assembly_stats {
    my $dbh=shift;
    my $unigene_build_id=shift;
    my (%unig_stats, %unig_stats_ext);

    unless (exists $unig_stats{'01_build_id'}) {
	$unig_stats{'01_build_id'}=$unigene_build_id;
    }
    unless (exists $unig_stats_ext{'01_build_id'}) {
	$unig_stats_ext{'01_build_id'}=$unigene_build_id;
    }
    my @organism_list;
    my $query1="SELECT organism_name FROM sgn.organism JOIN sgn.group_linkage ON sgn.organism.organism_id=sgn.group_linkage.member_id
                 JOIN sgn.unigene_build ON sgn.group_linkage.group_id=sgn.unigene_build.organism_group_id 
                 WHERE sgn.unigene_build.unigene_build_id=? ORDER BY sgn.organism.organism_id";
    my $sth1=$dbh->prepare($query1);
    $sth1->execute($unigene_build_id);
    while (my($organism_name)=$sth1->fetchrow_array() ) {
	my $organism_length=length($organism_name);
	my $organism_name_formated;
	if ($organism_length > 50) {
	    $organism_name_formated=substr $organism_name, 0, 47;
	    $organism_name_formated .= '...';
	} else {
	    $organism_name_formated=$organism_name;
	}
	push @organism_list, $organism_name_formated;
    }
    if (scalar(@organism_list) > 0) {
	unless (exists $unig_stats{'02_organism'}) {
	    $unig_stats{'02_organism'}=\@organism_list;
	}
    }
    my $query2="SELECT build_nr, build_date, status FROM sgn.unigene_build WHERE unigene_build_id=?";
    my $sth2=$dbh->prepare($query2);
    $sth2->execute($unigene_build_id);
    my ($build_nr, $build_date, $status)=$sth2->fetchrow_array();
    unless (!$build_nr && exists $unig_stats{'03_build_nr'}) {
	$unig_stats{'03_build_nr'}=$build_nr;
    } 
    unless (!$build_date && exists $unig_stats{'04_build_date'}) {
	$unig_stats{'04_build_date'}=$build_date;
    }
    unless (!$status && exists $unig_stats{'05_status'}) {
	if ($status eq 'C') {
	    $unig_stats{'05_status'}='current';
	} elsif ($status eq 'P') {
	    $unig_stats{'05_status'}='previous';
	} elsif ($status eq 'D') {
	    $unig_stats{'05_status'}='deprecated';
	}
    }
    my $query3="SELECT COUNT(DISTINCT unigene_id) AS unigene_nr, COUNT(DISTINCT est_id) AS member_id, 
                 COUNT(DISTINCT consensi_id) AS contig_nr FROM sgn.unigene JOIN sgn.unigene_member USING(unigene_id) 
                 WHERE unigene_build_id=?";
    my $sth3=$dbh->prepare($query3);
    $sth3->execute($unigene_build_id);
    my ($unigene_nr, $member_nr, $contig_nr)=$sth3->fetchrow_array();
    unless (exists $unig_stats{'06_unigene_nr'}) {
	$unig_stats{'06_unigene_nr'}=$unigene_nr;
    }
    unless (exists $unig_stats{'07_contig_nr'}) {
	$unig_stats{'07_contig_nr'}=$contig_nr;
    }
    unless (exists $unig_stats{'08_member_nr'}) {
	$unig_stats{'08_member_nr'}=$member_nr;
    }
    my $oper=$contig_nr * 100 / $unigene_nr;
    my $redun_calc=sprintf ('%.2f', $oper);
    $redun_calc .= ' %';
    unless (exists $unig_stats{'09_redundancy'}) {
	$unig_stats{'09_redundancy'}=$redun_calc;
    }
    my $query4="SELECT unigene_id FROM sgn.unigene WHERE unigene_build_id=? ORDER BY unigene_id ASC LIMIT 1";
    my $sth4=$dbh->prepare($query4);
    $sth4->execute($unigene_build_id);
    my ($first_unigene_id)=$sth4->fetchrow_array();
    my $first_unigene='SGN-U'.$first_unigene_id;
    unless (exists $unig_stats_ext{'10_first_unigene_id'}) {
	$unig_stats_ext{'10_first_unigene_id'}=$first_unigene;
    }
    my $query5="SELECT unigene_id FROM sgn.unigene WHERE unigene_build_id=? ORDER BY unigene_id DESC LIMIT 1";
    my $sth5=$dbh->prepare($query5);
    $sth5->execute($unigene_build_id);
    my ($last_unigene_id)=$sth5->fetchrow_array();
    my $last_unigene='SGN-U'.$last_unigene_id;
    unless (exists $unig_stats_ext{'11_last_unigene_id'}) {
	$unig_stats_ext{'11_last_unigene_id'}=$last_unigene;
    }
    return (\%unig_stats, \%unig_stats_ext);
}

=head2 get_unigene_annotation_stats

  Usage: my ($unig_blast_stats_href, $unig_prot_stats_href)=get_unigene_annotation_stats($dbh, $unigene_build_id);
  Desc: Get the unigene annotation stats (how many unigene has blast annotation, predicted protein...) in a href.
  Ret: Two $unig_stats_href, hash reference with key: annotation type and value: count (for blast annotations and protein annotations)
  Args: $dbh, a database conection object and $unigene_build_id, a unigene build id
  Side_Effects: the data is formated
  Example: my ($unig_annot_stats_blast_href, $unig_annot_stats_prot_href)=get_unigene_annotation_stats($dbh, $unigene_build_id);

=cut

sub get_unigene_annotation_stats {
    my $dbh=shift;
    my $unigene_build_id=shift;
    my (%unig_stats_blast, %unig_stats_prot);

    unless (exists $unig_stats_blast{'00_build_id'})  {
	$unig_stats_blast{'00_build_id'}=$unigene_build_id;
    }
    unless (exists $unig_stats_prot{'00_build_id'}) {
	$unig_stats_prot{'00_build_id'}=$unigene_build_id;
    }
    my $query0="SELECT db_name, blast_target_id FROM sgn.blast_targets ORDER BY blast_target_id ASC";
    my $sth0=$dbh->prepare($query0);
    $sth0->execute();
    while (my ($db_name, $blast_target_id)=$sth0->fetchrow_array()) {
	my $query1="SELECT COUNT(DISTINCT blast_annotation_id) FROM sgn.blast_annotations 
                     JOIN sgn.unigene ON sgn.blast_annotations.apply_id=sgn.unigene.unigene_id
                     JOIN sgn.blast_targets ON sgn.blast_annotations.blast_target_id=sgn.blast_targets.blast_target_id
                     WHERE sgn.unigene.unigene_build_id=? AND sgn.blast_annotations.blast_target_id=? AND n_hits > 0";
	my $sth1=$dbh->prepare($query1);
	$sth1->execute($unigene_build_id, $blast_target_id);
	my($blast_hits)=$sth1->fetchrow_array();
	unless (exists $unig_stats_blast{'01_blast_'.$db_name}) {
	    $unig_stats_blast{'01_blast_'.$db_name}=$blast_hits;
	}
    }
    my $n=1;
    my $query2="SELECT method FROM sgn.cds GROUP BY method ORDER BY method ASC";
    my $sth2=$dbh->prepare($query2);
    $sth2->execute();
    while (my ($method)=$sth2->fetchrow_array() ) {
	my $query3="SELECT COUNT(DISTINCT sgn.cds.cds_id), COUNT(DISTINCT sgn.domain_match.cds_id) FROM sgn.cds 
                    JOIN sgn.unigene USING (unigene_id) 
                    LEFT JOIN sgn.domain_match ON sgn.cds.cds_id=sgn.domain_match.cds_id 
                    LEFT JOIN sgn.domain ON sgn.domain_match.domain_id=sgn.domain.domain_id 
                    WHERE method=? AND unigene_build_id=?";
        my $sth3=$dbh->prepare($query3);
        $sth3->execute($method, $unigene_build_id);
        my ($prot_pred_hits, $prot_prod_with_domain)=$sth3->fetchrow_array();
	my $m=1;
	unless (exists $unig_stats_prot{$n.$m.'_protein_pred_by_'.$method}) {
	    $unig_stats_prot{$n.$m.'_protein_pred_by_'.$method}=$prot_pred_hits;
	    $m++;
	}
	unless (exists $unig_stats_prot{$n.$m.'_domain_matches_for_'.$method.'_method'}) {
	    $unig_stats_prot{$n.$m.'_domain_matches_for_'.$method.'_method'}=$prot_prod_with_domain;
	    $m++;
	}
	$n++;
    }
    return (\%unig_stats_blast, \%unig_stats_prot);
}

=head2 print_table

  Usage: print_table(\@data_href, $table_title);
  Desc: print a table, the argument should be a hash reference where keys=>headers and values=>data, if value is an array reference
        it will print the one array element per line
  Ret: none
  Args: \@data_href, an array reference with elements that are hash reference with keys=>headers and values=>data 
        and $table_title table title
  Side_Effects: change table title to uppercase
  Example: print_table(\@unig_stats_href, $table_title);

=cut

sub print_table {
    my $data_aref=shift;
    my $table_title=shift;
    my @data_href=@$data_aref;
    
    my @preheaders= sort keys %{$data_href[0]};
    my ($size_meass, $max_data)=(0,0);
    my %sizes;
    foreach my $data_href (@data_href) {
	my %data=%{$data_href};
	my @data_count=sort keys %data;
	my $data_count=scalar(@data_count);
	if ($data_count > $max_data) {
	    $max_data=$data_count;
	}
	my @preheaders= sort keys %data;
	foreach my $prehead (@preheaders) {
	    my $size_prehead=length($prehead);
	    $size_prehead -= 3;
	    unless (exists $sizes{$prehead}) {
		$sizes{$prehead}=$size_prehead;
	    }
	    my $value=$data{$prehead};
	    if (ref $value && $value != 0 ) {
		my @value_list=@$value;
		foreach my $single_value (@value_list) {
		    my $size_single_value=length($single_value);
		    if ($size_single_value > $sizes{$prehead}) {
			$sizes{$prehead}=$size_single_value;
		    }
		}
	    } else {
		my $size_value=length($value);
		if ($size_value > $sizes{$prehead}) {
		    $sizes{$prehead}=$size_value;
		}
	    }
	}
    }
    
    my @headers=sort keys %sizes;
    my $total_length;
    foreach my $line_z (@headers) { 
	my $z=$sizes{$line_z};
	$z += 2;
	$total_length += 1;
	my $x=0;
	while ($x < $z) {
	    $total_length += 1;
	    $x++;
	}
    }
    $total_length -= 1;
    print "+";
    my $g=0;
    while ($g < $total_length) {
	print "=";
	$g++;
    } 
    print "+\n";
    $total_length -= 2;
    my $title=sprintf '%-*s', $total_length, uc($table_title).':';
    print "| $title |\n";
    foreach my $line0 (@headers) { 
	my $s0=$sizes{$line0};
	$s0 += 2;
	print "+";
	my $a0=0;
	while ($a0 < $s0) {
	    print "=";
	    $a0++;
	}
    }
    print "+\n";
    foreach my $head (@headers) {
	my $headsize=$sizes{$head};
	my $headprint=$head;
	$headprint =~ s/^\d\d_//;
	my $print=sprintf '%-*s', $headsize, $headprint;
	print "| $print ";
    }
    print "|\n";
    foreach my $line1 (@headers) { 
	my $s1=$sizes{$line1};
	print "+";
	my $a1=0;
	$s1+=2;
	while ($a1 < $s1) {
	    print "=";
	    $a1++;
	}
    }
    print "+\n";
    my $lines_n=scalar(@data_href);
    my $l=0;
    my $white_l=0;
    $white_l -=1;
    my @arraydata;
    foreach my $line_href (@data_href) {
	my %line_data=%{$line_href};
	my @datatypes=sort keys %line_data;
	my $datatypes_n=scalar(@datatypes);
	my ($p, $size_run)=(0,0);
	$l++;
	
	foreach my $datatype (@datatypes) {
	    my $size=$sizes{$datatype};
	    my $value=$line_data{$datatype};
	    if (ref $value && $value != 0) {
		@arraydata=@$value;
		$white_l=$p;
		my $first_value=shift(@arraydata);
		my $line_fv=sprintf '%-*s', $size, $first_value;
		print "| $line_fv ";
		$size_run += $size;
		$size_run += 3;
	    } else {
		my $line_v=sprintf '%-*s', $size, $value;
		print "| $line_v ";
		$size_run += $size;
		$size_run += 3;
	    }
	    $p++;
	}
	if ($datatypes_n < $max_data) {
	    my $empty_fields_n=$datatypes_n;
	    my $diff=$total_length-$size_run;
	    my $empty_line=sprintf '%-*s', $diff, 'none';
	    print "| $empty_line ";
	}
		
	if ($white_l >= 0) {
	    my $arraydata_elements=scalar(@arraydata);
	    while ($arraydata_elements > 1) {
                $arraydata_elements=scalar(@arraydata);
		print "|\n";
		my $e=0;
		while ($e < $datatypes_n) {
		    if ($white_l == $e) {
			my $size_w=$sizes{$datatypes[$white_l]};
			my $value_wl=shift(@arraydata);
			my $line_wl=sprintf '%-*s', $size_w, $value_wl;
			print "| $line_wl ";
		    } else {
			my $empty_size=$sizes{$datatypes[$e]};
			my $empty_line=sprintf '%-*s', $empty_size, ' ';
			print "| $empty_line ";
		    }
		    $e++;
		}
	    }
	    $white_l -= 1;
	    
	}
	print "|\n";
	if ($l < $lines_n) {
	    foreach my $line0 (@headers) { 
	       my $s0=$sizes{$line0};
	       print "+";
	       $s0+=2;
	       my $t=0;
	       while ($t < $s0) {
	          print "-";
	          $t++;
	       }
            }
            print "+\n";
        }
    }
    foreach my $line0 (@headers) { 
	my $s0=$sizes{$line0};
	$s0 +=2;
	print "+";
	my $t=0;
	while ($t < $s0) {
	    print "=";
	    $t++;
	}
    }
    print "+\n";
}
