#!/usr/bin/perl

=head1 NAME

 unigene_versions_mapping_tool.pl
 A script to map different unigene versions.(version.1.1.).

=cut

=head1 SYPNOSIS

unigene_versions_mapping_tool.pl -H <dbhost> -D <dbname> -b <basename> [-o <organism_name> || -u <unigene_build>]
    
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
      


=item -h 

B<help>                  show the help

=back

=cut

=head1 DESCRIPTION

This script get the different versions of the unigenes and print a file with them

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

SGN_membership_mapping_tool.pl


=cut

use strict;
use Getopt::Std;
use CXGN::DB::InsertDBH;
use File::Basename;

our ($opt_H, $opt_D, $opt_b, $opt_o, $opt_u, $opt_h);
getopts("H:D:b:o:u:h");
if (!$opt_H && !$opt_D && !$opt_b && !$opt_o && !$opt_u && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

## Validating the input parameters

if ($opt_h) {
    help();
}

unless ($opt_H && $opt_D) {
    die ("Required argument -H <dbhost> AND -D <dbname> were not supplied.\n");
}
unless ($opt_b) {
    die ("Required argument -b <basename> was not supplied.\n");
}
unless ($opt_o) {
    die ("Required argument -o <organism_group_name> was not supplied.\n");
}

## Connecting with the database

our $dbh = CXGN::DB::InsertDBH::connect({	dbname => $opt_D ,
						dbhost => $opt_H ,
					  	dbargs => { RaiseError => 1, PrintError => 1 }
					})
	or die "Failed to connect to database ($DBI::errstr)";


## Creating the undef variables to define them later

my %unigene_groups;
my @mapping_unigenes;

## Selecting the group_id and the comment in the group table, using the -o option. To match with the data in the comment field it 
 ## use ILIKE, and print the options. I use the $group_query variable over the query because it don't work when it is used in
  ## execute function. Put the results into a hash for use later.

my $group_query = "'%".$opt_o."%'";
my $query = "SELECT group_id, comment FROM sgn.groups WHERE comment ILIKE $group_query AND type=1 ORDER BY group_id";
my $sth = $dbh->prepare($query);
$sth->execute();
print "The -o <organism_group_names> match with the follow groups:\n";
while (my ($group_id, $group_name) = $sth->fetchrow_array() ) {
    print "\tGroup_id=$group_id and group_name=$group_name\n";
    unless (exists $unigene_groups{$group_id}) {
	$unigene_groups{$group_id} = $group_name;
    }
}

## count how many groups are in the hash... if don't find any group... die, else continue the script.

my @groups = sort {$a <=> $b} keys %unigene_groups;
if (scalar(@groups) == 0) {
    print "\tnone\n";
    die "\nThere are not any group name that match with $opt_o.\n\n";
} else {
    foreach my $group (@groups) {

        ## Create and open a filehandles for the different groups. It means that it will store one version mapping per group

	my $organism_group_name = $unigene_groups{$group};
	$organism_group_name =~ s/\s/_/g;
	my $filename = $opt_b.'_'.$organism_group_name.'_unigene_version_mapping.tab';
	open my $fh, '>', $filename || die "Sorry, I can not oepn the file:$filename.\n\n";

        ## Getting the different unigene_build (order by unigene_build_id) and for each one, select the unigene_id and 
	 ## the members and put in a temp table. Finally put these temp tables names in an array. Also count to give a report 
          ##of how many unigenes and members get the script.

	my $unigene_build_query = "SELECT unigene_build_id FROM sgn.unigene_build WHERE organism_group_id=? ORDER BY unigene_build_id";
	my $sth1 = $dbh->prepare($unigene_build_query);
	$sth1->execute($group);
	print "\nThere are the following unigene builds associated to the group:$group:\n";
	while (my ($unigene_build_id) = $sth1->fetchrow_array() ) {
	    my $temp_table='temp_table_'.$unigene_build_id;
	    print "\t=> $unigene_build_id\t";
	    my (%unigenes, %unigenes_count);
	    my $query_unigene = "SELECT sgn.unigene.unigene_id, sgn.unigene_member.est_id INTO TEMP TABLE $temp_table FROM sgn.unigene 
                                JOIN sgn.unigene_member USING(unigene_id) WHERE unigene_build_id=?";
	    my $sth2 = $dbh->prepare($query_unigene);
	    $sth2->execute($unigene_build_id);
	    my $count1 = "SELECT COUNT(DISTINCT unigene_id), COUNT(DISTINCT est_id) FROM $temp_table";
	    $sth2 = $dbh->prepare($count1);
	    $sth2->execute();
	    my ($unigene_number, $unigene_member) = $sth2->fetchrow_array();
	    print "(unigene_nr:$unigene_number; unigene_member_nr:$unigene_member)\n";
	    push @mapping_unigenes, $temp_table;
	    
	}
	print "\n";

        ## Now the script get all the est_id for the organism that are members of the group_id and put them in a temp table with
         ## the library_shortname and the flags.

	my @est_id_list;
	my (%library_shortname, %flags, %status);
	my $query_organism = "SELECT member_id, organism_name FROM sgn.group_linkage JOIN sgn.organism 
                            ON sgn.group_linkage.member_id=sgn.organism.organism_id WHERE group_id=? ORDER BY organism_id";
	my $sth3=$dbh->prepare($query_organism);
	$sth3->execute($group);
	print "The group_id:$group has the following organism associated:\n";
	my $organism_list='(';
	while (my ($organism_id, $organism_name)=$sth3->fetchrow_array()) {
	    print "\t$organism_name\n";
	    $organism_list .= $organism_id.',';
	}
	$organism_list =~ s/,$/\)/;
	my $temp_table_est='temp_table_est';
        my $est_query = "SELECT sgn.library.library_shortname, sgn.est.est_id, sgn.est.flags INTO TEMP TABLE $temp_table_est 
                             FROM sgn.est JOIN sgn.seqread USING(read_id) JOIN sgn.clone USING(clone_id) 
                             JOIN sgn.library USING(library_id) JOIN sgn.organism USING(organism_id) WHERE organism_id IN $organism_list 
                             ORDER BY library_shortname, est_id";
        my $sth4 = $dbh->prepare($est_query);
        $sth4->execute();
        my $count_est = "SELECT COUNT(DISTINCT est_id) FROM $temp_table_est";
        $sth4 = $dbh->prepare($count_est);
        $sth4->execute();
        my ($est_nr)=$sth4->fetchrow_array();
        print "Total est_nr:$est_nr\n";
	 
	print "\n\n";

        ## Finally use a query search with the temp tables where it was stored the unigene member information. 
         ## If don't find any unigene_id for a concrete replace the variable for 'absent'. Also sustitute the flags for a tag using
	  ## a function. Print all the results into a file. 

	my $estnumber = 0;
	my $esttotal = scalar(@est_id_list);
        my $complexquery = "SELECT $temp_table_est.library_shortname, $temp_table_est.est_id, $temp_table_est.flags ";
        foreach my $temp_table_unigene (@mapping_unigenes) {
           my $varname=$temp_table_unigene.'_unigene_id';
	   $complexquery .= ", $temp_table_unigene.unigene_id AS $varname ";
        }
	$complexquery .= "FROM $temp_table_est ";
	foreach my $temp_table_unigene_other (@mapping_unigenes) {
	    $complexquery .= "LEFT JOIN $temp_table_unigene_other ON $temp_table_est.est_id=$temp_table_unigene_other.est_id ";
	}
	$complexquery .= "ORDER BY $temp_table_est.library_shortname, $temp_table_est.est_id";
	print "Doing the unigene query...:\n$complexquery\n";
	my $sth5 = $dbh->prepare($complexquery);
	$sth5->execute();
	while (my @data=$sth5->fetchrow_array) {
	    $estnumber++;
	    my @unigene_mapped;
	    my $est_library = shift(@data);
	    my $est_accession = "SGN-E".shift(@data);
 	    my $flags_tags = flags_equivalences(shift(@data));					
	    foreach my $unigene_map (@data) {
		unless (defined $unigene_map) {
		    push @unigene_mapped, 'absent';
		} else {
		    push @unigene_mapped, 'SGN-U'.$unigene_map;
		}
	    }
	    my $unigene_string=join("\t", @unigene_mapped);
	    my $string="$est_library\t$est_accession\t$flags_tags\t$unigene_string";

	    print STDERR "Processing the est_number:$estnumber ($est_accession)    \r";
	   
	    print $fh "$string\n";
	}	
    }
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
     This script gets the different versions for the unigene build for a organism
    Usage: 
     unigene_versions_mapping_tool.pl [-h] -H <dbhost> -D <dbname> -b <basename> [-o <organism_group_name>]    
    Flags:
      -H database hostname      for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name          sandbox or cxgn etc (mandatory)
      -b basename               basename for the output_file (mandatory)
      -o organism_group_name    organism_name for the search (mandatory)
      -h this help

EOF
exit (1);
}

=head2 flags_equivalences

  Usage: my $flag_tag=flags_equivalences($flag);
  Desc: Give a tag using a flag value according the equivalence hash
  Ret: $flag_tag, a scalar
  Args: $flag, a scalar, an integer
  Side_Effects: none
  Example: $flag_tag=flags_equivalences($flag);

=cut

sub flags_equivalences {
    my $flag = shift;
    my $tag;

    if ($flag != 0) {
	if ($flag >= 1024) {
	    $tag .= 'plastid contamination;';
	    $flag -= 1024;
	}
	if ($flag >= 512) {
	    $tag .= 'manually censored;';
	    $flag -= 512;
	}
	if ($flag >= 256) {
	    $tag .= 'possible chimera (unigene_assembly);';
	    $flag -= 256;
	}
	if ($flag >= 128) {
	    $tag .= 'possible chimera (arabidopsis screen);';
	    $flag -= 128;
	}
	if ($flag >= 64) {
	    $tag .= 'rRNA contamination;';
	    $flag -= 64;
	} 
	if ($flag >= 32) {
	    $tag .= 'cloning host contamination;';
	    $flag -= 32;
	}
	if ($flag >= 16) {
	    $tag .= 'low complexity;';
	    $flag -= 16;
	}
	if ($flag >= 8) {
	    $tag .= 'high expected base error;';
	    $flag -= 8;
	}
	if ($flag >= 4) {
	    $tag .= 'insert too short;';
	    $flag -= 4;
	}
	if ($flag >= 2) {
	    $tag .= 'possible chimera (vector screen);';
	    $flag -= 2;
	}
    } else {
	$tag = 'WITHOUT-FLAGS';
    }
    return $tag;
}


