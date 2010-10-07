#!/usr/bin/perl

=head1 NAME

 get_proteinseq_from_sgndb.pl.
 A script to get the protein predicted sequences from the unigenes of the SGN databases.(version.2.0.).

=cut

=head1 SYPNOSIS

get_proteinseq_from_sgndb.pl -H <dbhost> -D <dbname> -b <basename> [-o <organism_name> OR -u <unigene_build_id>] 
-p <protein prediction method> [-C OR/AND -P] [-M] [-A|-a<specify_annotation_parameter>][-h]

=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory)

=item -b

B<basename>             basename for the output fasta files (mandatory)
      
=item -o

B<organism name>        organism name to get the sequences, this option get only the protein predicted for the last unigene build (comandatory with -u).

=item -u

B<unigene build id>     unigene build id to get the protein predicted. If you do not specify anything (like SGN_pep-from-db.pl search all the unigene build ids and the organism of them and print (comandatory with -o).

=item -p

B<prot pred method>     protein prediction method used. (the most used are 'estscan','longest6frame' or 'preferred')(mandatory)

=item -C

B<get cds>              enable the option get cds (comandatory with -P)

=item -P

B<get protein>          enable the option get protein (comandatory with -C)

=item -M

B<get cds-unigene map>  enable the option map the cds_id against unigene_id and put in a file 

=item -U

B<use SGN-U as ids>     use the unigene ids as sequences ids

=item -A

b<ids with annotations> enable the option to get the sequences with annotation from unigene blast results

=item -a

B<specify annotation>   specify a) description length per database (example:'100') or b) database (example:'genbank')

=item -h

B<help>                 print help

=back 

=cut

=head1 DESCRIPTION

 This script get in fasta format the DNA sequences for the predicted cds OR/AND the predicted protein sequences for an unigene build.

 Note: If the cds/protein has cds_seq='NULL' and protein_seq='NULL', this script will skip the protein.
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

get_proteinseq_from_sgndb.pl


=cut


use strict;
use warnings;

use Getopt::Std;
use CXGN::DB::InsertDBH;
use Bio::SeqIO;
use Bio::Seq;
use File::Basename;

our ($opt_H, $opt_D, $opt_b, $opt_o, $opt_u, $opt_p, $opt_C, $opt_P, $opt_M, $opt_U, $opt_A, $opt_a, $opt_h);
getopts("H:D:b:o:u:p:CPMUAa:h");
if (!$opt_H && !$opt_D && !$opt_b && !$opt_o && !$opt_u && !$opt_p && !$opt_C && !$opt_P && !$opt_M && !$opt_A && 
    !$opt_a && !$opt_U && !$opt_h) {
               print "There are n\'t any tags. Print help\n\n";
               help();
}

my ($dbh, $unigene_build_id, $pp_method, $basename)= validate_input();
my ($cds_data_aref, $pep_data_aref, $cdsmap_file)=get_sequences($dbh, $unigene_build_id, $pp_method, $basename);

print "\n----------------------------------------\n";
print "Report:";
print "\n----------------------------------------\n";
print "Input options:\n";
print "\t* Unigene build id: $unigene_build_id\n";
print "\t* Protein prediction method: $pp_method\n";

if ($opt_A) {
    print "\t* Annotation parameters: Length = 250; Databases = All\n";
} elsif ($opt_a) {
    if ($opt_a =~ m/^\d+$/) {
	print "\t* Annotation parameters: Length = $opt_a; Databases = All\n";
    } else {
	if ($opt_a ne 'proteinpilot') {
	    print "\t* Annotation parameters: Length = 250; Databases = $opt_a\n";
	} else {
	     print "\t* Annotation parameters: Length = 250; Databases = All; Enabled ProteinPilot filter (add lcl| to ID)\n";
	 }
    }
} else {
    print "\t* Without any annotations\n";
}

print "Output data:\n";
if ($opt_C) {
    my @cds_data=@$cds_data_aref;
    print "\t* cds sequences:\n";
    print "\t\t- cds fasta file: $cds_data[0]\n";
    print "\t\t- number of sequences: $cds_data[1]\n";
    print "\t\t- total bases: $cds_data[2]\n";
    my $cds_average_length=$cds_data[2]/$cds_data[1];
    my $short_cds_average_length=int $cds_average_length;
    print "\t\t- average sequences length: $short_cds_average_length\n";
    print "\t\t- maximum sequence length: $cds_data[4]\n";
    print "\t\t- minimum sequence length: $cds_data[3]\n";    
}
if ($opt_P) {
    my @pep_data=@$pep_data_aref;
    print "\t* protein sequences:\n";
    print "\t\t- protein fasta file: $pep_data[0]\n";
    print "\t\t- number of sequences: $pep_data[1]\n";
    print "\t\t- total aminoacids: $pep_data[2]\n";
    my $pep_average_length=$pep_data[2]/$pep_data[1];
    my $short_pep_average_length=int $pep_average_length;
    print "\t\t- average sequences length: $short_pep_average_length\n";
    print "\t\t- maximum sequence length: $pep_data[4]\n";
    print "\t\t- minimum sequence length: $pep_data[3]\n";    
}
if ($opt_M) {
    print "\t* cds unigene mapping:\n";
    print "\t\t- cds-unigene map tab file: $cdsmap_file\n";
}
print "\n----------------------------------------\n";


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
      This script get the cds (-C) or protein (-P) sequences or the protein ids <=> unigene ids mapping (-M) 
      associated to an organism (-o) or an unigene build (-u) for a predicted protein method (-p). The sequences
      output file are in fasta format. The id can be protein ids (SGN-P) or unigene ids (SGN-U) (-U).
      
      This sequence in fasta format can contain the blast annotation associated to unigenes. There are two ways to
      get these annotations:
          1- Using -A parameter, it get the best blast score for each blast annotation database (genbank, arabidopsis
             and/or swissprot). The format will be the follow:
	                         
                                 >ID DESCRIPTION...........
	                         ...........sequence.......
	     where:
               -id is sequence id '(SGN-P#cds_id))' or '(SGN-U#unigene_id)' if you use -U option
               -description as '(*database_code) (e-value) defline;' for each blast database
	       -sequence, nucleotide (-C option) or protein (-p option).

             Note: The genbank defline for nr can contain more than one accession (nr has merged deflines for
                   the same sequences). This script only takes the first annotation (accession+description) for
                   this defline, and discard the rest.
		   
          2- Using -a <filter> you can filter the description using two ways:
		   a) the line length per database (by default is used 250 characters).
                   b) the database (by default get all, but you can spacify genbank, or swissprot)


    Note:
      If the cds/protein has cds_seq='NULL' and protein_seq='NULL', this script will skip the protein.
    
    Usage: 
      get_proteinseq_from_sgndb.pl -H <dbhost> -D <dbname> -b <basename> [-o <organism_name> OR -u <unigene_build_id>] 
      -p <protein prediction method> [-C OR/AND -P] [-M] [-A|-a <specific_annotation> [-U] [-h]
    
    Example:
      SGN_get_pep-from-db.pl -H localhost -D sandbox -b nt_protein_prediction -o 'Nicotiana tabacum'
      -p 'estscan' -C -P

    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -b basename             basename for the output file, the name will be basename
                              +protein prediction method+(cds or pep)+.fasta (mandatory)
      -o organism name        organism name for the unigene build. Get the lastest unigene (with status 'C') 
                              build of this organism. This option can be used only when the unigene build is
                              composed by one organism. For other cases, use argument -u. (comandatory with -u)
      -u unigene build id     unigene build id. If do not specify any, report all the possible options 
                              (write the differents unigene builds, by organism, date, status and unigene build id)
                              (comandatory with -o)
      -p prot. pred. method   method used in the prediction of the cds and the protein. The most commons are 'estscan' 
                              and 'longest6frame' or 'preferred' to get the preferred proteins
      -C get cds              get the predicted cds sequences in fasta format. The outfile will be basename+
                              cds.fasta (comandatory with -P)
      -P get protein          get the predicted protein sequences in fasta format. The outfile will be basename+
                              pep.fasta (comandatory with -C)
      -M map unigene-cds      get the cds_id (SGN-P------) and the unigene_id (SGN-U------) and print in a file
                              in a tab format
      -U use SGN-U ids        use the SGN-U------ as sequences ids (with unigene_id)
      -A get annotations      get complete annotations for the id.
			       (by default get 250 characters per database annotation).
      -a specify annotations  specify the length of the annotations (per database) or the database. Example can be
			      '100', '200' characters length, 'genbank', 'arabidopsis' or 'swissprot' or 'proteinpilot'
                              (special format compatible with some programs where the ID begin with lcl| and there is not more 
                               special characters)
      -h this help

EOF
exit (1);
}


=head2 validate_input

  Usage: my ($dbh, $unigene_build_id, $pp_method, $basename)=validate_input(); 
  Desc: Check the input arguments
  Ret: Four scalar, $dbh (conection database object), $unigene_build_id (an integer), $pp_method (protein prediction 
       method) and $basename.
  Args: none 
  Side_Effects: If something is wrong die. Also print information about unigene build id if -u is wrong.
  Example: my ($dbh, $unigene_build_id, $pp_method, $basename)=validate_input();

=cut

sub validate_input {
    if ($opt_h) {
	help();
    }
    unless ($opt_H) {
       die ("Required argument -H <database host> was not supplied.\n");
    }
    unless ($opt_D) {
       die ("Required argument -D <database name> was not supplied.\n");
    }
    unless ($opt_b) {
       die ("Required argument -b <basename for output files> was not supplied.\n");
    }
    unless ($opt_p) {
	die ("Required argument -p <protein prediction method> was not supplied.\n");
    }
    unless ($opt_o || $opt_u) {
	die ("Required arguments -o <organism name> OR -u <unigene build id> were not supplied.\n");
    }
    if ($opt_o && $opt_u) {
	die ("Arguments conflict, -o and -u are exclusive arguments. You only can use one, not both.\n");
    }
    unless ($opt_C || $opt_P) {
	die ("Requiered arguments -C <get cds sequences> AND/OR -P <get protein sequences> were not suppled.\n");
    }
    if ($opt_u) {     
        unless ($opt_u =~ m/\d+/ || $opt_u eq '?') {
	   die ("Invalid argument -u.It must be an integer.\nIf you do not know, and you want print the build ids use ?");
        }
    }

    our $dbh = CXGN::DB::InsertDBH::connect({	dbname=>$opt_D,
						dbhost=>$opt_H
	                                     })
	or die "Failed to connect to database ($DBI::errstr)";

    my $unigene_build_id;
    if ($opt_o) {
	my $query="SELECT COUNT(organism_group_id) FROM sgn.unigene_build 
                   JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id 
                   JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id 
                   WHERE sgn.unigene_build.status = 'C' and sgn.organism.organism_name = ?";
        my $sth=$dbh->prepare($query);
        $sth->execute($opt_o);
	my ($organism_group_id_count)=$sth->fetchrow_array();
	if ($organism_group_id_count == 0) {
	    die ("There are not any unigene build with status 'C' associated to $opt_o organism.\n");
	} elsif ($organism_group_id_count > 1) {
	    print "There are more than one unigene build with status 'C' associated to $opt_o organism.\n";
            my $only_member=0;
	    my $query2="SELECT organism_group_id FROM sgn.unigene_build 
                        JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id 
                        JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id 
                        WHERE sgn.unigene_build.status = 'C' and sgn.organism.organism_name = ?";
            $sth=$dbh->prepare($query2);
            $sth->execute($opt_o);
	    while (my ($organism_group_id)=$sth->fetchrow_array() ) {
		my $query3="SELECT COUNT(member_id) FROM sgn.group_linkage WHERE group_id = ?";
		my $sth2=$dbh->prepare($query3);
                $sth2->execute($organism_group_id);
		my ($member_count) = $sth2->fetchrow_array();
		if ($member_count == 1) {
		    $only_member=$organism_group_id;
		}
	    }
	    if ($only_member == 0) {
		print "There are not any unigene build with status 'C' associated ONLY to $opt_o organism."; 
                die (" Use the -u option\n");
	    } else {
		my $query4="SELECT unigene_build_id FROM sgn.unigene_build WHERE status = 'C' and organism_group_id = ?";
		$sth=$dbh->prepare($query4);
                $sth->execute($only_member);
		($unigene_build_id)=$sth->fetchrow_array();
	    }
	} else { 
	    my $query5="SELECT sgn.unigene_build.unigene_build_id FROM sgn.unigene_build 
                        JOIN sgn.group_linkage ON sgn.unigene_build.organism_group_id=sgn.group_linkage.group_id 
                        JOIN sgn.organism ON sgn.group_linkage.member_id=sgn.organism.organism_id 
                        WHERE sgn.unigene_build.status = 'C' and sgn.organism.organism_name = ?";
            my $sth=$dbh->prepare($query5);
            $sth->execute($opt_o);
	    ($unigene_build_id)=$sth->fetchrow_array();    
	}
	print "For the organism: $opt_o, the unigene build id with status 'C' is $unigene_build_id.\n";
    } elsif ($opt_u) {
        my $organism_group_id;
        if ($opt_u ne '?') {
            my $query6="SELECT organism_group_id FROM sgn.unigene_build 
                        WHERE unigene_build_id = ?";
            my $sth=$dbh->prepare($query6);
            $sth->execute($opt_u);
	    ($organism_group_id)=$sth->fetchrow_array();
	}
	if (!$organism_group_id) {
	    print "\nSorry, the unigene build specified is not in the sgn.unigene_build table.\n";
	    print "Rerun the script with one of the follow options:\n";
	    print "---------------------------------------------------\n";
	    my $query7="SELECT unigene_build_id, build_date, status, organism_group_id FROM sgn.unigene_build 
                        ORDER BY organism_group_id, build_date";
            my $sth=$dbh->prepare($query7);
            $sth->execute();
	    while (my ($unigene_build_id, $build_date, $status, $organism_group_id)=$sth->fetchrow_array() ) {
		print "\n=> Unigene_build_id: $unigene_build_id\t";
		print "Build_date: $build_date\t";
		print "Status: $status\t\n";
		print "\t Unigene build composed by the follow organism:\n";
		my $query8="SELECT organism_name FROM sgn.organism JOIN sgn.group_linkage 
                            ON sgn.organism.organism_id=sgn.group_linkage.member_id WHERE group_id = ?";
		my $sth2=$dbh->prepare($query8);
                $sth2->execute($organism_group_id);
		while (my ($organism_name)=$sth2->fetchrow_array() ) {
		    print "\t\t* $organism_name\n";
		}
	     }
	     die ("---------------------------------------------------\n\n");
	} else {
	    $unigene_build_id=$opt_u;
	}
    }
    my ($query8, $pp_method);
    if ($opt_p eq 'preferred') {	
	$query8="SELECT COUNT(cds_id) FROM sgn.cds JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id 
                 WHERE preferred = ? AND unigene_build_id = ?";
	$pp_method = 't';
    } else {
	$query8="SELECT COUNT(cds_id) FROM sgn.cds JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id 
                 WHERE method = ? AND unigene_build_id = ?";
	$pp_method = $opt_p;
    }
   
    my $sth=$dbh->prepare($query8);
    $sth->execute($pp_method, $unigene_build_id);
    my ($cds_count)=$sth->fetchrow_array();
    if ($cds_count < 1) {
        die ("There are not any cds associated to build_unigene_id=$unigene_build_id and prot. pred. method: $opt_p\n");
    }
    return ($dbh, $unigene_build_id, $opt_p, $opt_b);
}

=head2 get_sequences

  Usage: my ($cds_data_aref, $pep_data_aref, $cdsmap_file)=get_sequences($dbh, $unigene_build_id, $pp_method, $basename);
  Desc: Get the cds AND/OR protein sequences, print them in a file, count them, and calculate size information 
  Ret: Three scalar, $cds_data_aref, an array reference with the cds data (file name, count, 
       total bp of the file and min and max length), $pep_data_aref, an array reference with the prptein data
       and the name of the file with the map between unigene and cds ids.
  Args: $dbh (database conection object), $unigene_build_id (integer with the unigene build id), 
        $pp_method (protein prediction method) and $basename (output file basename)
  Side_Effects: if do not find a sequence, print ERROR and the ids in the screen
  Example: my($cds_data_aref, $pep_data_aref, $cdsmap_file)=get_sequences($dbh, $unigene_build_id, $pp_method, $basename);

=cut

sub get_sequences {
    my $dbh=shift;
    my $unigene_build_id=shift;
    my $pp_method=shift;
    my $basename=shift;

    my ($cds_file, $pep_file, $cdsmap_file, $cds_seqio, $pep_seqio, $cdsmap_fh);
    my ($cdscount, $pepcount, $total_cds_length, $total_pep_length, $min_cds_length, 
        $max_cds_length, $min_pep_length, $max_pep_length, $n)=(0,0,0,0,1000000000,0,1000000000,0,0);

    if ($opt_C) {

	print STDERR "\nGet cds sequences from the database...\n";
        $cds_file=$basename."_cds_predicted_by_".$pp_method.".fasta";
	
	$cds_seqio = Bio::SeqIO->new( -file => ">$cds_file", -format => 'fasta');
	
    }
    if ($opt_P) {

	print STDERR "\nGet protein sequences from the database...\n";
        $pep_file=$basename."_protein_predicted_by_".$pp_method.".fasta";

        $pep_seqio = Bio::SeqIO->new( -file => ">$pep_file", -format => 'fasta');
    }
    if ($opt_M) {
	print STDERR "\nGet mapping data (SGN-P and SGN-U) from the database...\n";
        $cdsmap_file=$basename."_unigene_cds_mapping_for_".$pp_method.".tab";
	open $cdsmap_fh, '>', "$cdsmap_file" || die "Cannot create the file: $cdsmap_file\n";
    }
    my $query;
    if ($pp_method eq 'preferred') {
	$query= "SELECT sgn.cds.cds_id, sgn.cds.cds_seq, sgn.cds.protein_seq, sgn.cds.unigene_id 
                FROM sgn.cds JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id 
                WHERE unigene_build_id = ? AND preferred = ?";
	$pp_method = 't';
    } else {
	$query= "SELECT sgn.cds.cds_id, sgn.cds.cds_seq, sgn.cds.protein_seq, sgn.cds.unigene_id 
                FROM sgn.cds JOIN sgn.unigene ON sgn.cds.unigene_id=sgn.unigene.unigene_id 
                WHERE unigene_build_id = ? AND method = ?";
    }
    my $sth=$dbh->prepare($query);
    $sth->execute($unigene_build_id, $pp_method);
    while (my ($cds_id, $cds_seq, $pep_seq, $unigene_id)=$sth->fetchrow_array() ) {
	my $id;
	$n++;
	print STDERR "\tprocessing data $n...\t\t\t\r";
	if ($opt_U) {
	    $id='SGN-U'.$unigene_id;
	} else {
	    $id='SGN-P'.$cds_id;
	}
	if ($opt_C) {
	    if (!$cds_seq) {
		print "ERROR: The $id have not cds_seq.\n";
	    }
	    elsif ($cds_seq eq 'NULL') {
		print "ERROR: The $id have NULL cds_seq.\n";
	    }
	    else {
		## Clean the sequence from not word char.
		$cds_seq =~ s/\\|_|\/|\s//g;
	
		my $cdsseq = Bio::Seq->new( -display_id => $id,
		                            -seq        => $cds_seq );		

		if ($opt_A || $opt_a) {
		    my $annotations=get_annotations($unigene_id,$dbh);
		    if ($annotations) {
			$cdsseq->desc($annotations);
		    }
		}
		$cds_seqio->write_seq($cdsseq);
		
                $cdscount++;
	        my $cds_length=length $cds_seq;
	        $total_cds_length += $cds_length;
	        if ($min_cds_length > $cds_length) {
		    $min_cds_length=$cds_length;
	        } elsif ($max_cds_length < $cds_length) {
		    $max_cds_length=$cds_length;
	        }
            }
	}
        if ($opt_P) {
	    if (!$pep_seq) {
		print "ERROR: The $id have not protein_seq.\n";
	    }
	    elsif ($pep_seq eq 'NULL') {
		print "ERROR: The $id have NULL protein_seq.\n";
	    }
	    else {
		$pep_seq =~ s/\\|_|\/|\s//g;
	
		my $pepseq = Bio::Seq->new( -display_id => $id,
					    -seq        => $pep_seq );

	        if ($opt_A || $opt_a) {
		    my $annotations=get_annotations($unigene_id,$dbh);
		    if ($annotations) {
			if (defined $opt_a && $opt_a eq 'proteinpilot') {
			    $id =~ s/\|/ /g;
			    $pepseq->display_id($id);
			}
			$pepseq->desc($annotations);
		    } else {
			if (defined $opt_a && $opt_a eq 'proteinpilot') {
			    $id =~ s/\|/ /g;
			    $pepseq->display_id($id);
			}
		    }
		}
		$pep_seqio->write_seq($pepseq);
                
		$pepcount++;
	        my $pep_length=length $pep_seq;
	        $total_pep_length += $pep_length;
	        if ($min_pep_length > $pep_length) {
		    $min_pep_length=$pep_length;
	        } elsif ($max_pep_length < $pep_length) {
		    $max_pep_length=$pep_length;
	        }
            }
	}
	if ($opt_M) {
	    print $cdsmap_fh "SGN-P$cds_id\tSGN-U$unigene_id\n";
	}
    }

    my @cds_data = ($cds_file, $cdscount, $total_cds_length, $min_cds_length, $max_cds_length);
    my @pep_data = ($pep_file, $pepcount, $total_pep_length, $min_pep_length, $max_pep_length); 
    return (\@cds_data, \@pep_data, $cdsmap_file);
}

=head2 get_annotations

  Usage: my $annotation=get_annotations($unigene_id, $dbh)
  Desc: get and parse annotations from the sgn.blast tables (get the annotation with highest score)
  Ret: $annotation, a scalar with the annotations
  Args: $unigene_id and $dbh (database conection)
  Side_Effects: none
  Example: my $annotation=get_annotation($unigene_id, $dbh)

=cut

sub get_annotations {
  my $unigene_id=shift;
  my $dbh=shift;

  my ($annotations, $desclength, $defline_not, $proteinpilot, $sth);
  if ($opt_A) {
      $desclength=250;
  } elsif ($opt_a) {
      if ($opt_a =~ m/^\d+$/) {
          $desclength=$opt_a;
      } else {
	  if ($opt_a ne 'proteinpilot') {
	      $defline_not="'%".$opt_a."%'";
	      $desclength=250;
	  } else {
	      $proteinpilot=1;
	      $desclength=250;
	  }
      }
  }
  unless (defined $defline_not) {
      my $query="SELECT blast_target_id, db_name FROM sgn.blast_targets 
                 JOIN sgn.blast_annotations USING(blast_target_id) 
                 JOIN sgn.unigene ON sgn.blast_annotations.apply_id=sgn.unigene.unigene_id 
                 WHERE unigene_id=? GROUP BY db_name, sgn.blast_targets.blast_target_id 
                 ORDER BY sgn.blast_targets.blast_target_id;";
      $sth=$dbh->prepare($query);
      $sth->execute($unigene_id);
  } else {
      my $query="SELECT blast_target_id, db_name FROM sgn.blast_targets 
                 JOIN sgn.blast_annotations USING(blast_target_id) 
                 JOIN sgn.unigene ON sgn.blast_annotations.apply_id=sgn.unigene.unigene_id 
                 WHERE unigene_id=? AND db_name ILIKE ? 
                 GROUP BY db_name, sgn.blast_targets.blast_target_id 
                 ORDER BY sgn.blast_targets.blast_target_id;";
      $sth=$dbh->prepare($query);
      $sth->execute($unigene_id, $defline_not);
  }
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
  if (defined $proteinpilot && defined $annotations) {
      $annotations =~ s/\|/-/g;
  }
  return $annotations;
}
