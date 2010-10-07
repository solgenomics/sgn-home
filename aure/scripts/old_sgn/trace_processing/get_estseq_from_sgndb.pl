#!/usr/bin/perl -w

=head1 NAME

 get_estseq_from_sgndb.pl.
 A script to get the database sequences in fasta format with differents ways to select them .(version.1.0.).

=cut

=head1 SYPNOSIS

get_estseq_from_sgndb.pl -H <dbhost> -D <dbname> -b <basename> [-o <organism_names> OR -l <library_names>] [-C] [-F] [-S] [-Q]
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory)

=item -b

B<basename>             basename for the output files.

=item -A                

B<get_all_organism>     get all the organism

=item -o

B<organism_name>        organism name (or organisms names) to get the sequences

=item -L 

B<get_by_libraries>     used only with organism, get the sequences and write in files organized by libraries

=item -l

B<library_names>        Library name (or libraries names) to get the sequences.   

=item -C

B<enable clean option>  Get the sequences according to the cleaning coordenates (hqi_start and hqi_length from qc_report)

=item -F

B<enable flags option>  Get the sequences that have sgn.est.flags equal to 0.

=item -f 

B<flags to exclude>     Exclude the sequences that have these flags

=item -S

B<enable status option> Get the sequences that have sgn.est.status value equal to 0

=item -s

B<status to exclude>    Exclude the sequences that have these status

=item -Q

B<enable .qual option>  Get the quality values asociated to a sequences in a fasta format.
 
=item -R

B<enable reverse option> There are 3' sequence. With this option you the reverse of the 3' sequences. 

=item -h 

B<help>                  show the help

=back

=cut

=head1 DESCRIPTION

This script is a modification of the query-unigene-input.pl script that add some features like get the sequences by organism name or library shortname, get the sequences without trim or get all the sequences independently of the flags or status.

 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

get_estseq_from_sgndb.pl


=cut

use strict;
use warnings;

use CXGN::DB::InsertDBH;
use Getopt::Std;

our ($opt_H, $opt_D, ,$opt_b, $opt_A, $opt_o, $opt_L, $opt_l, $opt_C, $opt_F, $opt_f, $opt_S, $opt_s, $opt_Q, $opt_R, $opt_h);

getopts ("H:D:b:Ao:Ll:CFf:Ss:QRh");

if (!$opt_H && !$opt_D && !$opt_b && !$opt_A && !$opt_o && !$opt_L && !$opt_l && !$opt_C && !$opt_F && !$opt_f && !$opt_S && !$opt_s && !$opt_Q && !$opt_R && !$opt_h) {
    print STDERR "There are n\'t any tags. Print help\n\n";
    help();
}

my ($basename, $dbh, $arrayref_lib_ids) = validate_input();

## After validate the input, print a resume of the paramaters choosed by the user
print STDERR "\n\nThe request options are:\n";
if ($opt_A) {
    print STDERR "\tDataset => All the current organisms\n";
} elsif ($opt_o) {
    print STDERR "\tDataset => Organisms:$opt_o\n";
} elsif ($opt_l) {
    print STDERR "\tDataset => Libraries:$opt_l\n";
}
if ($opt_L) {    
    print STDERR "\tOutput_sequence_files => $opt_b.library_shortname.seq\n";
    if ($opt_Q) {
	print STDERR "\tOutput_quality_files => $opt_b.library_shortname.qual\n";
    }
} else {
     print STDERR "\tOutput_sequence_file => $opt_b.seq\n";
    if ($opt_Q) {
	print STDERR "\tOutput_quality_file => $opt_b.qual\n";
    }
}
if ($opt_C) {
    print STDERR "\tGet sequences trimmed => Enable\n";
} else {
    print STDERR "\tGet sequences trimmed => Disabled\n";
}

if ($opt_F) {
    print STDERR "\tGet only sequences with flags equals to 0\n";
} elsif ($opt_f) {
    print STDERR "\tGet only sequences with flags different to: $opt_f\n";
} else {
    print STDERR "\tGet sequences with any flags\n";
}

if ($opt_S) {
    print STDERR "\tGet only sequences with status equals to 0\n";
} elsif ($opt_s) {
    print STDERR "\tGet only sequences with status different to: $opt_s\n";
} else {
    print STDERR "\tGet seuqneces with any status\n";
}
print STDERR "\n\n";

my $output_fname = $basename . ".seq";
my $qual_fname = $basename . ".qual";

my @libraries = @$arrayref_lib_ids;
print "".(join ", ", @libraries)."\n";
if (@libraries==0) {
  print STDERR "No libraries selected.\n";
  exit -1;
}

my $lq = $dbh->prepare("SELECT library_shortname from library where library_id = ?");
my $cursorname = "cursor";

my $lib_query="DECLARE $cursorname CURSOR FOR SELECT clone.clone_name, direction, est.est_id, seq, qscore, hqi_start, hqi_length, flags, status from clone inner join seqread using (clone_id) inner join est using (read_id) inner join qc_report using (est_id) where clone.library_id=?";

if ($opt_F) {
    $lib_query .= " AND sgn.est.flags=0";
}
if ($opt_f) {
    $lib_query .= " AND sgn.est.flags NOT IN $opt_f";
}
if ($opt_S) {
    $lib_query .= " AND sgn.est.status=0";
}
if ($opt_s) {
    $lib_query .= " AND sgn.est.status NOT IN $opt_s";
}
print "\nQUERY:\n$lib_query\n";
my $sq = $dbh->prepare ($lib_query);
my $sth = $dbh->prepare ("FETCH 100 FROM $cursorname");

if (!$opt_L) {
    if ( -f $output_fname . ".seq" ) {
        print STDERR "\nError!!! Output file \"$output_fname\" already exists\n";
        exit (1);
    }

    if ( -f $qual_fname ) {
        print STDERR "\nError!!! Output quality file $qual_fname already exists\n";
        exit (1);
    }

    open SF, ">$output_fname"
	or die "\nError!!! Can't open sequence output file \"$output_fname\" ($!)";

    if ($opt_Q) {
	open QF, ">$qual_fname"
	    or die "\nError!!! Can't open quality output file \"$qual_fname\" ($!)";
    }
}

my ($n, $total) = (0,0);
foreach my $library_id ( @libraries ) {
  $lq->execute($library_id);
  if ($lq->rows == 0) {
    print STDERR "No library found for library_id \"$library_id\"\n";
    next;
  }
  
  my ($shortname) = $lq->fetchrow_array();

  if ($opt_L) { 
      $output_fname = $basename.'.'.$shortname.'.seq';
      if ( -f $output_fname . ".seq" ) {
        print STDERR "\nError!!! Output file \"$output_fname\" already exists\n";
        exit (1);
      }
      open SF, ">$output_fname"
	  or die "\nError!!! Can't open sequence output file \"$output_fname\" ($!)";

      if ($opt_Q) {
	  $qual_fname = $basename. '.'.$shortname.'.qual';
          if ( -f $qual_fname ) {
              print STDERR "\nError!!! Output quality file $qual_fname already exists\n";
              exit (1);
          }
	  open QF, ">$qual_fname"
	      or die "\nError!!! Can't open quality output file \"$qual_fname\" ($!)";
      }
      $total=0;
  }
  
  my ($orgs_nr, $organism);
  if ($opt_o) {
      my @orgs=split(/,/, $opt_o);
      $orgs_nr=scalar(@orgs);
  }
  my $query_org="SELECT organism_name FROM sgn.organism JOIN sgn.library USING(organism_id) WHERE library_id=?";
  my $sth_org=$dbh->prepare($query_org);
  $sth_org->execute($library_id);
  $organism=$sth_org->fetchrow_array();
  my $m=0;

  print STDERR "Getting sequences for the organism:$organism and library_shortname=$shortname...\n";
  $sq->execute($library_id);
  while ($sth->execute) {
    if ($sth->rows == 0) {
      last;
    }
    while(my ($clone_name, $dir, $est_id, $seq, $qscore, $hqi_start, $hqi_length, $flags, $status) =
	  $sth->fetchrow_array()) {
      $m++;    
      my $header;
      if ($opt_A || $orgs_nr > 1) {
	  $header= ">SGN-E$est_id #$organism [$clone_name]";
      } else {
	  $header= ">SGN-E$est_id [$clone_name]";
      }
      my @qual = split/\s+/,$qscore;

      if ($opt_C) {
          $seq = substr($seq,$hqi_start,$hqi_length);
          @qual = splice(@qual,$hqi_start,$hqi_length);
      }

      if ($opt_R) {
          if ($dir) {
              if ($dir eq "3") {
	          $seq =~ tr/ACGT/TGCA/;
 	          $seq = join("",reverse split//,$seq);
 	          @qual = reverse @qual;
 	          $header .= "(-)";
              }
          }
      } else {
	  if ($dir) {
	      if ($dir eq "3") {
		  $header .= "(3')";
	      } elsif ($dir eq "5") {
		  $header .= "(5')";
	      }
	  }
      }
      if ($seq) {
          print SF "$header\n";
          print SF "$seq\n";
          if ($opt_Q) {
              print QF "$header\n";
              print QF join(" ",@qual),"\n";
          }
      } else {
          if ($flags eq 0 && $status eq 0) {
              print "\tThis est_id: $est_id with status:$status and flags:$flags have not sequence.\n";
          }
      }
      
      $total++;
    }
    
  }
  $n += $total;
  if ($opt_L) {
      close SF;
      if ($opt_Q) {
          close QF;
      }
      print STDERR "\t... done\n\t=> There are $total sequences in the file:$output_fname.\n";
  } else {
      print STDERR "\t... done (Added $m sequences (total:$total) to the file:$output_fname)\n";
  }
  $dbh->do ("CLOSE $cursorname"); 
}

$sq->finish();

if (!$opt_L) {
    close SF;
    if ($opt_Q) {
	close QF;
    }
    print STDERR "\n\t=> There are $total sequences in the file:$output_fname.\n";
}
print STDERR "\n $n sequences has been processed by this script.\n\n";



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
      This script is a modification of the query-unigene-input.pl script that add some features like get the sequences by organism name or library shortname, get the sequences without trim or get all the sequences independently of the flags or status.
    
    Usage: 
      get_estseq_from_sgndb.pl -H <dbhost> -D <dbname> -b <basename> [-A OR -o <organism_names> [-L] OR -l <library_names>] [-C] [-F] [-S] [-Q] [-f <flags excluded>] [-s <status excluded>]

    Example (for crude sequences):
      get_estseq_from_sgndb.pl -H localhost -D sandbox -b nicotiana_seqs -o 'Nicotiana tabacum','Nicotiana benthamiana','Nicotiana sylvestris'

    Example (for curated sequences and quality files):
      get_estseq_from_sgndb.pl -H localhost -D sandbox -b nicotiana_seqs -o 'Nicotiana_tabacum','Nicotiana_sylvestris' -C -F -S -Q

    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -b basename             basename for the output files (mandatory)
      -A get ALL              get all the organism
      -o organism/s name/s    organism name (or organism names) of the sequences (with single quotes) 
                              separated by commas (,). This option does not get the libraries associated to the 
			      selected organism that have the type with the comment "demethylated genomic 
                              sequence library"
      -L organize by library  organize the sequence by library as multiple files			      
      -l library/ies name/s   library name (or libraries names) of the sequences (with single quotes)
                              separated by commas (,)
      -C enable clean option  Get the sequences according to the cleaning coordenates (hqi_start and 
                              hqi_length from qc_report)
      -F enable flags option  Get the sequences that have sgn.est.flags equal to 0.
      -f exclude some flags   Exclude some flags when get the sequences (example -f (64,65) )
      -S enable status option Get the sequences that have sgn.est.status value equal to 0
      -s exclude some status  Exclude some status when get the sequences (example -s (64) )
      -Q enable .qual option  Get the quality values asociated to a sequences in a fasta format.
      -R get reverse 3' seq   With this option you get the reverse of the 3' sequences.
      -h this help

EOF
exit (1);
}

=head2 validate_input

  Usage: my ($basename, $dbh, $arrayref_libraries_ids)=validate_input();
  Desc: check the input variables and return the list of the libraries.
  Ret: Three scalars, $basename (the basename for output files), $dbh  (database conection object) and 
       $arrayref_libraries_ids (array reference of the arguments)
  Args: none
  Side_Effects: if there are something wrong, die.
  Example: my ($basename, $dbh, $arrayref_libraries_ids)=validate_input();

=cut

sub validate_input {
  print "Check the variables...\n";
  if ($opt_h) {
    help ();
  }

  unless ($opt_b) {
    die ("Required argument -b <basename for output files> was not supplied.\n");
  }
  unless ($opt_o || $opt_l || $opt_A) {
    die ("Required argument -A (get_ALL_organism), -o <organism name> or -l <library_name> were not supplied.\n");
  }
  if ($opt_A && $opt_o) {
      die "Sorry, the arguments -A (get_ALL_organism) and -o <get_concrete_organism> are incompatible.\n";
  } elsif ($opt_A && $opt_l) {
      die "Sorry, the argument -A (get_ALL_organism) and -l <get_concrete_library> are incompatible.\n";
  } elsif ($opt_o && $opt_l) {
      die "Sorry, the argument -o (get concrete organism) and -l (get concrete library) are incompatible.\n";
  }
  if ($opt_L && $opt_l) {
      die "Sorry, the -L (get_sequences_by_library) is incompatiple with the -l <library> option\n";
  }
  if ($opt_F && $opt_f) {
      die "The option -F (get only sequences with flags 0) is incompatible with option -f (do not get the sequences with these flags\n";
  }
  if ($opt_S && $opt_s) {
      die "The option -S (get only sequences with status 0) is incompatible with option -s (do not get the sequences with these status";
  }

  if (!$opt_D || !$opt_H) { 
      die "\nneed -D (database name) and -H (host) options";
  }
  my $dbh = CXGN::DB::InsertDBH->new( { dbname => $opt_D,
				       dbhost => $opt_H,
				     } );

  my (@org_ids, @lib_ids);
  my ($organism_name, $library_shortname);
  if ($opt_A) {
      my $query="SELECT library_id FROM sgn.library JOIN sgn.organism USING(organism_id) WHERE type NOT IN (SELECT type_id 
                 FROM sgn.types WHERE comment= 'demethylated genomic sequence library') ORDER BY organism_name, library_shortname";
      my $sth=$dbh->prepare($query);
      $sth->execute();
      while (my $library_id=$sth->fetchrow_array()) {
		  push @lib_ids, $library_id;
      }
  } elsif ($opt_o) {
      my @organism_names=split(/,/, $opt_o);
      foreach $organism_name (@organism_names) {
	  print "Checking organism: $organism_name\t...\t";
	  my $query="SELECT organism_id FROM sgn.organism WHERE organism_name=?";
	  my $sth=$dbh->prepare($query);
	  $sth->execute($organism_name);
	  my $organism_id=$sth->fetchrow_array();
	  if ($organism_id) {
	      print "Ok\t...Organism_id: $organism_id.\n";
	      push @org_ids, $organism_id;
	      my $query2="SELECT library_id FROM sgn.library JOIN sgn.organism USING(organism_id) WHERE type NOT IN 
                          (SELECT type_id FROM sgn.types WHERE comment = 'demethylated genomic sequence library') AND organism_id = ? 
                          ORDER BY organism_name, library_shortname";
	      $sth=$dbh->prepare($query2);
	      $sth->execute($organism_id);
	      while (my $library_id=$sth->fetchrow_array()) {
		  push @lib_ids, $library_id;
	      }
	  } else {
	      print "WRONG!!!\t...this organism have not organism_id.\n";
	  }
      }
  }
  if ($opt_l) {
      my @library_shortnames=split(/,/, $opt_l);
      foreach $library_shortname (@library_shortnames) {
	  print "Checking library: $library_shortname\t...\t";
	  my $query="SELECT library_id FROM sgn.library WHERE library_shortname=?";
	  my $sth=$dbh->prepare($query);
	  $sth->execute($library_shortname);
	  my $library_id=$sth->fetchrow_array();
	  if ($library_id) {
	      print "Ok\t...Library_id: $library_id.\n";
	      push @lib_ids, $library_id;
	  } else {
	      print "WRONG!!!\t...this library have not library_id.\n";
	  }
      }
  }
  return ($opt_b, $dbh, \@lib_ids);
}

