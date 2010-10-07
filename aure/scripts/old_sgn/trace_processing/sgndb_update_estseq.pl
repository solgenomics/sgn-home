#!/usr/bin/perl

=head1 NAME

 sgndb_update_estseq.pl
 A script to reload (update) sequences and qscores for est's.

=head1 SYNOPSIS

 sgndb_update_estseq.pl [-h] -H <database_host> -D <database_name> -s <sequence_file> -q <qscore_file> -i <identifier_type> 

 where:
   -h help
   -H database_host (mandatory)
   -D database_name (mandatory)
   -s sequence_file in fasta format
   -q qscore_file in qual format
   -i identifier_type (mandatory)

=head1 DESCRIPTION

 This script reload the sequences and the qscores for the EST sequences based in the sequence identifier.

 It don't load the sequences that don't map with the identifier of the database.

 Identifiers can be:
  - Genbank_accessions, for example: BQ448782,
  - Sgn_est_accessions, for example: SGN-E1312559 
  - Clone_names, for example: CAEST134 
 

=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=cut 


use strict;
use warnings;

use Getopt::Std;
use File::Basename;
use CXGN::DB::InsertDBH;
use Bio::SeqIO;


our ($opt_D, $opt_H, $opt_s, $opt_q, $opt_i, $opt_h);
getopts("D:H:s:q:i:h");

## Print help if none argument or -h argument is supplied

if (!$opt_D && !$opt_H && !$opt_s && !$opt_q && !$opt_i && !$opt_h) {
    print "\nNone argument were supplied. Printing help\n\n";
    help();
}
if ($opt_h) {
   help();
}

## Check that all the mandatory arguments were supplied.

if (!$opt_D || !$opt_H) { 
    die "\nDatabase arguments were not supplied: -D <database_name> and/or -H <database_host>.\n";
}

if (!$opt_i) {
    die "\nIdentifier argument was not supplied: -i <identifier_type> (sgn_est_accession, genbank_accession or clone_name).\n";
}

if (!$opt_s && !$opt_q) {
    die "\nSequence or Qscore file arguments weren't supplied: please supply at least one of them, -s <sequence_file> -q <qscore_file>\n";
}

## Now open a database connection

my $dbh =  CXGN::DB::InsertDBH->new( 
                                     { 
                                       dbname => $opt_D, 
                                       dbhost => $opt_H, 
                                     }
                                   );



## Now we are going to open the sequence or the qscore file. 
## and put them in two hashes

my (%sequence, %qscore);

if (defined $opt_s) {
    my $input_sequence  = Bio::SeqIO->new(  
                                            -file   =>      $opt_s,
                                            -format =>      'fasta'
                                         );

    my $n = 0;
    my $N = `grep -c '>' $opt_s`;
    chomp($N);

    while (my $seq_obj = $input_sequence->next_seq() ) {
          
           my $id = $seq_obj->display_id;
           my $seq = $seq_obj->seq;
           
           unless (exists $sequence{$id}) {
               $id =~ s/SGN-E//;
               $sequence{$id} = $seq;
           }
           else {
               print STDERR "\nThe id=$id is redundant. Appears more than one time in the file=$opt_s.\n";
           }
           $n++;
           print STDERR "Processing the sequence:$id ($n of $N)                     \r";
    }
}

print STDERR "\n... done...\n";

if (defined $opt_q) {
    my $input_qscore  = Bio::SeqIO->new(  
                                            -file   =>      $opt_q,
                                            -format =>      'qual'
                                        );

    my $m = 0;
    my $M = `grep -c '>' $opt_q`;
    chomp($M);

    while (my $seq_obj = $input_qscore->next_seq() ) {
          
           my $id = $seq_obj->display_id;
           my $qual_aref = $seq_obj->qual;
           
           unless (exists $qscore{$id}) {
               $id =~ s/SGN-E//;
               $qscore{$id} = join(' ', @{$qual_aref});
           }
           else {
               print STDERR "\nThe id=$id is redundant. Appears more than one time in the file=$opt_q.\n";
           }
           $m++;
           print STDERR "Processing the sequence:$id ($m of $M)                        \r";
    }
}

print STDERR "\n... done...\n";

## Now it will update the sequences and the quality.

my @ids = ();

if (defined $opt_s) {
    @ids = keys %sequence;
}
elsif (defined $opt_q) {
    @ids = keys %qscore;
}

my ($c, $f, $u) = (0, 0, 0);
my $C = scalar(@ids);

foreach my $upd_id (@ids) {       
    
    $c++;
    my $query;
    if ($opt_i eq 'sgn_est_accession') {
        $upd_id =~ s/SGN-E//;

        $query = "SELECT sgn.est.est_id FROM sgn.est WHERE est_id = ?";    
    }
    elsif ($opt_i eq 'genbank_accession') {
        $query = "SELECT sgn.est.est_id FROM sgn.est JOIN sgn.est_dbxref USING(est_id) 
                                                        JOIN public.dbxref USING(dbxref_id) 
                                                        WHERE accession = ?";
    }
    elsif ($opt_i eq 'clone_name') {
        $query = "SELECT sgn.est.est_id FROM sgn.est JOIN sgn.seqread USING(read_id) 
                                                        JOIN sgn.clone USING(clone_id) 
                                                        WHERE clone_name = ?";
    }
    my ($est_id) = $dbh->selectrow_array($query, undef, $upd_id);

    unless (defined $est_id) {
        $f++;
        print STDERR "\nThe id=$upd_id is not in the database.\n";
    }
    else {
	my $update;
        if (defined $sequence{$upd_id} && defined $qscore{$upd_id} ) {
            $update = "UPDATE sgn.est SET seq = '$sequence{$upd_id}', qscore = '$qscore{$upd_id}'  WHERE est_id = $est_id";
            $u++;
        }
        elsif (defined $sequence{$upd_id}) {
            $update = "UPDATE sgn.est SET seq = '$sequence{$upd_id}' WHERE est_id = $est_id";
            $u++;   
        } 
        elsif (defined $qscore{$upd_id}) {
            $update = "UPDATE sgn.est SET qscore = '$qscore{$upd_id}' WHERE est_id = $est_id";
            $u++;   
        } 
        else {
            $f++;
            print STDERR "\nThe id=$upd_id have not defined any sequence or qscore. Check the files.\n";
        }
	$dbh->do($update);
        print STDERR "Updating the id=$upd_id\t($c of $C)\t(ids updated=$u and ids failed=$f)                  \r";
    }
}

print STDERR "\n\n...done...\n\n";

print STDERR "Commit?\n(yes|no, default no)> ";
if (<STDIN> =~ m/^y(es)/i) {
    print STDERR "Committing...";
    $dbh->commit;
    print STDERR "okay.\n";
} 
else {
    print STDERR "Rolling back...";
    $dbh->rollback;
    print STDERR "done.\n";
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
      This script reload the sequences and the qscores for the EST sequences 
    based in the sequence identifier.

      It do not load the sequences that do not map with the identifier of 
    the database.

      Identifiers can be:
        - Genbank_accessions, for example: BQ448782,
        - Sgn_est_accessions, for example: SGN-E1312559 
        - Clone_names, for example: CAEST134       
    
    Usage:      
      sgndb_update_estseq.pl [-h] -H <database_host> -D <database_name> 
                                  -s <sequence_file> -q <qscore_file> -i <identifier_type>      
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -s sequence_file        sequence file in fasta format
      -q qscore_file          quality file in qual format
      -i identifier_type      identifier type (mandatory), as sgn_est_accession,
                              or genbank_accession or clone_name

EOF
exit (1);
}
