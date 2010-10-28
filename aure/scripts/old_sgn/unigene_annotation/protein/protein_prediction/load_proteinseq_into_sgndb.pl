#!/usr/bin/perl

=head1 NAME

 load_proteinseq_into_sgndb.pl
 Load protein sequences into the SGN database

 (old version name: load_protein_data.pl)

=head1 USAGE

load_proteinseq_into_sgndb.pl -H <dbhost> -D <dbname> -c <cds file in fasta format> -p <protein file in fasta format> [-e | -l ]

=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory)

=item -c

B<cds file>             cds file in fasta format (It must have the same number and order in the unigene ids than -p) 
                        (mandatory)

=item -p

B<protein file>         protein file in fasta format (It must have the same number and order in the unigene ids than -c)
                        (mandatory)

=item -e

B<method='estscan'>     the protein prediction method was 'estscan'

=item -l

B<method='longest6frame'> the protein prediction method was 'longest6frame'

=item -T

B<enable test mode>     enable the test mode (and do not store the data in the database)

=item -h

B<help>                 print help


=head1 DESCRIPTION

This script load protein data for SGN unigenes into the sgn.cds table of the SGN database. Proteins can either be predicted with ESTScan or with the longest 6-frame translation. 

The ids in the cds and protein files need to correspond to SGN unigenes identifiers.

When loading ESTScan data, use the -e option, and ESTScan specific data will be parsed from the description line.

When loading longest six frame translations, use the -l option to parse the direction information that is appended to the identifier.

The corresponding sequences in the protein and cds files need to be in the file in the same order.

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=head1 MODIFIED BY

Aureliano Bombarely <ab782@cornell.edu>

=head1 METHODS

load_proteinseq_into_sgndb.pl

=cut

use strict;

use Getopt::Std;
use CXGN::DB::InsertDBH;
use CXGN::Transcript::Unigene;
use CXGN::Transcript::CDS;
use Bio::SeqIO;

our ($opt_H, $opt_D, $opt_c, $opt_p, $opt_e, $opt_l, $opt_T, $opt_h);
getopts('H:D:c:p:elTh');

if (!$opt_H && !$opt_D && !$opt_c && !$opt_p && !$opt_e && !$opt_l && !$opt_T && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my ($dbh, $cds_fasta, $protein_fasta)=validate_input();

my $last_run_id=get_last_run_id($dbh);
$last_run_id++;
my $current_run_id=$last_run_id;

my $last_sgncds_before_storage=get_last_number('sgn.cds_cds_id_seq');

my $io = Bio::SeqIO->new( -format=>"fasta", -file=>$protein_fasta);
my $cds_io = Bio::SeqIO->new( -format=>"fasta", -file=>$cds_fasta);

my $old_unigene;

my $entries_count=0;
my $store_count=0;
my $skipped_count = 0;
while (my $seq = $io->next_seq())  {
    $entries_count++;
    my $id = $seq->id();
    my $prot = $seq->seq();

    my $cds_obj = $cds_io ->next_seq();
    my $cds_id = $cds_obj->id();
    my $cds_seq = $cds_obj ->seq();
    my $cds_desc = $cds_obj ->desc();

    if ($cds_id ne $id ) { 
	print STDERR "Warning! CDS and protein ids don't match ($cds_id vs $id)\n";
    }

    if ($id !~ /SGN-U\d+/i) { die "Need a unigene id for each sequence\n"; }
#    my $cds = CXGN::Transcript::CDS->new_with_unigene_id($dbh, $unigene_id);
    my $cds = CXGN::Transcript::CDS->new($dbh);


    my $unigene_id;
    my $frame;
    my $score;
    my $seq_edits;

    my $forward_reverse = "F"; 
    

    if ($opt_e) { 

	if ($id=~/SGN-U(\d+)/i) { 
	    $unigene_id=$1;
	}
	else { 
	    print STDERR "Skipping $id -- illegal id format [SGN-U\d+].\n";
	    next();
	}

	my $unigene = CXGN::Transcript::Unigene->new($dbh, $unigene_id);
	if (!$unigene) { print STDERR "Skipping $unigene_id... Not found in database.\n"; 
			 $skipped_count++;
			 next();
		     }
	my $seq_text = $cds_seq;
	my $seq_edits = $cds_seq; 
	$seq_edits =~ s/[a-z]//g; # remove lower case characters- these represented deleted nucleotides.
	
	if ($cds_desc=~/minus strand/i) { 
	    $forward_reverse = "R";
	}

	my ($score, $start, $end, @rest) = split /\s+/, $cds_desc;
#	print STDERR "SCORE: $score, START: $start, END: $end\n";
	$cds ->set_unigene_id($unigene_id);
	$cds ->set_seq_text($seq_text);
	$cds ->set_seq_edits($seq_edits);
	$cds ->set_cds_seq($seq_edits);
	$cds ->set_direction($forward_reverse);
	$cds ->set_begin($start);
	$cds ->set_end($end);
	$cds ->set_score($score);
	$cds ->set_protein_seq($prot);
	$cds ->set_method('estscan');
	$cds ->set_run_id($current_run_id);
    }

    elsif ($opt_l) { 
	if ($id =~ /(\d+)([-+]\d)/) { 
	    $unigene_id = $1;
	    $frame = $2;
	}
	if (!$unigene_id) { die "Need a unigene id"; }

	my $direction = "F";
	$cds->set_protein_seq($prot);
	$cds->set_unigene_id($unigene_id);
	$cds->set_method("longest6frame");
	$cds->set_cds_seq($cds_seq);
	if ($frame <0) { 
	    $direction = "R";
	}
	$cds->set_direction($direction);
	$cds->set_frame($frame);
        $cds->set_run_id($current_run_id);
    }

    

    if (!$opt_T) { 
	print STDERR "Storing cds for $id...";
	my $new_id = $cds->store();
	#print STDERR "Done\n";
	print "[$new_id]              \n";
	if ($new_id) { 
	    $store_count++;
	}
	else { 
	    die "An error occurred during store() for seq $id.";
	}
    }
    else { print STDERR "test mode, not storing... SGN-U$unigene_id\n"; }
    $cds = undef;
    $old_unigene = $unigene_id;
}

my $last_sgncds_after_storage=get_last_number('sgn.cds_cds_id_seq');
my $cds_store=$last_sgncds_after_storage-$last_sgncds_before_storage;

if ($cds_store ne $store_count) {
    print "WARNING: The number of new cds in the sgn.table ($cds_store) is not the same than the protein storage count ($store_count).\n";
}

print "\n----------------------------------------\n";
print "REPORT:";
print "\n----------------------------------------\n";
print STDERR "Number of sequences for the input files: $entries_count\n";
print STDERR "Store $store_count proteins.\nSkipped $skipped_count proteins.\nRun_id= $current_run_id\n";
print "\n----------------------------------------\n\n";


if (!$opt_T) {
    print STDERR "Commit?\n";
    my $yes = <>;
    if ($yes !~ /n/i) { 
        print STDERR "Committing load...\n";
        $dbh->commit();
    } 
    else { 
        print STDERR "Rolling back...\n";
        $dbh->rollback();
    }
} else {
    $dbh->rollback();
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
      This script load the protein prediction sequences in the sgn.cds table. 
    
    Usage: 
      load_proteinseq_into_sgndb.pl -H <dbhost> -D <dbname> -c <cds_file> -p <protein_file> (-e OR -l) 
    
    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -c cds input file       cds input file in fasta format (mandatory)
      -p protein input file   protein input file in fasta format (mandatory)
      -e load 'estscan'       protein prediction method 'estscan'. This argument is incompatible with -l.
      -l load 'longest6frame' protein prediction method 'longest6frame'. This argument is incompatible with -e.
      -T test mode            the script DO NOT LOAD the sequences in the database.
      -h this help

EOF
exit (1);
}


=head2 validate_input

  Usage: my ($dbh, $cds_file, $protein_file)=validate_input();
  Desc: check if all the arguments are right
  Ret: $dbh (database conection object), $cds_file (cds file name) and $protein_file (protein file name)
  Args: none
  Side_Effects: die if there are something wrong. Connect with the database.
  Example: my ($dbh, $cds_file, $protein_file)=validate_input();

=cut
 
sub validate_input {
    if ($opt_h) {
	help();
    }
    print "Checking arguments...\n";
    unless ($opt_c) {
	die "Required argument -c <cds_input_file> was not supplied.\n";
    }
    unless ($opt_p) {
	die "Required argument -p <protein_input_file> was not supplied.\n";
    }
    unless ($opt_e || $opt_l) {
	die "Required argument -e <input files from 'estscan'> OR -l <input files from 'longest6frame'> was not supplied.\n";
    }
    if ($opt_e && $opt_l) {
	die "The use of the argument -e <input files from 'estscan'> exclude the use of -l <input files from 'longest6frame' and vice-versa.\n";
    }
    unless ($opt_H && $opt_D) {
	die "Database required arguments (-H <dbhost> and -D <dbname>) were not supplied.\n";
    }
    unless (-s $opt_c) {
	die "The input file (-c) $opt_c, do not exists or if it exists has zero size.\n ";
    }
    unless (-s $opt_p) {
	die "The input file (-p) $opt_p, do not exists or if it exists has zero size.\n";
    }
    my $cds_entries=`grep '>' $opt_c | wc -l`;
    chomp($cds_entries);
    my $protein_entries=`grep '>' $opt_p | wc -l`;
    chomp($protein_entries);

    if ($cds_entries ne $protein_entries) {
	die "cds file ($opt_c) and protein file ($opt_p) have a different entries number ($cds_entries and $protein_entries respectively.\n";
    }

    my $dbh = CXGN::DB::InsertDBH::connect({dbhost=>$opt_H, dbname=>$opt_D} );

    print "\nArguments checked, they are right\n";
    return ($dbh, $opt_c, $opt_p);
}

=head2 get_last_run_id 

  Usage: my $last_run_id=get_last_run_id($dbh);
  Desc: get the last run_id from the sgn.cds (the run_id with the highest value)
  Ret: A scalar, an integer ($last_run_id)
  Args: $dbh (database object connection)
  Side_Effects: none
  Example: my $last_run_id=get_last_run_id($dbh);

=cut

sub get_last_run_id {
    my $dbh=shift;
    my $last_run=0;

    my $query="SELECT DISTINCT(run_id) FROM sgn.cds";
    my $sth=$dbh->prepare($query);
    $sth->execute();
    while (my ($run_id)=$sth->fetchrow_array() ) {
	if ($run_id > $last_run) {
	    $last_run=$run_id;
	}
    }
    return $last_run;
}

=head2 get_last_number

  Usage: my $last_number=get_last_number($sequence_number_id);
  Desc: get the number of the last entry in a database table.
  Ret: a scalar, the last_number
  Args: $dbseq (the name of the sequence that control the counts in the table).
  Side_Effects: none.
  Example: my $preload_last_lib_id = get_last_number('sgn.library_library_id_seq');

=cut
 
sub get_last_number {
    my $dbseq = shift;
    my $query1 = "SELECT last_value FROM $dbseq";
    my $sth=$dbh->prepare($query1);
    $sth->execute();
    my ($last_count)=$sth->fetchrow_array();
    return $last_count;
}
