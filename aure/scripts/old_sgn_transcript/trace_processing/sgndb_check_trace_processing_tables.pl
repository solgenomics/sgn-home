#!/usr/bin/perl

=head1 NAME

 sgndb_check_trace_processing_tables.pl.
 A script to check the information store in the trace tables (version.1.0.).

=head1 SYPNOSIS
  
 sgndb_check_trace_processing_tables.pl [-h] -H <dbhost> -D <dbname> [-o <organism_name> OR -l <libraries_shortname>] 
 [-t <trace path to the chromatogram dir>]

=head1 DESCRIPTION

This script check different data in the trace tables. The trace tables are:
 
=over

=item - sgn.organism 

 This table store data about the organism (taxonomy lineage, genomic like chromosome number or polyploid...) 
 
=item - sgn.library

 It store data about the biological source of the sequences like tissue, development, or cultivar. Also store another information like how the library was made (kits used, adaptors, host bacter strain...) and authors. 
 The libraries are associated to an organism.

=item - sgn.clone

 It store the clone name, so means the name for group of identical cells that store the DNA sequenced. Sometimes you can find the chromatogram name of the DNA sequencing, the locus or the GenBank accession name of the sequence. 
 The clones are associated to a library (it have the biological information of the clone). One clone can have more than one read (it could have been sequenced more than one time).

=item - sgn.seqread

 It store the information about the sequencing process of the clone, so means that it store data like the submiter, the sequencing batch, the primer used, the direction of the sequencing, the trace name (the name of the read or if exists of the chromatogram file) and trace location (if exists chromatogram... where is it?). The reads are associated to a clone and have associated one sequence of sgn.est. 

=item - sgn.est

 It store the sequence, quality data (15 for each nucleotide by default), or the call positions of this sequence in the chromatogram (if exists). Also store quality tags of the sequences in flags and status, for example flags=1 means that the sequence is a vector contamination and is not a sequence of the genetic material that the library said.
 Sometime a est_id can have not a sequence (the data field can be empty... because the phred program could not calculate a sequence from the chromatogram).
 The est_id are associated to a read_id and have associated one qc_report. Also they can have associated a dbxref in the sgn.est_dbxref table.

=item - sgn.qc_report

 It store data about the cleaning process of the sequence (for example data produced by phred program when process the chromatogram like entropy or expected error of the base). Also store the coordenates for the sequence clean (sequence without adaptor, vector, poliA...). The qc_ids are associated to a sequence (est_id).

=item - sgn.est_dbxref and public.dbxref

 This table relate the table with the sequences (sgn.est) and the table with the references of external database. In the sequence case store the accession and the version of GenBank database for a concrete sequence.


=over

=head1 AUTHOR

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=head1 METHODS
 
sgndb_check_trace_processing_tables.pl


=cut

use strict;

use File::Basename;
use CXGN::DB::InsertDBH;
use Getopt::Std;


our ($opt_D, $opt_H, $opt_o, $opt_l, $opt_t, $opt_d, $opt_h);
getopts("D:H:o:l:t:d:h");

if (!$opt_D && !$opt_H && !$opt_o && !$opt_l && !$opt_t && !$opt_d && !$opt_h) {
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



my @lib_ids;
my ($lib_id, $clone_id, $read_id, $est_id, $qc_id);
our $outputdir;

if ($opt_d) {
    $outputdir=$opt_o;
} else {
    my $currentdir=`pwd`;
    chomp($currentdir);
    $outputdir=$currentdir."/trace_processing_check";
}
mkdir "$outputdir", 0755 or warn "Cannot make $outputdir directory: $!";

if ($opt_o && $opt_l) {
    die "Sorry you only can specify organism OR libraries but not both.\n";
} elsif ($opt_o && !$opt_l) {
    @lib_ids=check_organism();
} elsif ($opt_l && !$opt_o) {
    @lib_ids=check_library();
} elsif (!$opt_o && !$opt_l) {
    die "Required argument -o <organism name> or -l <library name> was not supplied.\n";
}

my $count_lib=scalar @lib_ids;

my $count_clone_ids;
my $count_read_ids;
my $count_est_ids;
my $count_qc_ids;
our ($count_problems_library, $count_problems1_clone, $count_problems2_clone, $count_problems1_read, 
     $count_problems2_read, $count_problems3_read, $count_problems4_read, $count_problems1_est, 
     $count_problems2_est, $count_problems3_est, $count_problems4_est, $count_problems5_est, $count_problems6_est, 
     $count_problems7_est, $count_problems_qcreport)=(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0);

my $outputfile_l1=$outputdir."/library_problem_type1";
open our $outputfile_handlel1, '>', "$outputfile_l1" or die "Cannot create the file $outputfile_l1: $!"; 
my $outputfile_c1=$outputdir."/clone_problem_type1";
open our $outputfile_handlec1, '>', "$outputfile_c1" or die "Cannot create the file $outputfile_c1: $!"; 
my $outputfile_cw=$outputdir."/clone_warning";
open our $outputfile_handlecw, '>', "$outputfile_cw" or die "Cannot create the file $outputfile_cw: $!"; 
my $outputfile_rw=$outputdir."/read_warning";
open our $outputfile_handlerw, '>', "$outputfile_rw" or die "Cannot create the file $outputfile_rw: $!";; 
my $outputfile_r3=$outputdir."/read_problem_type3";
open our $outputfile_handler3, '>', "$outputfile_r3" or die "Cannot create the file $outputfile_r3: $!"; 
my $outputfile_r1=$outputdir."/read_problem_type1";
open our $outputfile_handler1, '>', "$outputfile_r1" or die "Cannot create the file $outputfile_r1: $!";;  
my $outputfile_r2=$outputdir."/read_problem_type2";
open our $outputfile_handler2, '>', "$outputfile_r2" or die "Cannot create the file $outputfile_r2: $!";;
my $outputfile_e1=$outputdir."/est_problem_type1";
open our $outputfile_handlee1, '>', "$outputfile_e1" or die "Cannot create the file $outputfile_e1: $!"; 
my $outputfile_e2=$outputdir."/est_problem_type2";
open our $outputfile_handlee2, '>', "$outputfile_e2" or die "Cannot create the file $outputfile_e2: $!"; 
my $outputfile_e3=$outputdir."/est_problem_type3";
open our $outputfile_handlee3, '>', "$outputfile_e3" or die "Cannot create the file $outputfile_e3: $!"; 
my $outputfile_e4=$outputdir."/est_problem_type4";
open our $outputfile_handlee4, '>', "$outputfile_e4" or die "Cannot create the file $outputfile_e4: $!";; 
my $outputfile_e5=$outputdir."/est_problem_type5";
open our $outputfile_handlee5, '>', "$outputfile_e5" or die "Cannot create the file $outputfile_e5: $!"; 
my $outputfile_ew1=$outputdir."/est_warning_type1";
open our $outputfile_handleew1, '>', "$outputfile_ew1" or die "Cannot create the file $outputfile_ew1: $!";;  
my $outputfile_ew2=$outputdir."/est_warning_type2";
open our $outputfile_handleew2, '>', "$outputfile_ew2" or die "Cannot create the file $outputfile_ew2: $!";;
my $outputfile_q1=$outputdir."/qc_report_problem_type1";
open our $outputfile_handleq1, '>', "$outputfile_q1" or die "Cannot create the file $outputfile_q1: $!";;


foreach $lib_id (@lib_ids) {
   my @clone_ids=get_clone_id($lib_id);
   $count_clone_ids += scalar @clone_ids;
   foreach $clone_id (@clone_ids) {
      my @read_ids=get_seqread_id($clone_id);
      $count_read_ids += scalar @read_ids;
      foreach $read_id (@read_ids) {
         if ($opt_t) {
             check_trace_location();
         }
         my @est_ids=get_est_id($read_id);
         $count_est_ids += scalar @est_ids;
         foreach $est_id (@est_ids) {
            check_est_seq($est_id);
            my @qc_ids=get_qc_id ($est_id);
            $count_qc_ids += scalar (@qc_ids);
            foreach $qc_id (@qc_ids) {
               check_qc_report($qc_id);
            }
         }
      }
   }
}

print "\n------------------------------\n";
print "REPORT:";
print "\n------------------------------\n";
if ($opt_o) {
    print "\n\tOrganism: $opt_o\n";
} elsif ($opt_l) {
    print "\n\tLibraries: $opt_l\n";
}
print "\n\tLibraries analyzed: $count_lib\n";
print "\t\t=> major problem (type1, without clones): $count_problems_library\n";
print "\n\tClones analyzed: $count_clone_ids\n";
print "\t\t=> major problem (type1, without reads): $count_problems1_clone\n";
print "\t\t=> warning (with more than one read): $count_problems2_clone\n";
print "\n\tReads analyzed: $count_read_ids\n";
print "\t\t=> major problem (type1, without est): $count_problems3_read\n";
print "\t\t=> major problem (type2, with more than one est and flag and status equal 0): $count_problems4_read\n";
if ($opt_t) {    
    print "\t\t=> major problem (type3, trace location wrong): $count_problems2_read\n";
    print "\t\t=> warning (without trace name or trace location): $count_problems1_read\n";
}
print "\n\tEsts analyzed: $count_est_ids\n";
print "\t\t=> major problem (type1, without qc report): $count_problems6_est\n";
print "\t\t=> major problem (type2, with more than one qc report): $count_problems7_est\n";
print "\t\t=> major problem (type3, without sequence but flags and status equal 0): $count_problems1_est\n";
print "\t\t=> major problem (type4, the sequence length and the qscore length are different): $count_problems3_est\n";
print "\t\t=> major problem (type5, without flags or status): $count_problems5_est\n";
print "\t\t=> warning (type 1, sequence with length minor than 100 nt haven\'t flags and status not equal to 0): $count_problems4_est\n";
print "\t\t=> warning (type 2, without sequence and flags and status not equal to 0): $count_problems2_est\n";

print "\n\tQc reports analyzed: $count_qc_ids\n";
print "\t\t=> major problem (type1, without hqi_start or hqi_length): $count_problems_qcreport\n";

close $outputfile_handlel1; 
close $outputfile_handlec1; 
close $outputfile_handlecw; 
close $outputfile_handlerw; 
close $outputfile_handler3; 
close $outputfile_handler1;  
close $outputfile_handler2;
close $outputfile_handlee1; 
close $outputfile_handlee2; 
close $outputfile_handlee3; 
close $outputfile_handlee4; 
close $outputfile_handlee5; 
close $outputfile_handleew1;  
close $outputfile_handleew2;
close $outputfile_handleq1;

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
      This script check different data in the trace tables. The trace tables are:
       
         - sgn.organism, this table store data about the organism (taxonomy lineage, genomic 
           like chromosome number or polyploid...) 
        
         - sgn.library, it store data about the biological source of the sequences like tissue, 
           development, or cultivar. Also store another information like how the library was made 
           (kits used, adaptors, host bacter strain...) and authors. 
           The libraries are associated to an organism.

         - sgn.clone, it store the clone name, so means the name for group of identical cells that 
           store the DNA sequenced. Sometimes you can find the chromatogram name of the DNA sequencing, 
           the locus or the GenBank accession name of the sequence. The clones are associated to a library 
           (it have the biological information of the clone). One clone can have more than one read 
           (it could have been sequenced more than one time).

         - sgn.seqread, it store the information about the sequencing process of the clone, so means 
           that it store data like the submiter, the sequencing batch, the primer used, the direction 
           of the sequencing, the trace name (the name of the read or if exists of the chromatogram file)
           and trace location (if exists chromatogram... where is it?). The reads are associated to a clone
           and have associated one sequence of sgn.est. 

         - sgn.est, it store the sequence, quality data (15 for each nucleotide by default), or the call 
           positions of this sequence in the chromatogram (if exists). Also store quality tags of the 
           sequences in flags and status, for example flags=1 means that the sequence is a vector contamination 
           and is not a sequence of the genetic material that the library said.
           Sometime a est_id can have not a sequence (the data field can be empty... because the phred program
           could not calculate a sequence from the chromatogram).
           The est_id are associated to a read_id and have associated one qc_report. 
           Also they can have associated a dbxref in the sgn.est_dbxref table.

         - sgn.qc_report, it store data about the cleaning process of the sequence (for example data produced 
           by phred program when process the chromatogram like entropy or expected error of the base). Also 
           store the coordenates for the sequence clean (sequence without adaptor, vector, poliA...). 
           The qc_ids are associated to a sequence (est_id).

         - sgn.est_dbxref and public.dbxref, this table relate the table with the sequences (sgn.est) and the table
           with the references of external database. In the sequence case store the accession and the version of 
           GenBank database for a concrete sequence.

    Report:
      This script give the report in two different ways:
	   
         1- Screen message: 
	    - Libraries analyzed: NN
                => major problem (type1, without clones): XX
            - Clones analyzed: NN
                => major problem (type1, without reads): XX
                => warning (with more than one read): XX
            - Reads analyzed: NN
                => major problem (type1, without est): XX
                => major problem (type2, with more than one est and flag and status equal 0): XX
                => major problem (type3, trace location wrong): XX
                => warning (without trace name or trace location): XX
            - Ests analyzed: NN
                => major problem (type1, without qc report): XX
                => major problem (type2, with more than one qc report): XX
                => major problem (type3, without sequence but flags and status equal 0): XX
                => major problem (type4, the sequence length and the qscore length are different): XX
                => major problem (type5, without flags or status): XX
                => warning (type 1, sequence with length minor than 100 nt haven\'t flags and status not equal to 0): XX
                => warning (type 2, without sequence and flags and status not equal to 0): XX
            - Qc reports analyzed: NN
                => major problem (type1, without hqi_start or hqi_length): XX

         2- Output dir with the files that contains the list of library_shortname, SGN-C, SGN-T and 
            SGN-E with these problems. 

    Usage: 
      sgndb_check_trace_processing_tables.pl [-h] -H <dbhost> -D <dbname> [-o <organism_name> OR -l <libraries_shortname>] 
      [-t <trace path to the chromatogram dir>]
    
    Example:
       sgndb_check_trace_processing_tables.pl -H localhost -D sandbox -o 'Nicotiana tabacum'

    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -o organism name        the organism name inside single quotes ('xxxx')
      -l library name         the library (or libraries shortname, with commas, xxxx,wwww,yyyy)
      -o trace path           the complete path to the chromatograms directories.
      -h print this help

EOF
exit (1);
}

=head2 check_organism

  Usage: my @library_ids=check_organism();
  Desc: Check if exists a organism in the database and return the library id of the libraries with organism_id 
        equal to this organism.
  Ret: An array, with the library ids
  Args: None (use $opt_o that is our variable)
  Side_Effects: die if the organism do not exists in the database
  Example: my $lib_ids=check_organism();
 
=cut

sub check_organism {
    my $query="SELECT organism_id FROM sgn.organism WHERE organism_name = ?";
    my $sth=$dbh->prepare($query);
    $sth->execute($opt_o);
    my ($org_id) = $sth->fetchrow_array();

    if (!$org_id) {
	die "The organism specified ($opt_o) is not in the database ($opt_D).\n";
    }
    

    my @lib_ids;
    my $query2="SELECT library_id FROM sgn.library WHERE organism_id = ?";
    $sth=$dbh->prepare($query2);
    $sth->execute($org_id);
    while (my ($lib_id) = $sth->fetchrow_array()) {
	push @lib_ids, $lib_id;
    }

    return @lib_ids;
 }

=head2 check_library

  Usage: my @library_ids=check_library();
  Desc: Convert the -l argument into an array (by commas), check if exists the libraries and return an array 
        with their ids
  Ret: An array with the library ids
  Args: None (use $opt_l that is an our variable)
  Side_Effects: die if only there is one library and it is not in the database
  Example: my @lib_ids=check_library();
 
=cut
    
sub check_library {
    my @libraries=split(',', $opt_l);
    my $lib;
    my @lib_ids;

    foreach $lib (@libraries) {
	my $query="SELECT library_id FROM sgn.library WHERE library_shortname = ?";
	my $sth=$dbh->prepare($query);
	$sth->execute($lib);
	my $lib_id=$sth->fechrow_array();
	push @lib_ids, $lib_id;
    }
    my $count=scalar @lib_ids;
    if ($count < 1) {
	die "The library (or libraries) specified (@libraries) is (are) not in the database ($opt_D).\n";
    }
    return @lib_ids;
}

=head2 get_clone_id

  Usage: my @clone_ids=get_clone_id($library_id);
  Desc: Get the clone id for a concrete library id
  Ret: An array with the clone ids of the specified library
  Args: $library_id (the library id)
  Side_Effects: Print a message in the screen (and the library shortname in a file) if the library have not clone ids
                associated.
  Example: my @clone_ids=get_clone_id($lib_id);
 
=cut

sub get_clone_id {
    my $lib_id=shift;
    my @clone_ids;

    my $query="SELECT clone_id FROM sgn.clone WHERE library_id = ? ORDER BY clone_id";
    my $sth=$dbh->prepare($query);
    $sth->execute($lib_id);
    while (my ($clone_id) = $sth->fetchrow_array()) {
	push @clone_ids, $clone_id;
    }
    my $count=scalar @clone_ids;
    if ($count < 1) {
        my $query2="SELECT library_shortname FROM sgn.library WHERE library_id = ?";
        $sth=$dbh->prepare($query2);
        $sth->execute($lib_id);
	my ($library_shortname) = $sth->fetchrow_array();
        print $outputfile_handlel1 "$library_shortname\n";
	print "The library $library_shortname ($lib_id) have not any clone associated.\n";
        $count_problems_library++;
    }
    return @clone_ids;
}

=head2 get_seqread_id

  Usage: my @read_ids=get_seqread_id($clone_id);
  Desc: Get the read ids for the clone id specified
  Ret: An array with the read ids associated to a clone id
  Args: $clone_id (clone id)
  Side_Effects: Print a message in the screen (and the reference in a file) if the clone have not associated reads or 
                have more than one.
  Example: my @read_ids=get_seqread_id($clone_id);
 
=cut

sub get_seqread_id {
    my $clone_id=shift;
    my @read_ids;

    my $query="SELECT read_id FROM sgn.seqread WHERE clone_id = ? ORDER BY clone_id";
    my $sth=$dbh->prepare($query);
    $sth->execute($clone_id);
    while (my ($read_id) = $sth->fetchrow_array()) {
	push @read_ids, $read_id;
    }
    my $count=scalar @read_ids;
    my $query2="SELECT clone_name FROM sgn.clone WHERE clone_id = ?";
    $sth=$dbh->prepare($query2);
    $sth->execute($clone_id);
    my ($clone_name) = $sth->fetchrow_array();

    if ($count < 1) {
        print $outputfile_handlec1 "SGN-C$clone_id";
	print "The clone $clone_name (SGN-C$clone_id) have not any read associated.\n";
        $count_problems1_clone++;
    } elsif ($count > 1) {
        print $outputfile_handlecw "SGN-C$clone_id\t$count\n";
	print "The clone $clone_name (SGN-C$clone_id) have more than one read associated.\n";
        $count_problems2_clone++;
    }
    return @read_ids;
}

 # NOTE ABOUT check_trace_location:
 # --------------------------------
 # In the est.pl script the full_pathname use CXGN::VHost to get the complete pathway. 
 # It depends of the enviroment variables, but you can not find these configuration en all the computers.
 # I have use $opt_t to detail these path. For one hand if exists $opt_t run this subroutine, for other hand
 # is most easy to check where are the files if you don't know what configuration have in your computer, but
 # you know where is the chromatograms (/data/prod/public/chromatograms)

=head2 check_trace_location

  Usage: check_trace_location($read_id);
  Desc: check if the chromatogram file have a right trace location and if it exists.
  Ret: none
  Args: $read_id
  Side_Effects: if the check is wrong, increase the global error count in one and print the SGN id (SGN-T) in a file.
  Example: check_trace_location($read_id)
 
=cut

sub check_trace_location {
    my $read_id=shift;

    my $query="SELECT trace_name, trace_location FROM sgn.seqread WHERE read_id = ? ORDER BY read_id";
    my $sth=$dbh->prepare($query);
    $sth->execute($read_id);
    my ($trace_name, $trace_location)=$sth->fetchrow_array();

    if (!$trace_name || !$trace_location) {
        print $outputfile_handlerw "SGN-T$read_id\n";
        print "\tThe trace SGN-T$read_id have not trace name ($trace_name) or trace location ($trace_location).\n";
        $count_problems1_read++;
    } else {
        my $full_pathname=$opt_t.$trace_location.$trace_name;
        unless (-f $full_pathname) {
            print $outputfile_handler3 "SGN-T$read_id\n";
	    print "The trace SGN-T$read_id have trace name and trace location but there is not the file $trace_name in the location $trace_location.\n";
	    $count_problems2_read++;
	}
    }
}

 # NOTE ABOUT get_est_id:
 # ----------------------
 # The rule is that one read_id only can have associated one est_id with flags and status = 0, if
 # it find more show like an error.

=head2 get_est_id

  Usage: my @est_ids=get_est_id($read_id);
  Desc: get the est_id associate to an read id and check if it have more than one est ids. 
  Ret: An array with the est ids assocaited to a read
  Args: $read_id
  Side_Effects: if it find that a read have not est id associate or have more than one increase the error count by one and
                write the read id in a error file
  Example: my @est_ids=get_est_id($read_id);
 
=cut
 
sub get_est_id {
    my $read_id=shift;
    my (@est_ids, @status, @flags);
    my $status_zero_count=0;
    my $flags_zero_count=0;
 

    my $query="SELECT est_id, status, flags FROM sgn.est WHERE read_id = ? ORDER BY est_id";
    my $sth=$dbh->prepare($query);
    $sth->execute($read_id);
    while (my ($est_id, $status, $flags)=$sth->fetchrow_array()) {
	push @est_ids, $est_id;
	if ($status && $status eq 0) {
	    $status_zero_count++;
	}
	if ($flags && $flags eq 0) {
	    $flags_zero_count++;
	}
    }

    my $count=scalar @est_ids;
    if ($count < 1) {
        print $outputfile_handler1 "SGN-T$read_id\n";
	print "The read SGN-T$read_id have not any est associated.\n";
        $count_problems3_read++;
    } elsif ($count > 1 && $status_zero_count < 1 && $flags_zero_count < 1) {
        print $outputfile_handler2 "SGN-T$read_id\t$count\n";
	print "The read SGN-T$read_id have more than one est associated with flags and status equal 0.\n";
        $count_problems4_read++;
    }
    return @est_ids;
}

 # NOTE ABOUT check_est_seq:
 # ------------------------
 #   In a est entries it check if the est have sequence. Sometimes you don't find it because the phred program
 # have not find a good sequence in the chromatogram (with the quality values used). It is not wrong, but the 
 # flags and the status must be different than 0 (for some scripts it is better if find one N instead the empty
 # field, the SGN_phred script (version 2.0) do it). 
 #   The current cleaning criterion is that the sequences with less than 100 nucleotides are not good (for some 
 # process like unigene build, so it is better if put a flags and a status that reflect that), but these can change
 # in a future (for example if you use this script to check 454 sequences in the database).

=head2 check_est_seq

  Usage: check_est_seq($est_id);
  Desc: check different aspect of the est table like if the est have sequence or flag or status, 
        or the qscore have different length than the seq.
  Ret: none
  Args: $est_id
  Side_Effects: Increase the global error count if find one error and write the est id (SGN-E) in a file
  Example: check_est_seq($est_id);
 
=cut

sub check_est_seq {
    my $est_id=shift;
    
    my $query="SELECT seq, qscore, status, flags FROM sgn.est WHERE est_id = ? ORDER BY est_id";
    my $sth=$dbh->prepare($query);
    $sth->execute($est_id);
    my ($seq, $qscore, $status, $flags)=$sth->fetchrow_array();
   
    unless (defined ($flags && $status) ) {
        print $outputfile_handlee5 "SGN-E$est_id\n";
	print "The SGN-E$est_id have not flags or status.\n";
        $count_problems5_est++;
    } else {
        if (!$seq && $status eq 0 && $flags eq 0) {
            print $outputfile_handlee3 "SGN-E$est_id\n";
	    print "The SGN-E$est_id have not sequence but the flags and the status have 0 value.\n";
            $count_problems1_est++;
        } elsif (!$seq) {
            print $outputfile_handleew2 "SGN_E$est_id\n";
	    print "The SGN-E$est_id have not sequence (the flags and the status are rigth).\n";
            $count_problems2_est++;
        } else {
            my $lenseq=length($seq);
            my $prelength_qscore=$qscore;
            $prelength_qscore =~ s/\d+\s*/n/gi;
            my $lenqscore=length($prelength_qscore);

            if ($lenseq ne $lenqscore) {
                print $outputfile_handlee3 "SGN-E$est_id\n";
	        print "The SGN-E$est_id have a qscore with a different length than the sequence.\n";
                $count_problems3_est++;
            }
            if ($lenseq < 100 && $status eq 0 && $flags eq 0) {
                print $outputfile_handleew1 "SGN-E$est_id\n";
	        print "The SGN-E$est_id have a length sequence with less than 100 nt but the flags and the status have 0 value.\n";
            $count_problems4_est++;
            }
	}
    }
}

=head2 get_qc_id

  Usage: my @qc_ids=get_qc_id($est_id)
  Desc: get the qc_ids associate to est_id.
  Ret: An array with the qc_ids associated to est id
  Args: $est_id
  Side_Effects: If find that an est have not qc_id or have more than one, increase the error count by one and write the
                est id in a file. 
  Example: my @qc_ids=get_qc_id($est_id)
 
=cut

sub get_qc_id {
    my $est_id=shift;
    my @qc_ids;

    my $query="SELECT qc_id FROM sgn.qc_report WHERE est_id = ? ORDER BY est_id";
    my $sth=$dbh->prepare($query);
    $sth->execute($est_id);
    while (my ($qc_id) = $sth->fetchrow_array()) {
	push @qc_ids, $qc_id;
    }
    my $count=scalar @qc_ids;
    
    if ($count < 1) {
        print $outputfile_handlee1 "SGN-E$est_id\n";
	print "The est SGN-E$est_id have not any qc report associated.\n";
        $count_problems6_est++;
    } elsif ($count > 1) {
        print $outputfile_handlee2 "SGN-E$est_id\t$count\n";
	print "The est SGN-E$est_id have more than one qc report associated.\n";
        $count_problems7_est++;
    }
    return @qc_ids;
}

=head2 check_qc_report

  Usage: check_qc_report($qc_id)
  Desc: check if a qc report entry have not any hqi_start or hqi_length when sgn.est.status and sgn.est.flags are
        equal to zero.
  Ret: none 
  Args: $qc_id
  Side_Effects: If it do not find any hqi_start or hqi_length increase the error count by one and write the ids
                (SGN-E) in a file.
  Example: check_qc_report($qc_id)
 
=cut

sub check_qc_report {
    my $qc_id=shift;
 
    my $query="SELECT hqi_start, hqi_length, flags, status, sgn.est.est_id FROM sgn.qc_report JOIN sgn.est
               ON sgn.qc_report.est_id=sgn.est.est_id WHERE qc_id=? ORDER BY sgn.est.est_id";
    my $sth=$dbh->prepare($query);
    $sth->execute($qc_id);
    my ($hqi_start, $hqi_length, $flags, $status, $est_id)=$sth->fetchrow_array();

    unless (defined ($hqi_start & $hqi_length) ) {
	if ($flags eq 0 && $status eq 0) {
            print $outputfile_handleq1 "SGN-E$est_id\n";
	    print "The qc_report associated to SGN-E$est_id (with flags and status eq 0) have not hqi_start or hqi_length.\n";
            $count_problems_qcreport++;
	}
    }
}
