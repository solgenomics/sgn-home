#!/usr/bin/perl

=head1 NAME

generate_pathologic_file.pl
create pathologic files from blast reports and annotation info

=head1 SYNOPSIS

generate_pathologic_file.pl -a <aracyc_dump_file> 
                            -s <swissprot_dat_file> 
                            -x <arabidopsis_blast_m8_bestmatch_with_annotations>
                            -y <swissprot_blast_m8_bestmatch_with_annotations>
                            -z <genbank_blast_m8_bestmatch_with_annotations>
                            -e <evalue_cutoff>
                            [-V] [-U] [-l <dblink_name>]

=head2 I<Flags:>

=over

=item -a

B<aracyc_dump_file>        AraCyc pathway file (mandatory)

=item -s

B<swissprot_dat_file>      Swissprot Dat file (optional)

=item -x

B<arabidopsis_blast>       Arabidopsis blast report m8 + annot. (mandatory)

=item -y

B<swissprot_blast>         Swissprot blast report m8 + annot. (optional)

=item -z

B<genbank_blast>           Genbank blast report m8 + annot. (optional)

=item -e

B<evalue_cutoff>           evalue cutoff for blast result parsing

=item -l

B<dblink_name>             dblink name to add the DBLINK field (optional)

=item -U

B<uppercase_id>            print IDs with uppercase (optional)

=item -V

B<be_verbose>              print more messages about parsing status

=item -h

B<print_help>              print this help


=head1 DESCRIPTION

    This script generate a file in a pathologic format used by PathwayTools 
 (http://bioinformatics.ai.sri.com/ptools/tpal.pf).

    As base file to transfer the metabolic annotations uses AraCyc dump file 
 (ftp://ftp.plantcyc.org/Pathways/aracyc_pathways.*). This file is a tabular
 format with 4 columns containing Arabidopsis pathways annotations.
   1. pathway
   2. reaction
   3. enzyme name
   4. Arabidopsis locus code

    As input file uses three different blast report files in m8 (tab) format
  with an additional column with annotations.
  (you can use 'blast_result_handle_tool.pl' script to generate these files)

    i-   Arabidopsis.blast.bestmatch.m8 
    ii-  Swissprot.bestmatch.m8  
    iii- GenBank.bestmatch.m8

    The script parses these files and search the annotation for each of them
  in the order: Arabidopsis ID, Swissprot E.C. code and finally GenBank
  description. It takes the Swissprot E.C. from swissprot_dat_file 
  (ftp://ftp.expasy.org/databases/swiss-prot/release/uniprot_sprot.dat.gz). 
  Blast hits with evalues higher than -e (or 1e-30) will be
  discarded.
    
    Once these information have transfered to the sequence match, it will print
  as STDOUT one entry per query_id such as:
 
   ID            query_id
   NAME          query_id
   FUNCTION      annotation_transfered
   DBLINK        -l:query_id (optional by option -l)
   PRODUCT-TYPE  P
   GENE-COMMENT  comment (annotation origin) 
   EC            ec
   //

=head1 AUTHORS

Created by:
Lukas Mueller <lam87@cornell.edu>

Refactored by:
Aureliano Bombarely <ab782@cornell.edu>

=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;

our ($opt_a, $opt_s, $opt_x, $opt_y, $opt_z, $opt_e, $opt_l, $opt_U, $opt_V, 
     $opt_h);
getopts("a:s:x:y:z:e:l:UVh");
if (!$opt_a && !$opt_s && !$opt_x && !$opt_y && !$opt_z && !$opt_e && 
    !$opt_U && !$opt_V && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}


## Check input files and define default values.
## Mandatory files.

my $aracyc_dump = $opt_a || 
die("INPUT ERROR: No -a <aracyc_dump> file was supplied.\n");

my $blast_at = $opt_x ||
die("INPUT ERROR: No -x <arabidopsis_blast> file was supplied.\n");

## Optional files

my $swpdat = $opt_s; 
my $blast_swp = $opt_y;
my $blast_gb = $opt_z;
my $dblink = $opt_l;

## Default values

my $evalue_cutoff = $opt_e || 1e-30;

## Variable declaration

my %ids = (); 
my %agi_annot =();
my %agi_ec_number = ();
my %swissprot_ec_number = ();
my %swissprot_evalue = ();
my %swissprot = ();
my %swissprot_match = ();
my %genbank = ();
my $ec_number; 


## FIRST STEP ########################################################

my $date = `date`;
chomp($date);
print STDERR "\n\n===================================================\n";
print STDERR "Initializing Script [$date].\n";
print STDERR "===================================================\n\n";

print STDERR "1) Parsing AraCyc dump file\n";
open (my $D, '<', $aracyc_dump) || die "Can't open $aracyc_dump.\n";

my $dl = 0;

while (<$D>) {
 
    chomp($_);
    $dl++;
    my ($pathway, $reaction, $protein_name, $AGI) = split(/\t/, $_);
    
    $AGI = uc($AGI);

    if ($opt_V) {
        print STDERR "\tProcessing line: $dl (ID:$AGI)     \r";
    }
    
    if ($AGI ne "UNKNOWN") { 
	$agi_annot{$AGI}=$protein_name; 
	if ($reaction =~ /\d+\.\d+\.\d+\.\d+/) {
            if (exists $agi_ec_number{$AGI}) { 
	        push @{$agi_ec_number{$AGI}}, $reaction;
            } 
            else {
                $agi_ec_number{$AGI} = [$reaction];
            } 
	}
    }
}

close($D);
print STDERR "\n\n";


## SECOND STEP ##########################################################

if (defined $swpdat) {

    print STDERR "2) Parsing swissprot dat file.\n";

    open (my $E, '<', $swpdat) || die "Can't open swissprot_dat_file.\n";
    my $sl = 0;

    ## Define the catching variables
    my $ec_number = "";
    my $acc = "";
  
    ## It will count as many IDs has the file to print the status
    my $swpid = 0;
    my $swpids = 0;
    if ($opt_V) {
        $swpids = `grep -c '^ID' $swpdat`;
        chomp($swpids);
    }
    
    while (<$E>) { 
    
         chomp($_);
         $sl++;

         if (/^ID\s+/) {
             $swpid++;
             if ($opt_V) {
                 print STDERR "\tProcessing $swpid of $swpids  \r";
             }
         }

         ## Accession catching
         if (/^AC\s+(.*)\;/) { 
	     $acc = $1;
         }

         my @accs = split /\;\s*/, $acc;
    
         ## EC number catching
         if (/^DE/) { 
             if (/EC=(\d+\.\d+\.\d+\.\d+)/) {
	         $ec_number = $1;	
             }
         }

         ## End of the section, pass variables to hash and reset them
         if (/^\/\//) { 

	    if ($ec_number) { 
	        foreach my $a (@accs) { 
                    if (exists $swissprot_ec_number{$a}) {
                        push @{$swissprot_ec_number{$a}}, $ec_number;
                    }
                    else {
		        $swissprot_ec_number{$a} = [$ec_number];
                    } 
	        }
	    }

            ## Reset the variables
	    $acc = "";
	    $ec_number = "";
         }
     }
     close($E);
}
else {
     print STDERR "2) No swissprot dat file was supplied. Skipping step 2.\n";
}
	
print STDERR "\n\n";    




## THIRD STEP ########################################################

my $old_id ="";
if (defined $blast_gb) {

    print STDERR "3) Parsing Genbank Blast.\n";

    open (my $G, '<', $blast_gb) || die "Can't open $blast_gb\n";
    my $gl = 0;  

    while (<$G>) { 
    
           chomp;
           $gl++;
    
           my ($id, $match, $evalue, $annot) = (split /\t/)[0,1,10,12];
           $ids{$id}=1;

           if ($opt_V) {
               print STDERR "\tProcessing line:$gl (id=$id, GB=$match)     \r";
           }    

           if ($id ne $old_id) { 
	       if ($evalue < $evalue_cutoff) { $genbank{$id}=$annot; }
           }
           $old_id = $id;
     }
     close($G);
}
else {
    print STDERR "3) No genbank blast file was supplied. Skipping Step 3.\n";
}
print STDERR "\n\n";


## FORTH STEP #############################################################

if (defined $blast_swp) {

    print STDERR "4) Parsing SwissProt Blast.\n";
    open (my $Q, '<', $blast_swp) || die "Can't open file $blast_swp.\n";
    my $bsl = 0;

    while(<$Q>) {
    
        chomp;
        $bsl++;
        my ($id, $match, $evalue, $annot) = (split /\t/)[0,1,10,12];

        if ($opt_V) {
            print STDERR "\tProcessing line:$bsl (ID=$id, SWP=$match).    \r";
        }

        $ids{$id}=1;
        if ($id ne $old_id) { 
	    if ($evalue < $evalue_cutoff) { 

 	        $match =~ s/sp\|(.*?)\|.*/$1/i;
	        $swissprot_match{$id} = $match;
                $swissprot_evalue{$id} = $evalue;

	        $swissprot{$id} = $annot;
	        $swissprot{$id} =~ s/(.+?)\s+OS\=.*/$1/; 
	    }
	}    
        $old_id= $id;
    }
    close($Q);
}
else {
    print STDERR "4) No Swissprot blast file was supplied. Skipping step 4.\n";
}
print STDERR "\n\n";


## FIFTH STEP ############################################################


print STDERR "5) Parsing Arabidopsis blast.\n";

open(my $F, '<', $blast_at) || 
     die "Can't open Arabidopsis blast file:$blast_at\n";

my $bal = 0;

my $agi_count = 0;
my $agiec_count = 0;
my $swp_count = 0;
my $swpec_count = 0;
my $gb_count = 0;

$old_id="";

while (<$F>) { 

    my @ec_number = ();
    chomp;
    $bal++;

    my ($id, $AGI, $evalue, $blast_annot) = (split /\t/)[0,1,10,12];   
    $ids{$id}++;

    if ($opt_V) {
        print STDERR "\tProcessing line:$bal (ID:$id, AGI:$AGI)           \r";
    }

    my $comment = "";

    if (($id ne $old_id) && ($evalue < $evalue_cutoff)) { 

	$id =~ s/\|//g;
	
	$AGI=uc($AGI);
	$AGI=~s/(AT\dG\d+)\.\d+/$1/;
	my $annotation = "";
	if (exists($agi_annot{$AGI})) { 
           
            if (exists $agi_ec_number{$AGI}) {
	        @ec_number = @{$agi_ec_number{$AGI}};           
                $agiec_count++;
            } 
	    $annotation = $agi_annot{$AGI};
	    $agi_count ++;
	    $comment = "similar to Arabidopsis $AGI [evalue: $evalue]";
	}	
	elsif (exists($swissprot{$id}))  {
	    
	    $annotation = $swissprot{$id};
	    if (exists($swissprot_ec_number{$swissprot_match{$id}})) { 
		@ec_number = @{$swissprot_ec_number{$swissprot_match{$id}}};
                $swpec_count++;
	    }
	    
	    $comment = "based on swissprot match ($swissprot_match{$id}";
            $comment .= " [evalue: $swissprot_evalue{$id}])";
            $swp_count++;
	}
	else { 
	    $annotation = $genbank{$id};
	    $comment = "based on genbank match.";
            $gb_count++;
	}
  
        if ($opt_U) {
            $id = uc($id);
        }
        print STDOUT "ID\t$id\nNAME\t$id\n";
        if (defined $annotation) {
	    print STDOUT "FUNCTION\t$annotation\n";
        }
        else {
            print STDOUT "FUNCTION\tunknown\n";
        }

        if (defined $dblink) {
            print STDOUT "DBLINK\t$dblink:$id\n";
        }

        print STDOUT "PRODUCT-TYPE\tP\nGENE-COMMENT\t$comment\n"; 

	if (scalar(@ec_number) > 0) {
              foreach my $ec (@ec_number) { 
                      print STDOUT "EC\t$ec\n";
              } 
        }
	print STDOUT "\/\/\n";

        $old_id = $id;
    
    }
}
close($F);

print STDERR "\n\n\n===================================================";
print STDERR "\nDone. Transferred:\n";
print STDERR "\tAraCyc annotations:\t$agi_count (with EC:$agiec_count)\n"; 
print STDERR "\tSwissprot annotations:\t$swp_count (with EC:$swpec_count)\n";
print STDERR "\tGenBank annotations:\t$gb_count\n";
print STDERR "===================================================\n\n";







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
       
      This script generate a file in a pathologic format used by PathwayTools 
     (http://bioinformatics.ai.sri.com/ptools/tpal.pf).

      As base file to transfer the metabolic annotations uses AraCyc dump file 
     (ftp://ftp.plantcyc.org/Pathways/aracyc_pathways.*). This file is a tabular
     format with 4 columns containing Arabidopsis pathways annotations.
        1. pathway
        2. reaction
        3. enzyme name
        4. Arabidopsis locus code

      As input file uses three different blast report files in m8 (tab) format
     with an additional column with annotations.
     (you can use 'blast_result_handle_tool.pl' script to generate these files)

        i-   Arabidopsis.blast.bestmatch.m8 
        ii-  Swissprot.bestmatch.m8  
        iii- GenBank.bestmatch.m8

      The script parses these files and search the annotation for each of them
    in the order: Arabidopsis ID, Swissprot E.C. code and finally GenBank
    description. It takes the Swissprot E.C. from swissprot_dat_file 
    (ftp://ftp.expasy.org/databases/swiss-prot/release/uniprot_sprot.dat.gz). 
    Blast hits with evalues higher than -e (or 1e-30) will be
    discarded.
    
      Once these information have transfered to the sequence match, it will 
    print as STDOUT one entry per query_id such as:
 
      ID            query_id
      NAME          query_id
      FUNCTION      annotation_transfered
      DBLINK        -l:query_id (optional by option -l)
      PRODUCT-TYPE  P
      GENE-COMMENT  comment (annotation origin) 
      EC            ec
      //

    Usage:
      
       generate_pathologic_file.pl -a <aracyc_dump_file> 
                            -s <swissprot_dat_file> 
                            -x <arabidopsis_blast_m8_bestmatch_with_annotations>
                            -y <swissprot_blast_m8_bestmatch_with_annotations>
                            -z <genbank_blast_m8_bestmatch_with_annotations>
                            -e <evalue_cutoff>
                            [-V] [-U] [-l <dblink_name>]
      
    Flags:
 
     -a <aracyc_dump_file>      AraCyc pathway file (mandatory)
     -s <swissprot_dat_file>    Swissprot Dat file (optional)
     -x <arabidopsis_blast>     Arabidopsis blast report m8 + annot. (mandatory)
     -y <swissprot_blast>       Swissprot blast report m8 + annot. (optional)
     -z <genbank_blast>         Genbank blast report m8 + annot. (optional)
     -e <evalue_cutoff>         evalue cutoff for blast parsing (1e-30 default)
     -l <dblink_name>           dblink name to add the DBLINK field (optional)
     -U <uppercase_ids>         print uppercase IDs (optional)
     -V <be_verbose>            print more messages about parsing status
     -h <print_help>            print this help

EOF
exit (1);
}


    
####
1; #
####
