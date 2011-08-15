#!/usr/bin/perl

=head1 NAME

 embl2tab.pl
 Script to extract annotations from uniprot files.

=cut

=head1 SYPNOSIS

 uniprot2tab.pl [-h] -u <uniprot_file> -o <output_basename> -f <annot_fields>
                     [-i <id_list>] [-s <species_list>] [-A]
 
=head2 I<Flags:>

=over


=item -u

B<uniprot_file>           uniprot (swissprot) file (mandatory)

=item -o

B<output_basename>        output basename (by default 'out')

=item -f

B<annot_fields>           annotation fields list separated by commas (mandatory)

=item -i

B<id_list_file>           file with the list of ids to retrieve (optional)

=item -s

B<species_list>           specieslist to retrieve separated by commas (optional)

=item -A

B<use_accessions>         use accessions instead IDs

=item -V

B<be_verbose>             verbose (print status messages)

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This extract the annotations for a specific field, for example if 'Description'
 it will take all the description lines for each of the IDS.

 The fields can be choosen in the following way:
  * CATEGORY, two letters, if this is not followed by '=' it will take the
    complete line. More than one category can be specified separated by
    comma ','.
  * SUBCATEGORY, a word, it will take from Word= (or :) to ';' or end of line.
    more than one subcategory can be used under the same category separated
    by semicolon ';'.

   Example1:

   uniprot2tab.pl -u uniprot_sprot.dat -f DE

   Example2:

   uniprot2tab.pl -u uniprot_sprot.dat -f DE=AltName

   Example3:

   uniprot2tab.pl -u uniprot_sprot.dat -f 'DE=AltName;EC,DT'

 More than one file can be used.

   Example4: 

   uniprot2tab.pl -u uniprot_sprot.dat,uniprot_trembl.dat -f DE

 Specific IDs and/or species can be selected.

   Example5:

   uniprot2tab.pl -u uniprot_sprot.dat -f DE -i ids.txt 
                  -s 'Arabidopsis thaliana'
 

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 uniprot2tab.pl


=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::SeqIO;

our ($opt_u, $opt_o, $opt_f, $opt_i, $opt_s, $opt_A, $opt_V, $opt_h);
getopts("u:o:f:i:s:AVh");
if (!$opt_u && !$opt_o && !$opt_f && !$opt_i && !$opt_s && !$opt_A && !$opt_V 
    && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Check variables.

my $uniprot_files = $opt_u 
    || die("DATA ARGUMENT -u <uniprot_files> WAS NOT SUPPLIED.\n");
my @unipfiles = split(/,/, $uniprot_files);

my $annot_fields = $opt_f
    || die("DATA ARGUMENT -f <annot_fields> WAS NOT SUPPLIED.\n");
my @afields = split(/,/, $annot_fields);

my $outbase = $opt_o || 'out.uniprot2tab.txt';

my $idfile = $opt_i;

my @species = ();
if ($opt_s) {
    @species = split(/,/, $opt_s);
}


## Define the variable to store the entries

my %entries = ();

my $date = `date`;
chomp($date);

print STDERR "\n\n===========================================================";
print STDERR "\nuniprot2tab.pl initiation [$date]\n";
print STDERR "===========================================================\n\n";

print STDERR "0) Running parameters:\n\n";

## LOAD THE FIELDS INTO A HASH

my %fields = ();
foreach my $fd (@afields) {

    if ($fd =~ m/^(.+?)=(.+)$/) {
	my $cat = $1;
	my $subcat = $2;
	my @subcat = split(/;/, $subcat);
	
	foreach my $scat (@subcat) {
	    print STDERR "\tFIELD (SUBCATEGORY) SET: $cat\t$scat\n";
	}
	
	$fields{$cat} = \@subcat;
    }
    else { 
	print STDERR "\tFIELD (CATEGORY) SET: $fd\n";
	$fields{$fd} = [];
    }
}
print STDERR "\n";

## LOAD THE ORGANISM INTO A HASH ###############

my %species = ();
foreach my $sp (@species) {
    print STDERR "\tSPECIES FILTER ENABLED FOR: $sp\n";
    $species{$sp} = 1;
}


## PARSING ID FILE ###############

my %ids = ();

print STDERR "\n\n1) Parsing ID file.\n\n";

if (defined $idfile) {

    open my $ifh, '<', $idfile;
    my $iline = 0;

    while(<$ifh>) {
	chomp($_);
	$iline++;

	if ($opt_V) {
	    print STDERR "\tParsing line:$iline from -i $idfile.           \r";
	}

	if ($_ !~ m/^#/) {
	    my @data = split(/\t/, $_);
	    $ids{$data[0]} = 1;
	}
    }
}
else {
    print STDERR "\tNo -i <id_list_file> was supplied. Skipping step 1.";
}


## PARSING uniprot FILE #############

print STDERR "\n\n2) Parsing Uniprot files.\n\n";

foreach my $upfile (@unipfiles) {

    open my $ufh, '<', $upfile;
    my $l = 0;

    ## Define the unique variables
    
    my @selected = ();
    my ($id, $acc, $sps, $taxon);

    while(<$ufh>) {
	chomp($_);
	$l++;

	## Get the tagged lines
	if ($_ =~ m/^(\w\w)\s+(.+)$/) {
	    my $tag = $1;
	    my $line = $2;

	    if ($tag eq 'ID') {
		if ($line =~ m/(.+?)\s+/) {
		    $id = $1;
		}
	    }
	    elsif ($tag eq 'AC') {
		if ($line =~ m/(.+?);/) {
		    $acc = $1;
		}
	    }
	    elsif ($tag eq 'OS') {
		$sps .= $line;
	    }
	    elsif ($tag eq 'OC') {
		$taxon .= $line;
	    }

	    if (exists $fields{$tag}) {
		
		if (scalar(@{$fields{$tag}}) > 0) {
		    
		    ## Break the line using ';'
		    my @line = split(/;/, $line);

		    ## Take all the fields
		    my @fields = @{$fields{$tag}};
		    
		    foreach my $fd (@fields) {
			foreach my $subline (@line) {
	
			    if ($subline =~ m/$fd=(.+)$/) {
				my $sel = $1;
				$sel =~ s/\.$//;
				## Remove the last dot if exists
				push @selected, $sel;
			    }
			    elsif ($subline =~ m/$fd:(.+)$/) {
				my $sel = $1;
				$sel =~ s/\.$//;
				push @selected, $1;
			    }
			}
		    }		    
		}
		else {
		    $line =~ s/\.$//;
		    push @selected, $line;
		}
	    }
	}
	elsif ($_ =~ m/^\/\//) {

	    ## Process species and taxon.

	    my @taxon = split(/;\s*/, $taxon);
	    
	    ## Clean the species line
	    $sps =~ s/\(.+$//;
	    
	    push @taxon, $sps;
	    my @clean_taxon = ();
	    foreach my $tx (@taxon) {
		$tx =~ s/\s+$//;
		$tx =~ s/^\s+//;
		$tx =~ s/\.$//;
		push @clean_taxon, $tx;
	    }
	    @taxon = @clean_taxon;

	    my $test = join('|', @taxon);

	    ## Continue only if there are selected objects

	    if (scalar(@selected) > 0) {
		my $key = $id;
		if ($opt_A) {
		    $key = $acc;
		}
		
		## By default select every ID or AC, but if -i option is
		## enabled only will selected ID thate exists in the hash

		my $selected_data = 1;
		if (scalar(keys %ids) > 0) {
		    unless (exists $ids{$key}) {
			$selected_data = 0;
		    }
		}

		## Select taxon (if species is enabled, no select the data 
	        ## if there is not taxon)
		if (scalar(keys %species) > 0) {
		    $selected_data = 0;
		    foreach my $tx (@taxon) {
			if (exists $species{$tx}) {

			    $selected_data = 1;
			}
		    }
		}
		
		if ($selected_data == 1) {
		    my @cp_sel = @selected;
		    if (exists $entries{$key}) {
			push @{$entries{$key}}, @cp_sel;
		    }
		    else {
			$entries{$key} = \@cp_sel;
		    }
		}
	    }
	    ## Clean the ID, ACC and TAXON
	    if ($opt_V) {
		my $sl = scalar(@selected);
		my $et = scalar(keys %entries);
		my $counts = "Line:$l, Entries:$et, (ID:$id AC:$acc SEL:$sl)";
		print STDERR "\tParsing $counts         \r";
	    }
	    undef($id);
	    undef($acc);
	    undef($sps);
	    undef($taxon);
	    @selected = ();
	}
    }
    print STDERR "\n";
}

## PRINTING THE OUTFILE ######################

print STDERR "\n\n3) Printing the output file.\n";
open my $ofh, '>', $outbase;

my $o_cnt = 0;
foreach my $key (sort keys %entries) {
    
    ## Remove the redundancy
    my %nr = ();
    foreach my $data (@{$entries{$key}}) {
	$nr{$data} = 1;
    }

    foreach my $dt (sort keys %nr ) {
	print $ofh "$key\t$dt\n";
	$o_cnt++;
    }
}
print STDERR "\t$o_cnt lines have been printed for $outbase file.\n";



## SCRIPTS ENDS ##

$date = `date`;
chomp($date);

print STDERR "\n\n=======================================================";
print STDERR "\ngff2pathologic.pl ending [$date]\n";
print STDERR "=======================================================\n\n";


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
    
    This extract the annotations for a specific field, for example if 
    'Description' it will take all the description lines for each of the IDS.
 
    The fields can be choosen in the following way:
      * CATEGORY, two letters, if this is not followed by '=' it will take the
        complete line. More than one category can be specified separated by
        comma ','.
      * SUBCATEGORY, a word, it will take from ThisWord= to ';' or end of line.
        more than one subcategory can be used under the same category separated
        by semicolon ';'.

      Example1:

      uniprot2tab.pl -u uniprot_sprot.dat -f DE

      Example2:

      uniprot2tab.pl -u uniprot_sprot.dat -f DE=AltName

      Example3:

      uniprot2tab.pl -u uniprot_sprot.dat -f 'DE=AltName;EC,DT'

    More than one file can be used.

      Example4: 

      uniprot2tab.pl -u uniprot_sprot.dat,uniprot_trembl.dat -f Description

    Specific IDs and or species can be selected.

      Example5:

      uniprot2tab.pl -u uniprot_sprot.dat -f Description -i ids.txt 
                     -s 'Arabidopsis thaliana'
    
    Usage:
      
       uniprot2tab.pl [-h] -u <uniprot_file> -o <output_basename> 
                           -f <annot_fields>
                           [-i <id_list>] [-s <species_list>] [-A]
      
    Flags:
 
      -u <uniprot_file>       uniprot (swissprot) file (mandatory)
      -o <output_basename>    output basename (by default 'out')
      -f <annot_fields>       annotation fields separated by commas (mandatory)
      -i <id_list_file>       file with the list of ids to retrieve (optional)
      -s <species_list>       species to retrieve separated by commas (optional)
      -A <use_accessions>     use accessions instead IDs
      -V <be_verbose>         verbose (print status messages)
      -h <help>               print the help

EOF
exit (1);
}



####
1; #
####
