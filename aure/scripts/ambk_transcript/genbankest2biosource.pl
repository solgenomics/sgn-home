#!/usr/bin/perl

=head1 NAME

 genbankest2biosource.pl
 Script to extract library (sample) data from GenBank EST format (version.0.1).

=cut

=head1 SYPNOSIS

 genbankest2biosource.pl -e <genbank_est> [-F] [-L] 

=head2 I<Flags:>

=over


=item -e

B<genbank_est>                  genbank est format file (mandatory)

=item -F

B<fasta_file>                   produces a fasta file with all the est sequences

=item -L

B<library_file>                 produces a library file in a gff3 format

=item -h

B<help>                         print the help

=back

=cut

=head1 DESCRIPTION

 This script parse genbank.est files and produces 3 files: 
  + biosource file, to be loaded as samples in biosource schema.
  + fasta file, with the genbank accession as seq_id
  + tab column with lib\tseq

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 genbankest2biosource.pl

=cut

use strict;
use warnings;
use autodie;

use File::Basename;
use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;
use Math::BigFloat;
use CXGN::Biosource::GB::Est;

our ($opt_e, $opt_F, $opt_L, $opt_h);
getopts("e:FLh");
if (!$opt_e && !$opt_F && !$opt_L && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Get the arguments and check the mandatory ones

my $estfile = $opt_e
    || die("MANDATORY ARGUMENT ERROR: -e <est_file> was not supplied.\n");

## First, parse the file and get the data

print STDERR "\n\n1) PARSING EST FILE.\n\n";

open my $estfile_fh, '<', $estfile;

my %gb_fields = ( 
    'IDENTIFIERS' => [ 'EST name', 'dbEST Id', 'GenBank Acc', 'GenBank gi' ],
    'CLONE INFO'  => [ 'Clone Id', 'Plate', 'Source', 'Insert length', 
		       'Id as DNA', 'DNA type' ],
    'PRIMERS'     => [ 'PCR forward', 'PCR backward', 'Sequencing', 
		       'PolyA Tail' ],
    'SEQUENCE'    => [ 'Quality', 'Entry Created', 'Last Updated' ],
    'COMMENTS'    => [],
    'PUTATIVE ID' => [],
    'LIBRARY'     => [ 'dbEST lib id', 'Lib Name', 'Organism', 'Subspecies', 
		       'Cultivar', 'Strain', 'Organ', 'Tissue type', 
		       'Lab host', 'Vector', 'R. Site 1','R. Site 2',
		       'Develop. stage', 'Description'],
    'SUBMITTER'   => [ 'Name', 'Lab', 'Institution', 'Address', 'Tel', 'Fax', 
		       'E-mail' ],
    'CITATIONS'   => [ 'PubMed ID', 'Title', 'Authors', 'Year', 'Status', 
		       'Citation'],
    );


my $l = 0;

my $gb_est;
my %est = ();
my $seq = '';
my $com = '';
my $putid = '';
my $description = '';
my $enable_catchseq = 0;
my $enable_catchcom = 0;
my $enable_catchputid = 0;
my $enable_catchdescrip = 0;


while (<$estfile_fh>) {
    $l++;

    print STDERR "\tParsing line=$l for file=$estfile.\r";
    chomp($_);

    my $catcher = $enable_catchseq + 
	$enable_catchcom + 
	$enable_catchdescrip +
	$enable_catchputid;

    ## Ignore the empty lines
    
    unless ($_ =~ m/^$/) {

	## Every new entry starts with 'IDENTIFIERS', it create a new
	## CXGN::Biosource::GB::Est object
	
	if ($_ =~ m/^IDENTIFIERS$/) {
	    $gb_est = CXGN::Biosource::GB::Est->new();
	}

	## Catch the field and the data only if the catchers are disabled

	if ( $catcher == 0 ) {
	    
	    if ($_ =~ m/^(.+?):\s+(.+)$/) {
		my $field = $1;
		my $data = $2;
		$gb_est->add_data({ $field => $data});
	    }
	}

	## Catch the sequence from SEQUENCE to Entry Created field

	if ($_ =~ m/^SEQUENCE$/) {
	    $enable_catchseq = 1;
	}
	elsif ($_ =~ m/^Entry Created|^Quality/ && $enable_catchseq == 1) {
	    $enable_catchseq = 0;
	    $gb_est->add_data({ 'sequence' => $seq });
	    $seq = '';
	}
	elsif ($enable_catchseq == 1) {
	    my $seqfrag = $_;
	    $seqfrag =~ s/\s+//g;
	    $seq .= $seqfrag;
	}
	
	## Catch the comment from COMMENT entry

	if ($_ =~ m/^COMMENTS$/) {
	    $enable_catchcom = 1;
	}
	elsif ($_ =~ m/^LIBRARY/) {
	    $enable_catchcom = 0;
	    $gb_est->add_data({ 'comments' => $com });
	    $com = '';
	}
	elsif ($enable_catchcom == 1) {
	    my $comfrag = $_;
	    $comfrag =~ s/^\s+//g;
	    $com .= ' ' . $comfrag;
	}

	## Catch the multiline PUTATIVE ID

	if ($_ =~ m/^PUTATIVE ID/) {
	    $enable_catchputid = 1;
	}
	elsif ($_ =~ m/^LIBRARY/ && $enable_catchputid == 1) {
	    $enable_catchputid = 0;
	    $gb_est->add_data({ 'putative id' => $putid });
	    $putid = '';
	}
	elsif ($_ =~ m/^COMMENTS/ && $enable_catchputid == 1) {
	    $enable_catchputid = 0;
	    $gb_est->add_data({ 'putative id' => $putid });
	    $putid = '';
	}
	elsif ($enable_catchputid == 1) {
	    my $putidfrag = $_;
	    $putidfrag =~ s/^\s+//g;
	    $putid .= ' ' . $putidfrag;
	}

	## Catch description for the library (multiline)
	
	if ($_ =~ m/^Description:\s+(.+)$/) {
	    $enable_catchdescrip = 1;
	    $description .= $1;
	}
	elsif ($_ =~ m/^SUBMITTER/) {
	    $enable_catchdescrip = 0;
	    $gb_est->add_data({ 'Description' => $description });
	    $description = '';
	}
	elsif ($enable_catchdescrip == 1) {
	    my $descripfrag = $_;
	    $descripfrag =~ s/^\s+//;
	    $description .= $descripfrag;
	}

	## Finally if finish the entry (||) put the data into a hash

	if ($_ =~ m/^\|\|$/) {
	    my $gb_accession = $gb_est->get_data('GenBank Acc');
	    $est{$gb_accession} = $gb_est;
	    $gb_est = '';
	}
    }
}
print STDERR "\n\n";

print STDERR "2) EXTRACTING LIBRARY INFORMATION.\n\n";

my $n = 0;
my %libs = ();
foreach my $gb_acc (keys %est) {
    my $lib_name = $est{$gb_acc}->get_data('Lib Name');

    unless (exists $libs{$lib_name}) {
	
	print STDERR "\tNew library ($lib_name) added to the pool. $n total libraries.\r";
	my @lib_fields = @{$gb_fields{'LIBRARY'}};
	my @cit_fields = @{$gb_fields{'CITATIONS'}};
	push @lib_fields, @cit_fields;
	my %new_lib;
	foreach my $libfield (@lib_fields) {
	    my $data = $est{$gb_acc}->get_data($libfield);
	    $new_lib{$libfield} = $data;
	    $new_lib{'EST_N'} = 1;
	    $n++;
	}
	$libs{$lib_name} = \%new_lib;
    }
    else {
	$libs{$lib_name}->{'EST_N'}++;
    }
}
print STDERR "\n\n";

print STDERR "3) PRINTING LIBRARY INFORMATION.\n\n";

my $output_lib = $estfile . '.estlibrary.bs';
open my $lib_fh, '>', $output_lib;


my $ln = 0;
my $LT = scalar(keys % libs);

foreach my $lib (sort keys %libs) {
    $ln++;
    my %data_libs = %{$libs{$lib}};
    print $lib_fh "\n\n####################\n## Library $ln of $LT (ESTs=$data_libs{'EST_N'})\n####################\n\n\n";
    print $lib_fh "*DATA_TYPE:             [sample]\n";
    print $lib_fh "*SAMPLE_NAME:           [$data_libs{'Lib Name'}]\n";
    print $lib_fh "*SAMPLE_TYPE_NAME:      [cDNA library]\n";
    
    if (defined $data_libs{'Description'}) {
	print $lib_fh "*SAMPLE_DESCRIPTION:    [$data_libs{'Description'}]\n";
    }

    print $lib_fh "*ORGANISM_NAME:         [$data_libs{'Organism'}]\n";
    print $lib_fh "\n\n## Publication Link ##\n\n";
    print $lib_fh "*DATA_TYPE:             [pub]\n";
    print $lib_fh "*SAMPLE_NAME:           [$data_libs{'Lib Name'}]\n";
    print $lib_fh "*TITLE:                 [$data_libs{'Title'}]\n";
    print $lib_fh "\n\n## PO Link ##\n\n";
    print $lib_fh "*DATA_TYPE:             [dbxref]\n";
    print $lib_fh "*SAMPLE_NAME:           [$data_libs{'Lib Name'}]\n";
    print $lib_fh "*DBNAME:                [PO]\n";
    
    if (defined $data_libs{'Organ'}) {
	print $lib_fh "## Search manually PO terms for organ=$data_libs{'Organ'}\n";
    }

    if (defined $data_libs{'Tissue type'}) {
	print $lib_fh "## Search manually PO terms for tissue=$data_libs{'Tissue type'}\n";
    }
    
    if (defined $data_libs{'Develop. stage'}) {
	print $lib_fh "## Search manually PO terms for development=$data_libs{'Develop. stage'}\n";
    }
    print $lib_fh "## *ACCESSIONS:            []\n";
    print $lib_fh "\n\n\n## Other Data that can be used for Description/Protocols: ##\n";
    foreach my $libkey (keys %data_libs) {
	if (defined $data_libs{$libkey}) {
	    print $lib_fh "##\t$libkey:\t$data_libs{$libkey}\n";
	}
    }
    print $lib_fh "\n\n";
}

print STDERR "4) PRINTING SEQUENCES FASTA.\n\n";
if ($opt_F) {
    my $output_fasta = $estfile . '.est.fasta';
    
    my $seqio = Bio::SeqIO->new( -format => 'fasta', -file => ">$output_fasta");

    my $s = 0;
    my $S = scalar(keys %est);

    foreach my $gbacc (keys %est) {

	$s++;
	print STDERR "\twritting seq $s of $S\r"; 
	my $lib_name = $est{$gbacc}->get_data('Lib Name');
	my $clone_name = $est{$gbacc}->get_data('Clone Id') || 'undef';
	my $est_lib_N = $libs{$lib_name}->{'EST_N'};
	my $alt_lib_name = 'General_EST_library';
	if ($est_lib_N > 99) {
	    $alt_lib_name = $lib_name;
	}
	
	my $seq = $est{$gbacc}->get_data('sequence');
	
	my $seqdescription = "CLONE_ID='$clone_name';LIB_NAME1='$lib_name';LIB_NAME2='$alt_lib_name';";
	my $seqobj = Bio::Seq->new( -id => $gbacc, -description => $seqdescription, -seq => $seq );
	$seqio->write_seq($seqobj);
    }
    print STDERR "\n\n";

}
else {
    print STDERR "\tOption Disabled.\n\n";
}

print STDERR "5) PRINTING LIBRARY MEMBERSHIP.\n\n";
if ($opt_L) {
    my $output_tab = $estfile . '.librarymembership.tab';
    
    open my $tab_io, '>', $output_tab;

    my $t = 0;
    my $T = scalar(keys %est);

    foreach my $gbacc (keys %est) {

	$t++;
	print STDERR "\twritting seq member $t of $T\r"; 
	my $lib_name = $est{$gbacc}->get_data('Lib Name');
	my $clone_name = $est{$gbacc}->get_data('Clone Id') || 'undef';
	my $est_lib_N = $libs{$lib_name}->{'EST_N'};
	my $alt_lib_name = 'General_EST_library';
	if ($est_lib_N > 99) {
	    $alt_lib_name = $lib_name;
	}
	
	print $tab_io "$gbacc\t$clone_name\t$lib_name\t$alt_lib_name\n";
    }
}
else {
    print STDERR "\tOption Disabled.\n\n";
}
print STDERR "\n\n";


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

       This script parse genbank.est files and produces 3 files: 
         + biosource file, to be loaded as samples in biosource schema.
         + fasta file, with the genbank accession as seq_id
         + tab column with lib\tseq
         
    Usage:

      genbankest2biosource.pl -e <genbank_est> [-F] [-L] 
       
    Flags:
      
      -e <genbank_est>                  genbank est format file (mandatory)
      -F <fasta_file>                   produces a fasta file with all the est sequences
      -L <library_file>                 produces a library file in a gff3 format
      -h <help>                         print the help

      

EOF
exit (1);
}


