#!/usr/bin/perl

=head1 NAME

 gff2pathologic.pl
 Script to create a pathologic file from gff3 file.

=cut

=head1 SYPNOSIS

 gff2pathologic.pl [-h] -g <gff3_file> -t <soterm> -o <output_basename> 
                        [-N][-V][-S][-C][-I][-A]
                        [-a <annotation_files>]
                        [-f <function_files>]
                        [-e <ec_codes_files>]
                        [-s <synonyms_files>] 
                        [-d <dblink_files>]
 
=head2 I<Flags:>

=over


=item -g

B<gff3_file>              gff3 file (mandatory)

=item -o

B<output_basename>        output basename (by default 'out')

=item -t

B<soterm>                 sequence ontology term catched ('gene' by default)

=item -a 

B<annotation_file>        blast result file with an extra column with 
                          annotations (optional)

=item -f

B<function_file>          function file (tabular), X columns where 1st=ID (opt)

=item -e

B<ec_codes_file>          ec code file (tabular), X columns where 1st=ID (opt)

=item -s

B<synonym_file>           synonym file (tabular), X columns where 1st=ID (opt)

=item -d

B<dblink_file>            dblink file (tabular), X columns where 1st=ID (opt)

=item -S

B<use_synonyms>           use synonyms for mapping with -f and -e files

iterm -A

B<use_matched_accession>  use accession from annotat. file to map with -f or -e

=item -N

B<note2function>          parse note field to try to extract a function.

=item -I

B<extract_introns>        extract intron data information

=item -C

B<use_cds_coords>         use CDS coordinates as STARTBASE-ENDBASE

=item -V

B<be_verbose>             verbose (print status messages)

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script create one pathologic file per seqid (column 1 from gff). 

 GFF3 To Pathologic:
    - Column 1: "seqid"      => Output files (one per file).
    - Column 2: "source"     => Ignored.
    - Column 3: "type"       => Selected lines through -t parameter.
    - Column 4: "start"      => STARTBASE for each entry
    - Column 5: "end"        => ENDBASE for each entry
    - Column 6: "score"      => Ignored.
    - Column 7: "strand"     => Switch STARTBASE <-> ENDBASE.
    - Column 8: "phase"      => Ignored.
    - Column 9: "attributes" => ID    -> ID
                                Name  -> NAME
                                Alias -> SYNONYM
                                Note  -> GENE-COMMENT 

  New functions, ec numbers and synonyms can be loaded with extra imput files 
  using -f, -e and -s options.

  Annotation files could be used to extract the function from a blast hit 
  result. To do it, use -a argument as dbname=blastm8_annotatedfile.  
  The result will be collected in the following sections:
     
    FUNCTION      gene match description
    GENE-COMMENT  similar to ACCESSION (DB=DBNAME, EVALUE=evalue, %IDENTITY=%ID)

  Databases SWP (SWISSPROT), ATH (ARABIDOPSIS TAIR) and NR (GENBANK) will
  parse accession and function according some internal parameters.

  When more than one annotation file is used, it will take tag the annotations
  containing the following words as bad annotations:
   'hypothetical protein', 'unknown', 'predicted protein', 
   'unnammed protein product' and 'Protein of unknown function'

  If the option -A is used it will get the accession to retrieve functions or/
  and EC from the -f and -e files. 
    
  Stringent blast result filtering is advisable before use this annotation file.
  An useful script could be: blast_result_handle_tool.pl

  ============================================================================
  Suggested Working Pipelines:
  ============================================================================

  I- From annotated gff3 file. 
     (Note=Function, contains EC number.) 

        gff2pathologic.pl -g myfile.gff3 -o outfiles -N

  II- From gff3 file using Swissport/Trembl annotations.
      (GeneID are one of the GN subcategories for Swissprot file)

      Step1: Map geneIDs and Uniprot accessions

             uniprot2tab -u uniprot_sprot.dat -f 'GN=Name' -s 'My Species' 
                         -A -o swp_acc2gene_id.map.txt
             
             awk '{ print $2"\t"$1}' swp_acc2gene_id.map.txt > geneid2acc.txt


      Step2: Extract synonyms from Swissprot/Trembl files.

             uniprot2tab -u uniprot_sprot.dat -f 'GN=OrderedLocusNames,ORFNames'
                         -s 'My Species' -A -o swp_synonyms.txt

      Step3: Extract functions from Swissprot/Trembl files

              uniprot2tab -u uniprot_sprot.dat -f 'DE=Full' -s 'My Species'
                          -A -o swp_functions.txt

      Step4: Extract EC from Swissprot/Trembl files:
 
              uniprot2tab -u uniprot_sprot.dat -f 'DE=EC' -s 'My Species'
                          -A -o swp_ec.txt

      Step5: Extract DBLINKS from Swissprot/Trembl files:

             uniprot2tab -u uniprot_sprot.dat -f 'DR' -s 'My Species'
                         -A -o pre.swp_dblinks.txt

             grep 'GO;' pre.swp_dblinks.txt | sed -r 's/GO; /GO:/' | 
                  sed -r 's/; /\t/g' | cut -f1,2 > swp_dblinks.go.txt

             grep 'RefSeq;' pre.swp_dblinks.txt | sed -r 's/RefSeq; /RefSeq:/g'
                  | sed -r 's/; /\tRefSeq:/g' | sed -r 's/\.$//' > 
                  swp_dblinks.refseq.txt

             seq -r 's/\t/SWISSPROT:/' geneid2acc.txt > swp_dblinks.swp.txt

       Step6: Generate the pathologic file

              gff2pathologic.pl -g myfile.gff3 -o outfiles -S 
                                -s geneid2acc.txt,swp_synonyms.txt
                                -f swp_functions.txt
                                -e swp_ec.txt
                                -d swp_dblinks.swp.txt,swp_dblinks.go.txt,
                                   swp_dblinks.refseq.txt

 =============================================================================
 Mapping examples:
 =============================================================================

 Example1:

  seqid       Chr1                            Filename = Chr1.pf
  source      TAIR9                           ===================
  type        gene                           
  start       3631                           
  end         5899                            ID           AT1G01010
  score       .                               NAME         AT1G01010 
  strand      +                         ===>  STARTBASE    3631
  phase       .                               ENDBASE      5899
  attributes    ID=AT1G01010;                 FUNCTION     unknown
                Note=protein_coding_gene;     GENE-COMMENT protein_coding_gene
                Name=AT1G01010

 Example2:

  seqid       SL2.40ch01                      Filename = SL2.40ch01.pf
  source      ITAG_eugene                     ========================
  type        mRNA                            (using -t mRNA)
  start       839826  
  end         851722                          ID         mRNA:Solyc01g00622.2.1
  score       .                               NAME       Solyc01g006220.2.1
  strand      -                               STARTBASE  851722
  phase       .                               ENDBASE    839826
  attributes  ID=mRNA:Solyc01g006220.2.1;     FUNCTION   unknown
              Name=Solyc01g006220.2.1;        GENE-COMMENT Histonre-lysine...
              Note=Histone-lysine 
                   N-methyltransferase 

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 gff2 pathologic.pl


=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::FeatureIO;

our ($opt_g, $opt_t, $opt_a, $opt_f, $opt_e, $opt_s, $opt_o, $opt_d, 
     $opt_N, $opt_S, $opt_A, $opt_I, $opt_C, $opt_V, $opt_h);
getopts("g:t:a:f:e:s:o:d:NSAICVh");
if (!$opt_g && !$opt_t && !$opt_a && !$opt_f && !$opt_e && !$opt_s && !$opt_o 
    && !$opt_N && !$opt_d && !$opt_S && !$opt_I && !$opt_C && !$opt_V && 
    !$opt_A && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

my $date = `date`;
chomp($date);

print STDERR "\n\n===========================================================";
print STDERR "\ngff2pathologic.pl initiation [$date]\n";
print STDERR "===========================================================\n\n";

## Check variables.

my $gff_file = $opt_g 
    || die("DATA ARGUMENT -g <gff_file> WAS NOT SUPPLIED.\n");

my $req_soterm = $opt_t || 'gene';

my (@function_files, @ec_code_files, @syn_files, @dblink_files);
my %annot_files = ();

if ($opt_f) {
    @function_files = split(/,/, $opt_f);
}
if ($opt_e) {
    @ec_code_files = split(/,/, $opt_e);
}
if ($opt_s) {
    @syn_files = split(/,/, $opt_s);
}
if ($opt_d) {
    @dblink_files = split(/,/, $opt_d);
}
if ($opt_a) {
    my @annot_files = split(/,/, $opt_a);
    foreach my $annpair (@annot_files) {
	if ($annpair =~ m/^(.+?)=(.+)$/) {
	    $annot_files{$1} = $2;
	}
	else {
	    print STDERR "\nWARNING:-a $annpair doesn't have format db=file.\n";
	}
    }
}


my $outbase = $opt_o || 'out';

## Define the variable to store the entries

my %sel_entries = ();

## PRINT SOME OF THE OPTIONS

if ($opt_V) {

    print STDERR "\n\tGFF3 FILE:\t\t$gff_file\n";
    print STDERR "\tSO TERM:\t\t$req_soterm\n";
    print STDERR "\n\tOPTIONAL INPUT FILES:\n\n";
    print STDERR "\t\t(+) SYNONYMS FILES:\n";

    if ($opt_s) {
	foreach my $sfile (@syn_files) {
	    print STDERR "\t\t\t$sfile\n";
	}
    }
    else {
	print STDERR "\t\t\tNone\n";
    }

    print STDERR "\t\t(+) FUNCTION FILES:\n";    
    if ($opt_f) {
	foreach my $ffile (@function_files) {
	    print STDERR "\t\t\t$ffile\n";
	}
    }
    else {
	print STDERR "\t\t\tNone\n";
    }

    print STDERR "\t\t(+) EC NUMBER FILES:\n";    
    if ($opt_e) {
	foreach my $efile (@ec_code_files) {
	    print STDERR "\t\t\t$efile\n";
	}
    }
    else {
	print STDERR "\t\t\tNone\n";
    }

    print STDERR "\t\t(+) DBLINKS FILES:\n";
    if ($opt_d) {
	foreach my $dfile (@dblink_files) {
	    print STDERR "\t\t\t$dfile\n";
	}
    }
    else {
	print STDERR "\t\t\tNone\n";
    }

    print STDERR "\t\t(+) BLAST ANNOTATION FILES:\n";    
    if ($opt_a) {
	foreach my $adb (keys %annot_files) {
	    print STDERR "\t\t\t$adb = $annot_files{$adb}\n";
	}
    }
    else {
	print STDERR "\t\t\tNone\n";
    }

    print STDERR "\n\n\tUSE SYNONYMS\t\t\t\t\t[";
    if ($opt_S) { print STDERR "Enabled]\n"; }
    else { print STDERR "Disabled]\n"; }

    print STDERR "\n\tUSE BLAST ACCESSIONS MATCHED\t\t\t[";
    if ($opt_A) { print STDERR "Enabled]\n"; }
    else { print STDERR "Disabled]\n"; }

    print STDERR "\n\tEXTRACT FUNCTION FROM NOTE\t\t\t[";
    if ($opt_N) { print STDERR "Enabled]\n"; }
    else { print STDERR "Disabled]\n"; }

    print STDERR "\n\tREPLACE SELECTED COORDS. WITH CDS COORDS\t[";
    if ($opt_C) { print STDERR "Enabled]\n"; }
    else { print STDERR "Disabled]\n"; }

    print STDERR "\n\tEXTRACT INTRONS\t\t\t\t\t[";
    if ($opt_I) { print STDERR "Enabled]\n"; }
    else { print STDERR "Disabled]\n"; }
}

## Parse the gff file

my $feaio = Bio::FeatureIO->new( 
    -file    => $gff_file, 
    -format  => 'GFF',
    -version => 3,
    );

my $ecnt = 0;

my %geneids = ();
my %cds_endcoords = ();
my %cds_stacoords = ();
my %introns = ();

## PARSING GFF# FILE #############

print STDERR "\n\n1) Parsing gff3 file for SO term = $req_soterm.\n\n";

while ( my $feature = $feaio->next_feature() ) {
    
    $ecnt++;

    my $seqid = $feature->seq_id();
    my $soterm = $feature->type()->term()->name();
    my $strand = $feature->strand();
    my $st = $feature->start();
    my $en = $feature->end();

    if ($opt_V) {
	print STDERR "\tParsing entry ($ecnt):$seqid|$soterm|$st|$en        \r";
    }

    if ($soterm eq $req_soterm) {

	if (exists $sel_entries{$seqid}) {
	
	    push @{$sel_entries{$seqid}}, $feature;
	}
	else {
	    $sel_entries{$seqid} = [$feature];
	}
	
	## Store the gene id list for -C or -I
	if ($opt_C || $opt_I) {

	    my ($feat_id) = $feature->get_Annotations('ID');
	    if (defined $feat_id) {

		my $id = $feat_id->value();
		$id =~ s/.+?://;
		$geneids{$id} = 1;
	    }
	}
    }
}


## If -C it will catch the CDS coordinates for calculate start and end

if ($opt_C || $opt_I) { 

    print STDERR "\n\n\tOPTION -C or/and -I enabled.\n";
    print STDERR "\tReparsing gff3 to extract CDS or/and Intron data\n\n";

    my $c_feaio = Bio::FeatureIO->new( 
	-file    => $gff_file, 
	-format  => 'GFF',
	-version => 3,
	);

    my $c_ecnt = 0;

    while ( my $feature = $c_feaio->next_feature() ) {

	$c_ecnt++;
	my $m_msg = '';

	my $seqid = $feature->seq_id();
	my $soterm = $feature->type()->term()->name();
	my $st = $feature->start();
	my $en = $feature->end();
        
	my ($feat_id) = $feature->get_Annotations('ID');
	
	## But introns doesn't have ID, just Parent ID.
	if ($soterm eq 'intron') {
	    ($feat_id) = $feature->get_Annotations('Parent');
	}

	if (defined $feat_id) {

	    my $id = $feat_id->value();
	    $id =~ s/.+?://;

	    if (defined $id && $soterm =~ m/^(CDS|intron)$/) {

		## CDS id and Intron id could be different from gene id
		## to match with the gene id there are different options,
		## easier, store different variations too.
		## For example: F10G8.9b comes from F10G8.9
		##              Solyc00g005000.2.1.1 comes from Solyc00g005000.2
		
                ## It will scan if dont exist

		unless (exists $geneids{$id}) {

		    my $match = 0;
		    my $del = 0; ## Let say, try 8 characters

		    while($match == 0 && $del < 8 && $del < length($id)) {
			$id =~ s/.$//;
			if (exists $geneids{$id}) {
			    $match = 1;
			}
		    }
		}

		## Add the coords if it found something
	     
		if (exists $geneids{$id}) {

		    $m_msg = "Matching ID:$id";

		    if ($opt_C && $soterm eq 'CDS') {
			if (exists $cds_endcoords{$id}) {
			    push @{$cds_endcoords{$id}}, $en;
			}
			else {
			    $cds_endcoords{$id} = [$en];
			}
			if (exists $cds_stacoords{$id}) {
			    push @{$cds_stacoords{$id}}, $st;
			}
			else {
			    $cds_stacoords{$id} = [$st];
			}

			$m_msg .= ' +CDS';
		    }

		    if ($opt_I && $soterm eq 'intron') {
	
			if (exists $introns{$id}) {
			    push @{$introns{$id}}, $st . "-" . $en;
			}
			else {
			    $introns{$id} = [$st . "-" . $en];
			}
			$m_msg .= ' +Intron';
		    }
		}
		else {
		    $m_msg = "NO Matching ID:$id"
		}
	    }
	}
	if ($opt_V) {
	    print STDERR "\t\tParsing line ($c_ecnt-$soterm): $m_msg        \r";
	}
    }
}

##PARSING BLAST ANNOTATION FILES ###########################

print STDERR "\n\n2) Parsing blast annotation files.\n\n";

my %blast = ();

foreach my $bldb (sort keys %annot_files) {

    print STDERR "\t$bldb = $annot_files{$bldb} parsing.\n";

    my %bldata = ();

    open my $bfh, '<', $annot_files{$bldb};
    my $bline = 0;

    while (<$bfh>) {
	$bline++;
	chomp($_);
	
	my @data = split(/\t/, $_);
	my $id = $data[0];
	my $acc = $data[1];
	my $ident = $data[2];
	my $aln_l = $data[3];
	my $evalue = $data[10];
	my $desc = $data[12];

	my ($accession, $function);

	## Parse data according the specific database

	if ($bldb =~ m/^(SWP|SWISSPROT|TR|TREMBL)$/i) {
	    if ($acc =~ m/^(sp|tr)\|(.+?)\|.+/) {
		$accession = $2;
	    }
	    if ($desc =~ m/^(.+?)\s+OS=/) {
		$function = $1;
	    }
	}
	elsif ($bldb =~ m/^(ATH|TAIR|ARABIDOPSIS)$/i) {
	    $accession = $acc;
	    if ($desc =~ m/\|?\s+Symbols:.+?\|\s+(.+)\s+\|/) {
		$function = $1;
	    }
	}
	elsif ($bldb =~ m/^(NR|GENBANK)$/i) {
	    if ($acc =~ m/gi\|\d+\|\w+\|(.+?)\|/) {
		$accession = $1;
	    }
	    $function = $desc;
	    $function =~ s/\[.+//;
	}
	
	## Default parsing
	unless (defined $accession){
	    $accession = $data[1];
	}
	unless (defined $function) {
	    $function = $data[12] || 'Unknown';
	}

	my $genecmt = "similar to $accession (DB=$bldb, EVALUE=$data[10], ";
	$genecmt .= "PERC_IDENT=$data[2], ALIGN_LENGTH=$data[3])";

	if ($function =~ m/(unknown|hypothetical|predicted|unnammed)/i) {
	    $function = 'ORF';
	}

	## Take always the first

	unless (exists $bldata{$id}) {
	    $bldata{$id} = [$accession, $function, $genecmt];
	}
	if ($opt_V) {
	    print STDERR "\t\tParsing line:$bline  (ID:$id|ACC:$accession)  \r";
	}
    }
    close($bfh);
    $blast{$bldb} = \%bldata;
    if ($opt_V) {
	print STDERR "\n\n";
    }
}


## PARSING SYNONYM FILE ####################################

print STDERR "\n\n3) Parsing synonym files.\n\n";

my %synonyms = ();

if (scalar(@syn_files) > 0) {
    foreach my $syn_file (@syn_files) {
	%synonyms = parse_annotation_file($syn_file, \%synonyms);
	if ($opt_V) { 
	    print STDERR "\n";
	}
    }
}
else {
    print STDERR "\tNo synonym file was supplied. Skipping step.";
}

## PARSING EC CODE FILE #######################################

print STDERR "\n\n4) Parsing ec code file.\n\n";

my %ecs = ();

if (scalar(@ec_code_files) > 0) {
    foreach my $ec_code_file (@ec_code_files) {
	%ecs = parse_annotation_file($ec_code_file, \%ecs);
	if ($opt_V) { 
	    print STDERR "\n\n";
	}
    }
}
else {
    print STDERR "\tNo ec code file was supplied. Skipping step.";
}


## PARSING FUNCTION FILE ###########################################

print STDERR "\n\n5) Parsing function file.\n\n";

my %functions = ();

if (scalar(@function_files) > 0) {
    foreach my $function_file (@function_files) {
	%functions = parse_annotation_file($function_file, \%functions);
	if ($opt_V) { 
	    print STDERR "\n\n";
	}
    }
}
else {
    print STDERR "\tNo function file was supplied. Skipping step.";
}

## PARSING DBLINK FILE ###########################################

print STDERR "\n\n6) Parsing dblink file.\n\n";

my %dblinks = ();

if (scalar(@dblink_files) > 0) {
    foreach my $dblink_file (@dblink_files) {
	%dblinks = parse_annotation_file($dblink_file, \%dblinks);
	if ($opt_V) { 
	    print STDERR "\n\n";
	}
    }
}
else {
    print STDERR "\tNo dblink file was supplied. Skipping step.";
}

## CREATING OUTPUT ################################################

print STDERR "\n\n7) Creating the output files:\n";

foreach my $selseq (sort keys %sel_entries) {

    my $filename = $outbase . '.' . $selseq . '.pf'; 
    open my $outfh, '>', $filename;

    my $c = 0;

    foreach my $ft (@{$sel_entries{$selseq}}) {
    
	## Get the name:

	## ID
	my ($idobj) = $ft->get_Annotations('ID');

	if (defined $idobj) {
	    
	    my $id = $idobj->value();

	    $c++;
	    $id =~ s/.+?://;  ## Remove everything before the colon ':'

	    print $outfh "ID\t$id\n";

	    ## NAME (mandatory)
	    my $name;
	    my ($nameobj) = $ft->get_Annotations('Name');
	    if (defined $nameobj) {
		$name = $nameobj->value();
	    }
	    else {
		$name = $id;
	    }
	    print $outfh "NAME\t$name\n";
	    

	    ## SYNONYM
	    my @synonyms = (); 
	    my @synobjs = $ft->get_Annotations('Alias');
	    foreach my $synobj (@synobjs) {
		my $syn = $synobj->value();
		push @synonyms, $syn;
	    }
	    ## Add the synonyms parsed in the file
	    if (exists $synonyms{$id}) {
		push @synonyms, @{$synonyms{$id}};
		if ($opt_S) {
		    my @syn_from_syn = ();
		    foreach my $synon (@synonyms) {
			if (exists $synonyms{$synon}) {
			    push @syn_from_syn, @{$synonyms{$synon}};
			}
		    }
		    push @synonyms, @syn_from_syn;
		}
	    }
	    ## print the synonyms without redundancy
	    my %nr_synonyms = ();
	    foreach my $synterm (@synonyms) {
		unless (exists $nr_synonyms{$synterm}) {
		    print $outfh "SYNONYM\t$synterm\n";
		    $nr_synonyms{$synterm} = 1;
		}
	    }
	    @synonyms = sort keys %nr_synonyms;


	    ## STARTBASE and ENDBASE
	    if ($ft->strand() < 0) {
		if ($opt_C && $cds_endcoords{$id} && $cds_stacoords{$id}) {
		    print $outfh "STARTBASE\t" .$cds_endcoords{$id}->[-1]. "\n";
		    print $outfh "ENDBASE\t" . $cds_stacoords{$id}->[0] . "\n";
		}
		else {
		    print $outfh "STARTBASE\t" . $ft->end() . "\n";
		    print $outfh "ENDBASE\t" . $ft->start() . "\n";
		}
	    }
	    else {
		if ($opt_C && $cds_endcoords{$id} && $cds_stacoords{$id}) {
		     print $outfh "STARTBASE\t" .$cds_stacoords{$id}->[0]. "\n";
		    print $outfh "ENDBASE\t" . $cds_endcoords{$id}->[-1] . "\n";
		}
		else {
		    print $outfh "STARTBASE\t" . $ft->start() . "\n";
		    print $outfh "ENDBASE\t" . $ft->end() . "\n";
		}
	    }

	    ## If -I, print INTRONS

	    if (exists $introns{$id}) {
		my @introns = @{$introns{$id}};
		foreach my $intron (@introns) {
		    print $outfh "INTRON\t$intron\n";
		}
	    }
	    
	    ## FUNCTION
	    my @functions;
	    my @notes;
	    my @functionobjs = $ft->get_Annotations('Function');
	    my @noteobjs = $ft->get_Annotations('Note');

	    foreach my $noteobj (@noteobjs) {
		push @notes, $noteobj->value();
	    }
	    
	    ## Load gff3 functions into the array
	    foreach my $functionobj (@functionobjs) {
		push @functions, $functionobj->value();
	    }
	    ## Load file funtions into the array
	    if (exists $functions{$id}) {
		push @functions, @{$functions{$id}};
	    }
	    elsif (defined $name && exists $functions{$name}) {
		push @functions, @{$functions{$name}};
	    }
	    elsif ($opt_S) {
		foreach my $synon (@synonyms) {
		    if (exists $functions{$synon}) {
			push @functions, @{$functions{$synon}};
		    }
		}
	    }
	    elsif ($opt_N) {
		if (scalar(@notes) > 0) {
		    foreach my $note (@notes) {
			
			if ($note =~ m/^(unknown|hypothetical)\s+protein/i) {
			    push @functions, 'ORF';
			}
			else {
			    push @functions, $note;
			}		    
		    }
		}
	    }
	    
	    ## Add functions if annotation files are used. 

	    if ($opt_a) {  ## Get functions from blast accessions
		foreach my $bldb (sort keys %blast) {
		    my %bldata = %{$blast{$bldb}};

		    if (exists $bldata{$id}) {
		
			my $acc = $bldata{$id}->[0];
			my $fun = $bldata{$id}->[1];
			my $genecmt = $bldata{$id}->[2];
			if ($opt_A) {
			    if (defined $acc && exists $functions{$acc}) {
				push @functions, @{$functions{$acc}};
			    }
			    else {
				push @functions, $fun;
			    }
			}
			else {
			    push @functions, $fun;
			}
			push @notes, $genecmt;
		    }
		}
	    }

	    ## Add default and remove unknowns
	    if (scalar(@functions) == 0) {
		push @functions, 'ORF';
	    } 
	    else {
		
		my @clean_function = ();
		my $orf_match = 0;

		my @orfs = ('unknown', 'hypothetical', 'uncharacterized',
			    'unnamed', 'ORF');
		my $regexp_orf = join('|', @orfs);

		foreach my $fn (@functions) {
		    if ($fn =~ m/($regexp_orf)/i){			
			$orf_match++;
		    }
		    else {
			push @clean_function, $fn;
		    }
		}

		if (scalar(@clean_function) > 0) {
		    @functions = @clean_function;
		}
		else {
		    @functions = ('ORF');
		}
	    }

	    ## Now it will remove the redundancy based in the match 
	    ## with the shorter sequence

	    if (scalar(@functions) > 0) {
		my %newfunctions = ();
		my %length = ();

		foreach my $fnc (@functions) {
		    $length{$fnc} = length($fnc);
		}

		## order by length

		my @ord_func = sort {$length{$a} <=> $length{$b}} keys %length;
		my $last = $ord_func[-1];
		my $first = shift(@ord_func);

		## start the simplification, only if the first function
		## has more than 4 characters

		if (length($first) > 4) {

		    $newfunctions{$first} = 1;

		    while (scalar(@ord_func) > 0)  {
			
			$first =~ s/^\s+//;
			$first =~ s/\s+$//;
			    
			foreach my $ofnc (@ord_func) {

			    ## create two variables to contain the data and
			    ## format them

			    my $ref = $first;
			    my $tar = $ofnc;

			    ## Remove characters such as: '_', '-', 'spaces'
			    ## brackets, '/','\', '|', '+'. (Non-alphanum.)

			    $ref =~ s/\s+//g;
			    $tar =~ s/\s+//g;
			    $ref =~ s/\W//g;
			    $tar =~ s/\W//g;				

			    unless ($tar =~ m/\Q$ref\E/i) {

				## It is a new function, only if has not added
				unless (exists $newfunctions{$ofnc}) {
				    $newfunctions{$ofnc} = 1;
				}
			    }
			    else {
				$newfunctions{$ofnc} = 0;
			    }
			}
			$first = shift(@ord_func);
		    }

		    ## delete the not new functions

		    foreach my $selfnc (keys %newfunctions) {
			if ($newfunctions{$selfnc} == 0) {
			    delete($newfunctions{$selfnc}); 
			}
		    }
		    @functions = keys %newfunctions;
		}
	    }
	   

	    ## Print functions (no-redundant)
	    my %nr_functions = ();
	    foreach my $function (sort @functions) {
		unless (exists $nr_functions{$function}) {
		    print $outfh "FUNCTION\t$function\n";
		    $nr_functions{$function} = 1;
		}
	    }

	    ## PRODUCT-TYPE (always P)
	    print $outfh "PRODUCT-TYPE\tP\n";
	    	    
	    ## EC
	    my @ecs = ();
	    
	    ## Load ECs from gff3 
	    my @ecs_objs = $ft->get_Annotations('EC');
	    foreach my $ecobj (@ecs_objs) {
		push @ecs, $ecobj->value();
	    }
	    ## Load ECs from EC code file
	    if (exists $ecs{$id}) {
		push @ecs, @{$ecs{$id}};
	    }
	    elsif (defined $name && exists $ecs{$name}) {
		push @ecs, @{$ecs{name}};
	    }
	    if ($opt_S) {
		foreach my $synon (@synonyms) {
		    if (exists $ecs{$synon}) {
			push @ecs, @{$ecs{$synon}};
		    }
		}
	    }
	    if ($opt_A) {
		foreach my $bldb (sort keys %blast) {
		    my %bldata = %{$blast{$bldb}};

		    if (exists $bldata{$id}) {
			my $acc = $bldata{$id}->[0];
			if (exists $ecs{$acc}) {
			    push @ecs, @{$ecs{$acc}};
			}
		    }
		}
	    }

	    ## Print ECs (no redundant)
	    my %nr_ecs = ();
	    foreach my $ec (@ecs) {
		unless (exists $nr_ecs{$ec}) {
		    print $outfh "EC\t$ec\n";
		    $nr_ecs{$ec} = 1;
		}
	    }
	    
	    my @dblinks = ();
	    ## GO as DBLINKS (no redundant)
	    my @GOterms = $ft->get_Annotations('Ontology_term');
	    foreach my $goobj (@GOterms) {
		push @dblinks, $goobj->identifier();
	    }
	    ## Load the dblink from file
	    if (exists $dblinks{$id}) {
		push @dblinks, @{$dblinks{$id}};
	    }
	    elsif (defined $name && exists $dblinks{$name}) {
		push @dblinks, @{$dblinks{name}};
	    }
	    ## Add the dblink for the synonyms if exists
	    if ($opt_S) {
		foreach my $synon (@synonyms) {
		    if (exists $dblinks{$synon}) {
			push @dblinks, @{$dblinks{$synon}};
		    }
		}
	    }
	    ## print the dblinks
	    my %nr_dblinks = ();
	    foreach my $dblink (@dblinks) {
		unless (exists $nr_dblinks{$dblink}) {
		    print $outfh "DBLINK\t$dblink\n";
		    $nr_dblinks{$dblink} = 1;
		}
	    }

	    ## GENE-COMMENT
	    foreach my $note (@notes) {
		print $outfh "GENE-COMMENT\t$note\n";
	    }

	    print $outfh "//\n";
	}
    }

    print STDERR "\t$filename file has been created with $c entries\n";
}

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
       
      This script create one pathologic file per seqid (column 1 from gff). 

      GFF3 To Pathologic:
	  - Column 1: "seqid"      => Output files (one per file).
	  - Column 2: "source"     => Ignored.
	  - Column 3: "type"       => Selected lines through -t parameter.
	  - Column 4: "start"      => STARTBASE for each entry
	  - Column 5: "end"        => ENDBASE for each entry
	  - Column 6: "score"      => Ignored.
	  - Column 7: "strand"     => Switch STARTBASE <-> ENDBASE.
	  - Column 8: "phase"      => Ignored.
	  - Column 9: "attributes" => ID    -> ID
	                              Name  -> NAME
                                      Alias -> SYNONYM
                                      Note  -> GENE-COMMENT 

    New functions, ec numbers and synonyms can be loaded with extra imput files 
    using -f, -e and -s options.

    Annotation files could be used to extract the function from a blast hit 
    result. To do it, use -a argument as dbname=blastm8_annotatedfile.
    The result will be collected in the following sections:
     
      FUNCTION      gene match description
      GENE-COMMENT  similar to ACCESSION (DB=DBNAME, EVALUE=evalue,...)

    Databases SWP (SWISSPROT), ATH (ARABIDOPSIS TAIR) and NR (GENBANK) will
    parse accession and function according some internal parameters.

    When more than one annotation file is used, it will take tag the annotations
    containing the following words as bad annotations:
     'hypothetical protein', 'unknown', 'predicted protein', 
     'unnammed protein product' and 'Protein of unknown function'

    If the option -A is used it will get the accession to retrieve functions or
    and EC from the -f and -e files.
    
    Stringent blast result filtering is advisable before use this annotation 
    file. An useful script could be: blast_result_handle_tool.pl

  ============================================================================
  Suggested Working Pipelines:
  ============================================================================

  I- From annotated gff3 file. 
     (Note=Function, contains EC number.) 

        gff2pathologic.pl -g myfile.gff3 -o outfiles -N

  II- From gff3 file using Swissport/Trembl annotations.
      (GeneID are one of the GN subcategories for Swissprot file)

      Step1: Map geneIDs and Uniprot accessions

             uniprot2tab -u uniprot_sprot.dat -f 'GN=Name' -s 'My Species' 
                         -A -o swp_acc2gene_id.map.txt
             
             awk '{ print $2"\t"$1}' swp_acc2gene_id.map.txt > geneid2acc.txt


      Step2: Extract synonyms from Swissprot/Trembl files.

             uniprot2tab -u uniprot_sprot.dat -f 'GN=OrderedLocusNames,ORFNames'
                         -s 'My Species' -A -o swp_synonyms.txt

      Step3: Extract functions from Swissprot/Trembl files

              uniprot2tab -u uniprot_sprot.dat -f 'DE=Full' -s 'My Species'
                          -A -o swp_functions.txt

      Step4: Extract EC from Swissprot/Trembl files:
 
              uniprot2tab -u uniprot_sprot.dat -f 'DE=EC' -s 'My Species'
                          -A -o swp_ec.txt

      Step5: Extract DBLINKS from Swissprot/Trembl files:

             uniprot2tab -u uniprot_sprot.dat -f 'DR' -s 'My Species'
                         -A -o pre.swp_dblinks.txt

             grep 'GO;' pre.swp_dblinks.txt | sed -r 's/GO; /GO:/' | 
                  sed -r 's/; /\t/g' | cut -f1,2 > swp_dblinks.go.txt

             grep 'RefSeq;' pre.swp_dblinks.txt | sed -r 's/RefSeq; /RefSeq:/g'
                  | sed -r 's/; /\tRefSeq:/g' | sed -r 's/\.$//' > 
                  swp_dblinks.refseq.txt

             seq -r 's/\t/SWISSPROT:/' geneid2acc.txt > swp_dblinks.swp.txt

       Step6: Generate the pathologic file

              gff2pathologic.pl -g myfile.gff3 -o outfiles -S 
                                -s geneid2acc.txt,swp_synonyms.txt
                                -f swp_functions.txt
                                -e swp_ec.txt
                                -d swp_dblinks.swp.txt,swp_dblinks.go.txt,
                                   swp_dblinks.refseq.txt

   =============================================================================
   Mapping examples:
   =============================================================================

   Example1:

    seqid       Chr1                            Filename = Chr1.pf
    source      TAIR9                           ===================
    type        gene                           
    start       3631                           
    end         5899                            ID           AT1G01010
    score       .                               NAME         AT1G01010 
    strand      +                         ===>  STARTBASE    3631
    phase       .                               ENDBASE      5899
    attributes    ID=AT1G01010;                 FUNCTION     unknown
                  Note=protein_coding_gene;     GENE-COMMENT protein_coding_gene
                  Name=AT1G01010

   Example2:

    seqid       SL2.40ch01                      Filename = SL2.40ch01.pf
    source      ITAG_eugene                     ========================
    type        mRNA                            (using -t mRNA)
    start       839826  
    end         851722                          ID         mRNA:Solyc01g00622.2
    score       .                               NAME       Solyc01g006220.2.1
    strand      -                               STARTBASE  851722
    phase       .                               ENDBASE    839826
    attributes  ID=mRNA:Solyc01g006220.2.1;     FUNCTION   unknown
                Name=Solyc01g006220.2.1;        GENE-COMMENT Histonre-lysine...
                Note=Histone-lysine 
                     N-methyltransferase 
    
    Usage:
      
      gff2pathologic.pl [-h] -g <gff3_file> -t <soterm> -o <output_basename> 
                        [-N][-V][-S][-C][-I][-A]
                        [-a <annotation_files>]
                        [-f <function_files>]
                        [-e <ec_codes_files>]
                        [-s <synonyms_files>] 
                        [-d <dblink_files>]
      
    Flags:
 
      -g <gff3_file>              gff3 file (mandatory)
      -o <output_basename>        output basename (by default 'out')
      -t <soterm>                 seqontology term catched ('gene' by default)
      -a <annotation_file>        blast result file with an extra column with 
                                  annotations (optional)
      -f <function_file>          function file (tabular), X columns 
                                  where 1st=ID (opt)
      -e <ec_codes_file>          ec code file (tabular), X columns 
                                  where 1st=ID (opt)
      -s <synonym_file>           synonym file (tabular), X columns 
                                  where 1st=ID (opt)
      -d <dblink_file>            dblink file (tabular), X columns 
                                  where 1st=ID (opt)
      -S <use_synonyms>           use synonyms for mapping with -f and -e files
      -A <use_matched_accession>  use accession from annotat. file to map 
                                  with -f or -e
      -N <note2function>          parse note field to try to extract a function.
      -I <extract_introns>        extract intron data information
      -C <use_cds_coords>         use CDS coordinates as STARTBASE-ENDBASE 
      -V <be_verbose>             verbose (print status messages)
      -h <help>                   print the help

EOF
exit (1);
}

=head2 parse_annotation_file

  Usage: my %annotation = parse_annotation_file($file);

  Desc: Parse an annotation file with X columns, 1st for ID.

  Ret: %annotation, a hash with key=id, value=arrayref of elements (one per col)

  Args: $file, filename to parse

  Side_Effects: print status if $opt_V (verbose)

=cut

sub parse_annotation_file {

    my $filename = shift ||
	die("ERROR: No filename was supplied to parse_annotation_file()");
    my $annot_href = shift;

    my %annot = ();
    if (defined $annot_href && ref($annot_href) eq 'HASH') {
	%annot = %{$annot_href};
    }

    open my $afh, '<', $filename;
    my $line = 0;

    while(<$afh>) {
	chomp($_);
	$line++;

	my @data = split(/\t/, $_);
	my $id = shift(@data);

	if ($opt_V) {
	    print STDERR "\tParsing $filename line: $line (ID=$id)          \r";
	}

	if (defined $id && scalar(@data) > 0) {
	    if (exists $annot{$id}) {
		push @{$annot{$id}}, @data;
	    }
	    else {
		$annot{$id} = \@data;
	    }
	}
    }
    
    close($afh);
    return %annot;
}



####
1; #
####
