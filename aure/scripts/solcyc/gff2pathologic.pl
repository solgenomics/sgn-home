#!/usr/bin/perl

=head1 NAME

 gff2pathologic.pl
 Script to create a pathologic file from gff3 file.

=cut

=head1 SYPNOSIS

 gff2pathologic.pl [-h] -g <gff3_file> -t <soterm> -o <output_basename> [-N][-V]
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

=item -N

B<note2function>          parse note field to try to extract a function.

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

our ($opt_g, $opt_t, $opt_f, $opt_e, $opt_s, $opt_o, $opt_d, $opt_N, $opt_S, 
     $opt_V, $opt_h);
getopts("g:t:f:e:s:o:d:NSVh");
if (!$opt_g && !$opt_t && !$opt_f && !$opt_e && !$opt_s && !$opt_o && !$opt_N 
    && !$opt_d && !$opt_S && !$opt_V && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}

## Check variables.

my $gff_file = $opt_g 
    || die("DATA ARGUMENT -g <gff_file> WAS NOT SUPPLIED.\n");

my $req_soterm = $opt_t || 'gene';

my (@function_files, @ec_code_files, @syn_files, @dblink_files);
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

my $outbase = $opt_o || 'out';

## Define the variable to store the entries

my %sel_entries = ();

my $date = `date`;
chomp($date);

print STDERR "\n\n===========================================================";
print STDERR "\ngff2pathologic.pl initiation [$date]\n";
print STDERR "===========================================================\n\n";

## Parse the gff file

my $feaio = Bio::FeatureIO->new( 
    -file    => $gff_file, 
    -format  => 'GFF',
    -version => 3,
    );

my $ecnt = 0;

## PARSING GFF# FILE #############

print STDERR "\n\n1) Parsing gff3 file for SO term = $req_soterm.\n\n";

while ( my $feature = $feaio->next_feature() ) {
    
    $ecnt++;

    my $seqid = $feature->seq_id();
    my $name = $feature->name();
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
    }
}


## PARSING SYNONYM FILE ####################################

print STDERR "\n\n2) Parsing synonym file.\n\n";

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

print STDERR "\n\n3) Parsing ec code file.\n\n";

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

print STDERR "\n\n4) Parsing function file.\n\n";

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

print STDERR "\n\n5) Parsing dblink file.\n\n";

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

print STDERR "\n\n6) Creating the output files:\n";

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
		print $outfh "STARTBASE\t" . $ft->end() . "\n";
		print $outfh "ENDBASE\t" . $ft->start() . "\n";
	    }
	    else {
		print $outfh "STARTBASE\t" . $ft->start() . "\n";
		print $outfh "ENDBASE\t" . $ft->end() . "\n";
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

	    ## Add default
	    if (scalar(@functions) == 0) {
		push @functions, 'ORF';
	    } 

	    ## Print functions (no-redundant)
	    my %nr_functions = ();
	    foreach my $function (@functions) {
		unless (exists $nr_functions{$function}) {
		    if ($function =~ m/putative uncharacterized protein/i) {
			print $outfh "FUNCTION\tORF\n";
		    }
		    elsif ($function =~ m/unknown protein/i) {
			print $outfh "FUNCTION\tORF\n";
		    }
		    elsif ($function =~ m/hypothetical protein/i) {
			print $outfh "FUNCTION\tORF\n";
		    }
		    else {
			print $outfh "FUNCTION\t$function\n";
		    }
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
	    elsif ($opt_S) {
		foreach my $synon (@synonyms) {
		    if (exists $ecs{$synon}) {
			push @ecs, @{$ecs{$synon}};
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
	    elsif ($opt_S) {
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

        New functions, ec numbers and synonyms can be loaded with extra input 
      files (separated by commas) using -f, -e, -s and -d options.

    Example1:

     seqid       Chr1                           Filename = Chr1.pf
     source      TAIR9                          ===================
     type        gene                           
     start       3631                           
     end         5899                           ID           AT1G01010
     score       .                              NAME         AT1G01010 
     strand      +                        ===>  STARTBASE    3631
     phase       .                              ENDBASE      5899
     attributes   ID=AT1G01010;                 FUNCTION     unknown
                  Note=protein_coding_gene;     GENE-COMMENT protein_coding_gene
                  Name=AT1G01010

    Example2:

     seqid       SL2.40ch01                    Filename = SL2.40ch01.pf
     source      ITAG_eugene                   ========================
     type        mRNA                          (using -t mRNA)
     start       839826  
     end         851722                          ID         Solyc01g006220.2.1
     score       .                               NAME       Solyc01g006220.2.1
     strand      -                               STARTBASE  851722
     phase       .                               ENDBASE    839826
     attributes  ID=mRNA:Solyc01g006220.2.1;     FUNCTION   unknown
                 Name=Solyc01g006220.2.1;        GENE-COMMENT Histonre-lysine...
                 Note=Histone-lysine 
                      N-methyltransferase 
    
    Usage:
      
      gff2pathologic.pl [-h] -g <gff3_file> -t <soterm> -o <output_basename> 
                        [-f <function_files>] [-e <ec_code_files>] 
                        [-s <synonym_files>]  [-d <dblink_files>] [-N] [-S]
      
    Flags:
 
      -g <gff3_file>              gff3 file (mandatory)
      -o <output_basename>        output basename (by default 'out')
      -t <soterm>                 sequence ontology term catched 
                                  ('gene' by default)
      -f <function_file>          function file (tabular), X columns where
                                  1st=ID and Xth=function (optional)
      -e <ec_codes_file>          ec code file (tabular), X columns where
                                  1st=ID and Xth=ec_code (optional)
      -s <synonym_file>           synonym file (tabular), X columns where 
                                  1st=ID and Xth=synonyms (optional)
      -d <dblink_files>           dblink files (tabular), X columns where
                                  1st=ID and Xth=DBCODE:ACCESSSION (optional)
      -S <use_synonyms>           use synonyms for mapping with -f and -e files
      -N <note2function>          parse note field to try to extract a function.
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
	    print STDERR "\tParsing line: $line (ID=$id)                 \r";
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
