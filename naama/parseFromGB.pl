
=head1 NAME

parseCDSfromGB

=head1 DESCRIPTION

Reads sequences from a GenBank file, parses the CDS and 
prints them to a FASTA format outfile
This program user Bio::SeqIO

=head1 SYNOPSIS

    %parseCDSfromGB.pl -gbfile <genbank format file> -outfile <file to write output> [options]

=head1 COMMAND-LINE OPTIONS

 -gbfile          The file containing GenBank output 
 -organism        The organism for the data

 -a               print all sequences to outfile (CDS + sequences with no CDS) 

 -x               Dry run- rollback the transaction
 -o               output file
 -h               display this help text

=head1 AUTHOR

Naama Menda -  <nm249@cornell.edu>

=cut

use Bio::SeqIO;

use Getopt::Long;
use Bio::DB::GenBank;

my ($ORGANISM, $GBFILE, $dry_run, $out, $help);

GetOptions(
    'organism=s'  => \$ORGANISM,
    'gbfile|i=s'    => \$GBFILE,
    'outfile|o=s' =>\$out,
    'help|h'      =>\$help,
    ) or ( system( 'pod2text', $0 ), exit -1 );

#if ($help || !$gbfile ) {   ;}

###my $in  = Bio::SeqIO->new(-file => $GBFILE , '-format' => 'genbank');
open (OUTFILE, ">$out") || die "cannot open \"$out\": $!";
open FILE, $GBFILE or die $!;

my $print;
my $found=0;

 SEQ: while (my $accession =  <FILE> ) {
     my $cds=undef;
     
     my $gb = new Bio::DB::GenBank;
     my $seq = $gb->get_Seq_by_acc($accession); # Accession Number
     print STDOUT "Looking at accession $accession..\n";
     
     my $id= $seq->display_id();
     
     my $desc= $seq->desc();          #dbxref.description
     my $accession=$seq->accession(); #feature.name
     my $gi = $seq->primary_id();     #dbxref.accession, feature.dbxref_id
     

     my $keywords = $seq->keywords();
     
     my $length=$seq->length();       #feature.seqlen
     my $sequence= $seq->seq();       #feature.residues
     my $mol_type=$seq->molecule();   #feaature.type_id
     my @secondary_accessions=$seq->get_secondary_accessions(); # feature_dbxref.dbxref_id ??
     my $seq_version= $seq->seq_version();
     my $species_obj=$seq->species();  #check if the species matches $ORGANISM
     my $species=$species_obj->species();
     #print "$species\n";
     my $common_name=$species_obj->common_name();
     
     ###foreach my $feat ( $seq_obj->top_SeqFeatures ) {
     for my $feat ($seq->get_SeqFeatures) {
	 if ( ( $feat->primary_tag eq 'CDS' ) || ( $feat->primary_tag eq 'mRNA') ) {
	     my $fname = $feat->primary_tag;
	     my $location = $feat->start . ":" . $feat->end ;
	     my $f_obj = $feat->spliced_seq;
	     print STDOUT "found $fname (location:$location) for feature $gi\n";
	     print OUTFILE  ">gi|$gi|$fname|$location|$id|$desc\n",$f_obj->seq,"\n\n";
	     #}else {print "skipping species $species (gi=$gi)...\n"; }
	     $cds=1;
	 }
     }
     if (!$cds) {  
	 print STDOUT "*No CDS/mRNA for feature $gi. \n";
	 ##print OUTFILE  ">gi|$gi|$id|$desc\n",$sequence,"\n\n"; 
     }
     
}

close OUTFILE;
