#!/usr/bin/perl 

=head1 NAME

get_longest_protein.pl -- a script that takes a DNA fasta file and outputs the longest translated sequence of the 6 possible frames. 

=head1 DESCRIPTION

get_longest_protein.pl <sequence_file.fasta>

=head2 Options:

=over 18

=item --proteinpilot

produces proteinpilot compatible files

=item  --verbose 

produces verbose output.

=back

=head2 Output: 

The output goes to two fasta files, named [inputfile].protein and [inputfile].cds, containing the protein sequences and the corresponding cds sequences, respectively.

=head2 Notes

proteinpilot is an obscure program for identifying proteins using protemoics techniques and requires the fasta file db to be in the format:
 
 >lcl|identifier description without special chars
 PROTEINSEQUENCE

any derivation from this will cause the program to crash (I think).



a fasta file with protein sequences that represent the longest possible protein sequence from all 6 translation frames. 

The corresponding frame is added to the id of the sequence (+1..+3, -1..-3). 

The longest possible sequence is defined as the longest sequence occurring between stop codons (which are represented as asterisks in the sequence).

=head1 VERSION

 0.2, October 1, 2007. Added cds file generation and altered the output to go to files.

=head1 AUTHOR

Lukas Mueller (lam87@cornell.edu)

=cut

use strict;
use Bio::Seq;
use Bio::SeqIO;
use Getopt::Long;

my $proteinpilot = "";
my $verbose      = "";
GetOptions(
    proteinpilot => \$proteinpilot,
    verbose      => \$verbose
);



if ($proteinpilot) { print STDERR "Generating proteinpilot output...\n"; }

my $filename = shift || die("\n\nERROR: None fasta file was supplied to the script. For more information use perldoc.\n\n");

my $file = Bio::SeqIO->new( '-format' => 'Fasta', -file => $filename );
my $proteinfile = Bio::SeqIO->new(
    '-format' => 'Fasta',
    -file     => ">" . $filename . ".protein"
);
my $cdsfile =
  Bio::SeqIO->new( '-format' => 'Fasta', -file => ">" . $filename . ".cds" );

while ( my $seqobj = $file->next_seq() ) {

    # Get line and chop it up in little pieces...

    my $accno = $seqobj->accession_number();
    my $desc2 = $seqobj->desc();
    my $seq   = $seqobj->seq();
    my $id    = uc( $seqobj->id() );

    if ($verbose) { print STDERR "Processing $id...\n$seq\n"; }

    my ( $peptides_ref, $cds_ref ) = get_longest_peptide($seqobj);

    my @peptides = @$peptides_ref;
    my @cds      = @$cds_ref;

    my $id_count = 1;
    for ( my $n == 0 ; $n < (@peptides) ; $n++ ) {
        if ($proteinpilot) {
            my $new_id = "lcl|" . ( $peptides[$n]->id );
            $peptides[$n]->id($new_id);
            my $desc = $desc2;
            $desc =~ s/\|/ /g;    # removed pipes from description, these
                                  # cause problems with proteinpilot
            $peptides[$n]->desc($desc);
        }
        $proteinfile->write_seq( $peptides[$n] );

        if ($verbose) {
            print STDERR "WRITING: "
              . $cds[$n]->id() . "\n"
              . $cds[$n]->seq() . "\n";
        }
        $cdsfile->write_seq( $cds[$n] );
        $id_count++;
    }

}

$proteinfile->close();
$cdsfile->close();
$file->close();

sub get_longest_peptide {

    # takes a DNA Bio::Seq object as a parameter and
    # returns the longest peptides from a six frame translation
    # as Bio::Seq protein objects.
    #
    my $seqobj = shift;
    my $id     = $seqobj->id();

    my %frame = ();

    # generate objects that start at frame 1, 2 and 3
    #
    $frame{"+1"} = $seqobj;
    $frame{"+1"}->id( $id . "+1" );
    $frame{"+2"} = Bio::Seq->new();
    $frame{"+2"}->seq( $seqobj->subseq( 2, $seqobj->length() ) );
    $frame{"+2"}->id( $id . "+2" );
    $frame{"+3"} = Bio::Seq->new();
    $frame{"+3"}->seq( $seqobj->subseq( 3, $seqobj->length() ) );
    $frame{"+3"}->id( $id . "+3" );

    # generate the reverse sequences and reverse frames
    #
    $frame{"-1"} = $seqobj->revcom();
    $frame{"-1"}->id( $id . "-1" );
    $frame{"-2"} = Bio::Seq->new();
    $frame{"-2"}->seq( $frame{"-1"}->subseq( 2, $frame{"-1"}->length() ) );
    $frame{"-2"}->id( $id . "-2" );
    $frame{"-3"} = Bio::Seq->new();
    $frame{"-3"}->seq( $frame{"-1"}->subseq( 3, $frame{"-1"}->length() ) );
    $frame{"-3"}->id( $id . "-3" );

    # translate the forward and reverse frames
    #
    my @proteins;
    my @cds;
    foreach my $r ( "+1", "+2", "+3", "-1", "-2", "-3" ) {
        my $p = $frame{$r}->translate();
        if ($verbose) {
            print STDERR "translated seq: " . $p->id() . ":" . $p->seq() . "\n";
            print STDERR "coding sequence: " . $frame{$r}->seq() . "\n";
        }
        if ( $p->seq() ) {
            push @proteins, $p;
            push @cds,      $frame{$r};
        }
    }

    # generate all the peptide fragments between stop codons (*)
    # and identify the longest one.
    #
    my @total_frags     = ();
    my @total_cds_frags = ();

    for ( my $i = 0 ; $i < @proteins ; $i++ ) {

        my @frag_objs     = ();
        my @cds_frag_objs = ();

        my @frags = split /\*/, $proteins[$i]->seq();
        my @cds_frags = split_coding( $cds[$i]->seq() );

        if ($verbose) {
            print STDERR "Number of protein fragments: "
              . scalar(@frags) . "\n";
            print STDERR "Number of cds fragments : "
              . scalar(@cds_frags) . "\n";
        }

        for ( my $n = 0 ; $n < @frags ; $n++ ) {
            my $frag_obj = Bio::Seq->new();
            $frag_obj->id( $proteins[$i]->id() );

            my $cds_obj = Bio::Seq->new();
            $cds_obj->id( $proteins[$i]->id() );

            if ( length( $frags[$n] ) > 0 ) {

                if ($verbose) {
                    print STDERR "PROTEIN: ["
                      . length( $frags[$n] )
                      . "] $frags[$n]\nCDS: ["
                      . length( $cds_frags[$n] )
                      . "] $cds_frags[$n]\n\n";
                }
                if ( ( length( $frags[$n] ) * 3 ) != length( $cds_frags[$n] ) )
                {
                    print STDERR "WARNING for "
                      . $frag_obj->id()
                      . ": protein ["
                      . length( $frags[$n] ) . " ~ "
                      . ( 3 * length( $frags[$n] ) )
                      . "] and cds ["
                      . length( $cds_frags[$n] )
                      . "] fragment length don\'t match!!!!!!!\n";
                }
                $frag_obj->seq( $frags[$n] );
                $cds_obj->seq( $cds_frags[$n] );

                push @frag_objs,     $frag_obj;
                push @cds_frag_objs, $cds_obj;
            }

        }
        @total_frags     = ( @total_frags,     @frag_objs );
        @total_cds_frags = ( @total_cds_frags, @cds_frag_objs );

    }
    my $longest_length = 0;
    foreach my $f (@total_frags) {
        if ( $f->length() > $longest_length ) {
            $longest_length = $f->length();
        }
    }
    my @return_set     = ();
    my @cds_return_set = ();
    for ( my $n = 0 ; $n < (@total_frags) ; $n++ ) {
        if ( $total_frags[$n]->length() == $longest_length ) {

            push @return_set,     $total_frags[$n];
            push @cds_return_set, $total_cds_frags[$n];
        }
    }
    return \@return_set, \@cds_return_set;
}

sub split_coding {
    my $sequence = shift;

    for my $n ( 1 .. ( 3 - length($sequence) % 3 ) )
    {    # add Ns as buffers at the end of the sequence
        $sequence .= "N";
    }

 # still,sometimes the last fragment won't match the length of the last peptide,
 # because bioperl is too smart...

    my $previous_end = 0;
    my @sub_seqs     = ();
    for (
        my $codon_start = 0 ;
        $codon_start < length($sequence) ;
        $codon_start += 3
      )
    {

        my $codon_seq = substr( $sequence, $codon_start, 3 );

        if ( $codon_seq =~ /TAA|TGA|TAG/i ) {
            my $seq =
              substr( $sequence, $previous_end,
                ( $codon_start - $previous_end ) );
            push @sub_seqs, $seq;
            if ($verbose) {
                print STDERR "CDS SEQ: $seq\n";
            }
            $previous_end = $codon_start + 3;
        }

    }

    # push last
    my $sequence_end =
      length($sequence) -
      ( length($sequence) % 3 );    # round to the closest coden

    #print STDERR "SEQUENCE END: $sequence_end.\n";
    my $seq = substr( $sequence, $previous_end, $sequence_end - $previous_end );
    if ($seq) { push @sub_seqs, $seq; }

    return @sub_seqs;
}

