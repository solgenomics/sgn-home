#!/usr/bin/env perl
use strict;
use warnings;
use 5.10.0;

use lib '/home/rob/dev/bioperl/Bio-FeatureIO/lib';

use IO::File;

use List::AllUtils qw( part );

use Bio::Chado::Schema;
use Bio::FeatureIO;

my $bcs_dsn = shift;

my $feature_in = Bio::FeatureIO->new(
    -format  => 'gff',
    -version => 3,
    -file    => shift,
    );

open my $fasta_in, '<', shift
    or die "cannot open genome fasta file (3rd argument) for reading";

# convert the gff3 to tbl
my %handlers = (
    gene => \&handle_gene,
);
while( my @f = $feature_in->next_feature_group ) {
    for my $feature ( @f ) {
        my $fh = get_tbl_output_fh( $feature )
            or next;

        no strict 'refs';

        my $type = $feature->primary_tag;
        my $handler = $handlers{$type}
            or die "no handler defined for top-level $type features";
        $handler->( $fh, $feature );
    }
}

# write corresponding .fsa files with the fasta sequences of each chromosome
{
    my $fh;
    while ( my $line = <$fasta_in> ) {
        if ( my ($chr_name) = $line =~ /^>\s*(\S+)/ ) {
            if ( my $accession = accession_for( $chr_name ) ) {
                print "\nwriting $accession.fsa\n";
                $fh = IO::File->new( "$accession.fsa", '>' ) or die "$! writing $accession.fsa";
                $fh->print( ">gb|$accession|\n" );
            } else {
                warn "skipping fasta output for $chr_name, I have no genbank accession for it";
                $fh = undef;
            }
        } elsif ( $fh ) {
            $fh->print( $line );
        }
    }
}

exit;

################ subroutines ################

sub handle_gene {
    my ( $fh, $gene ) = @_;

    $gene->primary_tag eq 'gene'
        or die "not a gene?";

    $fh->print_tab(
        [ $gene->start, $gene->end, 'gene' ],
        format_gene_attributes( $gene ),
    );

    for my $mrna ( $gene->get_SeqFeatures ) {
        $mrna->primary_tag eq 'mRNA'
            or die "second-level feature not an mRNA?";

        my @exon_tabs =
            map [ $_->start, $_->end, undef ],
            grep $_->primary_tag eq 'exon',
            $mrna->get_SeqFeatures;
        $exon_tabs[0][2] = 'mRNA';

        push @exon_tabs,
            format_mrna_attributes( $mrna );

        $fh->print_tab( @exon_tabs );

        my @cds_tabs =
            map [ $_->start, $_->end, undef ],
            grep $_->primary_tag eq 'CDS',
            $mrna->get_SeqFeatures;
        $cds_tabs[0][2] = 'CDS';

        push @cds_tabs,
            format_cds_attributes( $mrna );

        $fh->print_tab( @cds_tabs );
    }
}


sub format_mrna_attributes {
    my ( $mrna ) = @_;
    my ( $locus_tag ) = $mrna->get_tag_values( 'Name' );
    $locus_tag =~ s/(\.\d+)/.x/;

    my $prefix = "gnl|WGS:AEKE|";
    my @attributes = (
        [ protein_id    => $prefix.$locus_tag ],
        [ transcript_id => $prefix."mRNA:$locus_tag" ],
    );

    if( my ( $desc ) = eval { $mrna->get_tag_values( 'Note' ) }) {
        if( my ( $product ) = $desc =~ /^(.+)(?=\(AHRD)/ ) {
            $product =~ s/^\s+|\s+$//g;
            $product =~ s/(\w+)/munge_product_word($1)/eg;
            push @attributes, [ product => $product ];
        }

        if( my @interpro = $desc =~ /\b(IPR\d+)\b/g ) {
            push @attributes, map [ db_xref => "InterPro:$_" ], @interpro;
            push @attributes, map [ inference => "protein motif:InterPro:$_" ], @interpro;
        }

        if( my ( $ev_desc, $similar ) = $desc =~ m! \(AHRD \s+ V\d+ \s+ ([-\*]+) \s+ ([^\s\)]+) \) !x ) {
            # currently not using the evidence desc
            push @attributes, [ note => "similar to $similar" ];
        }
    }

    if( my @terms = eval { $mrna->get_tag_values( 'Ontology_term' ) }) {
        my ( $goterms, $other ) = part { /^GO:/ ? 0 : 1 } @terms;
        push @attributes, map format_go_attributes( $_ ), @{ $goterms || [] };
        push @attributes, map [ db_xref => $_ ], @{ $other || [] };
    }

    return format_attributes( \@attributes );
}

{ my $schema;
  sub bcs {
      $schema ||= Bio::Chado::Schema->connect( $bcs_dsn );
  }

  my $go_db; sub go_db { $go_db ||= bcs->resultset('General::Db')->find({ name => 'GO' }) }
}

sub munge_product_word {
    my $w = shift;
    # lowercase non-acronym words that are 3 chars or more
    $w = lc $w if $w =~ m![A-Z][a-z]{2,}!;
    # replace 'unknown' with 'hypothetical'
    return 'hypothetical' if $w eq 'unknown';
    return 'putative' if $w =~ / ^ ( probable | possible | predicted ) $ /x;
    return $w;
}

sub format_go_attributes {
    my ( $go_term ) = @_;
    my ( $go_num ) = $go_term =~ /(\d+)/;
    my $dbxref = go_db()->search_related('dbxrefs', { accession => $go_num } );

    my @matching_terms = do {
        my %seen_term;

        grep !$seen_term{$_->cvterm_id}++, (

            $dbxref->search_related('cvterm', undef, { prefetch => 'cv' } ),

            $dbxref->search_related( 'cvterm_dbxrefs' )
                   ->search_related( 'cvterm', undef, { prefetch => 'cv' } ),

        )
    };

    die "no cvterm found for $go_term" unless @matching_terms;

    return map {
        my $cvterm = $_;

        my $go_qualifier = {
            biological_process => 'go_process',
            molecular_function => 'go_function',
            cellular_component => 'go_component',
        }->{ $cvterm->cv->name } || die "unknown cv ".$cvterm->cv->name;

        my $name = $cvterm->name;
        $name =~ s/\(obsolete.+//;

        [ $go_qualifier => join( '|', (
            $name,
            $go_num,
            '',
            'ISS',
          )),
        ]
    } @matching_terms;
}


sub format_cds_attributes {
    my ( $mrna ) = @_;
    return format_mrna_attributes( $mrna );
}
sub format_gene_attributes {
    my ( $gene ) = @_;
    return format_attributes([
        map {
            ( my $tag = $_ ) =~ s/Solyc/Solyc_/;
            [ locus_tag => $tag ]
        } $gene->get_tag_values('Alias')
    ]);
}

sub format_attributes {
    my ( $attributes ) = @_;

    return
        map [ undef, undef, undef, @$_ ],
        sort { $a->[0] cmp $b->[0] }
        @$attributes;
}


use Data::Dumper;

my %filehandle_cache;
sub get_tbl_output_fh {
    my ( $feat ) = @_;
    my $seq_id = $feat->seq_id or die "no seq_id on feature!";

    return $filehandle_cache{$seq_id} //= do {

        my $accession = accession_for( $seq_id );

        unless( $accession ) {
            warn "I don't know the genbank accession for $seq_id, skipping all features on that sequence";
            0
        } else {
            my $fh = mytabfile->new( "$accession.tbl", '>' )
                or die "$! writing $accession.tbl";

            $fh->print( ">Feature gb|$accession|\n" );
            $fh
        }
    };
}

sub accession_for {
    my ( $id ) = @_;

     my $accession = {qw(
        SL2.40ch01	CM001064.1
        SL2.40ch02	CM001065.1
        SL2.40ch03	CM001066.1
        SL2.40ch04	CM001067.1
        SL2.40ch05	CM001068.1
        SL2.40ch06	CM001069.1
        SL2.40ch07	CM001070.1
        SL2.40ch08	CM001071.1
        SL2.40ch09	CM001072.1
        SL2.40ch10	CM001073.1
        SL2.40ch11	CM001074.1
        SL2.40ch12	CM001075.1
    )}->{ $id };

    return $accession;
}


package mytabfile;
use base 'IO::File';

sub print_tab {
    my $self = shift;
    no warnings 'uninitialized';
    $self->print( join( "\t", @$_ ), "\n" ) for @_;
}

