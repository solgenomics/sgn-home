#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use autodie ':all';

use File::Temp;
use Bio::SeqIO;
use Bio::GFF3::LowLevel 'gff3_parse_feature';

my $regions = <<'';
chr01 31453743 33253742 063 W PGSC0003DMB000000023 1 1800000 0
chr04 61427953 62187453 139 W PGSC0003DMB000000315 1 0759501 0
chr04 62237454 63783430 141 W PGSC0003DMB000000129 1 1545977 0
chr06 46908540 48443546 107 W PGSC0003DMB000000133 1 1535007 0
chr08 23937447 24327299 045 W PGSC0003DMB000000496 1 0389853 0
chr10 41275334 42771946 079 W PGSC0003DMB000000138 1 1496613 0
chr11 07501774 08921477 017 W PGSC0003DMB000000155 1 1419704 0
chr06 31762806 32246757 069 W PGSC0003DMB000000440 1 0483952 0

my %regions;
while( $regions =~ /(chr\d\d)\s(\d+)\s(\d+)/g ) {
     push @{$regions{$1}}, [$2,$3];
}

for my $fasta_file ( glob '{repeats,seq}/*.fasta' ) {
     my $in  = Bio::SeqIO->new( -file => $fasta_file, -format => 'fasta' );
     my $tempfile = File::Temp->new; $tempfile->close;
     my $width = $fasta_file =~ /repeats/ ? 50: 60;
     my $out = Bio::SeqIO->new( -file => ">$tempfile", -format => 'fasta', -width => $width );
     my $changed = 0;
     while( my $seq = $in->next_seq ) {
          for my $r ( @{ $regions{$seq->id} || [] } ) {
               my $start  = $r->[0] - 1;
               my $length = $r->[1] - $start;
               $seq->seq( do{ my $s = $seq->seq; substr( $s, $start, $length, 'N' x $length ); $s } );
               $changed = 1;
          }
          $out->write_seq( $seq );
     }
     if( $changed ) {
          system mv => $tempfile => $fasta_file;
     }
}

# now drop features from gff3 files that are in those regions
GFF3:
for my $gff3_file ( glob '*/*.gff3' ) {
     open my $g, '<', $gff3_file;
     my $temp = File::Temp->new;
     my $changed = 0;
     while( my $l = <$g> ) {
          if( $l !~ /^#/ && $l =~ /^(\S+)\t\S+\t\S+\t(\d+)\t(\d+)/ && (my $r = $regions{$1}) ) {
               for ( @$r ) {
                    if( $2 >= $_->[0] && $2 <= $_->[1] || $3 >= $_->[0] && $3 <= $_->[1] ) {
                         $l = '';
                         $changed = 1;
                    }
               }
          }
          $temp->print($l);
     }
     if( $changed ) {
          system mv => $temp => $gff3_file;
     }
}
