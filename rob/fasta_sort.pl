#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;

use File::Temp qw/tempfile/;

use Bio::Index::Fasta;
use Bio::Index::Qual;
use Bio::SeqIO;

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options] fasta_file fasta_file ...

  Sort a FASTA-format sequence or qual file alphabetically by identifier name.

  Options:

    -q treat inputs as qual files (defaults to false, meaning input is
       treated as nucleotide or protein sequence)

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('q',\%opt) or usage();

#read the identifiers from the files
my @idents;
foreach my $file (@ARGV) {
  open my $f, "grep '>' $file |" or die "$! grepping '$file'";

  while( my $line = <$f> ) {
    chomp $line;
    my ($ident) = $line =~ /^\s*>\s*(\S+)/
      or die "invalid identifier line '$line'";
    push @idents, $ident;
  }
}

my (undef,$index_file) = tempfile(File::Spec->catfile( File::Spec->tmpdir, 'fasta-sort-index-XXXXXX'), UNLINK => 1);

my $inx_class = $opt{q} ? 'Bio::Index::Qual' : 'Bio::Index::Fasta';
my $inx = $inx_class->new(-filename => $index_file,
				 -write_flag => 1);
$inx->make_index(@ARGV);

@idents = sort { $a cmp $b} @idents;

my $out = Bio::SeqIO->new( -format => ($opt{q} ? 'qual' : 'fasta'), -fh => \*STDOUT );
foreach my $ident (@idents) {
  my $s = $inx->fetch($ident)
    or die "sanity check failed, could not find seq '$ident' in index!!!";
  $out->write_seq( $opt{q} ? (-source => $s) : ($s) );
}

$out->close;


