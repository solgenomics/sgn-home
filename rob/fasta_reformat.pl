#!/usr/bin/env perl
use strict;
use warnings;
use English;
use Carp;
use FindBin;
use Getopt::Std;

#use Data::Dumper;
use Bio::SeqIO;

########### CONFIG/DEFAULTS ##########

our $print_width = 60;

######################################

sub usage {
  my $message = shift || '';
  $message = "Error: $message\n" if $message;
  die <<EOU;
$message
Usage:
  $FindBin::Script [options]  seq_file seq_file

  For each sequence or qual file given, reformats it by opening it
  with Bio::SeqIO and writing it back out to STDOUT.  This makes the
  sequence lines nicely wrapped, etc.  Also, using the options below,
  you can do other stuff like change the identifiers, and the
  deflines, and reverse-complement the sequence.

  If you give it multiple sequence files, it's sort of like a
  reformatting 'cat' command.

  Options:

    -i <pattern>
      run the given perl operation on the sequence identifiers

    -l <pattern>
      run the given perl operation on the sequence deflines

    -s <pattern>
      run the given perl operation on the sequences themselves

    -r reverse complement all the sequences

    -w <width>
      printing width for output
      Default: $print_width

    -q force treating files as qual

   Example:
     $FindBin::Script -l 's/foo/bar/; s/bar/baz/;' myfile.seq > myfile2.seq

EOU
}
sub HELP_MESSAGE {usage()}

our %opt;
getopts('i:l:s:rw:q',\%opt) or usage();
@ARGV or usage();

#check whether we got qual files
my $qual_mode = $opt{q} || do {
  open my $f, $ARGV[0] or die "$! opening $ARGV[0]";
  my @twolines;
  while(@twolines <2 and my $l = <$f>) {
    push @twolines,$l if $l =~ /\S/;
  }

  $twolines[0] =~ /^>/ && $twolines[1] =~ /^[ \d]+\n$/;
};
#warn "qual_mode is $qual_mode\n";

#compile each of the options into a subroutine
foreach my $x (qw/ i l s/) {
  next unless $opt{$x};

  my $s = eval "sub { no strict;  \$_ = shift; $opt{$x}; return \$_ }"
    or die "error compiling -$x '$opt{$x}' ($EVAL_ERROR)";

  $opt{$x} = $s;
}

$print_width /= 3 if $qual_mode;
$opt{w} ||= $print_width;

my $out = Bio::SeqIO->new( -fh => \*STDOUT, -format => $qual_mode ? 'qual' : 'fasta');
$out->width($opt{w}) unless $qual_mode;
foreach my $file (@ARGV) {
  my $in = Bio::SeqIO->new( -file => $file, -format => $qual_mode ? 'qual' : 'fasta');
  while(my $s = $in->next_seq) {
#     use Data::Dumper;
#     die Dumper $s;
    if($opt{i}) {
      munge( $s, display_id => $opt{i} );
    }

    if($opt{l}) {
      munge( $s, description => $opt{l} );
    }

    if($opt{s}) {
      if($qual_mode) {
	qualmunge( $s, $opt{s} );
      } else {
	munge( $s, seq => $opt{s} );
      }
    }

    $s = $s->revcom if $opt{r};

    $out->write_seq($s);
  }
}

sub qualmunge {
  my ($s,$op) = @_;
  my $result = $op->(join ' ', @{$s->qual});
  $s->qual($result);
}
sub munge {
  my ($s,$prop,$op) = @_;
  my $result = $op->($s->$prop());
  #warn "munge ".$s->$prop()." into $result\n";
  $s->$prop($result);
}
