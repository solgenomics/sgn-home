#!/usr/bin/perl


=head1 NAME

fasta_grep.pl - a utility to quickly extract stuff from a fasta file.

=head1 SYNOPSIS

fasta_grep "string" fastafile

=head1 VERSION



=head2 Version history



=head1 DESCRIPTION


=cut


use Getopt::Std;
use Bio::Seq;
use Bio::SeqIO;

use strict;

use vars qw($opt_h $opt_f $opt_v $opt_i $opt_r);

my $nooptions = !scalar(@ARGV);

getopts('hf:vir');

if ($opt_h) { 
    print <<HELP;

fasta_grep.pl [-hi] (-f patternfile | pattern) fastafile

-h: help
-i: match precisly on id only
-f: use patternfile as pattern source
-r: for repeat masker like ids, only match id part (before #)

HELP

   exit(0);
}

my @patterns;
if ($opt_f) { 
    # load patternfile
    if ($opt_v) { print STDERR "Reading pattern file...\n"; }
    open(F, "<$opt_f") || die "Can't open patternfile $opt_f"; 
    while (<F>) { 
	chomp;
	if ($opt_v)  { print STDERR "$_\n"; }
	push @patterns, $_;
    }
    close (f);
}
else { 
    $patterns[0]=shift;
    if ($opt_v) { print STDERR "Using pattern: $patterns[0]\n"; }
}

my $filename = shift;

my $file;
my $accno;
my $desc;
my $seq;
my $id;
my $seqobj;


my $seqcount=0;
my $match_count =0;

my $output_file = Bio::SeqIO -> new('-format'=>'Fasta');

$file = Bio::SeqIO -> new ('-format' => 'Fasta', -file => $filename);
while ($seqobj = $file -> next_seq() ) {
  # Get line and chop it up in little pieces...
  
  $accno = $seqobj -> accession_number();
  $desc  = $seqobj -> desc(); 
  $seq  = $seqobj -> seq();
  $id   = uc($seqobj -> id());


  $seqcount++;
  if (($seqcount % 500) == 0) { print STDERR "Reading seq $seqcount\r";}

  my $match_flag=0;
  foreach my $search_pattern (@patterns) { 
      $search_pattern =~ s/\///g;
      $id=~s/\///g;
      $search_pattern =~s/\|//g;
      $id=~s/\|//g;
      if ($opt_r) { 
	  $search_pattern =~ s/(.*)\#.*/$1/g;
	  $id=~s/(.*)\#.*/$1/g;
      }
#      print STDERR "Comparing patter $search_pattern to id $id\n";
      if ($opt_i) {
	  if ($opt_v) { 
	      print STDERR "checking pattern $search_pattern against $id\n";
	  }
	  if ($search_pattern && ($id =~ /^$search_pattern$/i)) { 
	      $output_file->write_seq($seqobj);
	      $match_flag=1;
	  }
      }
      elsif (($id =~ /$search_pattern/i) || ($desc=~/$search_pattern/i)) {
	  $output_file -> write_seq($seqobj);
	  $match_flag=1;
      }
  }
  if ($match_flag) { $match_count++; }
}

$file -> close();
$output_file -> close();

print STDERR "Matched $match_count sequences out of $seqcount.\n\n";
