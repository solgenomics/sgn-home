#!/usr/bin/perl -w
use strict;
use CXGN::BioTools::FastaParser;

if (!defined($ARGV[0]) || $ARGV[0] eq "help") {
    print STDERR <<EOF;

     Quick program to take the result of a CAP3 unigene build (all the separate
     files) and concatenate the results into a single unigene FASTA file

     Expects name of directory containing input to CAP3 build as FASTA
     files, typically the output of the preclustering step in tmp/cap3/

     Usage: <cap3 working directory>

     Outputs files unigene.seq and unigene.qual in current directory


     If singletons.seq and singletons.seq.qual are present (these are the
     singletons from preclustering), they will be added to the output files.
EOF
    exit 0;
}

my $dir = shift;
if (! -d "$dir" ) {
  die "\"$dir\" is not a directory or is inaccessible\n";
}

if ($dir !~ m/\/$/) {
  $dir .= "/";
}




open UNIGENESEQ, ">unigene.seq";
open UNIGENEQUAL, ">unigene.qual";

print STDERR "Adding the contig data...\n";

open FINDPIPE, "find ${dir} -name \"cluster-*.seq.cap.contigs\" |"
    or die "Can't open pipe to 'find' program ($!)";

my ($cluster_count, $singletons_count, $precluster_singletons_count) ;
my ($cluster_no);
while(<FINDPIPE>) {
    chomp;
    ($cluster_no) = m/cluster\-([0-9]+)\.seq/;

#    print STDERR "Using file $_ and qual file $_.qual\n";
    my $fh = CXGN::BioTools::FastaParser->load_file($_,$_.".qual");

    foreach ( $fh->get_seqids() ) {
	my @qual = ();
	my ($contig_no) = m/^Contig([0-9]+)/i;
	
	$cluster_count++;

	my $seq = $fh->get_sequence($_,\@qual);
	print UNIGENESEQ ">$cluster_no-$contig_no\n",$seq,"\n";
	print UNIGENEQUAL ">$cluster_no-$contig_no\n",join(" ",@qual),"\n";
    }

    undef $fh;
}

close FINDPIPE;

print STDERR "Adding the cap3 singletons...\n";

# Load assembly program generated singletons
open FINDPIPE, "find ${dir} -name \"cluster-*.seq.cap.singlets\" |"
    or die "Can't open pipe to 'find' program ($!)";

my ($cluster_file, $seq, $qual);

my (@qual);
while(<FINDPIPE>) {
    chomp;
    ($cluster_no) = m/cluster\-([0-9]+)\.seq/;
    ($cluster_file) = m/(.+?cluster\-[0-9]+.seq)/;
    my $qual_file = $cluster_file;
    $qual_file =~ s/seq$/qual/;

    # We have to open $_ (.singlets) because CAP3 doesn't put out the 
    # quality file. Thus, we read the seqids from S, but use $fh to get the
    # sequence and quality.
    my $fh = CXGN::BioTools::FastaParser->load_file($cluster_file, $qual_file); #"$cluster_file.qual");

    open S, "<$_"
      or die "Failed to open singlets file for cluster \"$cluster_file\" ($!)";

    while(<S>) {
      chomp;

      if (m/^>([^\n\t\ ]+)/) {

	my @qual = ();
	my $seq = $fh->get_sequence($1, \@qual);
	$singletons_count++;
	print UNIGENESEQ ">$cluster_no-S-$1\n", $seq, "\n";
	print UNIGENEQUAL ">$cluster_no-S-$1\n" . join(" ",@qual) . "\n";
      }

    }

    close S;
    undef $fh;
}

print STDERR "Adding the precluster singletons...\n";

if ( -f "${dir}singletons.seq" && -f "${dir}singletons.qual" ) {
  my $fh = CXGN::BioTools::FastaParser->load_file("${dir}singletons.seq",
				   "${dir}singletons.qual");
  if ($fh) {
    print STDERR "Appending preclustering singlets...\n";
    foreach my $seqid ( $fh->get_seqids() ) {
      my @qual;
      my $seq = $fh->get_sequence($seqid, \@qual);
      $precluster_singletons_count++;
      print UNIGENESEQ ">X-S-$seqid\n$seq\n";
      print UNIGENEQUAL ">X-S-$seqid\n" . join(" ",@qual) . "\n";
    }
    undef $fh;
  } else {
    die "No preclustering singlets were found. Aborting script. Files should be named singletons.seq and singletons.qual.\n";
  }

}


close UNIGENESEQ;
close UNIGENEQUAL;


print STDERR "Found $cluster_count sequences in clusters. $singletons_count singletons in clusters and $precluster_singletons_count from the precluster.\n";
print STDERR "Done.\n";

