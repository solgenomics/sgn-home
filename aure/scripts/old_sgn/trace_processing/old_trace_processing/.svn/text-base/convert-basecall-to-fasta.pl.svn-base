#!/usr/bin/perl -w
use strict;
use DBI;
use cmdline_args;

if (!$ARGV[0] || $ARGV[0] eq "help") {
  print <<EOF;

  Program convereted from load-sgn-from-phd.pl, to use the output of 
  parallel-phred.pl, which creates a directory of .phd files, to fasta
  files for processing without loading the sgn database first.

  Usage: -d phd_dir -f output_fname

EOF
  exit -1;
}

my ($output_fname, $phd_dir) = ("","");
my %args;
$args{-f} = [ \$output_fname, 1 ];
$args{-d} = [ \$phd_dir,      1 ];
cmdline_args::get_options(\%args, \@ARGV);

if ( ! $output_fname ) {
  print STDERR "You must specify an output filename with -f \n";
  exit -1;
}
$output_fname =~ s/\.seq$//;

if ( ! $phd_dir ) {
  print STDERR "You must specify the .phd directory with -d \n";
  exit -1;
}

my @trace_names;
while(<STDIN>) {
  chomp;
  s/^.*\///;
  s/\.gz$//;
  push @trace_names, $_;
}

my @fs_notfound;
foreach ( @trace_names ) {
  # This extra qname business is to get rid of chromatogram extension names
  # names like .esd, .scf, and .ab1, because they are typically not entered
  # into the database by scripts that preceed this. Since this is not handled
  # consistently either by us (SGN) or by submitters and facilities, some
  # work arounds are needed pretty much anywhere that handles the files.
  my $qname = $_;

  if (! -f "$phd_dir/$_.phd.1") {
    push @fs_notfound, $_;
  }
}

if ( @fs_notfound>0) {
  print STDERR "Following trace names do not have .phd.1 files:\n";
  print STDERR join("\n",@fs_notfound),"\n";
  print STDERR scalar(@fs_notfound)," not found on filesystem.\n";
  print STDERR "ERRORS DETECTED, FASTA Files not created.\n";
  exit -1;
} else {
  print STDERR "All .phd files found\n";
}

open SF, ">$output_fname.seq"
  or die "Can't open output file \"$output_fname.seq\" ($!)";

open QF, ">$output_fname.qual"
  or die "Can't open output file \"$output_fname.qual\" ($!)";

foreach my $id ( @trace_names ) {
  my $phd_file = "$id.phd.1";

  if (open PHDFILE, "<$phd_dir/$phd_file") {
    my @seq = ();
    my @qual = ();
    my @pos = ();
    my ($phred_version, $version, $clone_id, $read_id);

    while(<PHDFILE>) {
      chomp;
      if (m/^BEGIN_SEQUENCE\ (.+)/) {
	my $sequence_name = $1;
	if ($sequence_name ne $id) {
	  print STDERR "Warning: Sequence name in $phd_file \"$sequence_name\" does ".
	    "not match file name for $id.\n";
	}
      }
      if (m/^PHRED_VERSION: (.+)/) {
	$phred_version = $1;
      }
      if (m/^BEGIN_DNA/) {
	while(<PHDFILE>) {
	  chomp;
	  last if m/^END_DNA/;
	  my ($base,$qual,$pos) = split;
	  push @seq, uc $base;
	  push @qual, $qual;
	  push @pos, $pos;
	}
      }
      last if m/^END_DNA/;
    }

    close PHDFILE;

    my $header = ">$id\n";
    print SF $header;
    print SF join("",@seq),"\n";
    print QF $header;
    print QF join(" ",@qual),"\n";
	
  } else {
    print STDERR "Can't open PHD file \"$phd_file\" ($!)\n";
  }

}

close SF;
close QF;
