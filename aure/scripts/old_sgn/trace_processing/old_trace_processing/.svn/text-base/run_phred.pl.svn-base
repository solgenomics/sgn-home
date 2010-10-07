#!/usr/bin/perl -w
use strict;
#use job_dispatch;
#use fasta_parser;
use IO::File;

if (!$ARGV[0] || $ARGV[0] eq "help") {
  print <<EOF;

    Script for basecalling of chromatograms. Expects a list
    of chromatograms on STDIN. (use find)

    Usage: <list of input directories><output directory><temp file name> <additional argument for phred> <phred arg> <phred arg> ...
    
    The input and output file should be either complete path, or path from the directoty where the script run

EOF
exit -1;
}

my ($input_list, $outputdir, $temp, @phredargs) = @ARGV;

if ( -e $outputdir ) {
  print "File or directory $outputdir already exists\n";
} else {
  system("mkdir -p $outputdir");
}

#print "ENV: $ENV{PHRED_PARAMETER_FILE}\n";

# Used to measure the elapsed time of this script
my ($start_time, $stop_watch);
$start_time = time();


my $current_directory = `pwd`;
chomp $current_directory;

# Ensure that directory names have trailing slashes, for consistency
$outputdir =~ m/\/$/ or $outputdir .= "/";
$current_directory =~ m/\/$/ or $current_directory .= "/";

# Convert a relative path to an absolute path
if ($outputdir !~ m/^\//) {
  $outputdir = $current_directory . $outputdir;
}

my @traces = ();
open INLIST, $input_list or die "Couldn't find list of input files!";
while (<INLIST>){
  chomp;
  push @traces, $_;
}
close INLIST;

my $stop = 0;
for (my $i=0;$i<@traces;$i++) {
  my $id_file = "$outputdir/job-$i.id";;

  open ID_FILE, ">$id_file"
    or die "Failed to open id batch file $id_file ($!)";

 
    print ID_FILE "$traces[$i]\n";
    if ( ! -f $traces[$i] ) {
      print STDERR "Can't find trace $traces[$i] ($!)\n";
      $stop = 1;
    }
  

  close ID_FILE;

  my $morephredargs = join (' ',@phredargs);
  #print "$traces[$i]\n";
  my $command = "phred -id $traces[$i]  -pd $outputdir  -exit_nomatch -tags $morephredargs 2>$outputdir/phred-job-$i.out";

  #print "COMMAND: $command\n";
  system($command);
}

exit -1 if ($stop);


# We have noticed, over time, that chromatograms are mysteriously not processed
# for reasons unexplained. Check the output dir for .phd output files.
#foreach (@traces){
#  my $command = "find $_ -type f > '$temp'";
#  #system($command);
#  open LIST, $temp;
#  while (<LIST>){
#    chomp;
#    my ($in_file) = $_ =~ m/\/(.+\.ab1)$/;

#    if ( ! -f "$outputdir/$in_file.phd.1") {
#      print STDERR "No .phd file for $in_file\n";
#      # so there is a non-zero exit code if this script is ever automated
#      # with another
#      $stop = 1;
#    }
#  }
#}

exit -1 if ($stop);
showtime("Total Basecalling time", time() - $start_time);


sub showtime {
    my ($label, $seconds) = @_;
    my ($hours, $minutes);

    $hours = $seconds / 3600;
    $seconds %= 3600;
    $minutes = $seconds / 60;
    $seconds %= 60;

    printf STDERR "$label:\t %02d:%02d:%02d\n",$hours,$minutes,$seconds;
}


