#!/usr/bin/perl -w
use strict;
use CXGN::Garbage::job_dispatch;
use CXGN::BioTools::FastaParser;
use IO::File;

##### CONFIGURATION VARIABLES ######
my $phred_paramfile = $ENV{PHRED_PARAMETER_FILE} || '/usr/local/etc/phredpar.dat';
####################################

if (!$ARGV[0] || $ARGV[0] eq "help") {
  print <<EOF;

    Script for parallelizing basecalling of chromatograms. Expects a list
    of chromatograms on STDIN. (use find)

    Usage: <output directory> <additional argument for phred> <phred arg> <phred arg> ...

EOF
exit -1;
}

my ($outputdir,@phredargs) = @ARGV;

if ( -e $outputdir ) {
  print "File or directory $outputdir already exists\n";
} else {
  system("mkdir -p $outputdir");
}

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

my @traces = <STDIN>;
chomp @traces;

my $batch_size = int(@traces/28 + 1);
my @jobs = ();
my $stop = 0;
my ($i, $job) = (0,0);
while($i<@traces) {
  my $id_file = "$outputdir/job-$job.id";;

  open ID_FILE, ">$id_file"
    or die "Failed to open id batch file $id_file ($!)";

  for(my $j=0;$i<@traces && $j<$batch_size;$j++,$i++) {
    print ID_FILE "$traces[$i]\n";
    if ( ! -f $traces[$i] ) {
      print STDERR "Can't find trace $traces[$i] ($!)\n";
      $stop = 1;
    }
  }

  close ID_FILE;

  my $morephredargs = join (' ',@phredargs);
  my $exportparamfile = "export PHRED_PARAMETER_FILE=$phred_paramfile;";
  my $command = "$exportparamfile ( phred -if $id_file -zt . -pd $outputdir -exit_nomatch -tags $morephredargs 2>&1 > $outputdir/phred-job-$job.out )";

  push @jobs, [ $current_directory, $command, "phred job $job"];
  $job++;
}

exit -1 if ($stop);

my $failed = CXGN::Garbage::job_dispatch::dispatch_jobs(\@jobs);
if ($failed) {
  print STDERR "$failed base calling jobs failed. aborting script.\n";
  exit -1;
}

# We have noticed, over time, that chromatograms are mysteriously not processed
# for reasons unexplained. Check the output dir for .phd output files.
foreach ( @traces ) {
  my $id = substr $_, (rindex $_, "/")+1;
  $id =~ s/\.gz$//;

  if ( ! -f "$outputdir/$id.phd.1") {
    print STDERR "No .phd file for $id\n";
    # so there is a non-zero exit code if this script is ever automated
    # with another
    $stop = 1;
  }
}

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

my $dev_urandom = "";

sub tempname {

  if (! $dev_urandom) {
    $dev_urandom = new IO::File "</dev/urandom";
  }

  my $rand_string = "";

  $dev_urandom->read($rand_string,16);

  my @bytes = unpack("C16",$rand_string);

  $rand_string = "";
  foreach ( @bytes ) {
    $_ %= 62;
    if ($_<26) {
      $rand_string .= chr(65+$_);
    }
    elsif ($_<52) {
      $rand_string .= chr(97+($_-26));
    }
    else {
      $rand_string .= chr(48+($_-52));
    }
  }

  return $rand_string;
}

