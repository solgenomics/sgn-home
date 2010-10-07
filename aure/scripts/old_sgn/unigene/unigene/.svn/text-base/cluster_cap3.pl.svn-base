#!/usr/bin/perl

use strict;
use CXGN::Tools::Run;
use File::Spec;


my $dir = shift;
my $pattern = shift;


if (!$dir) { 
    die "Need a directory. Stop.\n";
}

if (!$pattern) { 
    $pattern = "*cluster-*.seq";
}

my $CAP3_OPTIONS = "-e 5000 -p 90 -d 10000 -b 60";

print "CAP3 Assembly Options: $CAP3_OPTIONS\n";

my $path = `pwd`;
chomp($path);

$path = File::Spec->catfile($path, $dir);

print STDERR "PATH: $path\n";

my @files = `find -L $path -iname "$pattern"`;

# run cap3 on each file.
#

my $jobs = ();
foreach my $f (@files) { 
    chomp($f);


    my $jobid;
    my $command = "/data/shared/bin64/cap3 $f $CAP3_OPTIONS";

    print STDERR "Processing command $command...\n";

    my $q = $f;
    $q=~s/seq$/qual/;

    if ( -e "$f.cap.contigs") { 
	print STDERR "$f already processed. Skipping...\n";
	next();
    }
    
    if (! -e "$q") { 
	print STDERR "Warning! $q qual file not found for file $f.\n";
    }
    
    CXGN::Tools::Run->run ("qsub", { in_file =>  $command, 
				     err_file => $f.".error",
				     out_file => $f.".stdout" });
    sleep (2);
    chomp $jobid;
    push @$jobs, $jobid;
}


close (F);
  
spin_a_while($jobs);

print STDERR "Done.\n";



sub spin_a_while {
  my $jobs = shift;
  # Don't flip out over this use of goto.  We want to get
  # out of the qstat stream as soon as a record matches 
  # a job id.
 QSTAT: while (1) {
    open (my $qstat,"qstat -f | grep 'Job Id:' |")
      or die ("failed to run qstat: $!");
    while (my $job = <$qstat>) {
      chomp $job;
      $job =~ s/Job.Id:\s+//; #
      if (grep { $_ eq $job; } @$jobs) {
	#print STDERR "found job $job\n";
	close ($qstat);
	sleep (60);
	goto QSTAT;
      }
    }
    # If we got here, there were no matching jobs.
    print STDERR "No matching jobs; task complete.\n";
    last QSTAT;
  }
}

 
