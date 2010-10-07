#!/usr/bin/perl

use strict;
use CXGN::Tools::Run;

my $dir = shift;

if (!$dir) { 
    die "Need a directory. Stop.\n";
}

my $CAP3_OPTIONS = "-e 5000 -p 90 -d 10000 -b 60";

print "CAP3 Assembly Options: $CAP3_OPTIONS\n";

my $path = `pwd`;
chomp($path);
print STDERR "PATH: $path\n";

my @files = `find -L . -iname "*cluster-*.seq"`;

# run cap3 on each file.
#

my $jobs = ();
foreach my $f (@files) { 
    chomp($f);
    print STDERR "Processing $path/$f...\n";
    my $jobid;
    my $command = "/data/shared/bin64/cap3 $path/$f $CAP3_OPTIONS";

    my $q = $f;
    $q=~s/seq$/qual/;

    if ( -e "$f.cap.contigs") { 
	print STDERR "$f already processed. Skipping...\n";
	next();
    }
    
    if (! -e "$q") { 
	print STDERR "Warning! $q qual file not found for file $f.\n";
    }
    
    system "$command";

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

 
