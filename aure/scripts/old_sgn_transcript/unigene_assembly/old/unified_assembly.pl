#!/usr/bin/perl -w
use strict;

use CXGN::Tools::Run;

use CXGN::Garbage::cmdline_args;
#use CXGN::Garbage::job_dispatch;
use CXGN::BioTools::FastaParser;

my $CAP3_OPTIONS = "-e 5000 -p 90 -d 10000 -b 60";

print "CAP3 Assembly Options: $CAP3_OPTIONS\n";

# Command line options: required for this script, or we bail reporting
# expected usage to the user
my ($chromadir, $outputdir, $vector_filename, $adaptor_sequence) = (0,0,0,0);

my %opts;
$opts{-o} = [ \$outputdir, 1 ];

usage() if (! defined($ARGV[0]) or $ARGV[0] eq "help");
usage() if CXGN::Garbage::cmdline_args::get_options(\%opts, \@ARGV);
 
$outputdir or usage();

# Used to measure the elapsed time of this script
my ($start_time, $stop_watch);
$start_time = time();

# In the Konix school of software, there are 4 /variables/ denoting
# directories that get swapped around, in addition to the process's
# cwd state.  Nonsense!  ALL DIRECTORY NAMES SHALL BE DENOTED BY THIS
# VARIABLE, WHICH IS DYNAMICALLY SCOPED.
local $main::wd = 
    canonicalize_directory_name(ensure_absolute_directory_name($outputdir));

system("mkdir -p ${main::wd}tmp");

# Randomize order of sequences in file for homology scanning
#precluster_sequences();

# Format homology scan database and lookup files

# Execute parallel homology scans

# Collect homology scan outputs into a single file

# Convert to binary adjanceny matrix

# Run clustering to scan for articulations

# Scan articulations for evidence of chimeras

# Re-run clustering censoring chimeras

# make FASTA cluster files from clustering output

# Run CAP3 in parallel to produce contigs

make_contigs();
# Collect singlets

# Collect contigs

# Collect ACE files, deposit into ACE directory

# Construct composite unigene file and quality file

# Done
showtime("Total Unigene Build time", time() - $start_time);
my ($unigenes, $contigs, $singlets);
$unigenes = qx'cat ${main::wd}unigene/unigene.seq | egrep -e "^>" | wc -l';
$contigs = qx'cat ${main::wd}unigene/unigene.seq | egrep -e "^>[0-9]+\-[0-9]+([\ \n\t]|$)" | wc -l';
chomp $unigenes;
chomp $contigs;
$singlets = $unigenes - $contigs;

print STDERR "$unigenes unigenes, in $contigs contigs and $singlets singlets\n";

### Directory name munging stuff
sub ensure_absolute_directory_name {
    my $dirname = shift;
    unless ($dirname =~ m|^/|) {
	my $d = `pwd`;
	chomp $d;
	warn ("Directory name $dirname is not absolute.  Merging with $d.");
	return ($d."/".$dirname);
    }
    return ($dirname);
}
sub canonicalize_directory_name {
    my $dirname = shift;
    unless ($dirname =~ m|^/|) {
	die ("Relative directory name $dirname passed to canoncicalize_directory_name");
    }
    # First normalize any repeated slashes.
    $dirname =~ s|/+|/|;
    # Then get rid of dotdots and dots
    my @dirbits = split ("/", $dirname);
    my @realbits = ();
    for my $bit (@dirbits) {
	#print STDERR "bit>> ".$bit."\n";
	if ($bit eq "") {
	    ; # skip it.
	} elsif ($bit eq ".") {
	    ; # skip it.
	} elsif ($bit eq "..") {
	    pop (@realbits);
	} else {
	    push (@realbits, $bit);
	}
    }
    my $ret = "/".(join "/", @realbits)."/";
    return ($ret);
}
sub enqueue_job {
  my ($command) = @_;
  my $jobid;
  print STDERR "Enqueueing: $command\n";
  CXGN::Tools::Run->run ("qsub", { in_file => \$command, out_file => \$jobid });
  sleep (2);
  chomp $jobid;
  #print ">>$jobid\n";
  return ($jobid);
}

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

sub precluster_sequences {
  my $stop_watch = time();
  my ($i, $j, $temp, $seq_id, $seq);
  # Randomize database ordering. This is because some libraries
  # may be much more "self matching" than others, and outputs of
  # basecalling is usually ordered by library. Thus, splitting
  # the homology scan into chunks can leave some chunks with much more
  # work to do than others, even if they are equally sized. Randomizing
  # the order of sequences before splitting into chunks avoids this
  my $fh = CXGN::BioTools::FastaParser->load_file("${main::wd}basecall/qc-passed.seq", undef);
  my @seqids = $fh->get_seqids();
  srand;
  for ($i=0;$i<@seqids;$i++) {
    $j = int(rand($#seqids+1));
    next if ($i == $j);
    $temp = $seqids[$i];
    $seqids[$i] = $seqids[$j];
    $seqids[$j] = $temp;
  }

  {
    local $main::wd = canonicalize_directory_name("${main::wd}tmp/precluster");
    print STDERR "wd>> $main::wd\n";
    system("mkdir -p $main::wd");

    {
      open (my $output_seq, ">${main::wd}qc-passed-rand.seq")
	or die "Failed opening randomized sequence file ($!)";

      foreach $seq_id ( @seqids ) {
	$seq = $fh->get_sequence($seq_id, undef);

	print $output_seq ">$seq_id\n$seq\n";
      }

      close ($output_seq);
    }

    system("format_seqdata --seqfile=${main::wd}qc-passed-rand.seq --basename=${main::wd}kp_cluster") 
      and die "Failed formatting sequence data for preclustering ($!)";
    system("format_lookup --database=${main::wd}kp_cluster --memsize=8")
      and die "Failed formatting lookup files for preclustering ($!)";

    {
      open (my $findpipe, "find $main::wd -name \"kp_cluster.lt.*\" |")
	or die ("Failed opening pipe to \"find\" ($!)");
      my $badcount;
      while (<$findpipe>) {
	chomp;
	unless ( -s $_ ) {
	  warn ("file $_ is empty");
	  $badcount++;
	}
      }
      if ($badcount) {
	die ("$badcount lookup tables are empty.  Quitting");
      }
      close ($findpipe);
    }

    {
      my $command;
      my @jobs = ();
      open (my $findpipe, "find $main::wd -name \"kp_cluster.lt.*\" |")
	or die "Failed opening pipe to \"find\" ($!)";
      while (<$findpipe>) {
	chomp;
	$command = "(/home/kreuter/sgn-tools/unigene/Kp/scan_sequences --seqfile=${main::wd}kp_cluster --lookupfile=$_ > $_-result)";
	my $jobid = enqueue_job ($command);
	chomp $jobid;
	push @jobs, $jobid;
      }
      close ($findpipe);
      unless (@jobs) {
	die ("No jobs enqueued.");
      }
      # All jobs have been enqueued.
      # Loop until all jobs are done.
      spin_a_while (\@jobs);
    }

    # {
    #     #unlink ("${main::wd}all-scan-results.txt");
    #       while (1) {
    # 	my $len = `wc -l ${main::wd}*-result | tail -1 | cut -d' ' -f2`;
    # 	chomp $len;
    # 	#$len =~ s/\s.*//;
    # 	open (my $findpipe, "find $main::wd -name \"kp_cluster.lt.*-result\" |")
    # 	  or die "Failed opening pipe to \"find\" ($!)";
    # 	while (<$findpipe>) {
    # 	  chomp;
    # 	  system("cat $_ >> ${main::wd}all-scan-results.txt")
    # 	    and die "Failed to concatenate sequence homology scan results ($!)";
    # 	}
    # 	my $lines = `wc -l ${main::wd}all-scan-results.txt`;
    # 	$lines =~ s/\s.*//;
    # 	chomp $lines;
    # 	print "len: $len\tlines: $lines\n";
    # 	if ($lines == $len) {
    # 	  last;
    # 	}
    # 	close ($findpipe);
    #       }
    #     }
    system ("find ${main::wd} -name \"kp_cluster.lt.*-result\" -exec cat {} \\; >${main::wd}all-scan-results.txt");
    my $cmd;

    $cmd = "cat ${main::wd}all-scan-results.txt | build_adjlist.pl > ${main::wd}adj_list.bin";
    system($cmd) and die ("Command failed: $cmd");

    # Note: dfs_cluster really does need to be run in a specific directory
    # in order to work properly.  Bogus.
    $cmd = "cd ${main::wd}; dfs_cluster --database=${main::wd}kp_cluster < ${main::wd}adj_list.bin > ${main::wd}clusters.txt";
    system($cmd) and die ("Command failed: $cmd");

    #$cmd = "cat ${main::wd}all-scan-results.txt | kp_chimera.pl ${main::wd}articulations.txt > ${main::wd}chimeric.txt";
    #system($cmd) and die ("Command failed: $cmd");

    #$cmd = "cat ${main::wd}adj_list.bin | dfs_cluster --database=${main::wd}kp_cluster --chimera=${main::wd}chimeric.txt > ${main::wd}clusters.txt";
    #system($cmd) and die ("Command failed: $cmd");
  }


  {
    local $main::wd = canonicalize_directory_name("${main::wd}tmp/cap3");
    print STDERR "wd>> $main::wd\n";
    system("mkdir -p $main::wd");
    system("cd $main::wd; make_clusters.pl ../../basecall/qc-passed.seq ../../basecall/qc-passed.qual < ../precluster/clusters.txt");
  }

  showtime("Preclustering Time", time() - $stop_watch);
}

sub make_contigs {
  my $stop_watch = time();
  {
    {
      local $main::wd = canonicalize_directory_name("${main::wd}tmp/cap3");
      print STDERR "wd>> $main::wd\n";

      open (my $findpipe, "find $main::wd -name \"cluster-*.seq\" |")
	or die ("Failed opening pipe to \"find\" ($!)");

      my @jobs = ();

      while (<$findpipe>) {
	chomp;
	my $file = $_;
	my $command = "(/home/kreuter/bin/cap3 $file $CAP3_OPTIONS &> $file.cap.out )";
	my $jobid = enqueue_job ($command);
	chomp $jobid;
	push @jobs, $jobid;
      }
      close ($findpipe);
      unless (@jobs) {
	die ("No jobs enqueued.");
      }
      # All jobs have been enqueued.  Spin until they complete.
      spin_a_while (\@jobs);
    }
    print STDERR ">>$main::wd\n";
    # Move CAP3 .ace file output to ACE directory
    {
      print "find $main::wd -type f -name \"cluster-*.cap.ace\" |\n";
      open (my $findpipe, "find $main::wd -type f -name \"cluster-*.cap.ace\" |")
	or die ("Failed opening pipe to \"find\" ($!)");

      local $main::wd = canonicalize_directory_name("${main::wd}ACE");
      print STDERR "wd>> $main::wd\n";
      system("mkdir -p ${main::wd}");

      while (<$findpipe>) {
	chomp;
	my $bn = $_;
	$bn =~ s|.*/||;
	unless ($bn) { # There had better be a basename.
	  die ("Find output $_ is not a file name with a non-directory component.");
	}
	#print STDERR "renaming $_ to ${main::wd}$bn\n";
	rename ($_, "${main::wd}$bn");
      }
      close ($findpipe);
    }


    # Find all the contigs and move them to unigene/ directory
    # Find all the singlets and move them to unigene/ directory
    # Move over preclustering singlets
    system("unigene_fasta.pl ${main::wd}");
  }
  system("mkdir -p ${main::wd}unigene");
  system("mv unigene.seq unigene.qual ${outputdir}unigene");

  showtime("CAP3 Contig Assembly time", time() - $stop_watch);
}


sub showtime {
    my ($label, $seconds) = @_;
    my ($hours, $minutes);

    $hours = $seconds / 3600;
    $seconds %= 3600;
    $minutes = $seconds / 60;
    $seconds %= 60;

    printf STDERR "$label:\t %02d:%02d:%02d\n",$hours,$minutes,$seconds;
}

sub usage {
    print STDERR <<EOF;
   
    
    Script for producing unigene set straight from source chromatograms
    
    Inputs/Options:
	-o <output directory tree> (required)

    Output:
	If execution is complete and successful, <output directory>
        will contain:

        ACE/  (directory containing CAP3 ACE file output)
	discard/ (rejected sequences and quality scores)
	    - qc_failed.seq 
		 quality control failed sequence
            - chimeric.seq
                 chimeric sequences
            NB: discard/ directory has not actually been implemented
	unigene/ (unigene sequence files and quality scores)
	    - contigs.seq
		contig consensus sequences
	    - unigene.seq
		composite file of singlets and contigs
	basecall/ (Outputs of basecalling and trimming process)
	    - raw sequence
	    - vector masked
	    - trimmed sequence
	tmp/ (Working directory of this script - may be removed)


EOF
exit 1;
}
