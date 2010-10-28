#!/usr/bin/env perl

# adjacency.pl: process the output of Koni's custom blast program
# into a bizarre binary format, for feeding into some badly 
# documented program that processes files in this format.

# The input file is a whitespace-delimited table, whose first, second, 
# and fifth columns represent two vertices of a graph and the distance
# between the vertices.  Vertices are identified by integers represented 
# as longs. We need to construct a binary file whose format is this:

# 0 The number of vertices,
# 1 For each vertex (in order, from zero up), the number of neighbors
# 2 For each vertex that has neighbors, a sequence consisting of
#   the pairs whose first element is id of the neighbor and whose 
#   second element is the distance; and the sequence must be sorted
#   by neighbor ids.

# We'll build up these three structures incrementally, in three files:

# 0count.bin
# 1index.bin
# 2adjacencies.bin

# The previous version of this program read the entire input file into 
# memory.  This didn't scale.  Now, we pre-process the input in a 
# two-stage pipeline.  Stage 1 extracts the two vertices v1 and v2, and
# the distance d, writes two records to an output file (v1, v2, d), 
# (v2, v1, d).  Stage 2 sorts the output file by the vertex id in the left
# column.  This represents *all* permutations of edges in the graph.

# Next, we loop through the all-permutations file, and for each record we
# write out the record to a temporary file, so long as we see the same first
# column (i.e., the temporary file will contain all edges for a single vertex).
# Once we've seen all the edges for the current vertex, we sort the
# temporary file by neighbor ids, and then write out the (neighbor, distance) pairs
# at the end of 2adjacencies.bin.  When we're done with the temp file, we write the 
# number of records we processed to the end of 1index.bin.

# Finally, we write the number of edges we processed to 0count.bin.

#system qq[awk '{ printf("%s\t%s\t%s\n", $1, $2, $5); printf("%s\t%s\t%s\n", $2, $1, $5); }' | sort -n -T. > all-scan-results-permuted.txt];

use strict;

my $permuted_file = "all-scan-results-permuted.txt";
my $count_file = "0count.bin";
my $idx_file = "1index.bin";
my $adj_file = "2adjacencies.bin";

# First, permute the columns into a tempfile.
# Note: the following pipeline ran about twice as fast as doing
# this in Perl, when I tried it.
# awk '{ printf("%s\t%s\t%s\n", $1, $2, $5); printf("%s\t%s\t%s\n", $2, $1, $5); }' | sort -n -T. > all-scan-results-permuted.txt
open (my $permutefh, "|sort -n --stable -T. > all-scan-results-permuted.txt") or die ("$0: $!");
while (<STDIN>) {
  chomp;
  split;
  my $vertex = $_[0];
  my $neighbor = $_[1];
  my $distance = $_[4];
  printf $permutefh "%s\t%s\t%s\n", $vertex, $neighbor, $distance;
  printf $permutefh "%s\t%s\t%s\n", $neighbor, $vertex, $distance;
}
close ($permutefh);

# Now we're going to look through the permuted file, each record of 
# which looks like this:

# <vertex>\t<neighbor>\t<distance>

# We can only write out the data for a given vertex /after/ we've 
# seen all the records for that vertex.  So as we go, we hold onto 
# the last vertex id we saw, and compare it to the vertex of the 
# current record.  As soon as they differ, we (1) close the file 
# handle to which we've been writing vertex data, (2) run the
# routine that appends the adjacency list and the neighbor count to
# the respective files, (3) open a new file for storing the current
# vertex's neighbors, and (4) write the current record to that new 
# file.
my $last_id = -1; # A sentinel value for the beginning of the loop.
my $last_id_fname; # the file name to which we dribble all input; 
                   # it changes as we go along
my $last_id_fh;    # the file handle to which we dribble all input;
                   # it changes as we go along

open (my $idxfh, ">$idx_file") or die ("$0: $!");
open (my $adjfh, ">$adj_file") or die ("$0: $!");
open (my $infh, "<$permuted_file") or die ("$0: $!");
while (<$infh>) {
  my $line = $_;
  chomp;
  split;
  my $curr_id = $_[0]; # the current vertex id
  if ($last_id != $curr_id) { # we're looking at a new id,
    if ($last_id > -1) { # we've seen an id before
                         # (i.e., we're not at start of the file)
      close ($last_id_fh); 
      flush_last_id ($last_id_fname);
    }
    my $diff = $curr_id - $last_id - 1;
    if ($diff > 0) { # this means there were some vertices between the last
                     # seen vertex and the current one with zero neighbors.
      for (my $i = 0; $i < $diff; $i++) {
	print $idxfh pack("l", 0);
      }
    }
    $last_id = $curr_id;
    $last_id_fname = "$last_id.tmp";
    open ($last_id_fh, ">$last_id_fname") or die ("$0: $!");
  }
  print $last_id_fh $line;
}
# At the end of this loop, we've still got 1 vertex to flush,
# the last one in the file.
flush_last_id ($last_id_fname);

# Write the count file.
open (my $cntfh, ">$count_file") or die ("$0: $!");
my $count = $last_id+1; # This calculation is suspect: if the ids run from 0 upto N
                        # then there were N+1 vertices; maybe I should count these
                        # in the loop?
print $cntfh pack ("l", $count);
close ($cntfh);
system ("cat $count_file $idx_file $adj_file");
print STDERR "$count ids processed.\n";

sub flush_last_id {
  my $last_id_fname = shift;
#  print STDERR "Writing out data for $last_id... ";
  my $count = 0;
  my $last_neighbor = -1;
  my $last_conf = undef;
  open (my $sortfh, "sort --stable -n -k 2 $last_id_fname |") or die ("$0: $!");
  while (<$sortfh>) {
    chomp;
    split;
    my $neighbor = $_[1];
    my $conf = $_[2];
    if (($last_neighbor > -1) && ($last_neighbor != $neighbor)) {
      print $adjfh pack("l", $last_neighbor), pack("l", $last_conf);
      $count++;
    }
    $last_neighbor = $neighbor;
    $last_conf = $conf;
  }
  if ($last_neighbor) {
    print $adjfh pack("l", $last_neighbor), pack("l", $last_conf);
    $count++;
  }
  close ($sortfh);
  print $idxfh pack("l", $count);
  unlink ($last_id_fname);
#  print STDERR "($count rows).\n";
}
