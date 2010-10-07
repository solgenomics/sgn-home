#!/usr/bin/perl
use strict;
use warnings;

use FindBin;


sub usage {

  die <<EOU;
$FindBin::Script <unigene assembly base directory>

Calls acefile_membership.pl, parses its output and reports the sizes
of each of the clusters.
EOU
}

@ARGV == 1 or usage;
-d $ARGV[0] or usage;

`which acefile_membership.pl`
  or die "acefile_membership.pl not found in path";

open my $memberpipe, "acefile_membership.pl $ARGV[0] |"
  or die "Could not run acefile_membership.pl: $!";
#open my $memberpipe, "$ARGV[0]";

while(<$memberpipe>) {
  if(/^Contig:|^Singlet:/) {
    my(undef,undef,undef,$size) = split;
    print "$size\n";
  }
}

