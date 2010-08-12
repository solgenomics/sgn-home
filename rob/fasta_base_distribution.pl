#!/usr/bin/env perl
use strict;
use warnings;

use Data::Dumper;

my %counts;
while(<>) {
    next if />/;
    chomp;
    $counts{$_}++ for split "",$_;
}

$counts{total} += $_ for values %counts;

print Dumper \%counts;
printf(q|%0.2f%% Ns|."\n", $counts{N}/$counts{total}*100 );
