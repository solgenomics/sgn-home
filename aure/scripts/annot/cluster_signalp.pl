#!/usr/bin/perl

use strict;
$| = 1;
use CXGN::Tools::Cluster::SignalP;

my $file_in = shift;
my $file_out = shift;

my $clusterer = CXGN::Tools::Cluster::SignalP->new
			({
					in => $file_in,
					out => $file_out,
					host => "solanine",
			});

$clusterer->submit();
$clusterer->spin(10);
$clusterer->concat();
print STDERR "\nDone!\n";
