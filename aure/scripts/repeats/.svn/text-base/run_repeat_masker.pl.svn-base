#!/usr/bin/perl

use strict;
use Getopt::Std;

use vars qw($opt_r);

getopts('r:');

my @files = @ARGV;

my $repeatmasker_bin = "/usr/local/bin/RepeatMasker";

my $repeat_lib_file = "/usr/local/blast_datasets/repeats.master.20051217";
if ($opt_r) { 
    $repeat_lib_file = $opt_r;
}
print STDERR "Using repeat library $repeat_lib_file...\n";

foreach my $f (@files) { 

    print STDERR "Now analyzing $f...\n";
    my $repeatmasker_outfile = $f.".repeats";
    my $repeatmasker_errfile = $f.".errors";

    system("$repeatmasker_bin -no_is -q -gff -arabidopsis -nolow -dir .  -lib $repeat_lib_file $f > $repeatmasker_outfile 2> repeatmasker_errfile");


}

print STDERR "Done.\n";
