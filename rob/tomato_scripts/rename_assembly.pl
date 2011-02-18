#!/usr/bin/env perl
use strict;
use warnings;
use 5.10.0;

use Data::Dumper;
use File::Spec;
use File::Basename;

use List::MoreUtils qw/ pairwise /;

use autodie ':all';

my %expand = ( ch => 'chromosomes',
               sc => 'scaffolds',
               ct => 'contigs',
              );

my $dest_dir = pop @ARGV;
my @files    = @ARGV;
my @targets = @files;

for( @targets ) {
    my ($bn,$dir) = fileparse($_);

    for ($bn) {
      s/^SL/S_lycopersicum_/;
      s/\.fasta/.fa/;
      while( my ($k, $v) = each %expand ) {
          s/(?<=[^a-z])$k(?=[^a-z])/$v/;
      }
      s/^(\D+)(\d\.\d+)([^\.]+)(\..+)$/$1$3.$2$4/;
      s/\.(fasta|agp)$/.$1.gz/;
    }

    $_ = $bn;
}

my %renames = pairwise {
    no warnings 'once';
    $a => $b
} @files, @targets;


#die Dumper \%renames;

mkdir $dest_dir;

while( my ($src, $dest) = each %renames ) {
    $dest = File::Spec->catfile( $dest_dir, $dest );
    if( $dest =~ /\.gz$/ && $src !~ /\.gz$/ ) {
        system "gzip -c '$src' > '$dest'";
    }
    else {
        system "cp '$src' '$dest'";
    }
}
