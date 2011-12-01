
=head1 NAME

detect_overlap_gff3.pl - detect overlapping features in two gff3 files

=head1 SYNOPSIS

detect_overlap_gff3.pl reference.gff mapped.gff

Output:

two tab separated columns to standard out:

 1) feature name
 2) feature description

Some stats are printed to STDERR at the end.

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut

use strict;
use warnings;

use File::Slurp;
use Getopt::Std;

our %args;

getopts('t:', \%args);


my $genome_gff = shift;
my $marker_gff = shift;

my @genome_gff = read_file($genome_gff);

my %features =();

my $count =0;
foreach my $g (@genome_gff) { 
    chomp($g);
    my ($ref, $source, $type, $start, $end, $dot, $dir, $frame, $desc) = split /\t/, $g;

    if ($type eq $args{t}) { 
	$desc =~ s/.*ID=[.*?:]*(.*?)\.\d+;.*/$1/;
	
    
	$features{$desc}->{start} = $start;
	$features{$desc}->{end} = $end;
	$features{$desc}->{ref} = $ref;
	$count++;
    }
    
}
print STDERR "Parsed $count genes.\n";
my @marker_gff = read_file($marker_gff);
my $total_covered =0;
my $feature_count = 0;

foreach my $m (@marker_gff) { 
    chomp($m);
    my ($ref, $source, $type, $start, $end, $dot, $dir, $frame, $desc) = split /\t/, $m;
    
    $desc =~ s/.*ID=(.*?);.*/$1/;

    foreach my $k (keys(%features)) { 
	if ($ref eq $features{$k}->{ref}) { # on same chr?
	    my $overlap = detect_overlap($start, $end, $features{$k}->{start}, $features{$k}->{end});
	    
	    if ($overlap) { 
		print STDERR "detected overlap: $k $desc               \r";
		print "$k\t$desc\n";
		$total_covered += $features{$k}->{end} - $features{$k}->{start};
		$feature_count++;
	    }
	}
    }
}
print STDERR "Done. A total of $feature_count features had overlap; the features totalled $total_covered basepairs.\n";

sub detect_overlap { 
    my $start = shift;
    my $end = shift;
    my $gstart = shift;
    my $gend = shift;

    my $overlap = 0;
    if ($start > $gstart && $start < $gend) { 
	$overlap =1;
    }
    if ($end < $gend && $end > $gstart) { 
	$overlap =1;
    }
    if ($end > $gend && $start < $gstart) { 
	$overlap =1;
    }
    return $overlap;
}

