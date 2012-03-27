
=head1 NAME

detect_overlap_gff3.pl - detect overlapping features in two gff3 files

=head1 SYNOPSIS

detect_overlap_gff3.pl -t [feature_type] reference.gff mapped.gff

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

getopts('t:v', \%args);

if (!defined($args{t})) { 
    die "Need -t for feature type.";
}

my $genome_gff = shift;
my $marker_gff = shift;

my @genome_gff = read_file($genome_gff);

my %features =();
my $feature_count =0;
my %features_by_ref;

print STDERR "Processing...\n";

foreach my $g (@genome_gff) { 
    chomp($g);
    my ($ref, $source, $type, $start, $end, $dot, $dir, $frame, $desc) = split /\t/, $g;

    if ($type eq $args{t}) { 
	$desc =~ s/.*ID=[.*?:]*(.*?)\.\d+;.*/$1/;
	$features{$desc}->{start} = $start;
	$features{$desc}->{end} = $end;
	$features{$desc}->{ref} = $ref;
	$feature_count++;
	$features_by_ref{$ref}++;
    }    
}

my @marker_gff = read_file($marker_gff);
my $total_covered =0;
my $total_feature_count =0;
my %total_unique_coverage;
my %total_unique_count;

foreach my $m (@marker_gff) { 
    chomp($m);
    my ($ref, $source, $type, $start, $end, $dot, $dir, $frame, $desc) = split /\t/, $m;

    if (defined($desc)) { 
	$desc =~ s/.*ID=(.*?);.*/$1/;
    }

    foreach my $k (keys(%features)) { 
	if ($ref eq $features{$k}->{ref}) { # on same ref (chr)?
	    my $overlap = detect_overlap($start, $end, $features{$k}->{start}, $features{$k}->{end});
	    
	    if ($overlap) { 
		#print STDERR "detected overlap: $k $desc               \r";
		if ($args{v}) { print "$k\t$desc\n"; }
		my $coverage = $features{$k}->{end} - $features{$k}->{start};
		$total_covered += $coverage;
		$total_feature_count++;
		$total_unique_coverage{$ref}->{$k} += $coverage;
		$total_unique_count{$ref}->{$k}++;
	    }
	}
    }
}

print STDERR qq { Parsed $feature_count features of type $args{t} in file $genome_gff.\n };

print STDERR "$total_feature_count features covered $total_covered bp\n";
#printf STDERR "10%s" REFERENCE, "10%s" "#FEATURES", "10%s" "BASEPAIRS\n";

foreach my $ref (sort(keys(%total_unique_count))) { 
    my $total_unique_features = scalar(keys(%{$total_unique_count{$ref}}));
    my $total_unique_coverage = 0;
    foreach my $k (keys(%{$total_unique_coverage{$ref}})) { 
	if (defined($total_unique_coverage{$ref}->{$k}) ) { 
	    $total_unique_coverage += $total_unique_coverage{$ref}->{$k};
	}
    }

    printf STDERR "%10s %10d %10d %10d", $ref, $features_by_ref{$ref}, $total_unique_features, $total_unique_coverage;
    print STDERR "\n";


}

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

