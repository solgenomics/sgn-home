#!/usr/bin/perl

=head1 NAME

tpf_calculate.pl - calculates BAC alignments used for the integration into the tpf files.

=head1 SYNOPSIS

tpf_calculate.pl -C config_file -f <gff|tab> -c

 -C provides the config file
 -f provides output data
 -c calculate the alignments. If omitted, just format the present files.

Requires the mummer package.

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut

use strict;
use warnings;

use Config::JFDI;
use File::Slurp;
use File::Temp qw/tempdir/;
use Getopt::Std;

our %args;

getopts('i:f:cC:q:s:', \%args);

# i: input file
# f: output format (gff3)
# c: compute alignments (query file -q, subject file -s)
# s: subject file for -c
# q: query file for -c


if ($args{c} && $args{q} && $args{s}) { 
    calculate_alignments($args{s}, $args{q});
}

if (!$args{C}) { die "Need a conf file (specify with -C)\n"; }

my $conf = Config::JFDI->new(file=>$args{C})->get();

if ($args{c}) { 
    print STDERR "Using conf file $args{C}...\n";

    foreach my $a (@{$conf->{query}}) { 
	foreach my $b (@{$conf->{subject}}) {
	    calculate_alignments($b, $a);
	}
    }
}

if (!$args{f}) { exit; }

foreach my $f (@{$conf->{subject}}) { 
    format_file($args{f}, $f.".tiling");
}

sub format_file { 
    my $format = shift;
    my $file = shift;

    my @match_coords = read_file($file); # .tiling file from nucmer/show-tiling
#my @agp = read_file($args{a}); #current agp
#my @tpf = read_file($args{t}); # current tpf file

# find blocks of aligned sequences in coords
    my @list = ();
    my @group = ();
    my $ref;
    my $previous_end_coord = 0;
    for (my $i=0; $i<@match_coords; $i++) { 
	
	chomp($match_coords[$i]);
	if ($match_coords[$i]=~/^\>(.*?) .*/) { 
	    $ref = $1;
	    #print "CURRENT REF = $ref ($match_coords[$i])\n";
	    next();
	}
	my ($ref_start, $ref_end, $overlap, $length, $coverage, $identity, $orientation, $id) = split /\t/, $match_coords[$i];
	
	print STDERR "Read $ref_start, $ref_end, $id, $overlap\n";
	
	if ($ref_start < $previous_end_coord)  {  # overlap
	    push @group,  [ $ref, $id, $ref_start, $ref_end, $orientation ];
	}
	else { 
	    print STDERR "saving old group...\n";
	    push @list, [ @group ];
	    print STDERR "creating new group...\n";
	    @group = ( );
	    push @group, [$ref, $id, $ref_start, $ref_end, $orientation ];
	}
	
	$previous_end_coord = $ref_end;
    }
# deal with last group
    push @list, [ @group ];
    
    
    print STDERR "List: \n";
    foreach my $l (@list) { 
	if ($args{f} =~ /tab/) { output_tab($l); next; }
	if ($args{f} =~ /gff/) { output_gff($l); }
	if ($args{f} =~ /tpf/) { output_tpf($l); }	
    }
    print STDERR "Done!\n";
}

sub output_tab { 
    my $l = shift;
    print "GAP!\n";
    print "group coords:\t$l->[0]->[2]\t$l->[-1]->[3]\n";
    foreach my $g (@$l) { 
	print join "\t", @$g;
	print "\n";
    }
}

sub output_gff { 
    my $l = shift;
    
    print join "\t", ($l->[0]->[0], ".", "match", $l->[0]->[2], $l->[-1]->[3], ".", $l->[0]->[4], ".", "ID=$l->[0]->[1];Parent=$l->[0]->[0]" );
    print "\n";

    
}
    

sub calculate_alignments { 
    my $ref_fasta = shift;
    my $qry_fasta = shift;

    #my $dir = tempdir( CLEANUP=>0, DIR=>);
    #print STDERR "Note: Using temp dir $dir\n";

    #chdir($dir);

    `nucmer -g 500 -b 2000 --prefix=$ref_fasta $ref_fasta $qry_fasta`;

    `show-coords -rcl $ref_fasta.delta > $ref_fasta.coords`;

   ### `show-aligns ref_qry.delta refname qryname > ref_qry.aligns`;

    `show-tiling $ref_fasta.delta > $ref_fasta.tiling`;
}

sub parse_conf { 
    my $file = shift;

    my @conf = read_file($file);
    chomp(@conf);
    my $cref;
    foreach my $c (@conf) { 
	my ($command, $subject, $query) = split /.*\s.*/, $c;
	push @{$cref->{analyses}}, {command =>  $command,
				    subject =>  $subject,
				    query   =>  $query,
	};
    }
    return $cref;
}
