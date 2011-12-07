
use strict;

use Getopt::Std;
use File::Slurp;
use TPF;
use CXGN::BioTools::AGP qw | agp_parse |;

our(%args);
getopts('t:a:f:b:', \%args);

# -a : current agp file
# -t : current tpf file
# -f : file specifying the new records to be integrated into the tpf file (has to be generated with tpf_calculate, tab option)
# -b : file with two tab separeted columns: (1) local BAC IDs, (2) Genbank accession for BAC

print STDERR "Reading conversion file. ";
my @conversion = read_file($args{b});
my $gb;
foreach (@conversion) { 
    chomp;
    my ($local, $genbank) = split /\t/;
    $gb->{$local} = $genbank;
}
print STDERR "Done.\n";

print STDERR "Reading TPF file. ";
my $tpf = TPF->new(file=>$args{t});
$tpf->parse();

print STDERR scalar(@{$tpf->tpf})." lines read.\n";

print STDERR "Reading AGP file. ";
my $agp = agp_parse($args{a});
print STDERR scalar(@$agp)." lines read.\n";

print STDERR "Reading tab file. ";
my @tab = read_file($args{f});
my $tab = \@tab;
print STDERR "Done.\n";
open(ERR,">", "tpf_integrate.log") || die "can't open error log file";
my @group;
my $gstart;
my $gend;
my @overlaps;
foreach my $new (@tab) { 
    chomp($new);

    if ($new =~ /^group/i) { 
	@group = ();
	(undef, $gstart, $gend) = split /\t/, $new;
	#print  "Processing group ($gstart, $gend)\n";
	@overlaps = find_overlaps($agp, $gstart, $gend);
	#print  "Identified overlaps: ".(join "\t", @overlaps). "\n";
	
    }
    elsif ($new =~ /^GAP/i) { 
	if (@overlaps) { 
	    my $errflag = process_group($gb, $tpf, $overlaps[0], $overlaps[-1], \@group);
	    if ($errflag) { 
		print ERR "Failed to insert ".(join ", ",  @group)."\n";
	    }
	}
	@group = ();
	@overlaps = ();
	
	# ignore
    }   
    else { 
	my ($ref, $name, $start, $end) = split /\t/, $new;
	push @group, $name;
    }

}
#process last group
if (@overlaps) { 
    process_group($gb, $tpf, $overlaps[0], $overlaps[-1], \@group);
}
$tpf->render();
close(ERR);

sub process_group { 
    my $gb = shift;
    my $tpf = shift;
    my $start_id = shift;
    my $end_id = shift;
    my $group = shift;

    $start_id =~ s/(.*)\.\d+/$1/; # remove version
    $end_id =~ s/(.*)\.\d+/$1/;

    my @formatted;
    foreach my $g (@$group) { 
	push @formatted, [ $gb->{$g}, "?", $g ];
    }

    print STDERR "Removing the following ids: ".(join ", ", $tpf->interval_ids($start_id, $end_id))."\n";
    
    return $tpf->replace($start_id, $end_id, @formatted);

}

sub find_overlaps { 
    my $agp = shift;
    my $s = shift;
    my $e = shift;

    my @overlaps; 

    foreach my $a (@$agp) { 
	my $overlap = detect_overlap($s, $e, $a->{ostart}, $a->{oend});
	if ($overlap) { 
	    push @overlaps, $a->{ident};
	}
    }
    return @overlaps;
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
