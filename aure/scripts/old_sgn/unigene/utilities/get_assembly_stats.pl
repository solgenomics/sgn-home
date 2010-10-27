#!/usr/bin/perl

use strict;

my $assembly_dir = shift;

my @files = `ls $assembly_dir/`;

my %contigs = ();
my %contig_bac_count = ();
my %libs_bac_count =();
my %contig_length = ();
my %contig_total_bac_length = ();
my $total_bacs =0;
my $error_bacs =0;
my %singleton_bacs = ();
my $total_singleton_len = 0;
my $total_bac_len = 0;

print STDERR "Reading files...\n";

foreach my $f (@files) { 
    chomp ($f);

    if ($f=~/ace$/) { 
	print STDERR "Processing $f...\r";
	open (F, "<$assembly_dir/$f" ) || die "Can't open $f...\n";
	my $contig= "";
	while (<F>) { 
	    chomp;
	    if (/^CO (.*?)\s+(.*?)\s+/) { 
		$contig=$1;
		$contig_length{$contig}=$2;
	    }
	    if (/^RD (.*?)\s+(.*)\s+/) { 
		my $bac_name = $1;
		my $seq_len = $2;
		push @{$contigs{$contig}}, $bac_name;
		$contig_bac_count{$contig}++;
#		print STDERR "$bac_name Length=$seq_len\n";
		$contig_total_bac_length{$contig}+=$seq_len;
		$total_bacs++;
		$total_bac_len+=$seq_len;
	    }
	    if (/^BS (.*?)\s+(.*)?\s+/) { 
		$singleton_bacs{$1}++;
		$total_singleton_len += $2;
	    }
	}
    }
}

print STDERR "Processing results...\n";
# count the distinct library configurations
#
my %distinct_libs = ();
my %distinct_libs_contig_len = ();
my %distinct_libs_bac_len = ();
foreach my $k (keys %contigs) { 
    my %libs = ();
    foreach my $bac (@{$contigs{$k}}) { 
	my $lib = "";
	if ($bac=~/(Hba|Mbo|Eco)/i) { 
	    $lib = $1;
	    $libs{$lib}++;
	    $libs_bac_count{$lib}++;
	}
	else { 
	    $error_bacs++;
	    print STDERR "ERROR: $bac\n";
	}
    }
    my $distinct_libs_in_contig = join "-", sort(keys(%libs));
    $distinct_libs{$distinct_libs_in_contig}++;
    $distinct_libs_contig_len{$distinct_libs_in_contig} +=$contig_length{$k};
    $distinct_libs_bac_len{$distinct_libs_in_contig} +=$contig_total_bac_length{$k};
}

# print out results
#
foreach my $l (keys %libs_bac_count) { 
    print "Total for $l: $libs_bac_count{$l}\n";
}
foreach my $l (keys %distinct_libs) { 
    print "$l\tcount:\t$distinct_libs{$l}\tseq len \t$distinct_libs_contig_len{$l}\ttotal bac len:\t$distinct_libs_bac_len{$l}\n";
}

my $singleton_count = 0;
foreach my $l (keys %singleton_bacs) { 
    $singleton_count += $singleton_bacs{$l};
}


print "Total BACs in contigs processed: $total_bacs. Total length: $total_bac_len Singletons: $singleton_count. Lengt=$total_singleton_len. BACs with errors: $error_bacs\n";




