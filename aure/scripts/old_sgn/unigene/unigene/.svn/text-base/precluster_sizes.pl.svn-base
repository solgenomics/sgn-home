
use strict;

# given a directory, reports the number of sequences in each cluster*.seq file and sorts the output by size.

my $dir = shift;

print STDERR "Collecting the data from directory $dir.\n";
print STDERR "Finding the seq files...\n";
my @files = `find $dir -iname "cluster*.seq"`;

my @counts = ();

print STDERR "Getting the number of sequences in each file... (1 '.'=100 files, ".scalar(@files)." total)...\n";
my $file_count = 0;
foreach my $f (@files) { 
    chomp($f);
#    if ($file_count % 100 == 0) { print STDERR "."; }
    printf STDERR "%5.1f", ($file_count / scalar(@files) * 100);
    print STDERR "\% \r";
    my $count = `grep '>' $f | wc -l`;
    chomp($count);
#    print STDERR "Count: $count\n";

    if ($count > 100000) { 
	print STDERR "\n$count > 100000 in file $f. Limiting to 1000000.\n"; 
	$count = 100000; 
    }
    push @counts, $count;
    $file_count++;
}

sort @counts;

foreach my $c (@counts) { 
    print $c."\n";
}




print "Sequence Size Distribution:\n\n";

my @cumulative = ();

foreach my $c (@counts) { 
    $cumulative[$c]++;
}

my $i=0;
my @compressed =();
my %label = ();
my $divs = 20;
my $increment = @cumulative/$divs;

for (my $i=0; $i<@cumulative; $i++) {
    $compressed[$i/$increment]+=$cumulative[$i];
}
my $max = 0;
foreach my $i (@compressed) {
    if ($max < $i) { $max= $i; }
}
my $maxwidth = 60;
my $factor=$maxwidth/$max;

#  print "max=$max. factor=$factor\n";
for (my $i=0; $i<@compressed; $i++) {
    my $size=$i*$increment;
    printf "%6d",  ($size); print " ["; printf "%6d", "$compressed[$i]"; print "]";
    print " ";
    for (my $n=0; $n<($compressed[$i] * $factor); $n++) {
	print "*";
    }
    print "\n";
}


