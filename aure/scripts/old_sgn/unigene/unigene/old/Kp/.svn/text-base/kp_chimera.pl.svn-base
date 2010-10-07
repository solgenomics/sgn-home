#!/usr/bin/perl -w
use strict;

if (defined($ARGV[0]) && $ARGV[0] eq "help") {
  print <<EOF;

  Quick script to determine whether articulation points are chimeric clones
  or not.

  Expects name of file containing integer codes of articulation point sequences
  and complete results of Kp clustering on STDIN. Outputs potentially chimeric
  sequences on STDOUT

EOF
  exit 0;
}

my %interests ;

my $articulation_file = shift;
open ARTFILE, "<$articulation_file"
  or die "Failed opening articulation point file \"$articulation_file\" ($!)";

while(<ARTFILE>) {
  chomp;
  $interests{$_} = 1;
}

close ARTFILE;

my %hit_spans;
while(<>) {
  my ($start, $stop, $item);
  chomp;
  next if (m/RC$/);
  split;

  next if ($_[0] == $_[1]);

  if (defined($interests{$_[0]})) {
    $start = $_[7]+8;
    $stop = $_[8]-8;
    $item = $_[0];
    add_span($start, $stop, $item);
  }

  if (defined($interests{$_[1]})) {
    $start = $_[9]+8;
    $stop = $_[10]-8;
    $item = $_[1];
    add_span($start, $stop, $item);
  }

}

my ($seq_id);
foreach $seq_id ( keys %hit_spans ) {
    my ($span1, $span2, $span_merged, $i, $j, $c1, $c2);
    for($i=0;$i<@{$hit_spans{$seq_id}};$i++) {

	$span_merged = 0;
	for($j=$i+1;$j<@{$hit_spans{$seq_id}};$j++) {
	    $span1 = $hit_spans{$seq_id}->[$i];
	    $span2 = $hit_spans{$seq_id}->[$j];
	    $c1 = $span1->[1]-$span2->[0];
	    $c2 = $span2->[1]-$span1->[0];
	    if (($c1 < 0 && $c2 < 0) || ($c1 > 0 && $c2 > 0)) { 
		$span1->[0] = $span2->[0] if ($span2->[0] < $span1->[0]);
		$span1->[1] = $span2->[1] if ($span2->[0] > $span1->[0]);
		$span1->[2] += $span2->[2];
		$span_merged = 1;
		last;
	    }
	}
	if ($span_merged) {
	    # If we merged two spans, remove the second one, and start over
	    # checking the first one for spans which overlap with it
	    splice(@{$hit_spans{$seq_id}},$j,1);
	    $i--;
	}
    }
}
		
foreach $seq_id ( keys %hit_spans ) {
    my ($chimera, $span1, $span2, $i, $j);
    if (@{$hit_spans{$seq_id}}>1) {
	$chimera = 0;
	for($i=0;$i<@{$hit_spans{$seq_id}};$i++) {
	    $span1 = $hit_spans{$seq_id}->[$i];
	    next if ($span1->[2] <= 1); 
	    for($j=0;$j<@{$hit_spans{$seq_id}};$j++) {
		next if $i == $j;
		$span2 = $hit_spans{$seq_id}->[$j];
		next if $span2->[2] <= 1;

		if ($span1->[0] < $span2->[0]) {
		    if ($span2->[0] - $span1->[1] < 36) {
			$chimera = 1;
		    }
		} else {
		    if ($span1->[0] - $span2->[1] < 36) {
			$chimera = 1;
		    }
		}
	    }
	}
	if ($chimera) {
	    print ">$seq_id\n";
	    for($i=0;$i<@{$hit_spans{$seq_id}};$i++) {
		print "SPAN $i: ",join("\t",@{$hit_spans{$seq_id}->[$i]}),"\n";
	    }
	}
    }
}


sub add_span {
  my ($start, $end, $item) = @_;
  my ($new_span, $c1, $c2, $span);

  $new_span = 1;
  foreach $span ( @{$hit_spans{$item}} ) {
      $c1 = $span->[1] - $start;
      $c2 = $end - $span->[0];
      if (($c1 < 0 && $c2 < 0) || ($c1 > 0 && $c2 > 0)) {
#	  print "Merging to existing SPAN for item $item: [ $start, $end ] ->";
#	  print " [ $span->[0], $span->[1] ] ($span->[2])\n";
#	  $span->[0] = $start if ($start < $span->[0]);
	  $span->[1] = $end if ($end > $span->[1]);
	  $span->[2]++;
	  $new_span = 0;
	  last;
      }
  }

  if ($new_span == 1) {
#    print "Adding new span for item $item: [ $start, $end ]\n";
      push @{$hit_spans{$item}}, [ $start, $end, 1 ];
  }
}
