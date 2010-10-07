#!/usr/bin/perl -w
use strict;

my %hit_spans = ();
keys(%hit_spans) = 100000;

my ($c1, $c2, $span, $seq_id, $line);
while($line = <>) {
  my ($query, $sbjct, $ident, $length, $start, $end, $new_span, $n);
  next if ($line =~ m/^\#/);

  ($seq_id, $sbjct, $ident, $length, undef, undef, $start, $end) = 
      split/\s/,$line;
  next if ($seq_id eq $sbjct);
  next if ($length < 100);

  if ($start > $end) {
      my $t;
      # Reverse complement match
      $t = $end;
      $end = $start;
      $start = $t;
  }

  $new_span = 1;
  foreach $span ( @{$hit_spans{$seq_id}} ) {
      $c1 = $span->[1] - $start;
      $c2 = $end - $span->[0];
      if (($c1 < 0 && $c2 < 0) || ($c1 > 0 && $c2 > 0)) { 
	  $span->[0] = $start if ($start < $span->[0]);
	  $span->[1] = $end if ($end > $span->[1]);
	  $span->[2]++;
	  $new_span = 0;
	  last;
      }
  }
  if ($new_span == 1) {
      push @{$hit_spans{$seq_id}}, [ $start, $end, 1 ];
  }
}

foreach $seq_id ( keys %hit_spans ) {
    my ($span1, $span2, $span_merged, $i, $j);
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
		    if ($span2->[0] - $span1->[1] < 20) {
			$chimera = 1;
		    }
		} else {
		    if ($span1->[0] - $span2->[1] < 20) {
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
