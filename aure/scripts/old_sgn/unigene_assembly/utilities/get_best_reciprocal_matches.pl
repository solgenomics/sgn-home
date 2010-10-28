#!/usr/bin/perl 

=head1 NAME

get_best_reciprocal_matches.pl - a script to identify reciprocal best matches in a blast report in m8 format.

=head1 DESCRIPTION

Takes a m8 formatted blast report on STDIN and outputs on STDOUT.

=head1 AUTHOR

Lukas Mueller <lam87@conrell.edu>


=cut 

use strict;

my %q_best_match = ();
my %q_best_evalue = ();
my %s_best_match = ();
my %s_best_evalue = ();

while (<>) { 
    chomp;
    my ($qid, $sid, $evalue) = (split /\t/)[0,1,10];
    
    print STDERR "$qid, $sid, $evalue\n";

    if (! exists($q_best_match{$qid})) { 
	$q_best_match{$qid} = $sid;
	$q_best_evalue{$qid} = $evalue;
    }
    if (! exists($s_best_match{$sid})) { 
	$s_best_match{$sid} = $qid;
	$s_best_evalue{$sid}=$evalue;
    }


    if ($evalue < $q_best_evalue{$qid}) { 
	$q_best_match{$qid}= $sid;
	print STDERR "Best match for $qid is $sid\n";
    }

    if ($evalue < $s_best_evalue{$sid}) { 
	$s_best_match{$sid} = $qid;
	$s_best_evalue{$sid} = $evalue;
	print STDERR "Best match for $sid is $qid\n";
    }

}

foreach my $k (keys(%q_best_match)) { 
    
    print STDERR "checkin $k  against $q_best_match{$k}...\($s_best_match{$q_best_match{$k}}\)\n";
    if ( ($k eq $s_best_match{$q_best_match{$k}})) { 
	print "$k reciprocal best match with $q_best_match{$k}\n";
    }
}
print STDERR "Done.\n";
	
    
