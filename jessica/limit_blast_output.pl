
use strict;
use warnings;

my $blast_results = shift;
my $correct_length = shift;

###set as needed:
my $percent_allowed = 99.99;
my $missmatch_allowed = 5;
my $length_wiggle = 5000;
###

my %check_length;
open(my $current_value,"<",$correct_length) || die "Can't open $correct_length.";
while(<$current_value>){
    chomp;
    if($current_value=~/^\#/||!$current_value){next;} #skip comment lines
    my ($bac_id, $bac_length) = split/\t/;
    $check_length{$bac_id}=$bac_length;
}



open(my $line, "<", $blast_results) || die "can't open $blast_results.";
while(<$line>){
    if($line=~/^\#/||!$line){next;} #skip comment lines
    chomp;
    my ($q_id, $s_id, $per_id, $len, $mismatch, $gaps, $q_start, $q_end, $s_start, $s_end, $e_val, $bit) = split /\t/; #declare columns
    #check for desigered lenght match
    if(($len <= $check_length{$q_id}) && ($len >= ($check_length{$q_id} - $length_wiggle))){
	#check for EXACT same length, AND check for at least desigered accuracy, AND check maximum allowed mismatch number
	if(($per_id >= $percent_allowed) && ($mismatch <= $mismatch_allowed)){
	    print join("\t",($q_id, $s_id, $per_id, $len, $mismatch, $gaps, $q_start, $q_end, $s_start, $s_end, $e_val, $bit));
	    print "\n";
	}
    }
}
	    
