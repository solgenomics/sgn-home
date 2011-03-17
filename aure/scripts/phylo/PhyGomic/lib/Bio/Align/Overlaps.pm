
package Bio::Align::Overlaps;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck carp |;
use Try::Tiny;
use Math::BigFloat;

use Bio::Matrix::Generic;
use Bio::SimpleAlign;

###############
### PERLDOC ###
###############

=head1 NAME

Bio::Align::Overlaps.pm
A class to manipulate alignment overlaps

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  


=head1 DESCRIPTION

 Object to manipulate Bio::SimpleAlign according with the overlapping regions

=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 




=head2 calculate_overlaps

  Usage: my $mtx = calculate_overlaps($align);

  Desc: Calculate the overlaps for each pair of sequences in an alignment.
        It returns an Bio::Matrix::Generic with sequences ids as col/row
        names and a hashref. with start,end,length,identity keys.

  Ret: a Bio::Matrix::Generic object.

  Args: $align, an alignment object (Bio::SimpleAlign)
        
  Side_Effects: Die if the argument used is not a Bio::SimpleAlign object.

  Example: my $mtx = calculate_overlaps($align, $id);

=cut

sub calculate_overlaps {
    my $align = shift ||
	croak("ERROR: No align object was supplied to calculate_overlaps()");

    if (ref($align) ne 'Bio::SimpleAlign') {
	croak("ERROR: $align supplied to calculate_overlaps isnt SimpleAlign");
    }

    ## Get the id if no id arg. was supplied
    
    my $id = shift || $align->id();
    
    ## Define the mtx object

    my $mtx;

    ## Get members, member_ids and coordinates

    my @members = $align->each_seq();
    my %coords = ();
    foreach my $member (@members) {

	my $seqid = $member->display_id();
	my ($st, $en) = get_coordinates($member);

	$coords{$seqid} = { st => $st, en => $en, seq => $member };
    }
    my @member_ids = sort( keys %coords);

    ## Calculate overlaps, only if number of members are > 1

    if (scalar(@members) > 1) {

	## Define the default values:

	my $default_vals = { start => 0, end => 0, length => 0, identity => 0 };

	## Create a matrix object with the default values

	$mtx = Bio::Matrix::Generic->new( -rownames          => \@member_ids,
					  -colnames          => \@member_ids, 
					  -matrix_id         => $id,
					  -matrix_init_value => $default_vals,
	    );
	
	
	## It will compare coordinatdes between two pairs
	## and fill only the down part of the simetric matrix ($x > $y)

	my $x = 0;

	foreach my $id_a (@member_ids) {
	    $x++;
	    my $y = 0;

	    foreach my $id_b (@member_ids) {
		$y++;
	    
		## Skip selfoverlap (x = y) and the uppart of the
		## matrix (x < y)

		if ($x < $y) {

		    ## Define coords

		    my $a_st = $coords{$id_a}->{st};
		    my $a_en = $coords{$id_a}->{en};
		    my $b_st = $coords{$id_b}->{st};
		    my $b_en = $coords{$id_b}->{en};
		    
		    my ($st_ovl, $en_ovl) = ($b_st, $a_en);  ## Defaults

		    ## Compare coordinates and define cases

		    my @case = ('A', 'A', 'A', 'A');

		    if ($a_st >= $b_st) {
			$case[0] = 'B';
			$st_ovl = $a_st;
		    }
		    if ($a_en >= $b_en) {
			$case[1] = 'B';
			$en_ovl = $b_en;
		    }
		    if ($a_st >= $b_en) {
			$case[2] = 'B';
		    }
		    if ($a_en >= $b_st) {
			$case[3] = 'B';
		    }
		    
		    ## Assign cases to overlaps:
		    ##      +-------+   BBAB      +-------+      AAAB   
		    ##   +--------+                  +-------+
		    ##
		    ##      +--+        BAAB      +-------+      ABAB
		    ##   +--------+                  +--+
		    ##
		    ##        +---+     BBBB     +--+            AAAA
		    ##   +--+                         +---+

		    my $case = join('', @case);

		    if ($case =~ m/A/ && $case =~ m/B/) {  ## There is overlap
			
			## Start and end were calculated, and length is
			## $en_ovl - $st_ovl. Only identity remains.
			## To calculate it will create a new alignment with
			## this two sequences and it will slide with the
			## ovl region.

			my @pr = ($coords{$id_a}->{seq}, $coords{$id_b}->{seq});
			my $praln = Bio::SimpleAlign->new(-seqs => \@pr);
			my $ident = $praln->slice($st_ovl, $en_ovl)
			                  ->percentage_identity();

			my $vals = { start    => $st_ovl, 
				     end      => $en_ovl, 
				     length   => $en_ovl - $st_ovl + 1, 
				     identity => $ident };

			$mtx->entry($id_a, $id_b, $vals);
			$mtx->entry($id_b, $id_a, $vals);
		    }
		    else {

			$mtx->entry($id_a, $id_b, $default_vals);
			$mtx->entry($id_b, $id_a, $default_vals);
		    }
		}
	    }
	}
    }

    return $mtx
}


=head2 get_coordinates

  Usage: my ($start, $end) = get_coordinates($seq);

  Desc: Calculate where a sequence starts and ends in an alignment. It reads
        the sequence and start to count the sequence when it doesnt find
        a gap sign ('-' , '*' or '.').

  Ret: $start, the start coordinates
       $end, the end corrdinates

  Args: $seq, a sequence object.
        
  Side_Effects: Die if no argument is supplied.
                Die if the argument supplied is not a Bio::Seq object.

  Example: my ($start, $end) = get_coordinates($seq);

=cut

sub get_coordinates {
    my $seq = shift ||
	croak("ERROR: No seq object was supplied to get_coordinates.");

    if (ref($seq) !~ m/^Bio::(Seq::Meta|LocatableSeq)/) {
	croak("ERROR: Object $seq supplied to get_coordinates isnt Bio::Seq");
    }

    ## Declare coordinates (position, start and end)

    my ($po, $st, $en) = (0, 0, 0);

    ## Split the sequence into nt and scan it

    my $seqstr = $seq->seq();
    my @nts = split(//, $seqstr);

    foreach my $nt (@nts) {
	$po++;
	if ($nt !~ m/(-|\*|\.)/) {  ## If it is not a gap
	    if ($st == 0) {
		$st = $po;
	    }
	    $en = 0;
	}
	else {
	    if ($st > 0 && $en == 0) {
		$en = $po - 1;
	    }
	}
    }
    if ($en == 0) {
	$en = $po;
    }
    
    ## Return coordinates
    return ($st, $en);
}




####
1; #
####
