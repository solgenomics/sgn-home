
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

=head2 seed_list

  Usage: my @seeds_href = seed_list($mtx, $method, $filter_href);

  Desc: Calculate the seed list based in length, identity, or ovlscore

  Ret: An array with arrayrefs. as elements. Each arrayref will have two 
       elements (one per seed id)

  Args: $mtx, a Bio::Matrix::Generic object with the overlaps
        $method, method to order the seeds. Permited values are 'length', 
        'identity' or 'ovlscore'.
        $filter_href, a hashref with key=length,identity,start,end and
        value=integer. It will throw away all the overlaps with values lower
        than these values.
        
  Side_Effects: Die if no argument is used.
                Die if first argument is not Bio::Matrix::Generic object.
                If method is undef, or it doesnt match with length or identity
                ovlscore will be used as default. 

  Example: my @seeds_href = seed_list($mtx);
           my @seeds_href = seed_list($mtx, undef, { length   => 100, 
                                                     identity => 80  });

=cut

sub seed_list {
    my $mtx = shift ||
	croak("ERROR: No argument was supplied to seed_list.");

    if (ref($mtx) ne 'Bio::Matrix::Generic') {
	croak("ERROR: $mtx supplied to seed_list isnt a Bio::Matrix::Generic");
    }

    my $method = shift || 'ovlscore';

    if ($method !~ m/^(length|identity)$/) {
	$method = 'ovlscore';
    }

    ## Get filter
    
    my $filter_href = shift;
    my %filter;
    if (defined $filter_href) {
	if (ref($filter_href) ne 'HASH') {
	    croak("ERROR: Filter argument is not a hashref. for seed_list");
	}
	else {
	    %filter = %{$filter_href};
	}
    }

    ## Declare the hsh to store the pairs

    my %seedscoring = ();

    ## First, it will get the list ids from the matrix

    my @ids = $mtx->row_names();

    ## It will ignore the redundant entries.

    my $x = 0;
    
    foreach my $id_a (sort @ids) {
	
	$x++;
	my $y = 0;
	
	foreach my $id_b (sort @ids) {
	    
	    $y++;

	    if ($x < $y) {
		my $entry = $mtx->entry($id_a, $id_b);

		if ($entry->{length} > 0) {  ## ignore no-overlaps
		
		    my $pass = 1;
		    foreach my $param (keys %filter) {
			if (defined $entry->{$param}) {
			    if ($entry->{$param} < $filter{$param}) {
				$pass = 0;
			    }
			}
			else {
			    croak("ERROR: filter param. $param isnt permited.");
			}
		    }

		    if ($method eq 'length') {
			if ($pass == 1) {
			    $seedscoring{$id_a.':'.$id_b} = $entry->{length};
			}
		    }
		    elsif ($method eq 'identity') {
			if ($pass == 1) {
			    $seedscoring{$id_a.':'.$id_b} = $entry->{identity};
			}
		    }
		    else {
			my $idenfrac = $entry->{identity} / 100;
			my $ovlscore = $entry->{length} * $idenfrac * $idenfrac;
			if ($pass == 1) {
			    $seedscoring{$id_a . ':' . $id_b} = $ovlscore;
			}		
		    }
		}
	    }
	}
    }

    ## Finally it will order the pairs based in the score

    my @pairs = sort {$seedscoring{$b} <=> $seedscoring{$a}} keys %seedscoring;
    
    ## and to a new arrays as pairs
    
    my @seedlist = ();
    foreach my $pair (@pairs) {

	my @seed_ids = sort(split(/:/, $pair));
	push @seedlist, \@seed_ids;
    }
    
    return @seedlist;
}

=head2 calculate_overseeds

  Usage: my %overseeds = calculate_overseeds($align, $seed_aref, $start, $end);

  Desc: Calculate the overlaps for each sequence in an alignment with the
        seeds overlap region. 
        The identity is calculated as the average of the seeds and the new
        sequence

  Ret: a hash with key=sequence_id, value=hashref. 
         with keys=start,end,length and identity.

  Args: $align, an alignment object (Bio::SimpleAlign)
        $seed_aref, an array ref. with two seed ids.
        $start, start of the seed overlaping region
        $end, end of the seed overlapping region
        
  Side_Effects: Die if the argument used is not a Bio::SimpleAlign object.
                Die if one (or more) arguments are absent.

  Example: my %overseeds = calculate_overseeds($align, $seed_aref, 10, 100);

=cut

sub calculate_overseeds {
    my $align = shift ||
	croak("ERROR: No align object was supplied to calculate_overseeds()");

    if (ref($align) ne 'Bio::SimpleAlign') {
	croak("ERROR: $align supplied to calculate_overseeds isnt SimpleAlign");
    }

    my $seeds_aref = shift ||
	croak("ERROR: No seed aref. was supplied to calculate_overseeds()");

    my @seeds;
    if (ref($seeds_aref) ne 'ARRAY') {
	croak("ERROR: $seeds_aref supplied to calculate_overseeds isnt AREF.");
    }
    else {
	@seeds = sort @{$seeds_aref};
	if (scalar(@seeds) != 2) {
	    croak("ERROR: seed aref. supplied doesnt have 2 members");
	}
    }

    my $start = shift ||
	croak("ERROR: No start was supplied to calculate_overseeds()");

    if ($start !~ m/^\d+$/) {
	croak("ERROR: $start supplied to calculate_overseeds isnt an integer");
    }

    my $end = shift ||
	croak("ERROR: No end was supplied to calculate_overseeds()");

    if ($end !~ m/^\d+$/) {
	croak("ERROR: $end supplied to calculate_overseeds isnt an integer");
    }

    my $id = join(':', @seeds);
    
    ## now it will calculate the coordinates for each seq

    ## Define the overseed hash

    my %overseeds;

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
	
	## It will compare coordinates between the seed coordinate and 
	## the alignment member. It will skip the seed members

	foreach my $member_id (@member_ids) {
	    if ($member_id ne $seeds[0] && $member_id ne $seeds[1]) {
		
                ## Define coords
		
		my $m_st = $coords{$member_id}->{st};
		my $m_en = $coords{$member_id}->{en};
		    
		my ($st_ovl, $en_ovl) = ($m_st, $end);  ## Defaults

		## Compare coordinates and define cases
		
		my @case = ('A', 'A', 'A', 'A');
		
		if ($start >= $m_st) {
		    $case[0] = 'B';
		    $st_ovl = $start;
		}
		if ($end >= $m_en) {
		    $case[1] = 'B';
		    $en_ovl = $m_en;
		}
		if ($start >= $m_en) {
		    $case[2] = 'B';
		}
		if ($end >= $m_st) {
		    $case[3] = 'B';
		}
		
		my $case = join('', @case);

		if ($case =~ m/A/ && $case =~ m/B/) {

		    ## It will calculate the identity as an overall of
		    ## two seeds and the candidate

		    my $seed0 = $coords{$seeds[0]}->{seq};
		    my $seed1 = $coords{$seeds[1]}->{seq};
		    my $memb =  $coords{$member_id}->{seq};
		    my $seqlist = [$seed0,$seed1,$memb];

		    my $praln = Bio::SimpleAlign->new(-seqs => $seqlist);
		    my $ident = $praln->slice($st_ovl, $en_ovl)
			              ->percentage_identity();

		    my $vals = { start    => $st_ovl, 
				 end      => $en_ovl, 
				 length   => $en_ovl - $st_ovl + 1, 
				 identity => $ident };
		    
		    $overseeds{$member_id} = $vals; 
		}
	    }
	}	
    }
    return %overseeds;
}


=head2 extension_list

  Usage: my @extension = extension_list(\%overseeds, $method);

  Desc: Calculate the seed extension list based in length, identity, or ovlscore

  Ret: An array with sequence id elements.

  Args: $overseed_href, a hashref with member_id as keys and hashref with 
        overlap parameters as value.
        $method, method to order the seeds. Permited values are 'length', 
        'identity' or 'ovlscore'.
        
  Side_Effects: Die if no argument is used.
                Die if first argument is not an hashref.
                If method is undef, or it doesnt match with length or identity
                ovlscore will be used as default. 

  Example: my @extension = extension_list(\%overseeds, 'identity');

=cut

sub extension_list {
    my $ovlseed_href = shift ||
	croak("ERROR: No argument was supplied to extension_list.");

    if (ref($ovlseed_href) ne 'HASH') {
	croak("ERROR: $ovlseed_href supplied to extension_list isnt a HASHREF");
    }

    my $method = shift || 'ovlscore';

    if ($method !~ m/^(length|identity)$/) {
	$method = 'ovlscore';
    }

    ## Get filter
    
    my $filter_href = shift;
    my %filter;
    if (defined $filter_href) {
	if (ref($filter_href) ne 'HASH') {
	    croak("ERROR: Filter argumebt isnt a hashref. for extension_list");
	}
	else {
	    %filter = %{$filter_href};
	}
    }

    ## Declare the hash to store the list

    my %extscores = ();

    ## First, it will get the list ids from the matrix

    my @ids = keys %{$ovlseed_href};
    
    foreach my $id (sort @ids) {
	
	my $entry = $ovlseed_href->{$id};

	if ($entry->{length} > 0) {  ## ignore no-overlaps

	    my $pass = 1;            ## Check the filter options
	    foreach my $param (keys %filter) {
		if (defined $entry->{$param}) {
		    if ($entry->{$param} < $filter{$param}) {
				$pass = 0;
		    }
		}
		else {
		    croak("ERROR: filter param. $param isnt permited.");
		}
	    }
		
	    if ($method eq 'length') {
		if ($pass == 1) {
		    $extscores{$id} = $entry->{length};
		}
	    }
	    elsif ($method eq 'identity') {
		if ($pass == 1) {
		    $extscores{$id} = $entry->{identity};
		}
	    }
	    else {
		my $idenfrac = $entry->{identity} / 100;
		my $ovlscore = $entry->{length} * $idenfrac * $idenfrac;
		if ($pass == 1) {
		    $extscores{$id} = $ovlscore;
		}		
	    }
	}
    }

    ## Finally it will order the pairs based in the score

    my @extension = sort {$extscores{$b} <=> $extscores{$a}} keys %extscores;

    return @extension;
}

=head2 global_overlap

  Usage: my %ovldata = global_overlap($align, \@selectedmembers);

  Desc: Calculate the total overlap region.

  Ret: %ovldata, a hash ref with the following keys:
        start, an start coordinate
        end, an end coordinate.
        length, length of the overlap (0 if there are no overlap)

  Args: $align, a Bio::SimpleAlign object
        $selectedmembers_aref, an array reference with the selected members 
        
  Side_Effects: Die if no argument is used.
                Die if first argument is not an Bio::SimpleAlign object
                Die if second argument is not an array ref. and if it has
                less of 2 memberss

  Example: my %ovldata = global_overlap($align, \@selectedmemb);

=cut

sub global_overlap {
    my $align = shift ||
	croak("ERROR: No align object was supplied to global_ovelap()");
    
    if (ref($align) ne 'Bio::SimpleAlign') {
	croak("ERROR: $align supplied to global_overlap isnt SimpleAlign");
    }

    my $selmemb_aref = shift;

    my %selected;

    if (defined $selmemb_aref) {
	if (ref($selmemb_aref) ne 'ARRAY') {
	    croak("ERROR: $selmemb_aref supplied to global_overlap isnt AREF");
	}
	else {
	    if (scalar(@{$selmemb_aref}) < 2) {
		croak("ERROR: Less than 2 mmb were supplied to global_overlap");
	    }
	    foreach my $selid (@{$selmemb_aref}) {
		$selected{$selid} = 1;
	    }
	}
    }
    else {  ## By default it will use all the sequences in the alignment

	foreach my $seq ($align->each_seq()) {
	    $selected{$seq->display_id()} = 1;
	}
    }

    ## First get the coordinates for each member.

    my $max_en = 0;  ## Get also the max end
    
    my @members = $align->each_seq();
    foreach my $member (@members) {

	my $seqid = $member->display_id();

	if (exists $selected{$seqid}) {
	    my ($st, $en) = get_coordinates($member);
	    $selected{$seqid} = { st => $st, en => $en, seq => $member };
	    
	    if ($en >= $max_en) {   ## Set the max. end to the bigger end
		$max_en = $en;
	    }
	}
    }

    ## Check the selected member that don't exists in the alignment

    foreach my $sel_id ( keys %selected) {
	if (ref($selected{$sel_id}) ne 'HASH') {
	    croak("ERROR: $sel_id doesnt exist into the align=$align\n");
	}
    }

    ## Second recalculate the overlaps

    my ($max_st, $min_en) = (0, $max_en);
    
    foreach my $seq_id (sort keys %selected) {

	if ($selected{$seq_id}->{st} >= $max_st) {
	    $max_st = $selected{$seq_id}->{st};
	}
	
	if ($selected{$seq_id}->{en} <= $min_en) {
	    $min_en = $selected{$seq_id}->{en};
	}
    }
    
       my $length = $min_en - $max_st + 1;

    my %ovldata = ( start => $max_st, end => $min_en, length => $length ); 
 
    if ($length < 1) {
	 %ovldata =  ( start => 0, end => 0, length => 0 );
    }
   
    return %ovldata;
}




####
1; #
####
