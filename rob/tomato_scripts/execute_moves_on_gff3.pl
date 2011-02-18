#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use autodie ':all';

use File::Temp;

#use Smart::Comments;

use List::Util qw/ sum /;

my ( $gff3_file, $movesfile ) = @ARGV;
@ARGV == 2 or die "must give gff3_file, moves_file\n";

my %moves = read_moves( $movesfile );

use File::ReadBackwards;
my $gff3_in = File::ReadBackwards->new( $gff3_file );
my $temp = File::Temp->new;
while( my $line = $gff3_in->readline ) {
    my $delete_this_line = 0;
#open my $gff3_in, '<', $gff3_file;
#while( my $line = <$gff3_in> ) {
    if( $line =~ /^#/ ) {
        if( my ( $directive, $body ) =  $line =~ /^##(\S+)\s+(.+)/ ) {
            chomp $body;
            if( $directive eq 'sequence-region' ) {
                my ( $name, $start, $end ) = split /\s+/, $body;
                ### $name
                $end = calculate_new_length( $end, $moves{ normalize_name( $name ) } );
                $line = "##sequence-region $name $start $end\n";
            }
        }
	$temp->print( $line );
    } else {
        my ( $seq, $source, $type, $start, $end, $score, $strand, $phase, $attr ) =
            my @f =
                split /\t/, $line, 9;
        chomp for $attr, $f[-1];

        my $moves = $moves{ normalize_name( $seq ) };

        # apply the moves to this annotation
        if( $moves ) {

            for my $move ( @$moves ) {
                my ( $loc, $offset ) = @$move;

                # positive offsets are insertions of N, just move
                # coords forward if the annotation affects them
                # (i.e. if it is at a lower cooridnate)
                if( $offset > 0 ) {
		    if( $loc <= $end ) {
			$end += $offset;
			if( $loc <= $start ) {
			    $start += $offset;
			} else {
			    $temp->print( "# MOVES: in line above, $offset bp insertion into middle of feature at $loc.\n" );
			    invalidate_feature( $line );
			}
		    }
                }
                # negative offsets are deletions.
                elsif( $offset < 0 ) {
                    my $del_start = $loc + $offset;
                    my $del_end   = $loc - 1;
                    # 0 1 2 3 4 5 6    ( loc 5, offset -2 )
                    #       s   e
                    if( $del_end <= $start ) {
                        # deletion happens entirely before this
                        # feature, just shift the two coordinates
                        $_ += $offset for $start, $end;
                    }
                    elsif( $del_start > $end ) {
                        # deletion is entirely after this feature, do
                        # nothing
                    }
                    # now know:
                    #  $del_end   >  $start
                    #  $del_start <= $end
                    else {

                        my @coords = sort {
			    $a->[0] <=> $b->[0]
                        } ( [$del_start,'start del'], [$del_end,'end del'], [$start,'start feat'], [$end,'end feat'] );

                        # this is deleting part of this feature, need
                        # to delete the whole feature group
                        $temp->print( "# MOVES: in line above, $offset bp deletion (".join(', ', map "@$_", @coords).")\n" );
			invalidate_feature( $line );
                    }
                }
            }
        }

	scan_feature( $line );
	if( feature_is_valid( $line ) ) {
	    @f[3,4] = ( $start, $end );
	    $line = join( "\t", @f )."\n";
	    $temp->print( $line );
	} else {
	    $temp->print( "# MOVES: deleted: $line" );
	}
    }
}

$temp->close;
my $tfh = File::ReadBackwards->new( "$temp" ) or die;
while( my $line = $tfh->readline ) {
    unless( $line =~ /^#/ ) {
	$line = "# MOVES: deleted: $line" unless feature_is_valid( $line );
    }
    print $line;
}

exit;

############################3
my %features_with_children;
my %invalid_features;
sub invalidate_feature {
    my ( $line ) = @_;
    my $id = _get_attr( $line, 'ID' );
    if( $id && $features_with_children{$id} ) {
	return; # ignore if we have children
    }
    if( my $parent = _get_attr( $line, 'Parent' ) ) {
	invalidate_feature( $parent );
    }
    return $invalid_features{ $line } = 1;
}
sub feature_is_valid {
    my ( $line ) = @_;
    my $id = _get_attr( $line, 'ID' );
    my $parent = _get_attr( $line, 'Parent' );
    if( $invalid_features{$line} || $id && $invalid_features{$id} || $parent && $invalid_features{$parent} ) {
	invalidate_feature( $parent ) if $parent;
	return 0;
    }
    return 1;
}
sub scan_feature {
    my ( $line ) = @_;
    if( my $parent = _get_attr( $line, 'Parent' ) ) {
	$features_with_children{ $parent } = 1;
    }
}

sub _get_attr {
    my ( $line, $attrname ) = @_;
    $line =~ /$attrname=([^\s;]+)/
	or return;
    return $1;
}

sub read_moves {
    my ( $f ) = @_;
    open( my $fh, $f );
    my %m;
    while( <$fh> ) {
        chomp;
        my @f = split;
        my $move = shift @f;
        $move eq 'move' or die;
        my $id = normalize_name( shift @f );
        push @{$m{$id}}, \@f;
    }

    ### moves: %m

#     # sort the moves in *ascending* order by coordinate
#     for ( values %m ) {
#         $_ = [ sort { $a <=> $b } @$_ ];
#     }

    return %m;
}


sub normalize_name {
    ( my $n = shift ) =~ s/^SL2\.\d\d//;
    return $n;
}

sub calculate_new_length {
    my ( $current_length, $moves ) = @_;

    my $end_offset = sum(
          map $_->[1],
          grep { $_->[0] <= $current_length }
          @$moves
       ) || 0;
    ### $end_offset

    return $current_length + $end_offset;
}
