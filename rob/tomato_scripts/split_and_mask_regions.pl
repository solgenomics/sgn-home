#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';

use autodie ':all';
use Storable 'dclone';
use Getopt::Std;

use Smart::Comments;

use Data::Dump 'dump';
use List::MoreUtils 'indexes';

use Bio::SeqIO;

use CXGN::BioTools::AGP 'agp_parse', 'agp_write';

my %opt;
getopts('c',\%opt);

@ARGV == 4 or die "invalid arguments\n";
my ( $mask_region, $agp_file, $obj_seq, $comp_seq ) = @ARGV;

my ( $seqname, $start, $end ) = $mask_region =~ /^ ([^:]+) : (\d+) \.\. (\d+) $/x
    or die "invalid mask region '$mask_region'";
### $seqname
### $start
### $end

### parse the agp
my $agp = agp_parse( $agp_file );

my @line_nums = indexes { $_->{ident} && $_->{ident} eq $seqname } @$agp;
my $line_num = shift @line_nums;
# cannot handle more than one AGP location
die "mutiple locations for $seqname in AGP" if @line_nums;


### found seq at line: $line_num
die "component $seqname not found in AGP" if $line_num == -1;
die "cannot split first AGP line" if $line_num == 0;
die "cannot split last AGP line"  if $line_num == $#{$agp};

# split the AGP contig there
my $affected_obj_name;
my $obj_start;
my $obj_end;
{
    my $here = $agp->[$line_num];
    die "can only handle interior regions" if $start <= 1 || $end >= $here->{cend};
    $affected_obj_name = $here->{objname};

    # subref transforms component coordinates into object coordinates
    my $comp2obj = do {
        if( $here->{orient} eq '+' ) {
            my $component_offset = $here->{ostart} - $here->{cstart} + 1;
            ### $component_offset
            sub { shift() + $component_offset - 1 }
        } else {
            my $component_offset = $here->{oend} - $here->{cstart} + 1;
            ### $component_offset
            sub { $component_offset - shift() + 1 }
        }
    };

    my $frag1 = dclone( $here );
    my $gap   = dclone( $here );
    $gap->{gap_type} = 'fragment';
    $gap->{linkage}  = 'yes';
    my $gap_length = $end - $start + 1;


    my $frag2 = dclone( $here );

    if( $here->{orient} eq '+' ) {

        $frag1->{oend} = $comp2obj->( $start-1 );
        $frag1->{cend} = $start-1;
        $frag1->{ident} .= 'a';
        $frag1->{length} = $frag1->{oend} - $frag1->{ostart} + 1;

        $obj_start     = $gap->{ostart} = $comp2obj->( $start );
        $obj_end       = $gap->{oend}   = $comp2obj->( $end );
        $gap->{cstart} = $start;
        $gap->{cend}   = $end;
        $gap->{type}   = 'U';
        $gap->{is_gap} = 1;
        $gap->{length} = $gap_length;

        $frag2->{ostart} = $comp2obj->( $end+1 );
        $frag2->{cstart} = 1;
        $frag2->{cend}   -= $gap_length + $frag1->{length};
        $frag2->{ident} .= 'b';

    } else {

        $frag2->{ostart} = $comp2obj->( $start-1 );
        $frag2->{cend} = $start-1;
        $frag2->{ident} .= 'a';
        $frag2->{length} = $frag2->{oend} - $frag2->{ostart} + 1;

        $obj_start     = $gap->{ostart} = $comp2obj->( $end );
        $obj_end       = $gap->{oend}   = $comp2obj->( $start );
        $gap->{cstart} = $start;
        $gap->{cend}   = $end;
        $gap->{type}   = 'U';
        $gap->{is_gap} = 1;
        $gap->{length} = $end - $start + 1;

        $frag1->{oend}  = $comp2obj->( $end+1 );
        $frag1->{cend}  -= $gap_length + $frag2->{length};
        $frag1->{cstart} = 1;
        $frag1->{ident} .= 'b';
    }

    splice( @$agp, $line_num, 1, $frag1, $gap, $frag2 );
}

### $affected_obj_name
### $obj_start
### $obj_end

# split the component seq
my $temp_comp_seqs;
unless( $opt{c} ) {
    $temp_comp_seqs = filter_fasta_to_temp(
        $comp_seq,
        $seqname,
        sub {
            my ($orig) = @_;
            my $s = $orig->seq;
            my $len = $end - $start + 1;
            ### id: $orig->id
            ### length: $orig->length
            ### mask start: $start
            ### mask len: $len

            my $frag1 = Bio::PrimarySeq->new(
                -id => $orig->id.'a',
                -seq => substr( $s, 0, $start-1 ),
               );

            my $frag2 = Bio::PrimarySeq->new(
                -id => $orig->id.'b',
                -seq => substr( $s, $end ),
               );

            return $frag1, $frag2;
        }
       );
}

my $temp_obj_seqs = filter_fasta_to_temp(
    $obj_seq,
    $affected_obj_name,
    sub {
        my ($obj) = @_;
        my $s = $obj->seq;
        my $len = $obj_end - $obj_start + 1;
        substr( $s, $obj_start - 1, $len, 'N'x$len );
        $obj->seq( $s );
        return $obj;
    },
   );

system 'mv', "$temp_obj_seqs", $obj_seq;
system 'mv', "$temp_comp_seqs", $comp_seq unless $opt{c};
agp_write( $agp, $agp_file );

exit;

###############################

sub filter_fasta_to_temp {
    my ( $file, $seqname, $sub ) = @_;

    my $temp_seq = File::Temp->new;

    my $in = Bio::SeqIO->new(
        -format => 'fasta',
        -file   => $file,
       );
    my $out = Bio::SeqIO->new(
        -format => 'fasta',
        -fh     => $temp_seq,
        -width  => 80,
       );
    my $found = 0;
    while ( my $s = $in->next_seq ) {
        if( $s->id eq $seqname ) {
            $out->write_seq($_) for $sub->($s);
            $found = 1;
        } else {
            $out->write_seq($s);
        }

    }
    $found or die "'$seqname' not found in '$file'";

    return $temp_seq;
}


sub split_agp_line {
    my ( $agp, $line_num, $start, $end ) = @_;


    return ( $affected_obj_name, $obj_start, $obj_end );
}
