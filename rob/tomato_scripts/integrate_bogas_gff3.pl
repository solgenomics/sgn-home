#!/usr/bin/env perl
use strict;
use warnings FATAL => 'all';
use autodie ':all';

use feature 'state';

use Carp;
use File::Temp;

use Hash::Merge 'merge';
Hash::Merge::set_behavior('RETAINMENT_PRECEDENT');

use Path::Class;
use URI::Escape;

#use Smart::Comments;

die 'args are bogas_gff3, itag_gff3' unless @ARGV == 2;
my ( $bogas_gff3, $itag_gff3 ) = map file($_), @ARGV;

# algorithm:
#   loop over itag gff3:
#     harvest the attributes from the gene and mRNA features to be superseded
#     filter the superseded features out of the ITAG gff
#   loop over bogas features:
#     reformat the bogas features
#     add any missing attributes to the gene and mRNA bogas features
#     print out the bogas features at the end of the itag gff
#   re-sort the itag gff

my $temp_out = File::Temp->new;

my %itag_features;
my $itag_fh = $itag_gff3->openr;

LINE:
while( my $line = <$itag_fh> ) {

    unless( $line =~ /^#/ ) {
        my @f = split /\t/, $line, 9;
        chomp $f[-1];
        my ( $seq, $source, $type, $start, $end, $score, $strand, $phase, $attr ) = @f;

        if( $source =~ /_eugene$/ ) {

            if( $type eq 'mRNA' || $type eq 'gene' ) {
                $f[8] = $attr = parse_gff3_attrs( $attr );
                # remember these attributes, index them by ID
                $attr->{ID} or die "no ID attr on line $.: $line";
                $itag_features{ $attr->{ID}->[0] } = \@f;
                ### remembering attributes: $attr
            }

            # don't print it
            next LINE;
        }
    }
    $temp_out->print( $line );
}

my %increment_version; #< hashref of features we need to increment the version on
my $bogas_fh = $bogas_gff3->openr;
while( my $line = <$bogas_fh> ) {

    if( $line =~ /^#/ ) {
        $temp_out->print( $line );
    } else {
        my @f = split /\t/, $line, 9;
        my ( $ref, $source, $type, $start, $end, $score, $strand, $phase, $attr ) = 0..8;
        my $attrs = parse_gff3_attrs( pop @f);

        if( $f[$type] eq 'TE' ) {
            $f[$type] = 'transposable_element';
        }

        #### ad-hoc transformations

        $attrs->{from_BOGAS} = 1;
        $f[$source] = 'ITAG_eugene';
	# correct length attr, these are wrong in bogas dumps
	if( $attrs->{length} ) {
	    $attrs->{length} = $f[$end]-$f[$start]+1;
	}
	if( $attrs->{Nb_exon} ) {
	    $attrs->{nb_exon} = delete $attrs->{Nb_exon};
	}

        if(    $f[$type] eq 'gene'              ) {
            $attrs->{ID}[0] =~ s/\.\d+$//;
            $attrs->{Name}[0] = $attrs->{ID}[0];
            ( my $alias = $attrs->{Name}[0] ) =~ s/\.\d+$//;
            $attrs->{Alias} = [ $alias ];
        }
        elsif( $f[$type] eq 'mRNA'              ) {
            $attrs->{Parent} ||= [ $attrs->{ID}[0] ];
            $attrs->{Parent}[0] =~ s/\.\d+$//;
            $attrs->{Parent}[0] =~ s/^mRNA:/gene:/;
        }
        elsif( $f[$type] eq 'CDS'               ) {
        }
        elsif( $f[$type] eq 'sequence_assembly' ) {
            # filter out sequence_assembly features
            next;
        }
        elsif( $f[$type] =~ /^(exon|intron)$/   ) {
            # make the exons and introns belong to the mRNA
            $attrs->{Parent}[0] =~ s/^gene:/mRNA:/;
        }
	elsif( $f[$type] =~ /_UTR$/ ) {
	    #uniqify UTR ids
	    state %utr_ctr;
	    $attrs->{ID}[0] =~ s/(\d+)$//;
	    $attrs->{ID}[0] .= ++$utr_ctr{ $attrs->{ID}[0] }-1
	}

        # merge the attributes with the ones from the existing ITAG stuff

        if( $attrs->{ID} and my $itag_f = $itag_features{ $attrs->{ID}[0] } ) {
            my $itag_attrs = $itag_f->[8];
            ### merging: $attrs
            ### with: $itag_attrs
            delete $itag_attrs->{ID};
            delete $itag_attrs->{Parent};
            delete $itag_attrs->{Name} if $f[$type] eq 'gene';
            $attrs = merge( $attrs, $itag_attrs );

            #if the coordinates of the old and new are different,
            #note that on the next pass, we need to increment the version number of this gene model
            unless( $f[$start] == $itag_f->[$start] && $f[$end] == $itag_f->[$end] ) {
                my ( $solyc ) = $attrs->{ID}[0] =~ /(Solyc\d\d[rg]\d+)/ or die;
                $increment_version{ $solyc } = 1;
                ### feature changed from: $itag_f
                ### to: @f
            }

            ### result: $attrs
        }

        # remove mRNA: prefixes from Names and Aliases
        for my $aname ( 'Name','Alias' ) {
            $attrs->{$aname}[0] =~ s/^[^:]+:// if $attrs->{$aname};
        }

        $temp_out->print( format_gff3_line( \@f, $attrs ) );
    }
}

$temp_out->close;


my $temp_out_2 = File::Temp->new;
my $temp_out_read = file( "$temp_out")->openr;
while( my $line = <$temp_out_read> ) {
    if( $line =~ /^#/ ) {
        $temp_out_2->print( $line );
    } elsif( $line =~ /\t(gene|mRNA)\t/ ) {
        $temp_out_2->print( increment_version_for_line( $line ) );
    }
}
$temp_out_read->seek(0,0);
while( my $line = <$temp_out_read> ) {
    unless( $line =~ /^#/ ||  $line =~ /\t(gene|mRNA)\t/ ) {
        $temp_out_2->print( increment_version_for_line( $line ) );
    }
}
$temp_out_2->close;
system "gff3_reformat.pl -s $temp_out_2";

########## SUBROUTINES ##############

sub increment_version_for_line {
    my $line = shift;
    my @f = split /\t/, $line, 9;
    my ( $ref, $source, $type, $start, $end, $score, $strand, $phase, $attr ) = 0..8;
    my $attrs = parse_gff3_attrs( pop @f);

    for my $aname ('ID','Parent','Name') {
	if( $attrs->{$aname} && $attrs->{$aname}[0] ) {
	    my ( $solyc ) = $attrs->{$aname}[0] =~ /(Solyc\d\d[rg]\d+)/ or die;
	    if( $increment_version{ $solyc } ) {
		$attrs->{$aname}[0] =~ s/\.(\d+)/'.'.($1+1)/e;
	    }
	}
    }
    return format_gff3_line( \@f, $attrs );
}


sub format_gff3_line {
    my ( $f, $attr ) = @_;
    $attr ||= {};

    my $attr_string = join ';' => (
        map {
            my $key = $_;
            my $val = $attr->{$key};
            $val = [ $val ] unless ref $val;
            "$key=".join( ',', map gff3_escape($_), @$val );
        } sort keys %$attr
       );

    return join( "\t", @$f, $attr_string )."\n";
}

# take the attr column, returns a hashref of the parsed attrs
sub parse_gff3_attrs {
    my ( $attr_string ) = @_;
    chomp $attr_string;

    my %attrs;
    for my $a ( split /;/, $attr_string ) {
        next unless $a;
        my ( $name, $values ) = split /=/, $a, 2;
        defined $name   && length $name   or die;
        defined $values && length $values or die;

        push @{$attrs{$name}}, map uri_unescape($_), split /,/, $values;
    }

    return \%attrs;
}

sub gff3_escape {
    uri_escape( $_[0], '\n\r\t;=%&,\x00-\x1f\x7f-\xff' );
}

