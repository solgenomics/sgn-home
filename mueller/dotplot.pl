
=head1 NAME

dotplot.pl - a quick hack to generate a dotplot from blast m8 reports

=head1 DESCRIPTION

dotplot.pl -i blast_m8_file -N query_id -M subject_id [-e evalue] [-m minimum_match_len (default 100)] [-p miminum_percent_identity (default 90)] [-x image_width_in_pixels (default 400)] [-y image_height_in_pixels (default 400)] 

outputs a png image to STDOUT.

=head1 AUTHOR

Lukas Mueller <lam87@cornell.edu>

=cut

use strict;
use warnings;

use Getopt::Std;

use GD;

our %args;

getopts('i:x:y:e:N:M:p:m:', \%args);

my $x = $args{x} || 400;
my $y = $args{y} || 400;
my $evalue_cutoff = $args{e} || 1e-10;
my $percent_identity_cutoff = $args{p} || 90;
my $minimal_align_length = $args{m} || 30;


my $image = GD::Image->new($x, $y);

open(my $F, "<", $args{i}) || die "Can't open file $args{i}\n";

my $max_query = `grep $args{N} $args{i} | cut -f8 | sort -nr | head -1`;
chomp($max_query);
print STDERR "Max query = $max_query\n";

my $max_subject = `grep $args{M} $args{i} | cut -f10 | sort -nr | head -1`;
chomp($max_subject);
print STDERR "Max subject = $max_subject\n";

my $x_scale = $max_query / $x;
my $y_scale = $max_subject / $y;

print STDERR "x-scale: $x_scale. y-scale: $y_scale\n";
my $bg = $image->colorResolve(255, 255, 255);
my $color = $image->colorResolve(100, 100, 100);
my $red = $image->colorResolve(255, 0, 0);
my $blue = $image -> colorResolve(0, 0, 255);
my $brown = $image -> colorResolve(150, 150, 0);
my $green = $image-> colorResolve(0, 0, 150);

while (<$F>) { 
    chomp;
    my ($query, $subject, $identity, $length, $gaps, $what, $q_start, $q_end, $s_start, $s_end, $evalue, $score) = split /\t/;

#print STDERR "Processing $query - $subject ($evalue)...\n";
    
    if ( ($evalue <= $evalue_cutoff) && ($args{N} eq $query) && ($args{M} eq $subject) && ($percent_identity_cutoff > $identity) && ($minimal_align_length < $length) ) { 

	if ($identity > 70) { $color = $green; }
	if ($identity > 80) { $color = $blue; }
	if ($identity > 85) { $color = $red; }


	$image->line($q_start / $x_scale, $s_start / $y_scale, $q_end / $x_scale, $s_end / $y_scale, $color);
	
    }
}	

close($F);

print $image->png();	

    
