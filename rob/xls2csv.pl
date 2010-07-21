#!/usr/bin/env perl
use strict;
use warnings;

use Getopt::Std;
use Pod::Usage;

use Spreadsheet::ParseExcel;

my %opt;
getopts( 'is:', \%opt ) or pod2usage;
$opt{s} && $opt{s} =~ /^\d+$/ or pod2usage('must provide an integer -s arg');

my @files = @ARGV or pod2usage;
my $excelp = Spreadsheet::ParseExcel->new;

foreach my $file (@files) {
    my $book = $excelp->Parse($file);

    #     print "\n";
    #     print "Original Filename :", $book->{File} , "\n";
    #     print "Number of Sheets  :", $book->{SheetCount} , "\n";
    #     print "Author            :", $book->{Author} , "\n";
    #     print "\n";

    use Data::Dumper;
    $Data::Dumper::Maxdepth = 2;
    my $sheet = $book->{Worksheet}->[ $opt{s} - 1 ];
    defined $sheet or die "no sheet $opt{s} found in file $file\n";
    my $sheetname = $sheet->{Name};

    my $out_fh = do {
        if ( $opt{i} ) {
            my $f = $file;
            $f =~ s/\.xls$/.csv/i;
            open my $o, '>', $f or die "$! writing $f";
            $o;
        }
        else {
            \*STDOUT;
        }
    };

    my $cumulativeBlankLines = 0;

    my $minrow = $sheet->{MinRow};
    my $maxrow = $sheet->{MaxRow};
    my $mincol = $sheet->{MinCol};
    my $maxcol = $sheet->{MaxCol};

    #print "Minrow=$minrow Maxrow=$maxrow Mincol=$mincol Maxcol=$maxcol\n";

    for ( my $row = $minrow ; $row <= $maxrow ; $row++ ) {
        my $outputLine = "";

        for ( my $col = $mincol ; $col <= $maxcol ; $col++ ) {
            my $cell = $sheet->{Cells}[$row][$col];
            if ( defined($cell) ) {
                $_ = $cell->Value;    #{Val};

                # convert '#NUM!' strings to missing (empty) values
                s/#NUM!//;

                # escape double-quote characters in the data since
                # they are used as field delimiters
                s/\"/\\\"/g;
            }
            else {
                $_ = '';
            }

            $outputLine .= "\"" . $_ . "\"" if ( length($_) > 0 );

            # separate cells with commas
            $outputLine .= "," if ( $col != $maxcol );

        }

        #$outputLine =~ s/[, ]+$//g;  ## strip off trailing blanks and commas

        # skip blank/empty lines
        if ( $outputLine =~ /^[, ]*$/ ) {
            $cumulativeBlankLines++;
        }
        else {
            $out_fh->print("$outputLine\n");
        }
    }

    #print "  (Ignored $cumulativeBlankLines blank lines.)\n"
    #if ($cumulativeBlankLines);
    #print "\n";
}

=head1 NAME

xls2csv.pl

=head1 SYNOPSIS

xls2csv.pl [options] -s <sheetnum>  file file file ...

=head1 DESCRIPTION

Translate the Microsoft Excel spreadsheet file contained in <excel
file> into comma separated value format.  By default, dumps all CSV to
stdout.

Options:

  -i  in-place translation, make a new .csv file next to each .xls
      file containing CSV for the given sheet number
=cut
