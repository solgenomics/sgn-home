#!/usr/bin/perl

=head1

load_tab_delimited_file.pl

=head1 SYNOPSIS

    load_tab_delimited_file.pl -i [input tab-delimited file] -t

=head1 COMMAND-LINE OPTIONS

 -i  tab-delimited data file
 -n  no header in input file (default: expect header)
 -h  Help

=cut

use strict;
use Getopt::Std;

my $file_name;
my @example_header;
my @loaded_data;
my $use_header = 1;
my @data;
my $data_row;
my $cell;
our ($opt_i, $opt_n, $opt_h);

getopts('hni:');
if ($opt_h) {
  help();
  exit;
}
if ($opt_n) {
  $use_header = 0;
}
if (!$opt_i) {
  print STDERR "\nFilename required\n\n\n";
  help();
}

$file_name = $opt_i;
@example_header = qw(field_a field_b field_c);
if ($use_header == 1) {
  @data = load_file($file_name,@example_header);
}
else {
  @data = load_file($file_name);
}
foreach $data_row (@data) {
  my @row = @{$data_row};
  foreach $cell (@row) {
    print "$cell\t";
  }
  print "\n";
}

sub load_file {
  my $name_of_file_to_load = shift;
  my @expected_header = shift;
  my @header;
  my @loaded_data;
  my %bad_header_fields;
  my $array_counter = 0;
  my $current_line = 1;
  open FILE, "<", $name_of_file_to_load or die $!;
  @loaded_data = map { chomp; [split /\t/];} <FILE>;
  if ($expected_header[0]) {
    @header = @{shift @loaded_data};
    if (@header < @expected_header) {
      print STDERR @header;
      die "Header has too few fields\n";
    }
    while ($array_counter < @expected_header) {
      if ($header[$array_counter] ne $expected_header[$array_counter]) {
	$bad_header_fields{$expected_header[$array_counter]} = $header[$array_counter];
      }
      $array_counter++;
    }
    if (%bad_header_fields) {
      foreach my $bad_field (keys %bad_header_fields) {
	print STDERR "Header field ".$bad_header_fields{$bad_field}." is incorrect.  Expected ".$bad_field."\n";
      }
      die;
    }
  }
  #remove trailing blank lines
  while (@{$loaded_data[-1]} == 0) {
      pop (@loaded_data);
  }
  #check that each row has the same number of cells
  foreach my $data_line (@loaded_data) {
    if (@{$data_line} != @{$loaded_data[0]}) {
      die "Line $current_line has the wrong number of cells.  Expected ".@{$loaded_data[0]}." and found ".@{$data_line}."\n";
      $current_line++;
    }
  }
  return @loaded_data;
}

sub help {
  print STDERR <<EOF;
  $0:

    Description:

     This is a template script for loading data from a tab-delimited file.

    Usage:
     
      load_tab_delimited_file.pl -i <tab-delimited file name>

    Flags:

      -i <tab-delimited file name>  input data file (mandatory)
      -n                            no header included in input file
      -h <help>                     help

EOF
  exit (1);
}

=pod

  =back

  =head1 LICENSE

  Same as Perl.

  =head1 AUTHORS

  Jeremy D. Edwards <jde22@cornell.edu>

  =cut


