
#!/usr/bin/perl
=head1 NAME

 seq_coordinates_processing.pl.
 A script to process files with sequences coordenates (combine them in a global file).(version.2.0.).

=cut

=head1 SYPNOSIS

 seq_coordinates_procesing.pl [-h] -i <seq_coordenates_files> -t <coordenates_type> -o <output_file> -c <output_file_criterion> [-t <change_tag_format>]
  
    
=head2 I<Flags:>

=over

=item -i 

B<seq_coordenates_files>        sequence coordenates files, using single quotes and separated by commas (mandatory)

=item -f

B<qscore_in_fasta_format>       coordenates format, first base position equal 0 or 1 and L (lenght) or E (end) (begin position 0 and L by default)

=item -o

B<output_file>                  name of the output file (SeqCoordComb by default)

=item -c

B<output_file_criterion>       criterion (type) for the output file (begin=0 and length by default)

=item -t

B<tag_change_file>             tab file with two columns that contains the equivalences to change the seqclean clean tags into another (example: change 'UniVec' for "trim_lib1"). First column => old annotation, second column => new annotation

=item -r

B<message_for_id_removed>      in the seqclean format output, the message that can appear in the -f6 column (cleantag) when the id appears in the first files but not in the rest( by default: removed).

=item -h

B<help>                         print the help  

=back

=cut

=head1 DESCRIPTION

This script parse differents sequence coordenates file and combine them into a file.
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

seq_coordinates_processing.pl


=cut

use strict;
use warnings;

use File::Basename;
use Getopt::Std;

our ($opt_i, $opt_f, $opt_o, $opt_c, $opt_t, $opt_r, $opt_h);
getopts("i:f:o:c:t:r:h");
if (!$opt_i && !$opt_f && !$opt_o && !$opt_c && !$opt_t && !$opt_r && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

my ($aref_files, $aref_format, $outputfile, $outputformat)=check_input();
my @files=@$aref_files;
my @formats=@$aref_format;
my $n=0;
my (%global_n_perc,%global_coord5p,%global_coord3p,%global_raw_l,%global_clean,%global_trim);
my (%cleantags_changecode,%check_id);
if ($opt_t) {
    my $cleantags_changecode_href=parse_new_cleantags_file($opt_t);
    %cleantags_changecode=%$cleantags_changecode_href;
}
my @data_array=(\%global_n_perc,\%global_coord5p,\%global_coord3p,\%global_raw_l,\%global_clean,\%global_trim,
                 \%cleantags_changecode, \%check_id);
my $current_data_aref=\@data_array;

foreach my $file (@files) {
    my $nr=$n+1;
    my $new_data_aref=process_file($nr, $file, $formats[$n], $current_data_aref);
    $current_data_aref=$new_data_aref;
    $n++;
}

open my $out, '>', $outputfile || die "I can not open the output file $outputfile\n";
my $outputseqcleanfile=$outputfile;
$outputseqcleanfile =~ s/\.tab//;
$outputseqcleanfile .= '.cln';
my $outseqclean;

if ($outputformat eq 'S' || 's') {
    open $outseqclean, '>', $outputseqcleanfile || die "I can not open the output file $outputfile\n";
}

my @ids= sort keys %{$$current_data_aref[7]};
%global_n_perc = %{$$current_data_aref[0]};
%global_coord5p  = %{$$current_data_aref[1]};
%global_coord3p  = %{$$current_data_aref[2]};
%global_raw_l = %{$$current_data_aref[3]};
%global_clean  = %{$$current_data_aref[4]};
%global_trim  = %{$$current_data_aref[5]};

print "Printing data into the output file\n";
foreach my $id (@ids) {
    my $begin=$global_coord5p{$id};
    my $end=$global_coord3p{$id};    
    my ($first_coord, $second_coord)=transform_coord_to_format($outputformat, $begin, $end);
    print $out "$id\t$first_coord\t$second_coord\n";
    if ($outputformat eq 'S' || 's') {
	my $n_perc=$global_n_perc{$id}||'';
	my $raw_length=$global_raw_l{$id}||'';
	my $clean_tag=$global_clean{$id}||'';
	my $trim_tag=$global_trim{$id}||'';
	print $outseqclean "$id\t$n_perc\t$first_coord\t$second_coord\t$raw_length\t$clean_tag\t$trim_tag\n";
    }
}
print "The output file in $outputformat format is $outputfile\n";
if ($outputformat eq 'S' || 's') {
    print "The output file in seqclean format (seven columns) is $outputseqcleanfile\n";
    if (uc $formats[0] ne 'S') {
	print "\n\t<> Important message about the .cln file. <>\n\t<> The first file ($files[0]) has format ($formats[0]).<>\n\t<> This format has not sequence length so the sequence length into the output file:\n\t\t$outputseqcleanfile\n\t   is wrong. <>\n\t<> It should be something like length + X were X = unknown length from the files with format different of \'S\' <>\n\n";
    }
}


=head2 help

  Usage: help()
  Desc: print help of this script
  Ret: none
  Args: none
  Side_Effects: exit of the script
  Example: if (!@ARGV) {
               help();
           }

=cut

sub help {
  print STDERR <<EOF;
  $0: 

    Description: 
        This script parse different sequence coordenates file and combine them (calculate the final coordeantes). These files must have three columns: -f1=id, -f2=begin_coordenate (there are two format, first base is position 0 (format=0 by default)and first base is position 1 (format=1)), -f3=sequence_length (format=L, by default) or end_coordenate (format=E used as format=E in the -f option). The output file has the combination of the coordenates. The sequences coordeantes files and the formats must be in order (example: -i 'first,second,third' -f '0/L,0/L,1/E')

    Note about formats:
        SGN database has format 0/L
        Lucy coordenates output has format 1/E, or 0/0 for removed sequences (high expected error)
        Seqclean has format 1/E  
    
    Usage: 
       seq_coordinates_processing.pl [-h] -i <seq_coordenates_files> -f <seq_coordenates_type> -o <output_file> -c <output_coordeantes_criterion>
    
    Flags:

       -i <input_coordeantes_files>      input sequence coordeantes files (three columns, and they must be detailed by asc order, -i 'first,second,third...') (mandatory, at least two)

       -f <seq_coordenates_format>        the formats of the coordeantes. There are two formats, for the begin 0 or 1 (if the first base position is 0 or 1, 0 by default), and for the end (length (L) if is detailed as new sequence length or end if is detailed as end position in the original sequence (E), length(L) by default). The way to write as argument is -f '1/E,0/L,0/L'. Another format is 'S' that means seqclean. These means that use columns 1,3 and 4 with begin=1 and end=end_position.

       -o <output_file>                  the output file name (SeqCoordComb by default). 

       -c <output_file_criterion>        the criterion that is choosen for the output file (0/L by default) 

       -t <clean_tag_change_file>        It is a file tab with two columns to change the format of clean tag in the seqclean format (-f1=seqclean_original_tag -f2=new_tag, example -f1='UniVec' -f2='trim_lib1').

       -r <message for id removed>       Message that appears in the -f6 column (clean tag) in the seqclean format file when the id only appears in the first file (by default: removed) 
       -h <help>                         print the help    

EOF
exit (1);
}

=head2 check_input

  Usage: check_input()
  Desc: check if exists some mandatory inputs. Also check the compatibility of them (for example -r and -R).
  Ret: none
  Args: none
  Side_Effects: die if there are something wrong
  Example: check_input();

=cut

sub check_input {
  if ($opt_h) {
	help();
  }
  
  if (!$opt_i) {	
	die "Requiered argument -i <sequence_coordenates_files> was not supplied.\n";
  }
  my @inputfiles=split(/,/,$opt_i);
  my $n_files=scalar(@inputfiles);
  if ($n_files < 2) {
      print "Warning!!! the argument -i have one file (anyway you can use to convert coordenates formats)\n\n";
  }
  my $file;
  foreach $file (@inputfiles) {
     unless (-f $file) {
	 die "The file $file do not exists or it is a directory\n";
     }
     my $n_column_f1=`cut -f1 $file | wc -l`;
     chomp($n_column_f1);
     my $n_column_f2=`cut -f2 $file | wc -l`;
     chomp($n_column_f2);
     my $n_column_f3=`cut -f3 $file | wc -l`;
     chomp($n_column_f3);
     if ($n_column_f1 != $n_column_f2 || $n_column_f1 != $n_column_f3) {
	 die "The file $file have not homogeneous number of rows in each column \n\t-f1=$n_column_f1 \n\t-f2=$n_column_f2 \n\t-f3=$n_column_f3\n ";
     }
  }
  my $aref_files=\@inputfiles;
  my @filesformat;
  if ($opt_f) {
      @filesformat=split(/,/,$opt_f);
      my $n_formats=scalar(@filesformat);
      if ($n_files != $n_formats) {
	  die "Sorry, the number of files specified ($opt_i) is not the same than the number of formats specified for these files ($opt_f).\n";
      }
      my $formats;
      foreach $formats (@filesformat) {
	  unless ($formats eq '0/L'||'0/E'||'1/L'||'1/E'||'S' ) {
            die "The format $formats is not permited. You only can choose the values 0/L, 0/E, 1/L, 1/E or S. See help for more information\n";
	}
      }
  } else {
      foreach $file (@inputfiles) {
	  push @filesformat, '0/L';
      }
  }
  my $aref_formats=\@filesformat;
  my $outputfile=$opt_o||'SeqCoordComb.tab';
  my $outputformat=$opt_c||'0/L';

  return ($aref_files, $aref_formats, $outputfile, $outputformat);
}

=head2 process_file

  Usage: my $new_data_aref=process_file($nr, $file, $format, $data_aref);
  Desc: This subroutine parse a file, and process the data according the format
  Ret: An array reference of hash references with the data, $new_data_aref;
  Args: $nr, number of the file, $file, the name of the file, $format, the format of the file 
        and $data_aref, the array ref with the hash ref
  Side_Effects: use other subroutines and print processing status messages 
  Example: my $new_data_aref=process_file($file, $format, $data_aref);

=cut

sub process_file {
  my $nr=shift;
  my $file=shift;
  my $format=shift;
  my $data_aref=shift;

  open my $fh, '<', $file || die "I can not open the file $file\n";
  my $n_lines=`cut -f1 $file | wc -l`;
  chomp($n_lines);
  print "\n* Processing file ($nr):$file with format:$format and $n_lines lines\n";

  my %global_n_perc=%{$$data_aref[0]}; 
  my %global_coord5p=%{$$data_aref[1]};
  my %global_coord3p=%{$$data_aref[2]};
  my %global_raw_l=%{$$data_aref[3]};
  my %global_clean=%{$$data_aref[4]};
  my %global_trim=%{$$data_aref[5]};
  my %cleantags_changecode=%{$$data_aref[6]};
  my %check_id=%{$$data_aref[7]};

  my $count=0;
  while (<$fh>) {
      chomp ($_);
      my @data;
      my @predata=split(/\t/, $_);
      foreach my $dta (@predata) {
	  $dta =~ s/\s+//;
	  push @data, $dta;
      }
      my $id=$data[0];
      $count++;
      print "\t=> processing the line $count of $n_lines file lines\r";

      unless (exists $check_id{$id} ) {
	  $check_id{$id}='1';
      } else {
	  $check_id{$id}='+';
      }

      my ($new5coord, $new3coord);	  
      if (exists $global_coord5p{$id} && $global_coord5p{$id} ) {
	  my $old5coord=$global_coord5p{$id};
	  my $old3coord=$global_coord3p{$id};

	  if ($format eq 'S') {
	      my $n_percent=$data[1];
	      unless (exists $global_n_perc{$id}) {
		  $global_n_perc{$id}=$n_percent;
	      }
	      $new5coord=$data[2];
	      $new3coord=$data[3];
	      my $raw_length=$data[4];
	      unless (exists $global_raw_l{$id}) {
		  $global_raw_l{$id}=$raw_length;
	      }
	      my $clean_tag=$data[5];
	      if ($clean_tag) {
		  unless (exists $global_clean{$id}) {
		      if (!$opt_t) {
			  $global_clean{$id}=$clean_tag;
		      } else {
			  $global_clean{$id}=$cleantags_changecode{$clean_tag};
		      }
		  }
	      }
	      my $new_trim_tags=$data[6];
	      if ($new_trim_tags) {
		  my $old_trim_tags;
		  $new_trim_tags =~ s/;\s*$//;
		  $new_trim_tags =~ s/;\s+/;/g;
		  unless (exists $global_trim{$id}) {
		      $global_trim{$id}=combine_trim_tags($nr, $old_trim_tags, $new_trim_tags, $raw_length, 
							  $old5coord, $old3coord);
		  } else {
		      $old_trim_tags=$global_trim{$id};
		      $global_trim{$id}=combine_trim_tags($nr, $old_trim_tags, $new_trim_tags, $raw_length, 
							  $old5coord, $old3coord);
		  }
	      }
	  } else {
	      my $first_coord=$data[1];
	      my $second_coord=$data[2];
	      ($new5coord, $new3coord)=transform_coord_from_format($format, $first_coord, $second_coord);
	  }
	  my ($current_global5coord, $current_global3coord)=combine_coord($old5coord, $old3coord, $new5coord, $new3coord);
	  $global_coord5p{$id}=$current_global5coord;
	  $global_coord3p{$id}=$current_global3coord;
      } else {
	  if ($format eq 'S') {
	      $global_n_perc{$id}=$data[1];
	      $global_coord5p{$id}=$data[2];
	      $global_coord3p{$id}=$data[3];
	      $global_raw_l{$id}=$data[4];
	      if (!$opt_t) {
		  $global_clean{$id}=$data[5];
	      } else {
		  $global_clean{$id}=$cleantags_changecode{$data[5]};
	      }
	      my $first_trimtag_constructor;
	      if ($data[6]) {
		  $data[6] =~ s/;\s*$//;
		  $data[6] =~ s/;\s+/;/g;
		  my @first_trim_tags=split(/\;/, $data[6]);
		  foreach my $first_trim_tag (@first_trim_tags) {
		      if ($first_trim_tag =~ m/trimpoly/) {
			  my @trimpolytags=process_trimpoly_tag($nr,$first_trim_tag,$data[4]);
			  foreach my $trimpoly (@trimpolytags) {
			      $first_trimtag_constructor .= $trimpoly;
			      $first_trimtag_constructor .= '; ';
			  }
		      } else {
			  my $first_trim_tag_ntagged= '('.$nr.')_'.$first_trim_tag;
			  $first_trimtag_constructor .= $first_trim_tag_ntagged;
			  $first_trimtag_constructor .= '; ';
		      }
		  }
		  $first_trimtag_constructor =~ s/;\s+$//;
	      }
	      $global_trim{$id}=$first_trimtag_constructor;  
	  } else {
	      my $first_coord=$data[1];
	      my $second_coord=$data[2];
	      ($global_coord5p{$id}, $global_coord3p{$id})
		  =transform_coord_from_format($format, $first_coord, $second_coord);
	  }
      }
  }
  close $fh;

  my $added=0;
  my @ids=keys %check_id;
  foreach my $single_id (@ids) {
      my $check=$check_id{$single_id};
      if ($check eq '+') {
	  $check_id{$single_id}='-';
      } elsif ($check eq '1') {
	  $added++;
	  $check_id{$single_id}='-';
      } elsif ($check eq '-') {
	  unless (exists $global_clean{$single_id}) {
	      $global_clean{$single_id}=$opt_r||'removed';
	  }
	  $check_id{$single_id}='-';
      }
  }
  print "\n\t$added new ids were added to the processing hashes\n\n";
  my $new_data_aref=[\%global_n_perc, \%global_coord5p, \%global_coord3p, \%global_raw_l, \%global_clean, \%global_trim,
		     \%cleantags_changecode, \%check_id];
  return $new_data_aref;
}

=head2 transform_coord_from_format

  Usage: my ($begin, $end)= transfom_coord($format, $first, $second)
  Desc: change the coordenates acording with a format into the format that use BioPerl to trim the sequences
  Ret: Two scalars, the coordenates acording the format given
  Args: $format, a scalar with the format as (0/L,1/L,0/E or 1/E), and the coordenates ($first, $second)
  Side_Effects: none
  Example: my ($begin, $end)= transfom_coord($format, $first, $second)

=cut

sub transform_coord_from_format {
    my $format=shift;
    my $first=shift;
    my $second=shift;
    my ($begin,$end);

    if ($format eq 'S') {
	$format='1/E';
    }
    if ($format eq '0/L') {
	$begin=$first+1;
	$end=$second+$first-1;
    } elsif ($format eq '1/L') {
	$begin=$first;
	$end=$second+$first-1;
    } elsif ($format eq '0/E') {
	$begin=$first+1;
	$end=$second+1;
    } else {
	$begin=$first;
	$end=$second;
    }
    return ($begin,$end);
}

=head2 transform_coord_to_format

  Usage: my ($first, $second)= transfom_coord($format, $begin, $end)
  Desc: change the coordenates from BioPerl format (1/E) to the format $format 
  Ret: Two scalars, the coordenates acording the format given
  Args: $format, a scalar with the format as (0/L,1/L,0/E or 1/E), and the coordenates ($first, $second)
  Side_Effects: none
  Example: my ($first, $second)= transfom_coord($format, $begin, $end)

=cut

sub transform_coord_to_format {
    my $format=shift;
    my $begin=shift;
    my $end=shift;
    my ($first,$second);

    if ($format eq 'S') {
	$format='1/E';
    }
    if ($format eq '0/L') {
	$first=$begin-1;
	$second=$end+1-$begin;
    } elsif ($format eq '1/L') {
	$first=$begin;
	$second=$end+1-$begin;
    } elsif ($format eq '0/E') {
	$first=$begin-1;
	$second=$end-1;
    } else {
	$first=$begin;
	$second=$end;
    }
    return ($first,$second);
}

=head2 combine_coord

  Usage: my ($global5coord, $global3coord)=combine_coord($old5coord, $old3coord, $new5coord, $new3coord);
  Desc: combine coord in 1/E format
  Ret: Two scalars, the new global coordenates
  Args: The old coordenates and the new coordenates (example for file 1 and file 2);
  Side_Effects: none
  Example:my ($global5coord, $global3coord)=combine_coord($old5coord, $old3coord, $new5coord, $new3coord);

=cut

sub combine_coord {
  my $old_coord5p=shift;
  my $old_coord3p=shift;
  my $new_coord5p=shift;
  my $new_coord3p=shift;

  my $global_coord5p=$old_coord5p+$new_coord5p-1;
  my $global_coord3p=$old_coord5p+$new_coord3p-1;

  return ($global_coord5p, $global_coord3p);
}


=head2 combine_trim_tags

  Usage: my $global_trim_tags=combine_trim_tags($nr,$old_trim_tags,$new_trim_tags,$raw_length,$old_coord5p,$old_coord3p)
  Desc: Split the trim tag scalar into elements, calculate the new coordenates for poliA and other 
        elements and store into a new scalar 
  Ret: $global_trim_tags, a scalar with all the tags separated by ;
  Args: $nr, file number, $old_trim_tags and $new_trim_tags scalars, $raw_length and $old_coord5p 
        and $old_coord3p with the old coordenates
  Side_Effects: none
  Example: my $global_trim_tags=combine_trim_tags($nr, $old_trim_tags, $new_trim_tags, $old_coord5p, $old_coord3p)

=cut

sub combine_trim_tags {
    my $nr=shift;
    my $old_trim_tags=shift;
    my $new_trim_tags=shift;
    my $raw_length=shift;
    my $old_coord5p=shift;
    my $old_coord3p=shift;
    
    my @old_trim_tags_a;
    if ($old_trim_tags) {
	@old_trim_tags_a=split(/\;/, $old_trim_tags);
    }
    my @new_trim_tags_a=split(/\;/, $new_trim_tags);
    foreach my $trim_tag_element (@new_trim_tags_a) {
	if ($trim_tag_element =~ m/trimpoly/) {
	    my @new_trimpolytags=process_trimpoly_tag($nr, $trim_tag_element, $raw_length, $old_coord5p, $old_coord3p);
	    foreach my $new_trimpoly (@new_trimpolytags) {
		my $new_trimpoly_ntagged='('.$nr.')_'.$new_trimpoly;
		push @old_trim_tags_a, $new_trimpoly_ntagged;
	    }
	} elsif ($trim_tag_element =~ m/(.+)<(\d+)-(\d+)>/) {
	    my $name=$1;
            my $tag_coord;
	    my ($tag_coord5p,$tag_coord3p)=combine_coord($old_coord5p, $old_coord3p, $2, $3);
	    if ($name =~ m/^\((\d+)\)_/) {
	        my $old_nr=$1;
	        my $new_nr=$old_nr+$nr-1;
	        $name =~ s/\($old_nr\)/\($new_nr\)/;
                $tag_coord=$name.'<'.$tag_coord5p.'-'.$tag_coord3p.'>';
            } else {
	        $tag_coord='('.$nr.')_'.$name.'<'.$tag_coord5p.'-'.$tag_coord3p.'>';
            }
	    push @old_trim_tags_a, $tag_coord;
	} else {
	    my $trim_tag_element_ntagged='('.$nr.')_'.$trim_tag_element;
	    push @old_trim_tags_a, $trim_tag_element_ntagged;
	}
    }
    my $global_trim_tag;
    foreach my $trim_tag (@old_trim_tags_a) {
	$global_trim_tag .= $trim_tag.'; ';
    }
    $global_trim_tag =~ s/;\s+$//;
    return $global_trim_tag;
}

=head2 parse_new_cleantags_file

  Usage: my $cleantags_changecode_href=parse_new_cleantags_file($file);
  Desc: parse the clean code into a hash and return it
  Ret: A hash reference with the old clean tag like keys and the new like values
  Args: $file, file name of the clean code changes
  Side_Effects: die if there are more than one old clean tag with the same name
  Example: my $cleantags_changecode_href=parse_new_cleantags_file($file);

=cut 

sub parse_new_cleantags_file {
    my $file=shift;
    
    my %cleantags_changecode;
    open my $fh_cleantags, '<', $file || die "I can not open the new clean tags file: $file\n";
    while (<$fh_cleantags>) {
	chomp($_);
	my @input=split(/\t/, $_);
	unless (exists $cleantags_changecode{$input[0]} ) {
	    $cleantags_changecode{$input[0]}=$input[1];
	} else {
	    die "Sorry the old clean tags must be UNIQUES (There are more than one $cleantags_changecode{$input[0]})\n";
	}
    }
    close $fh_cleantags;
    return \%cleantags_changecode;    
}

=head2 process_trimpoly_tag

  Usage: my ($trimpoly5tag, $trimpoly3tag)=process_trimpoly_tag($nr,$trimpoly_tag,$raw_length,$old_coord5p,$old_coord3p);
  Desc: process the trimpoly tag into two variables with the sequences coordenates.
  Ret: Two scalars, $trimpoly5tag and $trimpoly3tag, new trimpoly tags with sequence coordenates
  Args: $nr, file number, $trimpoly_tag, old trimpoly tag, $raw_length, sequence length, 
        and $old_coord5p and $old_coord3p, seq coord.
  Side_Effects: none
  Example: my ($trimpoly5tag, $trimpoly3tag)=process_trimpoly_tag($nr,$trimpoly_tag,$raw_length,$old_coord5p,$old_coord3p);

=cut 

sub process_trimpoly_tag {
  my $nr=shift;
  my $trimpoly_tag=shift;
  my $raw_length=shift;
  my $old_coord5p=shift;
  my $old_coord3p=shift;
  if (!$old_coord5p && !$old_coord3p) {
      $old_coord5p=1;
      $old_coord3p=$raw_length;
  }
  $trimpoly_tag =~ m/trimpoly\[\+(\d+),\s*-(\d+)]/;
  my $five_polyA=$1;
  my $three_polyA=$2;
  my @new_trim_tags;

  if ($five_polyA > 0) {
      my $begin5=$old_coord5p;
      my $end5=$old_coord5p+$five_polyA-1;
      my $five_polyA_tag='('.$nr.')_putative_polyA<'.$begin5.'-'.$end5.'>';
      push @new_trim_tags, $five_polyA_tag;
  }
  if ($three_polyA > 0) {
      my $begin3=$old_coord3p-$three_polyA+1;
      my $end3=$old_coord3p;
      my $three_polyA_tag='('.$nr.')_putative_polyA<'.$begin3.'-'.$end3.'>';
      push @new_trim_tags, $three_polyA_tag;
  }
  return @new_trim_tags;
}
