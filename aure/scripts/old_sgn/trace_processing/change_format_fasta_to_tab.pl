#!/usr/bin/perl

=head1 NAME

 change_format_fasta_to_tab.pl.
 A script to change from fasta format to tab format, and add the quality values.(version.1.0.).

=cut

=head1 SYPNOSIS

change_format_fasta_to_tab.pl -H <dbhost> -D <dbname> -f <fastafile> -q <qual_file> [-r <int>] [-R]  
    
=head2 I<Flags:>

=over

=item -H 

B<database hostname>    for example localhost or db.sgn.cornell.edu (mandatory if -q is active)
      
=item -D 

B<database name>        sandbox or cxgn etc (mandatory if -q is active)
      
=item -f 

B<fasta_file>           sequence file in fasta format (mandatory)	
      
=item -q 

B<quality file>         quality values in fasta format
      
=item -r 

B<quality def value>    quality default value when there are not qual. file (15 by default)

=item -R

B<disable def. value>   disable the option that give a quality default value with the fasta file
      
=item -h 

B<help>                  show the help

=back

=cut

=head1 DESCRIPTION

This script change the format of a fasta file to tab file and add to it a default qual. value (15 by default, but you can specify with -r option). Another option is give a qual. file with the qual. value fof the sequences (-q option).
 
=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

change_format_fasta_to_tab.pl


=cut

use strict;
use Bio::SeqIO;
use File::Basename;
use Getopt::Std;
use CXGN::DB::InsertDBH;

our ($opt_H, $opt_D, $opt_f, $opt_q, $opt_r, $opt_R, $opt_h);
getopts("H:D:f:q:r:Rh");
if (!$opt_H && !$opt_D && !$opt_f && !$opt_q && !$opt_R && !$opt_r && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}
check_input();

print "\nProcessing ... $opt_f fasta file.\n\n";
my $fastafile_tab=fastaseq_to_tab();
if ($opt_q) {
        print "Processing ... $opt_q qual file.\n\n";
	my $qualfile_tab=qual_to_tab();
	print "Join the sequence tab and quality values tab files ...\n\n";
	my $output=load_db($fastafile_tab, $qualfile_tab);
	print "...Processed.\nThe output file is: $output\n\n";
} else {
	if ($opt_r && !$opt_R) {
		print "... Processed.\n\nThere is not .qual file. It script is using $opt_r by default.\n\n";
	} elsif ($opt_R) {
		print "... Processed without qual. values.\n\n";
	} else {
		print "... Processed.\n\nThere is not .qual file. It script use 15 by default.\n\n";
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
      SGN_fasta_to_tab is a script that get the fasta sequences from a file and put them in another file in tab format. If do not specify the .qual. value use 15 by default (but you can change this default value with -d tag). If you specify a .qual file, this script get the qual values and change to tab format. Also get both files and load in the database for combine them into a tab file with seq and qual values.
    
    Usage: 
      change_format_fasta_to_tab.pl -H <dbhost> -D <dbname> -f <fastafile> [-q <qual_file>] [-r <int>] [-R]

    Examples:
    1-With quality files:

      change_format_fasta_to_tab.pl -H localhost -D sandbox -f /home/aure/Nicotiana_tabacum/tgi_seq/SAL_KST.fasta 
      -q /home/aure/Nicotiana_tabacum/SAL_KST.qual
    
    2-Without quality files but the user wants a quality default value (15) in the tab file:

      change_format_fasta_to_tab.pl -f /home/aure/Sequences/test_seq.fasta

    3-The user only wants a tab file without quality values:

      change_format_fasta_to_tab.pl -f /home/aure/Sequences/test_seq.fasta -R

    Flags:
      -H database hostname    for example localhost or db.sgn.cornell.edu (mandatory)
      -D database name        sandbox or cxgn etc (mandatory)
      -f fasta file	      file with the sequences in fasta format.
      -q qual file	      file with the quality values in a fasta format.
      -r quality value def    quality value by default.
      -R disable feature      disable add qual values
      -h this help

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
  if (!$opt_f) {	
	die "Sorry, you do not specify a fasta file in the -f option.\n";
  }
  if ($opt_r =~ m/\D/) {
	die "Sorry, the quality default value must be a digit.\n";
  }
  if ($opt_q && $opt_r) {
	die "Sorry, if you give a qual. file you can not give a default value.\n";
  } elsif ($opt_R && $opt_r) {
	die "Sorry, you can not choose the option quality default value (-r) and not default value (-R).\n";
  }
  my $input_fastacount=`grep '>' $opt_f | wc -l`;
  if ($opt_q) {	
	my $input_qualcount=`grep '>' $opt_q | wc -l`;
        if ($input_fastacount ne $input_qualcount) {
	     print "Sorry, there are diferent number of id in the fasta file ($opt_f) and in the qual file ($opt_q).\n";
	}
        if (!$opt_H || !$opt_D) {
	     die "Sorry, you need specify database and localhost.\n";
        }
  }
}

=head2 fastaseq_to_tab

  Usage: my $fastaseq_tab=fastaseq_to_tab();
  Desc: Get the id and the sequences from a fasta file and write in an output file in a tab format. If there is not 
        -q or -R write by default 15 as quality value (to change this value -r option).   
  Ret: A scalar. The name of the tab format file.
  Args: none
  Side_Effects: die if can not open the input file.
  Example: my $fastaseq_tab=fastaseq_to_tab();

=cut

sub fastaseq_to_tab {
  my $fastaseq=$opt_f;
  my $outfile = $fastaseq.".tab";
  my $def_val;
  if ($opt_r) {
	$def_val = $opt_r;
  } else {
	$def_val = 15;
  }
  open my $out, '>', $outfile || die "I can not open the input file $outfile.\n";
  my $input = Bio::SeqIO->new(	-file 	=> 	$fastaseq,
  	              		-format => 	'fasta'
			);
  while ( my $seq_obj = $input->next_seq() ) {
	my $id=$seq_obj->display_id;
        my $seq=$seq_obj->seq;
	if (!$opt_q && !$opt_R) {
		my $qual_def=$seq;
		$qual_def =~ s/\w/$def_val /gi;
		print $out "$id\t$seq\t$qual_def\n";
	} elsif (!$opt_q && $opt_R) {
		print $out "$id\t$seq\n"; 
	} else {
        print $out "$id\t$seq\n";
	}
  }
  return $outfile;
}

=head2 qual_to_tab

  Usage: my $qual_tab=qual_to_tab();
  Desc: Get the id and the quality values from a fasta file and write in an output file in a tab format. 
  Ret: A scalar. The name of the tab format file.
  Args: none
  Side_Effects: die if can not open the input file.
  Example: my $qual_tab=qual_to_tab();

=cut

sub qual_to_tab {
  my $qual=$opt_q;
  my $outfile = $qual.".tab";
  open my $in, '<', $qual || die "I can not open the input file $qual.\n"; 
  open my $out, '>', $outfile || die "I can not open the input file $outfile.\n";
  my $n=1;
  while (<$in>) {
	if ($_ =~ m/^>(.+)$/gi) {
		my $id=$1;
                chomp($id);
                if ($n > 1) {
		    print $out "\n$id\t";
		} else {
		    print $out "$id\t";
		}
		$n++;
	} elsif ($_ =~ m/^((\d+\s)+)$/gi) {
		my $qual=$1;
        	chomp($qual);
                print $out "$qual ";
        }

  }
  return $outfile;
}

=head2 load_db

  Usage: my $outputfile=load_db($fastaseq_tab, $qual_tab);
  Desc: Copy the tab files (the sequences tab and the quality values tab) in the database and select id, 
        seq and quality values. Copy the selection in a output file.   
  Ret: A scalar. The name of the output tab format file.
  Args: none
  Side_Effects: die if can not open the output file.
  Example: my $output=load_db($fastaseq_tab, $qual_tab);

=cut

sub load_db {
  my $fastaseq_tab=shift;
  my $qual_tab=shift;
  my $fastaseq_temp = File::Basename::basename($fastaseq_tab);
  $fastaseq_temp =~ s/\W/_/gi;
  my $qual_temp = File::Basename::basename($qual_tab);
  $qual_temp =~ s/\W/_/gi;
  my $comp_temp=$fastaseq_temp.$qual_temp; 

  my $dbh = CXGN::DB::InsertDBH::connect ({ 	dbname => 'sandbox', 
						dbhost => 'localhost',
					  });

  $dbh->do("CREATE TEMP TABLE $fastaseq_temp (f_id varchar(250), seq text)");
  $dbh->do("COPY $fastaseq_temp FROM '$fastaseq_tab'");
  $dbh->do("CREATE TEMP TABLE $qual_temp (q_id varchar(250), qual text)");
  $dbh->do("COPY $qual_temp FROM '$qual_tab'");
  $dbh->do("SELECT f_id, seq, qual INTO TEMP TABLE $comp_temp FROM $fastaseq_temp JOIN $qual_temp ON f_id=q_id");
  my $completefile_tab=create_file('fastaqual.tab'); 
  $dbh->do("COPY $comp_temp TO '$completefile_tab'");
  return $completefile_tab;
}

=head2 create_dbloadfile

  Usage: my $newfilename=create_dbloadfile($sufix);
  Desc: make a new file and change the permissions to (-rw-rw-rw-).
  Ret: the name of the new file.
  Args: $sufix (sufix for the filename, for example if find 'est' the filename will be 'est.tab'.
  Side_Effects: none
  Example: my $output_filename=create_dbloadfile('out.tab');

=cut
 
sub create_file {
    my $sufix = shift;
    my $fastaseq_tab=$opt_f;
    
    my $dbload_dir=File::Basename::dirname($fastaseq_tab);
    my $outfile_name = $dbload_dir ."/". $sufix;
    print "Create $outfile_name file.\n";
    open my $outfile, '>', "$outfile_name" || die "Cannot create the file $outfile_name\n";
    chmod 0666, "$outfile_name";
    return $outfile_name;
}
