#!/usr/bin/perl

=head1 NAME

 fasta_htmltag.pl
 
 Print fasta sequence in fasta format with html tags (version.0.1).

=cut

=head1 SYPNOSIS

 fasta_htmltag.pl [-h] -f <fasta_input_file> 
                            -t <tag_file>
                            [-e <external_tag_file>]
                            -o <output_file>
                            [-R]

=head2 I<Flags:>

=over


=item -f

B<fasta_file>           fasta file (mandatory)

=item -t

B<tag_file>             snp file with 4 columns seq_id, start, end, 
                        html_tag and modification

=item -o

B<output_file>          output filename (outfile by default)

=item -e

B<external_tag_file>    print an external tag in a sequence region

=item -R

B<print_with_rule>      print the sequences with a nucleotide rule

=item -h

B<help>                   print the help

=back

=cut

=head1 DESCRIPTION

 This script replace the nucleotide in the specified region 
 by modification, adding a html tag before and after the modification
 region.

 For example:

 >example1
 ATGTGTGTGTGTGATGCGTGCATGGCTTTGCGATGCGTAGGACTGCCAGTTTGACTTGTG
 ACAACAACCAGCTGACGCAC

 example1 20 25 <span style="color: rgb(100, 0, 253);">catgg</span> 

 It will produce:

 >example1
 ATGTGTGTGTGTGATGCGTG<span style="color: rgb(100, 0, 253);">catgg</span>CTTTGCGATGCGTAGGACTGCCAGTTTGACTTGTG
 ACAACAACCAGCTGACGCAC

=cut

=head1 AUTHORS

  Aureliano Bombarely Gomez.
  (ab782@cornell.edu).

=cut

=head1 METHODS

 fasta_htmltag.pl

=cut

use strict;
use warnings;
use autodie;

use Getopt::Std;
use Bio::SeqIO;
use Math::BigFloat;


our ($opt_f, $opt_t, $opt_o, $opt_e, $opt_R, $opt_h);
getopts("f:t:o:e:Rh");
if (!$opt_f && !$opt_t && !$opt_o && !$opt_e && !$opt_R && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}
if ($opt_h) {
    help();
}


## Check arguments

my $fastafile = $opt_f 
    || die("DATA ARGUMENT ERROR: -f fasta file WAS NOT SUPPLIED.\n");

my $tagfile = $opt_t
    || die("DATA ARGUMENT ERROR: -t tag file WAS NOT SUPPLIED.\n");

my $outfile = $opt_o
    || 'outfile';

my $ext_tagfile = $opt_e;


## Parse tag file

my %tags = ();

## Use of autodie

open(my $tag_fh, '<', $tagfile); 

print STDERR "\nPARSING tag file:$tagfile\n\n";
my $l = 0;

while(<$tag_fh>) {
    chomp($_);
    $l++;

    print STDERR "\tParsing line $l from tag_file.\r";

    ## Ignore comments
    ## First column => seq_id
    ## Second column => start
    ## Third column => end
    ## Forth column => tag
	
    my @data = split(/\t/, $_);
    my $seq_id = $data[0];

    if (scalar(@data) >= 4) {
	if (exists $tags{$seq_id}) {
	    push @{$tags{$seq_id}}, [ $data[1], $data[2], $data[3] ];
	}
	else {
	    $tags{$seq_id} = [ [$data[1], $data[2], $data[3] ] ];
	}
    }
}   
print STDERR "\n\tDone.\n\n";

## Parse external tag file

my %extags = ();

if (defined $ext_tagfile) {
    
    open(my $extag_fh, '<', $ext_tagfile); 

    print STDERR "\nPARSING external tag file:$ext_tagfile\n\n";
    my $ll = 0;

    while(<$extag_fh>) {
	chomp($_);
	$ll++;

	print STDERR "\tParsing line $ll from external tag_file.\r";

	## Ignore comments
	## First column => seq_id
	## Second column => start
	## Third column => end
	## Forth column => tag
	
	my @data = split(/\t/, $_);
	my $seq_id = $data[0];

	if (scalar(@data) >= 4) {
	    if (exists $extags{$seq_id}) {
		push @{$extags{$seq_id}}, [ $data[1], $data[2], $data[3] ];
	    }
	    else {
		$extags{$seq_id} = [ [$data[1], $data[2], $data[3] ] ];
	    }
	}
    }   
    print STDERR "\n\tDone.\n\n";

}

## Now it will open the fasta file and it will do the nt changes

my $seqio = Bio::SeqIO->new( -file => $fastafile, -format => 'fasta' );
open (my $outfh, '>', $outfile );

print $outfh '<meta http-equiv="Content-Type" content="text/html; charset=ISO-8859-1">';
print $outfh "\n</head><body>\n";

print STDERR "Tagging file.\n\n";


while( my $seqobj = $seqio->next_seq() ) {
    
    ## It will print in the end of the script all the sequences from this parsing
    my @seqprint;

    my $id = $seqobj->id();
    my $seq = $seqobj->seq();
    my $length = $seqobj->length();
    my @nt_list = split(//, lc $seq);
    
    if (defined $tags{$id}) {

	my @tag_data = @{$tags{$id}};
	
	foreach my $tag_aref (@tag_data) {
	    my $st = $tag_aref->[0];
	    my $en = $tag_aref->[1];
	    my $region = $en-$st;
	    my $tag = $tag_aref->[2];
	    
	    if ($tag =~ m/^(<.+>)(.+?)(<.+>)$/) {
		my $init = $1;
		my $mod = $2;
		my $end = $3;
		
		my @mod_elements = split(//, $mod);
	    
		## It will count the length of modif
		my $mod_l = length($mod);
		if ($region == $mod_l) {
			
		    my $p = 0;
		    foreach my $el (@mod_elements) {
			$p++;
			if ($p == 1) {
			    
			    $nt_list[$st+$p-2] = $init.$el;
			    if ($mod_l == 1) {
				$nt_list[$st+$p-2] .= $end;
			    }
			    
			}
			elsif ($p == $mod_l) {
			    $nt_list[$st+$p-2] = $el.$end;
			}
			else {
			    $nt_list[$st+$p-2] = $el;
			}
		    }		
		}
		else {
		    print STDERR "MSG: Modification region have different length($st,$en <=> $mod) . Skipping tagging.\n";
		}
	    }
	    elsif ($tag =~ m/^(<.+>)(<.+>)$/) {
		
		my $init = $1;
		my $end = $2; 
		
		$nt_list[$st-1] = $init.$nt_list[$st-1];
		$nt_list[$en-1] = $nt_list[$en-1].$end;
	    }
	}
	
	## Now it will print the sequence in fragments of 60 nt.

	print $outfh "<pre>&gt;$id\n";

	my $ntcount = 0;
	my $total_count = 0;
	my (%ext_tag_start);

	foreach my $nt (@nt_list) {
	    print $outfh "$nt";
	    $ntcount++;
	    $total_count++;
	    if ($ntcount == 60) {
		
		if (defined $extags{$id}) {
		    my @ext_tags_aref = @{$extags{$id}};
		    foreach my $ext_tag_aref (@ext_tags_aref) {
			if ($total_count >= $ext_tag_aref->[0] && $total_count <= $ext_tag_aref->[1]) {
			    if (!$ext_tag_start{$ext_tag_aref->[2]}) {
				$ext_tag_start{$ext_tag_aref->[2]} = 1;
				print $outfh "\t= $ext_tag_aref->[2]"; 
			    }
			    else {
				print $outfh "\t= +"; 
			    }
			}
		    }		    
		}
		
		if ($opt_R) {
		   print $outfh "\t-\t$total_count"; 
		}
		print $outfh "\n";
		$ntcount = 0;
	    }
	}
    }
}
print $outfh "</pre>\n</body></html>\n";




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
      
      This script replace the nucleotide in the specified region 
     by modification, adding a html tag before and after the modification
     region.

     For example:

      >example1
      ATGTGTGTGTGTGATGCGTGCATGGCTTTGCGATGCGTAGGACTGCCAGTTTGACTTGTG
      ACAACAACCAGCTGACGCAC

      example1 20 25 <span style="color: rgb(100, 0, 253);">catgg</span> 

      It will produce:

      >example1
      ATGTGTGTGTGTGATGCGTG<span style="color: rgb(100, 0, 253);">catgg</span>CTTTGCGATGCGTAGGACTGCCAGTTTGACTTGTG
      ACAACAACCAGCTGACGCAC   
    
    Usage:
      fasta_htmltag.pl [-h] -f <fasta_input_file> -t <tag_file> 
                            -o <output_file> [-R]
      
    Example:
      fasta_htmltag.pl [-h] -f example.fasta -t tag.tab -o example_tagged.html
      
    Flags:

      -f <fasta_file>           fasta file (mandatory)
      -t <tag_file>             snp file with 4 columns seq_id, start, end, 
                                html_tag and modification
      -o <output_file>          output filename (outfile by default)
      -R <print_with_rule>      print the sequences with a nucleotide rule
      -h <help>                 print the help

     
EOF
exit (1);
}

