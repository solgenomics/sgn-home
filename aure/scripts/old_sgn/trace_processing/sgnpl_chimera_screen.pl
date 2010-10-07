#!/usr/bin/perl

=head1 NAME

sgnpl_chimera_screen.pl.
A script to search possible chimera sequences in a fasta file (version.2.0.).

=head1 SYPNOSIS
  
sgnpl_chimera_screen.pl -i <fasta_file> -d <dataset> [-b <blast_script>] [-p <blast_program>] [-s <selfblast file>] 
[-m <permuted_dataset>] [-e <e-value>]
    
=head2 I<Flags:>

=over

=item -i 

input fasta file (mandatory)

=item -b

script or program used to do the blast (blastall by default but you can use qsblast.pl)
      
=item -p 

program used for the blast (blastx by default)
      
=item -d 

dataset to compare in fasta file (mandatory) 
      
=item -s 

selfblast dataset results in m8 format
      
=item -m 

permuted dataset selfblast file
      
=item -e 

cuttoff e-value for the blast
      
=item -h 

print this help

=back

=head1 DESCRIPTION

This script search chimera sequences (sequences are made up by parts of diferent genes) in an input fasta file.
To do it, use a blast search over a protein sequence dataset. The files are:

=over

=item I.Input files:

=over

=item - 

Input sequences in fasta format where the script go to screen the possible chimeras.
                
=item - 

A file with a sequences dataset in a fasta format (-f <fasta> ).
Also it can use a selfblast results file in a tab format (m8) (-s <selfblast in m format>), 
or a permuted selfblast dataset file (-m <permuted selfblast dataset>).

=back
        
=item II. Output files:

=over
	       
=item - 

Output file with the ids of the sequences with possible chimeras in a tab format with three column: 

=over 

=item -f1 => 

id of the input sequences.
                                
=item -f2 =>

5'best blast result with the dataset file to compare,
                                
=item -f3 =>

3'best blast result with the dataset file to compare.   )
         
=back

=back

=back

=head1 AUTHOR

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=head1 METHODS
 
sgnpl_chimera_screen.pl


=cut

use strict;
use File::Basename;
use Getopt::Std;
use Bio::SeqIO;
use CXGN::DB::InsertDBH;

our ($opt_i, $opt_b, $opt_p, $opt_d, $opt_s, $opt_m, $opt_e, $opt_h);
getopts("i:b:p:d:s:m:e:h");

if (!$opt_i && !$opt_b && !$opt_p && !$opt_d && !$opt_s && !$opt_m && !$opt_e && !$opt_h) {
    print "There are n\'t any tags. Print help\n\n";
    help();
}

print "Checking arguments...\n";
my ($input_file, $dataset, $blast_script, $blast_program, $permuted_dataset, $evalue_cutoff)=validate_input();

our $chimera_dir=File::Basename::dirname($input_file)."/chimera_screening_files";
mkdir "$chimera_dir", 0755 or warn "Cannot make $chimera_dir directory: $!";
print "Create chimera screening directory ($chimera_dir)\n\n";

my $chimera_input=make_chimera_input($input_file);

my $chimera_blast=$chimera_dir."/".File::Basename::basename($input_file).".chimera_blast.m8";
$chimera_blast =~ s/\.fasta//gi;

print "Running $blast_script program ...\n";
print "$blast_script -p $blast_program -i $chimera_input -d $dataset -o $chimera_blast -e $evalue_cutoff -m 8\n\n";
system "$blast_script -p $blast_program -i $chimera_input -d $dataset -o $chimera_blast -e $evalue_cutoff -m 8";

my $chimera_analyze_file=analyze_chimera_blast($permuted_dataset, $chimera_blast);

print "\n----------------------------\n";
print "Report:";
print "\n----------------------------\n";
my $chimera_count=`cut -f1 $chimera_analyze_file | sort -u | wc -l`;
chomp ($chimera_count);
my $input_count=`grep '>' $input_file | wc -l`;
chomp ($input_count);
my $dataset_count=`grep '>' $dataset | wc -l`;
chomp ($dataset_count);
print "\n\tInput file:\t$input_file\n";
print "\tDataset:\t$dataset\n\n";
print "\tInput sequences:\t\t$input_count\n";
print "\tDataset sequences:\t\t$dataset_count\n";
print "\tChimera positive sequences:\t$chimera_count";
print "\n----------------------------\n";

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
      This script search chimera sequences (sequences are made up by parts of diferent genes) in an input fasta file.
      To do it, use a blast search over a protein sequence dataset. The files are:

	I.Input files:
		- Input sequences in fasta format where the script go to screen the possible chimeras.
                - A file with a sequences dataset in a fasta format (-f <fasta> ).
                  Also it can use a selfblast results file in a tab format (m8) (-s <selfblast in m format>), 
		  or a permuted selfblast dataset file (-m <permuted selfblast dataset>.								 

        II. Output files:
		- Output file with the ids of the sequences with possible chimeras in a tab format. 
                 (Three column: -f1=id of the input sequences,
                                -f2=5'best blast result with the dataset file to compare,
                                -f3=3'best blast result with the dataset file to compare.   )
         
    
    Usage: 
      sgnpl_chimera_screen.pl -i <fasta_file> -d <dataset> [-b <blast_script>] [-p <blast_program>] [-s <selfblast file>] 
      [-m <permuted_dataset>] [-e <e-value>]
    
    Example:
      sgnpl_chimera_screen.pl -i /home/aure/GenBank_sequences/Nicotiana_tabacum/tgi_sequences/SAL_KST/SAL_KST.fasta 
      -d /home/aure/screening/ATH1_pep -m /home/aure/screening/ath1_selfblastp.permuted

    Flags:
      -i input fasta file (mandatory)
      -b script used to do the blast (without script by default, so means... blastall)
      -p program used for the blast (blastx by default)
      -d dataset to compare in fasta file (mandatory) 
      -s selfblast dataset results in m8 format
      -m permuted dataset selfblast file
      -e cuttoff e-value for the blast
      -h print this help

EOF
exit (1);
}

=head2 validate_input

  Usage: my @input_arguments=validate_input();
  Desc: check the different input arguments and give defaults values when these do not exists
  Ret: Five scalars, the input filename, the input dataset filename, the blast program, the permuted dataset filename
       and the e-value cutoff for the blast searches. 
  Args: none
  Side_Effects: die if the check is wrong.
  Example: my ($input_file, $dataset, $program, $permuted_dataset, $evalue_cutoff)=validate_input();
               
=cut


sub validate_input {
    if ($opt_h) {
	help();
    }
    unless ($opt_i) {
	die("Required argument -i <input_fasta_file> was not supplied.\n");
    }
    unless ($opt_d) {
	die("Required argument -d <comparative dataset> was not supplied.\n");
    }
    if ($opt_p) {
	unless ($opt_p == ("blastn"|"blastx"|"tblastx") ) {
	    die "Sorry, the argument -b <blast_program> only can be \"blastn\", \"blastx\" or \"tblastx\"\n";
	}
    }

    my ($input_file, $dataset, $blast_script, $blast_program, $dataset_selfblast, $permuted_dataset, $dataset_dir_files, $evalue_cutoff);
    $input_file=$opt_i;
    $dataset=$opt_d;
  
    my ($dataset_in, $dataset_sq, $dataset_hr);
    my $dataset_name=File::Basename::basename($dataset);
    if (!$opt_p || $opt_p == "blastx") {
        $dataset_in=$dataset_name.".pin";
        $dataset_sq=$dataset_name.".psq";
        $dataset_hr=$dataset_name.".phr";
    } else {
	$dataset_in=$dataset_name.".nin";
        $dataset_sq=$dataset_name.".nsq";
        $dataset_hr=$dataset_name.".nhr";
    }
    my $check_ds_count=0;
    my $dataset_dir=File::Basename::dirname($opt_d);
    opendir DH, $dataset_dir or die "Cannot open $dataset_dir to check dataset format.\n";
    foreach $dataset_dir_files (readdir DH) {
	if ($dataset_dir_files eq $dataset_in) {
	    $check_ds_count++;
	} elsif ($dataset_dir_files eq $dataset_sq) {
	    $check_ds_count++;
	} elsif ($dataset_dir_files eq $dataset_hr) {
	    $check_ds_count++;
	}
    }

    if ($check_ds_count < 3) {
	print "The dataset is not formatted.\n";
	print "Formating the protein dataset ...\n";
        if (!$opt_p || $opt_p == "blastx") {
	    system "formatdb -i $opt_d -p T";
	} else {
	    system "formatdb -i $opt_d -p F";
	}
    }

    if (!$opt_e) {
	$evalue_cutoff='1e-6';
    } else {
	$evalue_cutoff=$opt_e;
    }

    if (!$opt_b) {
	$blast_script="blastall";
    } else {
	$blast_script=$opt_b;
    }
    if (!$opt_p) {
	$blast_program="blastx";
    } else {
	$blast_program=$opt_p;
    }

   
    if (!$opt_s && !$opt_m) {
	print "There are not dataset selfblast file (-s) or permuted dataset file (-m).\n";
	print "Running dataset selfblast ...\n";
        my $selfpblast;
        if (!$opt_p) {
	    $selfpblast="blastp";
        } elsif ($opt_p eq "blastx") {
            $selfpblast="blastp";
	} else {
	    $selfpblast="blastn";
	}
	$dataset_selfblast=File::Basename::dirname($opt_i)."/".File::Basename::basename($opt_d).".selfblast.m8";
        print "SELFBLAST:$blast_script -p $selfpblast -i $opt_d -d $opt_d -o $dataset_selfblast -e $evalue_cutoff -m 8\n";
	system "$blast_script -p $selfpblast -i $opt_d -d $opt_d -o $dataset_selfblast -e $evalue_cutoff -m 8";
	print "\t...make dataset selfblast file ($dataset_selfblast).\n";
	print "Making the permuted dataset file ...\n.";
	$permuted_dataset=make_permuted_ds_file($dataset_selfblast);
	print "\t...make permuted dataset file ($permuted_dataset).\n";
    } elsif ($opt_s && !$opt_m) {
	$dataset_selfblast=$opt_s;
        print "There are not permuted dataset file (-m).\n";
	print "Making the permuted dataset file ...\n.";
	$permuted_dataset=make_permuted_ds_file($dataset_selfblast);
	print "\t...make permuted dataset file ($permuted_dataset).\n";
    } elsif (!$opt_s && $opt_m) {
	$dataset_selfblast="It is not necessary";
	$permuted_dataset=$opt_m;
    } else {
        $dataset_selfblast=$opt_s;
	$permuted_dataset=$opt_m;
	print "Message: If you use -m <permuted dataset file> is not necessary use -s <dataset selfblast file>.\n";
    }
    
    my $input_seqcount=`grep '>' $opt_i | wc -l`;
    chomp($input_seqcount);
    if ($input_seqcount < 1) {
	die "Sorry, the input file have not sequences in fasta format.\n";
    } else {
	print "The input file has $input_seqcount sequences.\n";
    }
    
    return ($input_file, $dataset, $blast_script, $blast_program, $permuted_dataset, $evalue_cutoff);
}

=head2 make_permuted_ds_file

  Usage: my $permuted_dataset = make_permuted_ds_file($dataset_selfblast);
  Desc: this subroutine write in a file the possible permutations of the selfblast results file.
  Ret: A scalar, the name of the file with the results of the operation
  Args: $dataset_selfblast (the selfblast results file).
  Side_Effects: none
  Example: my $permuted_dataset = make_permuted_ds_file($dataset_selfblast);
               
=cut

sub make_permuted_ds_file {
    my $dataset_selfblast=shift;
    open my $ds_selfblast_handle, '<', "$dataset_selfblast" or die "cannot open the file $dataset_selfblast: $!";
    my $permuted_dataset=File::Basename::dirname($opt_i)."/".File::Basename::basename($opt_d).".permuted_ds";
    open my $permuted_ds_handle, '>', "$permuted_dataset";
    my %equivalence;
	  keys %equivalence = 32000;
	  while (<$ds_selfblast_handle>) {
	    my ($query, $subject) = split;
	    push @{$equivalence{$query}}, $subject 
	      unless (grep { $_ eq $subject } @{$equivalence{$query}});
	    push @{$equivalence{$subject}}, $query
	      unless (grep { $_ eq $query } @{$equivalence{$subject}});
	  }
	  foreach my $key (keys %equivalence) {
	    print  $permuted_ds_handle "$key\t", join ("\t", @{$equivalence{$key}}) . "\n";
	  }
    return $permuted_dataset;
}

=head2 make_chimera_input

  Usage: my $chimera_input_file=make_chimera_input($input_fasta);
  Desc: subroutine that get the the sequences of the input file, take two substring of 300 nt (one of the 5 end and 
        another of the 3 end).
  Ret: A scalar, the name of the fasta file with the 5 and 3 end sequences
  Args: $input_fasta (input fasta file).
  Side_Effects: none
  Example: my $chimera_input=make_chimera_input($input_file);
               
=cut

sub make_chimera_input {
    my $input_fasta=shift;
    my $in = Bio::SeqIO->new(	-file=>$input_fasta,
				-format=>"Fasta"
			 );

    my $chimera_input = $chimera_dir."/chimera_input_".File::Basename::basename($opt_i);  
    open my $chimera_input_handle, '>', "$chimera_input";
    while (my $chseq=$in->next_seq) {
	my $chseq_id = $chseq->id();
	my $chseq_length = $chseq->length();
	my $chseq_seq = $chseq->seq();
	if ($chseq_length > 300) {
		my $range = (( $chseq_length < 600 ) ? $chseq_length/2 : 300);
		my $foward = substr ($chseq_seq, 0, $range);
		my $reverse = substr ($chseq_seq, -1*$range);
		print $chimera_input_handle ">$chseq_id:5\n$foward\n>$chseq_id:3\n$reverse\n";
        }
     }
     print "Make the chimera_screen_input_fasta file ($chimera_input)\n\n";
     return $chimera_input;

}

=head2 analyze_chimera_blast

  Usage: my $analyze_chimera_blast_filename=analyze_chimera_blast($permuted_dataset, $chimera_blast);
  Desc: subroutine that calculate the possible chimera compare the blast results of the 5 and 3 sequences parts with 
        the permuted selfblast table.  
  Ret: A scalar, the name of the results file.
  Args: $permuted_dataset (file with the permuted selfblast results for the dataset) and $chimera_blast (blast results
        of the 5 and 3 sequences parts)
  Side_Effects: none.
  Example: my $results=analyze_chimera_blast($permuted_ds, $chimera_blastp);
               
=cut

sub analyze_chimera_blast {
    my $permuted_dataset=shift;
    my $chimera_blast=shift;

my %csa_equivalence;
	open my $permuted_ds_handle, '<', "$permuted_dataset" or die "cannot open the file $permuted_dataset: $!";
	while (<$permuted_ds_handle>) {
	  chomp;
	  my @record = split;
	  my $id = shift @record;
	  @{$csa_equivalence{$id}} = @record;
	}
	my %blast;
	keys %blast = 20000;

        open my $chimera_blast_handle, '<', "$chimera_blast" or die "cannot open the file $chimera_blast: $!";
        my $chimera_screen_analyze = $chimera_dir . "/OUT_chimera_screen_analyze.tab";
        open my $chimera_screen_analyze_handle, '>', "$chimera_screen_analyze";
        while (<$chimera_blast_handle>) {
	  my ($csa_query, $csa_subject) = split;
	  $blast{$csa_query} = $csa_subject unless (exists($blast{$csa_query}));
	}
	my %clones_already_seen;
	foreach my $csa_query (keys %blast) {
	  my $clone_id = $csa_query;
	  $clone_id =~ s/:[35]//;
	  next if ($clones_already_seen{$clone_id});
	  $clones_already_seen{$clone_id} = 1;
	  my $partner = $csa_query;
	  ($csa_query =~ m/3$/) ? $partner =~ s/.$/5/ : $partner =~ s/.$/3/;
	  if ($blast{$csa_query} && $blast{$partner}) {
	    if (($blast{$csa_query} ne $blast{$partner}) &&
	       (!(grep { $_ eq $blast{$csa_query} }
	       @{$csa_equivalence{$blast{$partner}}}))) {
	       print $chimera_screen_analyze_handle "$clone_id\t$blast{$csa_query}\t$blast{$partner}\n";
	       }
	    }
	}
       print "Create the chimera screen analyze file ($chimera_screen_analyze).\n";
       return $chimera_screen_analyze;
}
