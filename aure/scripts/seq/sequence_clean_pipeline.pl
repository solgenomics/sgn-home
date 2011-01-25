#!/usr/bin/perl

=head1 NAME
 
 sequence_clean_pipeline.pl.
 Perl script to use different sequence clean programs as seqclean, 
 nseg and trf.

=cut

=head1 SYNOPSIS

Set the environment variables for the programs.

export LUCY='lucy executable file path'
export SEQCLEAN='seqclean executable file path'
export NSEG='nseg executable file path'
export TRF='trg executable file path'

sequence_clean_pipeline.pl -i <fasta_file> -a <argument file> [-X] [-h]
    
=head2 I<Flags:>

=over
      
=item -i 

B<Input fasta>          input file for clean. (mandatory)

=item -q

B<Input qual>           qual file (mandatory when lucy program is used)
      
=item -a

B<Argument file>        argument file with all the argument for the 
                        programs seqclean, nseg and trf
   
=item -o

B<output basename>      output basename (by default: input fasta filename)

=item -X 

B<Print argument filr>  print an empty argument file
           
=item -h 

B<Help>                 print the help

=back

=cut

=head1 DESCRIPTION

  This script clean a fasta file. Use seqclean tool 
  (http://compbio.dfci.harvard.edu/tgi/software/) to do it, 
  but this script check the output file. If this output has in its report 
  a trimmed sequence number more than 0, rerun the tool over the clean 
  fasta file. 
  Finally take the all the cleanning coordenates files, calculate the final 
  clean coordenates and put in a tab file.
  This script make a cleanning files folder and moves all the files to it. 

  Also run over the sequences the programs nseg (low complexity) and trf 
  (tandem repeat finder)
      
  The output file (in the cleaning folder in the work dir) are:

=over

=item - 

  cleaning_1 folder: Is the folder where the seqclean program put the 
  cleaning slices files before ensambled all. 
  This folder are removed before run seqclean another time.
          
=item - 

  .cdix file: binary file with indexes for each run of the program.
          
=item - 

  .clean file: fasta file with the cleanned sequences for each run of the 
  program.
          
=item - 

  .cln file: report file with the clean coordenates for each run of the  
  program. 
          
=item - 

  seqcl<input_fasta_file>.log: report file of the process for each run of 
  the program.
          
=item - 

  err_seqcl<input_fasta_file>.log: report file with possible errors for each 
  run of the program
          
=item - 

  outparts_cln_sort: file with the all the search slices for the last run of 
  the program.

=item - 

  (tag -o): output file with the clean coordenates for the global process.
          
=item - 

  global_seqclean_report.txt: report of all seqclean process.
          
=item - 

  equiv_clean_codes: .tab file with the equivalences between files used in 
  the screening and the clean_code

Note:	  
	  
  The use of the argument -A in th seqclean program require from the use of 
  other tags like -N

=back

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez.
 (ab782@cornell.edu).

=cut

=head1 METHODS

sequence_clean_pipeline.pl


=cut

use strict;
use File::Basename;
use Getopt::Std;
use Bio::SeqIO;
use Bio::Seq::PrimaryQual;
use Bio::Tools::SeqStats;
use Math::BigFloat;
use Cwd;

our ($opt_i, $opt_q, $opt_a, $opt_o, $opt_X, $opt_h);
getopts("i:q:a:o:Xh");


if ($opt_X) {
    print_argument_file();
}
elsif (!$opt_i && !$opt_q && !$opt_a && !$opt_o && !$opt_X && !$opt_h) {
             print "There are n\'t any tags. Print help\n\n";
             help();
}

## FIRST VALIDATE INPUT (fasta file sequence and argument file)

my $input_fasta = $opt_i ||
    die("ARGUMENT ERROR: None -i <input_fasta_file> was used.\n");
my $input_qual = $opt_q;
my $arg_file = $opt_a ||
    die("ARGUMENT ERROR: None -a <argument file> was used.\n");

my $output_basename = $opt_o || $input_fasta;

## The processing sequence file will change 
## with different methods (they are concatenated)

my $currfasta = $input_fasta;
my $currqual = $input_qual;

## SECOND VALIDATE ENV VARIABLES

my $lucy_exec = $ENV{'LUCY'};
my $seqclean_exec = $ENV{'SEQCLEAN'};
my $nseg_exec = $ENV{'NSEG'};
my $trf_exec = $ENV{'TRF'};

if (!$lucy_exec && !$seqclean_exec && !$nseg_exec && !$trf_exec) {
    die("ENVIRONMENT VARIABLE: None of the environment variables with the program executables are set.\n")
}

## also it should return error if lucy is defined but there is not
## qual file

if($lucy_exec && !$input_qual) {
    die("ARGUMENT PROGRAM ERROR: -q <input_qual_file> was not defined for lucy.\n");
}

## THIRD PARSE THE ARGUMENT FILE

print STDERR "\n1) PARSE ARGUMENT FILE.\n";
my %args = parse_arg_file($arg_file);

## FORTH GET THE SEQ LENGTHS

my $seqio = Bio::SeqIO->new( -file => "$input_fasta", -format => 'fasta' );

print STDERR "\n2) EXTRACT SEQUENCE LENGTH FROM INPUT FILE TO CREATE BASE CLN FILE.\n\n";
my $cln_data_href = {};

my $s = 0;

while (my $seqobj = $seqio->next_seq) {

    $s++;
    print STDERR "\tProcessing sequence $s         \r";

    my $seqlength = $seqobj->length();
    my $seqid = $seqobj->primary_id();
    my $seq = $seqobj->seq();

    my $seqsts_href = Bio::Tools::SeqStats->new($seqobj)
	                                  ->count_monomers();

    my $nperc = ($seqsts_href->{'N'} / $seqlength) * 100;
    my $nperc_f = Math::BigFloat->new("$nperc")
	                      ->bfround(-2);

    $cln_data_href->{$seqid} = { 
	                          id       => $seqid,
				  n_perc   => $nperc_f,
				  s_start  => 1,
				  s_end    => $seqlength+1,
				  s_length => $seqlength,
				  rm_tag   => '',
				  trim_tag => ''
				  
                               };
}
print STDERR "\n";

## FORTH RUN LUCY

my $lucy_args_href = $args{'LUCY'};

print STDERR "\n3) RUNNING LUCY.\n";

if (defined $lucy_args_href && defined $lucy_exec) {
     
    print STDERR "\t[Enabled].\n\n";

    my ($out_seq, $out_qual, $out_cln) = run_lucy($currfasta, $currqual, \%args);
    
    ## Get the cln data from cln file

    my %new_cln = ();

    open my $cln_fh, '<', $out_cln 
	|| die("OPEN FILE ERROR: $out_cln file can not be openned (system error: $!).\n");

    while(<$cln_fh>) {
	chomp($_);

	my @data = split(/\t/, $_);

	$new_cln{$data[0]} = { 
	                          id       => $data[0],
				  n_perc   => $data[1],
				  s_start  => $data[2],
				  s_end    => $data[3],
				  s_length => $data[4],
				  rm_tag   => $data[5],
				  trim_tag => $data[6]				  
                               };
    }

    ## It will create a new seq file that probably other programs will use

    $currfasta = $out_seq;
    $currqual = $out_qual;

    ## Also the cln file should be added to the cln base file

    $cln_data_href = add_cln_data($cln_data_href, \%new_cln);
}
else {
    print STDERR "\t[disabled].\n";
}

## FIFTH RUN SEQCLEAN

my $seqclean_args_href = $args{'SEQCLEAN'};

print STDERR "\n4) RUNNING SEQCLEAN.\n";

if (defined $seqclean_args_href && defined $seqclean_exec) {

    print STDERR "\t[Enabled].\n\n";
    my ($outfile_cln, $outfile_clean) = run_seqclean($currfasta, \%args);
	
    ### Now it will parse the output
    
    my $fcln_data_href = {};

    my ($rm_count, $trim_count);
    ($fcln_data_href, $rm_count, $trim_count) = parse_seqclean_cln($outfile_cln);
    
    $cln_data_href = add_cln_data($cln_data_href, $fcln_data_href);

    my $r = 1;

    if (defined $args{'RERUN_SEQCLEAN'}) {

	if ($args{'RERUN_SEQCLEAN'}->{'-C'} =~ m/remove/i) {

	    while ($rm_count > 0) {

		print STDERR "\n\n==============================================================================================\n";
		print STDERR "* RERUN SEQCLEAN OPTION ENABLED: sequence removed are: $rm_count > 0... reruning seqclean\n";
		print STDERR "* (STEP: $r)";
		print STDERR "\n==============================================================================================\n\n";

		$r++;

		my $input_seqclean_file = $outfile_clean;

		($outfile_cln, $outfile_clean) = run_seqclean($input_seqclean_file, \%args); 
		my ($new_cln_data_href, $new_rm_count, $new_trim_count)  = parse_seqclean_cln($outfile_cln);
		
		$cln_data_href = add_cln_data($cln_data_href, $new_cln_data_href);
		$rm_count = $new_rm_count;
	    }
		
	}

	if ($args{'RERUN_SEQCLEAN'}->{'-C'} =~ m/trim/i) {

	    while ($trim_count > 0) {

		print STDERR "RERUN SEQCLEAN OPTION ENABLED: sequence trimmed are: $trim_count > 0... reruning seqclean ($r)\n\n";
		$r++;

		my $input_seqclean_file = $outfile_clean;

		($outfile_cln, $outfile_clean) = run_seqclean($input_seqclean_file, \%args); 
		my ($new_cln_data_href, $new_rm_count, $new_trim_count)  = parse_seqclean_cln($outfile_cln);
		
		$cln_data_href = add_cln_data($cln_data_href, $new_cln_data_href);
		$trim_count = $new_trim_count;
	    }		
	}		
    }
}
else {
    print STDERR "\t[disabled].\n";
}

## SIXTH RUN NSEG

my $nseg_args_href = $args{'NSEG'};

print STDERR "\n5) RUNNING NSEG.\n";

if (defined $nseg_args_href && defined $nseg_exec) {

    print STDERR "\t[Enabled].\n";
    my $outfile_nseg = run_nseg($currfasta, \%args);
    my $nseg_data_href = parse_masked_seq($outfile_nseg);

    my $nseg_lowcompl;
    if (defined $args{'FILTER_NSEG'}) {
	$nseg_lowcompl = $args{'FILTER_NSEG'}->{'locomp_perc'};
    }
    
    $cln_data_href = add_masked_data($cln_data_href, $nseg_data_href, 'nseg', 'low complexity', $nseg_lowcompl);
}
else {
    print STDERR "\t[Disabled].\n";
}

## SEVENTH RUN TRF

my $trf_args_href = $args{'TRF'};

print STDERR "\n6) RUNNING TRF.\n";

if (defined $trf_args_href && defined $trf_exec) {

    print STDERR "\t[Enabled].\n";

    my ($outfile_trf_masked, $outfile_trf_data) = run_trf($currfasta, \%args);
    my $trf_data_href = parse_masked_seq($outfile_trf_masked);

    my $trf_lowcompl;
    if (defined $args{'FILTER_TRF'}) {
	$trf_lowcompl = $args{'FILTER_TRF'}->{'locomp_perc'};
    }
    
    $cln_data_href = add_masked_data($cln_data_href, $trf_data_href, 'trf', 'simple tandem repeat', $trf_lowcompl);
}
else {
    print STDERR "\t[Disabled].\n";
}

## EIGHT PRINT THE OUTPUT FILES

## Now it will print into a file the data

print STDERR "\n\n7) PRINTING OUTPUT FILES:\n";
my $global_outfile = $output_basename . ".global.cln";

open my $og_fh, '>', $global_outfile ||
    die("OPEN FILE ERROR: File:$global_outfile could not be opened (system error:$!).\n");

my %clndata = %{ $cln_data_href};
foreach my $id (keys %clndata) {
    my $n_perc = $clndata{$id}->{'n_perc'};
    my $s_start = $clndata{$id}->{'s_start'};
    my $s_end = $clndata{$id}->{'s_end'};
    my $s_length = $clndata{$id}->{'s_length'};
    my $rm_tag = $clndata{$id}->{'rm_tag'};
    my $trim_tag = $clndata{$id}->{'trim_tag'};
    
    print $og_fh "$id\t$n_perc\t$s_start\t$s_end\t$s_length\t$rm_tag\t$trim_tag\n";
}

my $global_outfasta = $output_basename . ".global.fasta";
my $removed_outfasta = $output_basename . ".removed.fasta";

my $out_seqio = Bio::SeqIO->new( -file => ">$global_outfasta", -format => 'fasta');
my $rem_out_seqio = Bio::SeqIO->new( -file => ">$removed_outfasta", -format => 'fasta');

my $in_seqio = Bio::SeqIO->new( -file => "$input_fasta", -format => 'fasta');

while (my $seq_ob = $in_seqio->next_seq()) {
    my $s_id = $seq_ob->primary_id();
    my $s_rm_tag = $clndata{$s_id}->{'rm_tag'};

    unless ($s_rm_tag =~ m/\w+/) {
	my $trunc = $seq_ob->trunc($clndata{$s_id}->{'s_start'}, $clndata{$s_id}->{'s_end'});
	$out_seqio->write_seq($trunc);
    }
    else {
	$seq_ob->desc("[removed by: $s_rm_tag]");
	$rem_out_seqio->write_seq($seq_ob);
    }
}

if (defined $input_qual) {

    my $global_outqual = $output_basename . ".global.qual";

    my $out_qualio = Bio::SeqIO->new( -file => ">$global_outqual", -format => 'qual');

    my $in_qualio = Bio::SeqIO->new( -file => "$input_qual", -format => 'qual');

    while (my $qual_ob = $in_qualio->next_seq()) {
	my $q_id = $qual_ob->primary_id();
	my $q_rm_tag = $clndata{$q_id}->{'rm_tag'};

	unless ($q_rm_tag =~ m/\w+/) {
	    my $subqual_aref = $qual_ob->subqual($clndata{$q_id}->{'s_start'}, $clndata{$q_id}->{'s_end'});
	    my $newqual = Bio::Seq::PrimaryQual->new( -id => $q_id, -qual => join(' ', @{$subqual_aref}));
	    $out_qualio->write_seq($newqual);
	}
    }
}



print STDERR "\t...done.\n\n";


## NINTH, PRINT A REPORT

print STDERR "\n\n8) PRINTING SEQUENCE REMOVED REPORT:\n";

my %rm_tags_count;

foreach my $t_id (keys %clndata) {
    my $rm_data = $clndata{$t_id}->{'rm_tag'} || 'None';

    if (exists $rm_tags_count{$rm_data}) {
	$rm_tags_count{$rm_data}++;
    }
    else {
	$rm_tags_count{$rm_data} = 1;
    }
}

print STDERR "\n================================================================\n";
my $header = sprintf('%-20s','TAG:') . ' ' . sprintf('%20s', 'SEQ_COUNT:');
print STDERR "\t$header\n";
print STDERR "================================================================\n";

foreach my $tag (keys %rm_tags_count) {
    my $tag_f = sprintf('%-20s',$tag);
    my $count_f = sprintf('%20s',$rm_tags_count{$tag});

    print STDERR "\t$tag_f $count_f\n";
}
print STDERR "\n================================================================\n\n";



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
      This script clean a fasta file. Use seqclean tool 
      (http://compbio.dfci.harvard.edu/tgi/software/) to do it, with the 
      inprovement that check the output file. If this output has in its 
      report a trimmed sequence number more than 0, rerun the tool over the 
      clean fasta file. Finally take the all the cleanning coordenates files, 
      calculate the final clean coordenates and put in a tab file. This 
      script make a cleanning files folder and moves all the files to it. 

      The output file (in the cleaning folder in the work dir) are:
	  - cleaning_1 folder: Is the folder where the seqclean program put 
            the cleaning slices files before ensambled all. This folder are 
            removed before run seqclean another time.
          - .cdix file: binary file with indexes for each run of the program.
          - .clean file: fasta file with the cleanned sequences for each run 
            of the program.
          - .cln file: report file with the clean coordenates for each run of 
            the program. 
          - seqcl<input_fasta_file>.log: report file of the process for each 
            run of the program.
          - err_seqcl<input_fasta_file>.log: report file with possible errors 
            for each run of the program
          - outparts_cln_sort: file with the all the search slices for the 
            last run of the program.

          - (tag -o): output file with the clean coordenates for the global 
            process.
          - global_seqclean_report.txt: report of all seqclean process.
          - equiv_clean_codes: .tab file with the equivalences between files 
            used in the screening and the clean_code

    Usage: 
      Set the environment variables for the programs.

      export LUCY='lucy executable file path'
      export SEQCLEAN='seqclean executable file path'
      export NSEG='nseg executable file path'
      export TRF='trg executable file path'

      sequence_clean_pipeline.pl -i <fasta_file> -a <argument file> [-X] [-h]

    Flags:
      -i input fasta          input file for clean. (mandatory)
      -q input qual           quality file (mandatory when lucy option is used)
      -a argument file        argument file with all the arguments used in 
                              the programs (mandatory)
      -o output basename      output basename (by default: input fasta filename)
      -X print arg file       print argument file
      -h help                 print this help

EOF
exit (1);

}

=head2 run_lucy

 Usage: my ($outfile_seq, $outfile_qual, $outfile_cln) = run_lucy($input_fasta_file, $input_qual_file, $arg_file_href)
 Desc: Run the lucy command and return the output file names
 Ret: $outfile_seq, the fasta file cleanned
      $outfile_qual, the qual file cleanned
 Args: $input_file, a scalar, the filename of the input file in fasta format
       $arg_file_href, a hash reference with keys=programs and value=hash reference with
                       key=argument and  value=value
 Side_Effects: None
 Example: my ($outfile_seq, $outfile_qual, $outfile_cln) = run_lucy($input_file, $input_qual_file, $arg_file_href);

=cut 

sub run_lucy {
    my $input_fasta = shift ||
	die("ARGUMENT ERROR: None input fasta name was supplied to th internal function run_lucy.\n");
    my $input_qual = shift ||
	die("ARGUMENT ERROR: None input qual name was supplied to th internal function run_lucy.\n");
    my $args_file_href = shift ||
	die("ARGUMENT ERROR: Argument file hash reference was not supplied to the internal function run_seqclean.\n");

    my $lucy_args_href = $args_file_href->{'LUCY'};

    my %lucy_args = %{$lucy_args_href};
    
    my @args_command_line;

    foreach my $arg (keys %lucy_args) {
	my $value = $lucy_args{$arg};

	if (!$value) {
	    push @args_command_line, $arg;
	}
	elsif ($value =~ m/ENABLE/i) {
	    push @args_command_line, $arg;
	}
	else {
	    push @args_command_line, "$arg $value";
	}
    }

    ## Also it will define by default the output as:

    my $preout_fasta = 'lucyout_' . $input_fasta . '.lucyclean.fasta';
    my $preout_qual = 'lucyout_' . $input_fasta . '.lucyclean.qual';

    push @args_command_line, "-output $preout_fasta $preout_qual";
	
    my $arg_command_line = join(' ', @args_command_line);

    my $lucy_command = $ENV{'LUCY'} . " " . $arg_command_line . ' ' . $input_fasta . ' ' . $input_qual;
    
    print STDERR "Lucy command:\n$lucy_command\n\n";

    print STDERR "==== Lucy output messages ====\n\n";
    system($lucy_command);
    print STDERR "\n\n==== Lucy done ====\n\n";

    ## After the run, it will take both files to take the clean sequence and create the 
    ## .cln file

    ## First, get the clean coordinates

    my %cleancoord = ();

    my $out_fasta = $input_fasta . '.lucyclean.fasta';
    my $out_qual = $input_fasta . '.lucyclean.qual';
    my $out_cln = $input_fasta . '.lucyclean.cln';

    my $seqio = Bio::SeqIO->new( -file => "$preout_fasta", -format => 'fasta' );

    my $out_seqio = Bio::SeqIO->new( -file => ">$out_fasta", -format => 'fasta' );

    print STDERR "\n\n";
 
    while (my $seqobj = $seqio->next_seq()) {
	
	my $id = $seqobj->primary_id();
	my $descr = $seqobj->description();

	my @coord_data = split(/\s+/, $descr);
	my $start = $coord_data[3];
	my $end = $coord_data[4];

	print STDERR "\tExtracting sequence coordinates for lucy trim ($id, $start, $end).        \r";

	## Now it will trim the seqobject and create a new one from there
	
	my $subseq = $seqobj->subseq($start, $end);

	my $newseq = Bio::Seq->new( -seq => $subseq, -id => $id );
	$out_seqio->write_seq($newseq);

	$cleancoord{$id} = [$start, $end];
    }
    print STDERR "\n\n";

    ## Also it will trim the qual

    my $qualio = Bio::SeqIO->new( -file => "$preout_qual", -format => 'qual' );
    
    my $out_qualio = Bio::SeqIO->new( -file => ">$out_qual", -format => 'qual' );

    while (my $qualobj = $qualio->next_seq()) {

	my $qid = $qualobj->primary_id();
	my $qstart = $cleancoord{$qid}->[0];
	my $qend = $cleancoord{$qid}->[1];

	my $subqual_aref = $qualobj->subqual($qstart, $qend);
	my $subqual = join(' ', @{$subqual_aref});
	
	print STDERR "\tProcessing qual ($qid, $qstart, $qend).           \r";

	my $newqual = Bio::Seq::PrimaryQual->new( -id => $qid, -qual => $subqual );
	$out_qualio->write_seq($newqual);
    }
    print STDERR "\n\n";

    ## finally it will create the cln file. To do it it will open the 
    ## input fasta (the output fasta for lucy has not the sequences removed)

    my $in_seqio = Bio::SeqIO->new( -file => "$input_fasta", -format => 'fasta' );

    open my $outcln_fh, '>', $out_cln 
	|| die("OPEN FILE ERROR: $out_cln file can not be openned (system error: $!)\n");
    

    while (my $inseqobj = $in_seqio->next_seq() ) {
	my $id = $inseqobj->primary_id();
	my $seq = $inseqobj->seq();
	my $olength = length($seq);

	$seq =~ s/N//;
	
	my $flength = length($seq);
	my $nperc = 100 - ($flength/$olength)*100;

	## to take only two numbers 

	my $n = Math::BigFloat->new($nperc);
	my $f_nperc = $n->bfround(-2);

	## cln format will be (separated by tabs): id, n_perc, s_start, s_end, s_length
	## rm_tag, trim_tag

	my ($sstart, $send, $rm_tag, $trim_tag) = (1, $olength, '', '');

	if (exists $cleancoord{$id}) {

	    my @trimtags = ();

	    $sstart = $cleancoord{$id}->[0];
	    $send = $cleancoord{$id}->[1];

	    if ($sstart > 1) {
		push @trimtags, 'lucy_qualtrim:<1-' . $sstart . '>';
	    }
	    if ($send < $olength) {
		push @trimtags, 'lucy_qualtrim:<' . $send . '-' . $olength .'>';
	    }
	    $trim_tag = join('; ', @trimtags);
	}
	else {

	    $rm_tag = 'low_qual';

	    if (defined $lucy_args{'-minimum'}) {

		## In this case the sequence was removed by size
		
		if ($lucy_args{'-minimum'} > $olength) {	
		    $rm_tag = 'short';
		}		
	    }
	}
	print $outcln_fh "$id\t$f_nperc\t$sstart\t$send\t$olength\t$rm_tag\t$trim_tag\n";
	print STDERR "\tprinting cln file ($id, $f_nperc, $sstart, $send, $olength, $rm_tag, $trim_tag).    \r";
    }
    print STDERR "\n\n";

    return($out_fasta, $out_qual, $out_cln);
}



=head2 run_seqclean

 Usage: my ($outfile_cln, $outfile_clean) = run_seqclean($input_file, $arg_file_href)
 Desc: Run the seqclean command and return the output file names (report and sequence)
 Ret: $outfile_cln, a report file of the clean
      $outfile_clean, a clean sequence file in fasta format
 Args: $input_file, a scalar, the filename of the input file in fasta format
       $arg_file_href, a hash reference with keys=programs and value=hash reference with
                       key=argument and  value=value
 Side_Effects: None
 Example: my ($outfile_cln, $outfile_clean) = run_seqclean($input_file, $arg_file_href);

=cut 

sub run_seqclean {
    my $input_file = shift ||
	die("ARGUMENT ERROR: None input file name was supplied to th internal function run_seqclean.\n");
    my $args_file_href = shift ||
	die("ARGUMENT ERROR: Argument file hash reference was not supplied to the internal function run_seqclean.\n");

    my $seqclean_args_href = $args_file_href->{'SEQCLEAN'};

    my %seqclean_args = %{$seqclean_args_href};
    
    my @args_command_line;

    foreach my $arg (keys %seqclean_args) {
	my $value = $seqclean_args{$arg};

	if (!$value) {
	    push @args_command_line, $arg;
	}
	elsif ($value =~ m/ENABLE/i) {
	    push @args_command_line, $arg;
	}
	else {
	    push @args_command_line, "$arg $value";
	}
    }
    my $arg_command_line = join(' ', @args_command_line);

    my $log_file = $input_file . '.seqclean_analysis.log';

    my $seqclean_command = $ENV{'SEQCLEAN'} . " " . $input_file . " " . $arg_command_line . ' > ' . $log_file;

    print STDERR "==== Seqclean output messages ====\n\n";
    system($seqclean_command);
    print STDERR "\n\n==== Seqclean done ====\n\n";

    my $output_cln = $seqclean_args{'-r'} || $input_file . '.cln';
    my $output_clean = $seqclean_args{'-o'} || $input_file . '.clean';

    return ($output_cln, $output_clean);
}

=head2 parse_seqclean_cln

 Usage: my ($cln_data_href, $rm_count, $trim_count)  = parse_seqclean_cln($outfile_cln)
 Desc: Parse the output cln file for seqclean
 Ret: $cln_data_href, a hash reference with keys=id and value=hash reference with
                      keys=data_type and value=value
      $rm_count, number of sequences removed from the output sequence file
      $trim_count, number of sequences trimmed in the output sequence file
 Args: $outfile_cln, seqclean output cln file
 Side_Effects: None
 Example: my ($cln_data_href, $rm_count, $trim_count)  = parse_seqclean_cln($outfile_cln)

=cut 

sub parse_seqclean_cln {
    my $cln_file = shift ||
	die("ARGUMENT ERROR: None cln file name was supplied to th internal function parse_seqclean_cln.\n");
    
    my %clean_data;
    my ($rm_count, $trim_count) = (0, 0);

    open my $oc_fh, '<', $cln_file ||
	die("OPEN FILE ERROR: File:$cln_file can not be openned (system error: $!).\n");

    print STDERR "\n";
    my $l = 0;

    while(<$oc_fh>) {
	chomp($_);
	$l++;
     
	## There are many empty spaces in the parse. It will removed after split into elements
	
	my @predata = split(/\t/, $_);

	my @data;

	foreach my $predata (@predata) {
	    $predata =~ s/^\s+//;
	    $predata =~ s/\s+$//;
	    push @data, $predata;
	}

	## Load the data into an hash

	$clean_data{$data[0]} = { 
	                          id       => $data[0],
			   	  n_perc   => $data[1] || '',
				  s_start  => $data[2] || '',
				  s_end    => $data[3] || '',
				  s_length => $data[4] || '',
				  rm_tag   => $data[5] || '',
				  trim_tag => $data[6] || ''
	                        };

	if ($data[5] =~ m/\w+/) {
	    $rm_count++;
	}
	if ($data[6] =~ m/\w+/) {
	    $trim_count++;
	}
    }
    print STDERR "\n";
    return (\%clean_data, $rm_count, $trim_count);
}


=head2 add_cln_data

 Usage: $cln_data_href = add_cln_data($cln_data_href, $new_cln_data_href);
 Desc: add the cln data to an old cln data and calculate the new coordinates
 Ret: A hash reference with the new coordinates where keys=ids and values=hash reference
      with keys=data type and value=value
 Args: $cln_data_href and $new_cln_data_href, hash reference where keys=ids and 
       value=hash reference with keys=data type and value=value
 Side_Effects: none
 Example: $cln_data_href = add_cln_data($cln_data_href, $new_cln_data_href);

=cut

sub add_cln_data {
    my $cln_data_href = shift ||
	die("ARGUMENT ERROR: None cln hash reference was supplied to th internal function add_cln_data.\n");

    my $new_cln_data_href = shift ||
	die("ARGUMENT ERROR: None new cln hash reference was supplied to th internal function add_cln_data.\n");

    my %clndata_comb = ();

    my %clndata = %{$cln_data_href};
    my %newclndata = %{$new_cln_data_href};

    foreach my $id (keys %clndata) {
	
	## The data to treat are: id, n_perc, s_start, s_end, s_length, rm_tag and trim_tag.
	## Id will be the same but it the new cln probably loose some of the ids removed in the last cln
	## In that case it will preserve the old data
	
	if (defined $newclndata{$id}) {

	    $clndata_comb{$id} = { 'id' => $id };

	    ## The N percentage always should be refered to the first cln file
	    ## The trimming can change this value, but it will keep the data from the original sequence

	    $clndata_comb{$id}->{'n_perc'} = $clndata{$id}->{'n_perc'};

	    ## s_start always should be refered to the initial sequence

	    $clndata_comb{$id}->{'s_start'} = $clndata{$id}->{'s_start'} + $newclndata{$id}->{'s_start'} - 1;
	    $clndata_comb{$id}->{'s_end'} = $newclndata{$id}->{'s_end'} + $clndata{$id}->{'s_start'} - 1;

	    ## s_length always should be the initial 

	    $clndata_comb{$id}->{'s_length'} = $clndata{$id}->{'s_length'};

	    ## rm_tag always will be when it is not null	    	   

	    if ($clndata{$id}->{'rm_tag'} =~ m/\w+/) {
		$clndata_comb{$id}->{'rm_tag'} = $clndata{$id}->{'rm_tag'};
	    } 
	    else {
		$clndata_comb{$id}->{'rm_tag'} = $newclndata{$id}->{'rm_tag'};
	    }
	    
	    ## finally trim_tag need to parse the trim tag coordinates
	    ## There are two types: 
	    ##    + trimpoly[+N, -M] where N are the base removed from the 5' and M the bases removed from
	    ##      the 3'. In this case it will sum with the previous one
	    ##    + sequence_accession:<start-end> where start and end are the corrdinates for the match
	    
	    my (%trim_tags);

	    my $trim_tags = $clndata{$id}->{'trim_tag'};
	    if (defined $trim_tags) {
		
		my @trim_tags = split(/;/, $trim_tags);
		
		foreach my $trim (@trim_tags) {
		    $trim =~ s/^\s+//g;
		    $trim =~ s/\s+$//g;

		    if ($trim =~ m/trimpoly\[\+(\d+), -(\d+)\]/) {
			my $pa_five = $1;
			my $pa_three = $2;
			
			if (exists $trim_tags{'trimpoly'}) {
			    $trim_tags{'trimpoly'}->[0] += $pa_five;
			    $trim_tags{'trimpoly'}->[1] += $pa_three;			    
			}
			else {
			    $trim_tags{'trimpoly'} = [$pa_five, $pa_three];
			}
		    }
		    elsif ($trim =~ m/(.+):<(\d+)-(\d+)>/) {
			my $acc = $1;
			my $a_start = $2;
			my $a_end = $3;

			if (exists $trim_tags{$acc}) {
			    push @{$trim_tags{$acc}}, [$acc, $a_start, $a_end];
			}
			else {
			    $trim_tags{$acc} = [[$acc, $a_start, $a_end]];
			}	
		    }
		    else {
			$trim_tags{$trim} = [[$trim]];
		    }
		}
	    }
		
	    my $new_trim_tags = $newclndata{$id}->{'trim_tag'};
	    if (defined $new_trim_tags) {
	
		my @newtrim_tags = split(/;/, $new_trim_tags);
	    
		foreach my $newtrim (@newtrim_tags) {
		    $newtrim =~ s/^\s+//g;
		    $newtrim =~ s/\s+$//g;
		    
		    if ($newtrim =~ m/trimpoly\[\+(\d+), -(\d+)\]/) {
			my $pa_five = $1;
			my $pa_three = $2;
			
			if (exists $trim_tags{'trimpoly'}) {
			    $trim_tags{'trimpoly'}->[0] += $pa_five;
			    $trim_tags{'trimpoly'}->[1] += $pa_three;			    
			}
			else {
			    $trim_tags{'trimpoly'} = [$pa_five, $pa_three];
			}
		    }
		    elsif ($newtrim =~ m/(.+):<(\d+)-(\d+)>/) {
			my $acc = $1;
			my $a_start = $clndata{$id}->{'s_start'} + $2 - 1;
			my $a_end = $clndata{$id}->{'s_start'} - 1 + $3;

			## It wil need modify the coordinates as other start and ends
			if (exists $trim_tags{$acc}) {
			    push @{$trim_tags{$acc}}, [$acc, $a_start, $a_end];
			}
			else {
			    $trim_tags{$acc} = [[$acc, $a_start, $a_end]];
			}
		    }
		    else {
			$trim_tags{$newtrim} = [[$newtrim]];
		    }
		}	      
	    }
		
	    ## Now it will recompose the trim_tag
	    my @tag_list;

	    foreach my $trimtag (keys %trim_tags) {
		
		if ($trimtag =~ m/trimpoly/) {
		    my $f = $trim_tags{$trimtag}->[0];
		    my $t = $trim_tags{$trimtag}->[1];

		    push @tag_list, "trimpoly[+$f, -$t]";		    
		}
		else {

		    ## Trimtag will contain arrays of arrays

		    my @tags_arefs = @{$trim_tags{$trimtag}};

		    foreach my $tag_aref (@tags_arefs) {

			if (scalar(@{$tag_aref}) == 3) {
			    my $acc = $tag_aref->[0];
			    my $acc_start = $tag_aref->[1];
			    my $acc_end = $tag_aref->[2];
		    
			    push @tag_list, "$acc:<$acc_start-$acc_end>";
			}
			else {
			    push @tag_list, $tag_aref->[0];
			}
		    }
		}
	    
		$clndata_comb{$id}->{'trim_tag'} = join('; ', @tag_list);
	    } 
	}
	else {
	    $clndata_comb{$id} = $clndata{$id};
	}
    }

    return \%clndata_comb;
}

=head2 run_nseg

 Usage: my $outfile_nseg = run_nseg($input_file, $arg_file_href)
 Desc: Run the nseg command and return the output file name
 Ret: $outfile_nseg, the stdout for nseg
 Args: $input_file, a scalar, the filename of the input file in fasta format
       $arg_file_href, a hash reference with keys=programs and value=hash reference with
                       key=argument and  value=value
 Side_Effects: None
 Example: my ($outfile_nseg) = run_nseg($input_file, $arg_file_href);

=cut 

sub run_nseg {
    my $input_file = shift ||
	die("ARGUMENT ERROR: None input file name was supplied to th internal function run_nseg.\n");
    my $args_file_href = shift ||
	die("ARGUMENT ERROR: Argument file hash reference was not supplied to the internal function run_nseg.\n");


    my $outfile_nseg = $input_file . ".nseg_analysis.fasta";

    my $nseg_args_href = $args_file_href->{'NSEG'};

    my %nseg_args = %{$nseg_args_href};
    
    my @args_command_line;

    ## Add some deafult parameters as -x and -c 100 

    unless ($nseg_args{'-x'} =~ m/\w+/) {
	$nseg_args{'-x'} = 1;
    }
    unless ($nseg_args{'-c'} =~ m/\w+/) {
	$nseg_args{'-c'} = 100;
    }

    my @c_args = ('-window', '-locut', '-hicut');

    ## Window, locut and hicut args will be in order

    foreach my $c_arg (@c_args) {
	if (defined $nseg_args{$c_arg} ) {
	    push @args_command_line, $nseg_args{$c_arg};
	}
    }

    ## The rest of the argument do not need order

    foreach my $arg (keys %nseg_args) {

	my $value = $nseg_args{$arg};
	unless ($arg =~ m/window|locut|hicut/i) {
	    if ($arg =~ m/-c|-m|-t|-z/) {
		if ($value =~ m/\w+/) {
		    push @args_command_line, "$arg $value";
		}
	    }
	    else {
		if ($value =~ m/1|enable/i) {
		    push @args_command_line, "$arg";
		}
	    }
	}
    }
    my $arg_command_line = join(' ', @args_command_line);

    my $nseg_command = $ENV{'NSEG'} . " " . $input_file . " " . $arg_command_line . ' > ' . $outfile_nseg;

    print STDERR "\n==== Nseg output messages (if therer are any) ====\n\n";
    system($nseg_command);
    print STDERR "\n\n==== Nseg done ====\n\n";
 
    return $outfile_nseg;
}

=head2 parse_masked_seq

 Usage: my $data_href = parse_masked_seq($outfile_masked_fasta);
 Desc: Parse masked sequences in fasta format and return the percentage of masking
 Ret: $data_href, a hash reference with keys=id and value=masked_percentage                       
 Args: $outfile_fasta, the out filename
 Side_Effects: None
 Example: my $nseg_data_href = parse_nseg($outfile_nseg);

=cut

sub parse_masked_seq {
    my $out = shift ||
	die("ARGUMENT ERROR: None output nseg file name was supplied to th internal function parse_masked_seq.\n");

    my $seqio = Bio::SeqIO->new( -file => "$out", -format => 'fasta' );

    my %data = ();

    while (my $seqobj = $seqio->next_seq() ) {
	
	my $id = $seqobj->primary_id();
	my $seq = $seqobj->seq();

	my $t_length = length($seq);
	$seq =~ s/N//ig;

	my $length = length($seq);

	$data{$id} = 100 - ($length/$t_length)*100;
    }
    return \%data;
}

=head2 add_masked_data

 Usage: $cln_data_href = add_masked_data($cln_data_href, $masked_data_href, $rm_tag, $trim_tag, $mask_perc_cutoff);
 Desc: add the masked data to an old cln data
 Ret: A hash reference with the new coordinates where keys=ids and values=hash reference
      with keys=data type and value=value
 Args: $cln_data_href hash reference where keys=ids and 
                      value=hash reference with keys=data type and value=value
       $masked_data_href, hash reference with key=ids and
                      value=masked percentage
       $rm_tag, a scalar, rm_tag used in the masked
       $trim_tag, a scalar, trim_tag used for the cln file
       $mask_perc_cutoff, a scalar, an integer that will be used to decide what id
                      will be tag with nseg

 Side_Effects: none
 Example: $cln_data_href = add_masked_data($cln_data_href, $nseg_data_href, 'nseg', 'low_complexity');

=cut

sub add_masked_data {
   my $cln_data_href = shift ||
	die("ARGUMENT ERROR: None cln hash reference was supplied to th internal function add_masked_data.\n");

   my $mask_data_href = shift ||
	die("ARGUMENT ERROR: None mask hash reference was supplied to th internal function add_masked_data.\n");

   ## By default the values for rm_tag and trim_tag will be seq_masked and mask_percentage
   
   my $rm_tag = shift || 'seq_masked';
   
   my $trim_tag = shift || 'mask_precentage';

   my $mask_perc_cutoff = shift || 50; ## By default it will remove any seq with > 50% of low complexity

   my %clndata_comb = ();

   my %clndata = %{$cln_data_href};
   my %maskdata = %{$mask_data_href};
   
   foreach my $id (keys %clndata) {
	
       ## nseg only will modify rm_tag and trim_tag adding rm_tag=nseg and trim_tag=low_complexity(XX%)
       
       $clndata_comb{$id} = $clndata{$id}; ## By default it will take the same values

       if (defined $maskdata{$id}) {

	    my $mask_perc = $maskdata{$id};
	    my $simply_val;
	    if ($mask_perc =~ m/(\d+)\.(\d+)/) {
		$simply_val = $1;
	    }

	    if ($mask_perc > $mask_perc_cutoff) {

		if ($clndata{$id}->{'rm_tag'} =~ m/\w+/) {
		    $clndata_comb{$id}->{'rm_tag'} = $clndata{$id}->{'rm_tag'} . "; $rm_tag";
		} 
		else {
		    $clndata_comb{$id}->{'rm_tag'} = "$rm_tag";
		}
		
		if ($clndata{$id}->{'trim_tag'} =~ m/\w+/) {
		    $clndata_comb{$id}->{'trim_tag'} = $clndata{$id}->{'trim_tag'} . "; $trim_tag($simply_val%)";
		} 
		else {
		    $clndata_comb{$id}->{'trim_tag'} = "$trim_tag($simply_val%)";
		}
	    }	    
	}
    }

    return \%clndata_comb;
}

=head2 run_trf

 Usage: my $outfile_trf = run_trf($input_file, $arg_file_href)
 Desc: Run the trf command and return the output file names
 Ret: $outfile_nseg_mask, a fasta file with the sequences masked 
      $outfile_nseg_dat, a dat file
 Args: $input_file, a scalar, the filename of the input file in fasta format
       $arg_file_href, a hash reference with keys=programs and value=hash reference with
                       key=argument and  value=value
 Side_Effects: None
 Example: my ($outfile_trf_mask, $outfile_trf_dat) = run_nseg($input_file, $arg_file_href);

=cut 

sub run_trf {
    my $input_file = shift ||
	die("ARGUMENT ERROR: None input file name was supplied to th internal function run_trf.\n");
    my $args_file_href = shift ||
	die("ARGUMENT ERROR: Argument file hash reference was not supplied to the internal function run_trf.\n");

    ## First, define the outputfile names.

    my $trf_args_href = $args_file_href->{'TRF'};

    my %ta = %{$trf_args_href};

    ## Define the default variables.

    my %def_trfa = ( 
	             '-Match'     => 2,
		     '-Mismatch'  => 7, 
		     '-Delta'     => 7, 
		     '-PM'        => 80, 
		     '-PI'        => 10, 
		     '-Minscore'  => 50, 
		     '-MaxPeriod' => 500,
	             '-f'         => 1,
	             '-m'         => 1,
	             '-h'         => 1,
	           );

    ## Replace the trf arguments for the defaukt variables if they are not defined (without any value)

    foreach my $def_arg (keys %def_trfa) {

	if (defined $ta{$def_arg}) {
	    unless ($ta{$def_arg} =~ m/^\d+$/) {
		$ta{$def_arg} = $def_trfa{$def_arg};
	    }
	}
	else {
	    $ta{$def_arg} = $def_trfa{$def_arg};
	}
    }
    
    my @base = ($input_file, $ta{'-Match'}, $ta{'-Mismatch'}, $ta{'-Delta'}, $ta{'-PM'}, $ta{'-PI'}, $ta{'-Minscore'}, $ta{'-MaxPeriod'});
    my $base = join('.', @base);

    my $outfile_trf_mask = $base . '.mask';
    my $outfile_trf_dat = $base . '.dat';
    my $log_file = $base . '.log';

    my @args_command_line;

    ## The args for trf should be ordered

    my @trfas = ('-Match', '-Mismatch', '-Delta', '-PM', '-PI', '-Minscore', '-MaxPeriod' );
    
    foreach my $trfa (@trfas) {

	## All the arguments should be defined (if the user has not define them, they were defined with default values)
	
	push @args_command_line, $ta{$trfa};
    }
	
    ## Add the options when they are enabled or 1

    foreach my $trf_opt (keys %ta) {
	if ($trf_opt =~ m/-m|-f|-d|-h/) {
	    my $value = $ta{$trf_opt};
	    if ($value =~ m/^1|enable/i) {
		push @args_command_line, $trf_opt;
	    }
	}
    }
	   
    my $arg_command_line = join(' ', @args_command_line);

    my $trf_command = $ENV{'TRF'} . " " . $input_file . " " . $arg_command_line .' > ' . $log_file;

    print STDERR "\n==== Trf output messages (if there are any) ====\n\n";
    system($trf_command);
    print STDERR "\n\n==== Trf done ====\n\n";
 
    return($outfile_trf_mask, $outfile_trf_dat);
}



=head2 parse_arg_file

 Usage: my %arg_file = parse_arg_file($file)
 Desc: Parse the argument file and return an hash with the argumens 
 Ret: %arg_file, a hash with:
                 keys=program
                 value=hash reference with
                       keys=argument
                       value=value
                       
 Args: $file, a file with the arguments
 Side_Effects: None
 Example: my %arg_file = parse_arg_file($file)

=cut

sub parse_arg_file {
    my $file = shift ||
	die("ARGUMENT ERROR: None file was supplied to the parse argument function.\n");

    open my $afh, '<', $file ||
	die("OPEN FILE ERROR: File:$file can not be openned. (system error: $!).\n");

    my %args = ();

    while (<$afh>) {
	chomp($_);

	unless ($_ =~ m/^#/) {
	    if ($_ =~ m/^\*(\w+)-(\w+):\s+\[(.*)\]/) {
		my $prog = $1;
		my $arg = '-'.$2;
		my $val = $3;

		if (exists $args{$prog}) {
		    $args{$prog}->{$arg} = $val;
		}
		else {
		    $args{$prog} = {$arg => $val};
		}
	    }
	}
    }
    return %args;
}

=head2 print_argument_file

  Usage: print_argument_file()
  Desc: print an argument file
  Ret: none
  Args: none
  Side_Effects: create the argument file
  Example: print_argument_file()

=cut

sub print_argument_file {
    
    my $dir = getcwd;
    my $arg_file = $dir . '/seqclean_argument_file.log';

    open my $arg_fh, '>', $arg_file || die("OPEN FILE ERROR: File:$arg_file can not be openend (system error: $!).\n");

    print STDERR "PRINTING ARGUMENT FILE option...\n";

    my $info = '#

## Note1: To disable an option use # or delete the line.
## Note2: To enable options without values (for example: SEQCLEAN-A, write enable or 1). 
## Note3: About vector and adaptor removing. Use seqclean and keep lucy to trim
##        sequences based in qscore
## Note4: Some lucy options have been removed from this pipeline:
##        [-vector vector_sequence_file splice_site_file]
##        [-size vector_tag_size] 
##        [-threshold vector_cutoff]
##        [-output sequence_filename quality_filename] 
##        (it will create the files input_fasta.lucyclean.seq and input_fasta.lucyclean.qual)


#########################################################################################
# LUCY ARGUMENT #####################################################################
#########################################################################################


*LUCY-pass_along: [] min_value max_value med_value (separated by spaces)(default values: 0 0 0)
*LUCY-range:      [] area1 area2 area3 (default values: 40 60 100)
*LUCY-alignment:  [] area1 area2 area3 (default values: 8 12 16)
*LUCY-cdna:       [] [minimum_span maximum_error initial_search_range] (skip by default)
*LUCY-keep:       [] (skip bu default)
*LUCY-minimum:    [] good_sequence_length (100 by default)
*LUCY-debug:      [] [filename]
*LUCY-error:      [] max_avg_error max_error_at_ends (0.025 and 0.02 by default)
*LUCY-window:     [] window_size max_avg_error [window_size max_avg_error ...] (default values: 50 0.08 10 0.3)
*LUCY-bracket:    [] window_size max_avg_error (default values: 10 0.02)
*LUCY-quiet:      [] 
*LUCY-inform_me:  []
*LUCY-xtra:       [] cpu_threads


#########################################################################################
# SEQCLEAN ARGUMENT #####################################################################
#########################################################################################

*SEQCLEAN-c: []  use the specified number of CPUs on local machine (default 1) or a list of PVM nodes in <PVM_nodefile>
*SEQCLEAN-n: []  number of sequences taken at once in each search slice (default 2000)
*SEQCLEAN-v: []  comma delimited list of sequence files to use for end-trimming of <seqfile> sequences  (usually vector sequences)
*SEQCLEAN-l: []  during cleaning, consider invalid the sequences sorter than <minlen> (default 100)
*SEQCLEAN-s: []  comma delimited list of sequence files to use for screening <seqfile> sequences for contamination (mito or other sps)
*SEQCLEAN-r: []  write the cleaning report into file <reportfile> (default: <seqfile>.cln)
*SEQCLEAN-o: []  output the "cleaned" sequences to file <outfasta> (default: <seqfile>.clean)
*SEQCLEAN-x: []  minimum percent identity for an alignemnt with a contaminant (default 96)
*SEQCLEAN-y: []  minimum length of a terminal vector hit to be considered (>11, default 11)
*SEQCLEAN-N: []  disable trimming of ends rich in Ns (undetermined bases)
*SEQCLEAN-M: []  disable trashing of low quality sequences
*SEQCLEAN-A: []  disable trimming of polyA/T tails
*SEQCLEAN-L: []  disable low-complexity screening (dust)
*SEQCLEAN-I: []  do not rebuild the cdb index file
*SEQCLEAN-m: []  send e-mail notifications to <e-mail>
*SEQCLEAN-t: []  write the log file into the file <logfile>
*SEQCLEAN-e: []  write the error file into the file <errorfile>

*RERUN_SEQCLEAN-C: []  repeat the clean until does not return more clan tags

#####################################################################################
# NSEG ARGUMENT #####################################################################
#####################################################################################

*NSEG-window: []  OPTIONAL window size (default 21) 
*NSEG-locut:  []  OPTIONAL low (trigger) complexity (default 1.4) 
*NSEG-hicut:  []  OPTIONAL high (extension) complexity (default locut + 0.2)  
*NSEG-x:      []  each input seq is represented by a single output seq with low-complexity regions replaced by strings of x char 
*NSEG-c:      []  number of sequence characters/line (default 60)
*NSEG-m:      []  minimum length for high-complexity segment (default 0). Shorter segments are merged with adj. low-complexity segments 
*NSEG-l:      []  show only low-complexity segments (fasta format) 
*NSEG-h:      []  show only high-complexity segments (fasta format) 
*NSEG-a:      []  show all segments (fasta format) 
*NSEG-n:      []  do not add complexity information to the header line 
*NSEG-o:      []  show overlapping low-complexity segments (default merge) 
*NSEG-t:      []  maximum trimming of raw segment (default 100) 
*NSEG-p:      []  prettyprint each segmented sequence (tree format) 
*NSEG-q:      []  prettyprint each segmented sequence (block format) 
*NSEG-z:      []  period

*FILTER_NSEG-locomp_perc: [] low complexity percentage cutoff. 

#####################################################################################
# TRF ARGUMENT #####################################################################
#####################################################################################

*TRF-Match:     []  matching weight
*TRF-Mismatch:  []  mismatching penalty
*TRF-Delta:     []  indel penalty
*TRF-PM:        []  match probability (whole number)
*TRF-PI:        []  indel probability (whole number)
*TRF-Minscore:  []  minimum alignment score to report
*TRF-MaxPeriod: []  maximum period size to report
*TRF-m:         []  masked sequence file
*TRF-f:         []  flanking sequence
*TRF-d:         []  data file
*TRF-h:         []  suppress html output
*TRF-r:         []  no redundancy elimination

*FILTER_TRF-locomp_perc: [] tandem simple repeat percentage cutoff.

//';

    print $arg_fh "$info";
    print STDERR "...done (printed argument file: $arg_file)\n\n";
    exit (1);
}




####
1;##
####
