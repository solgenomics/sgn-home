
=head1 NAME 

Bio::Tools::Run::Phylo::Phylip::Dnaml 
- Wrapper for the phylip program dnaml

=head1 SYNOPSIS

 %runargs = ( 
    quiet       => 1,
    OUTGROUP    => 3,   
    );

 my $factory = Bio::Tools::Run::Phylo::Phylip::Dnaml->new(%runargs);
 my @tree = $factory->run($filename);


=head1 DESCRIPTION

Wrapper for dnaml Joseph Felsentein from a multiple alignment file or a
L<Bio::SimpleAlign> object and returns a L<Bio::Tree::Tree> object;

VERSION Support

This wrapper currently supports v3.5 of phylip. 
There is also support for v3.6.

=head1 PARAMETERS FOR DNAML COMPUTATION

=head2 BEST_TREE

Title		: BEST_TREE
Description	: (optional)

                 Allows the program to search for the best tree, and the 
                 User tree setting, which reads a tree or trees ("user trees") 
                 from the input tree file and evaluates them. The input tree 
                 file's default name is intree. In many cases the programs will
                 also tolerate having the trees be preceded by a line giving 
                 the number of trees

=head2 TRANS_RATIO

Title		: TRANS_RATIO
Description	: (optional)

                 Transition/transversion option, a real number greater than 
                 0.0, as the expected ratio of transitions to transversions. 
                 Note that this is not the ratio of the first to the second 
                 kinds of events, but the resulting expected ratio of 
                 transitions to transversions. The exact relationship between 
                 these two quantities depends on the frequencies in the base 
                 pools. The default value of the T parameter if you do not use 
                 the T option is 2.0.

=head2 BASEFREQ

Title		: BASEFREQ
Description	: (optional)

                 Use the empirical frequencies of the bases, observed in the 
                 input sequences, as the base frequencies, you simply use the 
                 default setting of the F option. These empirical frequencies 
                 are not really the maximum likelihood estimates of the base 
                 frequencies, but they will often be close to those values 
                 (they are the maximum likelihood estimates under a "star" or 
                 "explosion" phylogeny). If you change the setting of the 
                 BASEFREQ (F) option you will be prompted for the frequencies 
                 of the four bases. These must add to 1 and are to be typed on 
                 one line separated by blanks, not commas.

=head2 ONECATSITES

Title		: ONECATSITES
Description	: (optional)

                 The ONECATSITES (C) option allows user-defined rate categories.
                 The user is prompted for the number of user-defined rates, and
                 for the rates themselves, which cannot be negative but can be
                 zero. These numbers, which must be nonnegative (some could be 
                 0), are defined relative to each other, so that if rates for 
                 three categories are set to 1 : 3 : 2.5 this would have the 
                 same meaning as setting them to 2 : 6 : 5. The assignment of 
                 rates to sites is then made by reading a file whose default 
                 name is "categories". It should contain a string of digits 1 
                 through 9. A new line or a blank can occur after any character
                 in this string. Thus the categories file might look like this:
                    122231111122411155
                    1155333333444

=head2 RATESITES

Title		: RATESITES
Description	: (optional)

                The R (Hidden Markov Model rates) option allows the user to 
                approximate a Gamma distribution of rates among sites, or a 
                Gamma distribution plus a class of invariant sites, or to 
                specify how many categories of substitution rates there will 
                be in a Hidden Markov Model of rate variation, and what are the
                rates and probabilities for each. By repeatedly selecting the 
                RATESITES (R) option one toggles among no rate variation, the 
                Gamma, Gamma+I, and general HMM possibilities. constant rate
                bt default. If you choose Gamma or Gamma+I the program will 
                ask how many rate categories you want. In that case the format
                should be "Gamma\nX" where X is the number of categories.

=head2 SITESWEIGHT

Title		: SITESWEIGHT
Description	: (optional)

                This signals the program that, in addition to the data set, 
                you want to read in a series of weights that tell how many 
                times each character is to be counted. If the weight for a 
                character is zero (0) then that character is in effect to be 
                omitted when the tree is evaluated. If it is (1) the character
                is to be counted once. Some programs allow weights greater than
                1 as well.

=head2 SPEEDIER

Title		: SPEEDIER
Description	: (optional)

                Speedier but rougher analysis. Yes by default

=head2 GLOBALREAG

Title		: GLOBALREAG
Description	: (optional)

                After the last species is added to the tree, each possible 
                group to be removed and re-added. This improves the result, 
                since the position of every species is reconsidered. It 
                approximately triples the run-time of the program.

=head2 JUMBLE

  Title        : JUMBLE 
  Description  : (optional)

                 This enables you to tell the program to use a random
                 number generator to choose the input order of
                 species.  seed: an integer between 1 and 32767 and of
                 the form 4n+1 which means that it must give a
                 remainder of 1 when divided by 4.  Each different
                 seed leads to a different sequence of addition of
                 species.  By simply changing the random number seed
                 and re-running programs one can look for other, and
                 better trees.

  Usage        : @params = ('jumble'=>17); where 17 is the random seed.
                 Also JUMBLE => [1001, 3] where 1001 are seed number and 
                 3 the number of times.
                 Defaults to no jumble
                

=head2 OUTGROUP

Title		: OUTGROUP
Description	: (optional)

                  This option selects the species to be used as the outgroup
                  Usage @params = ('outgroup'=>'X'); 
                  where X is an positive integer not more than the 
                  number of sequences 
                  Defaults to 1

=head2 MULTIPLE

Title    : MULTIPLE
Description: (optional)

          This allows multiple distance matrices to be generated from multiple
          MSA.
          Usage: @params = ('MULTIPLE'=>100) where the value specifyies the 
          number of aligments given.




=head1 FEEDBACK

=head2 Mailing Lists

User feedback is an integral part of the evolution of this and other
Bioperl modules. Send your comments and suggestions preferably to one
of the Bioperl mailing lists.  Your participation is much appreciated.

  bioperl-l@bioperl.org                  - General discussion
  http://bioperl.org/wiki/Mailing_lists  - About the mailing lists

=head2 Support 
 
Please direct usage questions or support issues to the mailing list:
  
L<bioperl-l@bioperl.org>
  
rather than to the module maintainer directly. Many experienced and 
reponsive experts will be able look at the problem and quickly 
address it. Please include a thorough description of the problem 
with code and data examples if at all possible.

=head2 Reporting Bugs

Report bugs to the Bioperl bug tracking system to help us keep track
the bugs and their resolution.  Bug reports can be submitted via the
web:

  http://bugzilla.open-bio.org/

=head1 AUTHOR - Aureliano Bombarely
       based in Bio::Tools::Run::Phylo::Phylip modules created by Shawn Hoon

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut


	
package Bio::Tools::Run::Phylo::Phylip::Dnaml;

use vars qw($AUTOLOAD @ISA $PROGRAM $PROGRAMDIR $PROGRAMNAME
	    @DNAML_PARAMS @OTHER_SWITCHES
	    %OK_FIELD);

use strict;
use warnings;
use autodie;

use Cwd;

use Bio::SimpleAlign;
use Bio::AlignIO;
use Bio::TreeIO;
use Bio::Tools::Run::Phylo::Phylip::Base;


# inherit from Phylip::Base which has some methods for dealing with
# Phylip specifics

@ISA = qw(Bio::Tools::Run::Phylo::Phylip::Base);

# You will need to enable the protdist program. This
# can be done in (at least) 3 ways:
#
# 1. define an environmental variable PHYLIPDIR:
# export PHYLIPDIR=/home/shawnh/PHYLIP/bin
#
# 2. include a definition of an environmental variable CLUSTALDIR in
# every script that will use Clustal.pm.
# $ENV{PHYLIPDIR} = '/home/shawnh/PHYLIP/bin';
#
# 3. You can set the path to the program through doing:
# my @params('program'=>'/usr/local/bin/protdist');
# my $protdist_factory = Bio::Tools::Run::Phylo::Phylip::ProtDist->new(@params);
# 


BEGIN {
	@DNAML_PARAMS = ('BEST_TREE', 'TRANS_RATIO', 'BASEFREQ', 
			 'ONECATSITES', 'RATESITES',  'SITESWEIGHT',
			 'SPEEDIER', 'GLOBALREAG', 'JUMBLE',
			 'OUTGROUP', 'MULTIPLE'
	    );
	@OTHER_SWITCHES = qw(QUIET);
	foreach my $attr(@DNAML_PARAMS,@OTHER_SWITCHES) {
		$OK_FIELD{$attr}++;
	}
}

=head2 program_name

 Title   : program_name
 Usage   : Bio::Tools::Run::Phylo::Phylip::Dnaml::program_name()
 Function: holds the program name
 Returns : string with the program name (dnaml)
 Args    : None

=cut

sub program_name {
  return 'dnaml';
}

=head2 program_dir

 Title   : program_dir
 Usage   : Bio::Tools::Run::Phylo::Phylip::Dnaml::program_dir()
 Function: returns the program directory, obtained from ENV variable.
 Returns : string with the program dir
 Args    : none

=cut

sub program_dir {
  return Bio::Root::IO->catfile($ENV{PHYLIPDIR}) if $ENV{PHYLIPDIR};
}

=head2 new

 Title   : new
 Usage   : Bio::Tools::Run::Phylo::Phylip::Dnaml->new(@args)
 Function: returns the program directory, obtained from ENV variable.
 Returns : A Bio::Tools::Run::Phylo::Phylip::Dnaml object
 Args    : @args, an array with the following arguments.
           quiet       => 1 or 0
           BEST_TREE   => Yes or $tree or $treefile_newick, 
           TRANS_RATIO => $int, an integer, 
           BASEFREQ    => $string, an string with 'Nt{space}{freq}\n' each nt 
	   ONECATSITES => $string, an string with site categ. like '1:2.5:2'
           RATESITES   => $strain, an string with the rate method and categ.  
           SITESWEIGHT => $strain, an string with character weights
	   SPEEDIER    => Yes or No, 
           GLOBALREAG  => No or $strain, with the reorganization 
           JUMBLE      => A array reference with [$seed_odd_num, $jumble_times]
	   OUTGROUP    => No or $int, an integer with the out.sequence position 
           MULTIPLE    => No or "Yes\n$number_alignments"

=cut

sub new {
    my ($class,@args) = @_;
    my $self = $class->SUPER::new(@args);    
    my ($attr, $value);
    while (@args)  {
	$attr =   shift @args;
	$value =  shift @args;
	next if( $attr =~ /^-/ ); # don't want named parameters
	if ($attr =~/PROGRAM/i) {
	    $self->executable($value);
	    next;
	}
	if ($attr =~ /IDLENGTH/i){
	    $self->idlength($value);
	    next;
	}
	$self->$attr($value);	
    }
    return $self;
}


sub AUTOLOAD {
    my $self = shift;
    my $attr = $AUTOLOAD;
    $attr =~ s/.*:://;
    $attr = uc $attr;
    $self->throw("Unallowed parameter: $attr !") unless $OK_FIELD{$attr};
    $self->{$attr} = shift if @_;
    return $self->{$attr};
}

=head2 idlength 

 Title   : idlength 
 Usage   : $obj->idlength ($newval)
 Function: 
 Returns : value of idlength 
 Args    : newvalue (optional)

=cut

sub idlength{
   my $self = shift;
   if( @_ ) {
      my $value = shift;
      $self->{'idlength'} = $value;
    }
    return $self->{'idlength'};

}


=head2  run 

 Title   : run 
 Usage   : $inputfilename = 't/data/prot.phy';
           $tree= $dnamlfactory->run($inputfilename);
           or
           $seq_array_ref = \@seq_array; @seq_array is array of Seq objs
           $aln = $factory->align($seq_array_ref);
           $tree = $dnamlfactory->run($aln);

 Function: Create a distance matrix from a SimpleAlign object or a 
           multiple alignment file
 Returns : L<Bio::Tree::Tree>
 Args    : Name of a file containing a multiple alignment in Phylip format
           or an SimpleAlign object 
 Side Efects: Throws an exception if argument is not either a string (eg a
              filename) or a Bio::SimpleAlign object. If
              argument is string, throws exception if file corresponding to 
              string name can not be found. 

=cut

sub run{

    my ($self, $input, $running_mode) = @_;
    my ($infilename);

    # Create input file pointer
    $infilename = $self->_setinput($input);
    
    if (defined $running_mode) {
	if ($running_mode =~ m/STOP_WRONGINPUT/i) {
	    if (!$infilename) {
		$self->thrown("Problems setting up for dnaml. Check $input");
	    }
	}
	elsif ($running_mode !~ m/WARN_WRONGINPUT/i) {
	    if (!$infilename) {
		$self->warn("Problems setting up for dnaml. Check $input");
	    }
	}
    }
    
    if (defined $infilename) {
	# Create parameter string to pass to protdist program
	my $param_string = $self->_setparams();
    
	# run protdist
	my @mat = $self->_run($infilename, $param_string);
	return  wantarray ? @mat:\@mat;
    }
}

#################################################

=head2  _run

 Title   :  _run
 Usage   :  Internal function, not to be called directly	
 Function:  makes actual system call to dnaml program
 Example :
 Returns : Bio::Tree object
 Args    : Name of a file containing a set of multiple alignments in 
           Phylip format and a parameter string to be passed to protdist

=cut

sub _run {
    my ($self, $infile, $param_string) = @_;
    my $instring;
    my $curpath = cwd; 
   
    unless( File::Spec->file_name_is_absolute($infile) ) {
	$infile = $self->io->catfile($curpath, $infile);
    }
    $instring =  $infile."\n$param_string";
    
    $self->debug( "Program ".$self->executable." $instring\n");
    
    chdir($self->tempdir);

    #open a pipe to run dnaml to bypass interactive menus
    my $exe = $self->executable();
    
    if ($self->quiet() || $self->verbose() < 0) {
	$exe .= " > /dev/null";
    }

    my $fail = 0;
    open(my $DNAML, "|" . $exe);
    print $DNAML $instring;
    
    close($DNAML);    
    chdir($curpath);
   
    # get the results
    my $outfile = $self->io->catfile($self->tempdir, $self->outfile);
    my $treefile = $self->io->catfile($self->tempdir, $self->treefile);
    
    $self->throw("dnaml did not create tree correctly (expected $treefile)")
	unless (-e $treefile);

    # Get the tree
    my $treeio = Bio::TreeIO->new(-file => $treefile, -format  => 'newick');
    my @tree;
    while (my $tree = $treeio->next_tree){
        push @tree, $tree;
    }

    # Clean up the temporary files created along the way...
    unless ( $self->save_tempfiles ) { 
        unlink $outfile;
        unlink $treefile;
    }
    return @tree; 
}


=head2  _setinput()

 Title   : _setinput
 Usage   : Internal function, not to be called directly	
 Function: Create input file for dnaml program
 Example :
 Returns : name of file containing a multiple alignment in Phylip format 
 Args    : SimpleAlign object reference or input file name


=cut

sub _setinput {
    my ($self, $input) = @_;
    my ($alnfilename,$tfh);

    # suffix is used to distinguish alignment files from an align object
    # If $input is not a reference it better be the name of a file with the 
    ## sequence

    # a phy formatted alignment file 
    unless (ref $input) {
        
	# check that file exists or throw
        $alnfilename= $input;
        unless (-e $input) {
	    return 0;
	}
	return $alnfilename;
    }

    #  $input may be a SimpleAlign Object
    if ($input->isa("Bio::Align::AlignI")) {

	## First check how many sequences has the alignment
	## With 2 or less it can give problems

	my $members_num = $input->num_sequences();
	
	if ($members_num > 2) { 

	    my ($tfh, $filename) = $self->io->tempfile( -dir => $self->tempdir);
	
	    my $alnIO = Bio::AlignIO->new(-fh       => $tfh, 
					  -format   => 'phylip',
					  -idlength => $self->idlength());
	    $alnIO->write_aln($input);
	    close($tfh);
	    
	    $alnfilename = $filename;
	    return $alnfilename;
	}
	else {
	    return 0;
	}
    }

    ## Get the names of the alignments to be used as outgroup

    if (defined $alnfilename) {
	my %names;

	my $alignio = Bio::AlignIO->new( -file   => $alnfilename,
					 -format => 'phylip',
	    );
	foreach my $align ($alignio->next_aln()) {
	    my $i = 1;
	    foreach my $seq ($align->each_seq() ) {
		my $seqid = $seq->id();
		$names{$seqid} = $i;
		$i++;
	    }
	}
	$self->names(\%names);
    }
}

sub _input_nbr {
    my ($self,$val) = @_;
    if($val){
        $self->{'_input_nbr'} = $val;
    }
    return $self->{'_input_nbr'};
}

=head2  _setparams()

 Title   : _setparams
 Usage   : Internal function, not to be called directly	
 Function: Create parameter inputs for protdist program
 Example :
 Returns : parameter string to be passed to protdist
 Args    : name of calling object

=cut

sub _setparams {
    my ($attr, $value, $self);

    #do nothing for now
    $self = shift;
    my $param_string = "";
    my $cat = 0;

    ## This menu arguments are not in Bio::Tools::Run::Phylo::Phylip::PhylipConf
    ## It will create them here.

    my %menu = ('BEST_TREE'   => "U\n", 
		'TRANS_RATIO' => "T\n", 
		'BASEFREQ'    => "F\n", 
		'ONECATSITES' => "C\n", 
		'RATESITES'   => "R\n",
		'SITESWEIGHT' => "W\n",
		'SPEEDIER'    => "S\n", 
		'GLOBALREAG'  => "G\n", 
		'JUMBLE'      => "J\n",
		'OUTGROUP'    => "O\n",
		'MULTIPLE'    => "M\n",
		'SUBMIT'      => "Y\n",
	);

    ## Now it will create the string to run the program

    foreach  my $attr (@DNAML_PARAMS) {
    	$value = $self->$attr();
	
    	next unless (defined $value);
    	
	if ($attr =~/BEST_TREE/i){
	    if ($value =~ /Yes/i){
		$param_string .= $menu{'BEST_TREE'} . "$value\n";
	    }
	    else {  ## It will use or a file or a Bio::Tree::Tree object
		
		my $treenewick;
		if (ref($value) ne 'Bio::Tree::Tree') { ## It will be a file
		    $param_string .= $menu{'BEST_TREE'} . "$value\n";
		}
		else {  ## it will create a temp file
		    my ($tfh, $treefile) = $self->io
                                                ->tempfile( 
			                             -dir => $self->tempdir );
    
		    my $treeIO = Bio::TreeIO->new( -fh       => $tfh, 
			  			   -format   => 'newick',
						 );

		    $treeIO->write_tree($value);
		    $treeIO->close();
		    $param_string .= $menu{'BEST_TREE'} . "$treefile\n";
		}
	    }
	}
	if ($attr =~/TRANS_RATIO/i){
	    if ($value =~ /\d+/i){
		$param_string .= $menu{'TRANS_RATIO'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for TRANS_RATIO (digit)");
	    }
	}
	if ($attr =~/BASEFREQ/i){
	    if ($value =~ /\w+/i){
		$param_string .= $menu{'BASEFREQ'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for BASEFREQ");
	    }
	}
	if ($attr =~/ONECATSITES/i){
	    if ($value =~ /\w+/i){
		$param_string .= $menu{'ONECATSITES'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for ONECATSITES");
	    }
	}
	if ($attr =~/RATESITES/i){
	    if ($value =~ /constant rate/i){  ## For now it will use only that
		$param_string .= $menu{'RATESITES'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for RATESITES (constant rate)");
	    }
	}
	if ($attr =~/SITESWEIGHT/i){
	    if ($value =~ /Yes|No/i){
		$param_string .= $menu{'SITESWEIGHT'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for SITESWEIGHT (Yes|No)");
	    }
	}
	if ($attr =~/SPEEDIER/i){
	    if ($value =~ /Yes|No/i){
		$param_string .= $menu{'SPEEDIER'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for SPEEDIER (Yes|No)");
	    }
	}
	if ($attr =~/GLOBALREAG/i){
	    if ($value =~ /Yes|No/i){
		$param_string .= $menu{'GLOBALREAG'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for GLOBALREAG (Yes|No)");
	    }
	}
	if ($attr =~ /JUMBLE/i){
	    if (ref($value) ne 'ARRAY') {
		$self->throw("Unallowed value for random seed, need odd number")
		    unless ($value =~ /\d+/ && ($value % 2 == 1));
		if ($value =~ m/\d+\\n\d+/) {
		    $param_string .= $menu{'JUMBLE'}."$value\n";
		}
		else {
		    $param_string .= $menu{'JUMBLE'}."$value\n1\n";
		}
	    }
	    else {
		my @jumvals = @{$value};
		foreach my $jumval (@jumvals) {
		    unless ($jumval =~ m/^\d+$/ && ($jumval % 2 == 1)) {
			$self->throw("Unallowed random seed/iter. Only odd n.");
		    }
		}
		$param_string .= $menu{'JUMBLE'}."$jumvals[0]\n$jumvals[1]\n";
	    }
	}

	if ($attr =~/RANDOMINPUT/i){
	    if ($value =~ /Yes|No/i){
		$param_string .= $menu{'RANDOMINPUT'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for RANDOMINPUT (Yes|No)");
	    }
	}
	if ($attr =~/OUTGROUP/i){
	    if($value !~/^\d+$/){ # is a name so find the rank 
		my %names = %{$self->_alignseq_names};
		$names{$value} || $self->throw("Outgroup $value not found");
		$value = $names{$value};
		$param_string .= $menu{'OUTGROUP'} . "$value\n";
	    }
	    else {
		$param_string .= $menu{'OUTGROUP'} . "$value\n";
	    }
	}
	if ($attr =~/MULTIPLE/i){
	    if ($value =~ /Yes|No/i){
		$param_string .= $menu{'MULTIPLE'} . "$value\n";
	    }
	    else {
		$self->throw("Unallowed value for MULTIPLE (Yes|No)");
	    }
	}
    }

    #set multiple option is not set and there are more than one sequence
    
    if ( ($param_string !~ $menu{'MULTIPLE'}) && 
	 (defined ($self->_input_nbr)         &&
	  ($self->_input_nbr > 1))) 
    {
	$param_string.=$menu{'MULTIPLE'}.$self->_input_nbr."\n";
    }
    $param_string .= $menu{'SUBMIT'};
    return $param_string;
}

=head2  names()

 Title   :  names
 Usage   :  $tree->names(\%names)
 Function:  get/set for a hash ref for storing names in matrix
            with rank as values.
 Example :
 Returns : hash reference 
 Args    : hash reference 

=cut

sub names {
    my ($self,$name) = @_;
    if($name){
        $self->{'_names'} = $name;
    }
    return $self->{'_names'};
}

=head1 Bio::Tools::Run::Wrapper methods

=cut

=head2 outfile

 Title   : outfile
 Usage   : $obj->outfile($newval)
 Function: Get/Set default PHYLIP outfile name ('outfile' usually)
 Returns : value of outfile
 Args    : newvalue (optional)


=cut


=head2 treefile

 Title   : treefile
 Usage   : $obj->treefile($newval)
 Function: Get/Set the default PHYLIP treefile name ('treefile' usually)
 Returns : value of treefile
 Args    : newvalue (optional)


=cut

=head2 no_param_checks

 Title   : no_param_checks
 Usage   : $obj->no_param_checks($newval)
 Function: Boolean flag as to whether or not we should
           trust the sanity checks for parameter values  
 Returns : value of no_param_checks
 Args    : newvalue (optional)


=cut

=head2 save_tempfiles

 Title   : save_tempfiles
 Usage   : $obj->save_tempfiles($newval)
 Function: 
 Returns : value of save_tempfiles
 Args    : newvalue (optional)


=cut

=head2 outfile_name

 Title   : outfile_name
 Usage   : my $outfile = $protdist->outfile_name();
 Function: Get/Set the name of the output file for this run
           (if you wanted to do something special)
 Returns : string
 Args    : [optional] string to set value to


=cut


=head2 tempdir

 Title   : tempdir
 Usage   : my $tmpdir = $self->tempdir();
 Function: Retrieve a temporary directory name (which is created)
 Returns : string which is the name of the temporary directory
 Args    : none


=cut

=head2 cleanup

 Title   : cleanup
 Usage   : $codeml->cleanup();
 Function: Will cleanup the tempdir directory after a ProtDist run
 Returns : none
 Args    : none


=cut

=head2 io

 Title   : io
 Usage   : $obj->io($newval)
 Function:  Gets a L<Bio::Root::IO> object
 Returns : L<Bio::Root::IO>
 Args    : none


=cut

1; # Needed to keep compiler happy
