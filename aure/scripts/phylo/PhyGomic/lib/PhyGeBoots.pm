
package PhyGeBoots;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Math::BigFloat;

use Bio::Seq::Meta;
use Bio::LocatableSeq;
use Bio::Cluster::SequenceFamily;

use Bio::AlignIO;
use Bio::Align::DNAStatistics;

use Bio::Tools::Run::Phylo::Phylip::SeqBoot;
use Bio::Tools::Run::Phylo::Phylip::Neighbor;
use Bio::Tools::Run::Phylo::Phyml;
use Bio::Tools::Run::Phylo::Phylip::Consense;

use Bio::Matrix::IO;
use Bio::Matrix::Generic;

###############
### PERLDOC ###
###############

=head1 NAME

PhyGeBoots.pm
A class to bootstrap data.

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use PhyGeBoots;

  ## To get each cluster from a blast file

  my $phygeboots = PhyGeBoots->new({ seqfam        => $seqfam_obj, 
                                     run_bootstrap => $args_href1,
                                     run_distances => $args_href2,
                                     run_njtrees   => $args_href3,
                                     run_consensus => $args_href4,
                                   });


  ## Accessors

  my $seqfam = $phygeboots->get_seqfam();
  $phygeboots->set_seqfam($seqfam);

  my @bootsaligns = @{$phygeboots->get_aligns()};
  $phygeboots->set_aligns(\@bootsaligns);

  my @bootsdists = @{$phygeboots->get_dists()};
  $phygeboots->set_dists(\@bootsdists);

  my @bootstrees = @{$phygeboots->get_trees()};
  $phygeboots->set_trees(\@bootstrees);

  my $consensus = $phygeboots->get_consensus();
  $phygeboots->set_consensus($consensus);


  ## To run different tools over the phygeboots object.

  my @aligns = $phygeboots->run_bootstrap($args_href);
  my @dists = $phygeboots->run_distances($args_href);
  my @trees = $phygeboots->run_njtrees($args_href);
  my @trees = $phygeboots->run_mltrees($args_href);
  my $consensus = $phygeboots->run_consensus($args_href);


=head1 DESCRIPTION

 Object to manipulate bootstrapping of a Bio::Cluster::SequenceFamily object


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $phygeboots = PhyGeBoots->new($arguments_href);

  Desc: Create a phygeboots object and run some of the functions if they are
        specified.

  Ret: a PhyGeBoots.pm object

  Args: A hash reference with the following key-value pairs: 
         + seqfam => a Bio::Cluster::SequenceFamily object, 
         + run_bootstrap => a hash ref. with args. (see run_bootstrap function)
         + run_distances => a hash ref. with args. (see run_distances function)
         + run_trees     => a hash ref. with args. (see run_trees function)
         + run_consensus => a hash ref. with args. (see run_consensus function)

  Side_Effects: Die if the argument used is not a hash or there are argument
                incompatibility (for example run_trees without run_distances).

  Example: my $phygeboots = PhyGeBoots->new();
           my $phygeboots = PhyGeBoots->new({ seqfam => $seqfam });

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my $seqfam = '';
    my @aligns = ();
    my @dists = ();
    my @trees = ();
    my $consensus = '';

    ## Check argument compatibility

    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ARGUMENT ERROR: $args_href used for new() is not HASH REF.");
	}
    }

    my $err = "ARGUMENT INCOMPATIBILITY for new() function: ";
    unless (defined $args_href->{'seqfam'}) {

	if (defined $args_href->{'run_bootstrap'}) {
	    	$err .= "run_bootstrap arg. can not be used without seqfam";
		croak($err);
	}
    }
    else {
	$seqfam = $args_href->{'seqfam'};
    }
    unless (defined $args_href->{'run_bootstrap'}) {
	if (defined $args_href->{'run_distances'}) {
	    $err .= "run_distances arg. can not be used without ";
	    $err .= "run_bootstrap";
	    croak($err);	    
	}
    }
    unless (defined $args_href->{'run_distances'}) {
	if (defined $args_href->{'run_trees'}) {
	    $err .= "run_trees arg. can not be used without ";
	    $err .= "run_distances";
	    croak($err);	    
	}
    }
    unless (defined $args_href->{'run_trees'}) {
	if (defined $args_href->{'run_consensus'}) {
	    $err .= "run_consensus arg. can not be used without ";
	    $err .= "run_trees";
	    croak($err);	    
	}
    }

    ## Run the different options sequencially

    if (defined $seqfam) {
	$self->set_seqfam($seqfam);

	if (defined $args_href->{'run_bootstrap'}) {
	    @aligns = $self->run_bootstrap($args_href->{'run_bootstrap'});
	    $self->set_aligns(\@aligns);
	
	    if (defined $args_href->{'run_distances'}) {
		@dists = $self->run_distances($args_href->{'run_distances'});
		$self->set_dists(\@dists);
	    
		if (defined $args_href->{'run_trees'}) {
		    @trees = $self->run_trees($args_href->{'run_trees'});
		    $self->set_trees(\@trees);
		
		    if (defined $args_href->{'run_consensus'}) {
			my $cons_href = $args_href->{'run_consensus'};
			$consensus = $self->run_consensus($cons_href);
			$self->set_consensus($consensus);		    
		    }
		}
	    } 
	}
    }
    else {  ## create an empty object.
	
	$self->set_seqfam($seqfam);
	$self->set_aligns(\@aligns);
	$self->set_dists(\@dists);
	$self->set_trees(\@trees);
	$self->set_consensus($consensus);
    }
    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get_seqfam

  Usage: my $seqfam = $phygeboots->get_seqfam(); 

  Desc: Get the Bio::Cluster::SequenceFamily object from PhyGeBoots object.

  Ret: a  Bio::Cluster::SequenceFamily object

  Args: None

  Side_Effects: None

  Example: my $seqfam = $phygeboots->get_seqfam(); 

=cut

sub get_seqfam {
    my $self = shift;

    return $self->{seqfam};
}

=head2 set_seqfam

  Usage: $phygeboots->set_seqfam($seqfam);

  Desc: Get the Bio::Cluster::SequenceFamily object from PhyGeBoots object.

  Ret: None

  Args: a Bio::Cluster::SequenceFamily object

  Side_Effects: Die if no argument is used or if the object is not a 
                Bio::Cluster::SequenceFamily object.

  Example: $phygeboots->set_seqfam($seqfam);

=cut

sub set_seqfam {
    my $self = shift;
    my $seqfam = shift;
   
    unless (defined $seqfam) {
	croak("ARG. ERROR: No arg. was used for set_seqfam function");
    }
    else {
	if ($seqfam =~ m/\w+/) {
	    my $ob = 'Bio::Cluster::SequenceFamily';
	    unless (ref($seqfam) eq $ob) {
		croak("ARG. ERROR: Arg=$seqfam used for set_seqfam() isnt $ob");
	    } 
	}
	$self->{seqfam} = $seqfam;    
    }
}


=head2 get_aligns

  Usage: my @aligns = $self->get_aligns();

  Desc: Get an array of Bio::SimpleAlign object from PhyGeBoots object

  Ret: An array of Bio::SimpleAlign objects

  Args: None

  Side_Effects: None

  Example: my @aligns = $self->get_aligns();

=cut

sub get_aligns {
    my $self = shift;
    return @{$self->{aligns}};
}

=head2 set_aligns

  Usage: $phygeboots->set_aligns(\@aligns);

  Desc: Set an array of Bio::SimpleAlign object in the PhyGeBoots object

  Ret: None

  Args: An array of Bio::SimpleAlign objects.

  Side_Effects: Die if no argument is specified, if this argument is not an
                array reference or if at least one of its elements is not a 
                Bio::Cluster::SimpleAlign object

  Example: $phygeboots->set_aligns(\@aligns);

=cut

sub set_aligns {
    my $self = shift;
    my $aref = shift ||
	croak("ARG. ERROR: No arg. was supplied to set_aligns function");

    my $obj = 'Bio::SimpleAlign';
    unless (ref($aref) eq 'ARRAY') {
	croak("ARG. ERROR: Arg=$aref supplied to set_aligns isn't a ARRAYREF");
    }
    else {
	foreach my $elobj (@{$aref}) {
	    unless (ref($elobj) eq $obj) {
		croak("ARG. ERROR: Element=$elobj for set_aligns isnt a $obj");
	    }
	}
    }
    $self->{aligns} = $aref;
}

=head2 get_dists

  Usage: my @dists = $self->get_dists();

  Desc: Get an array of Bio::Matrix::PhylipDist object from PhyGeBoots object

  Ret: An array of Bio::Matrix::PhylipDist objects

  Args: None

  Side_Effects: None

  Example: my @dists = $self->get_dists();

=cut

sub get_dists {
    my $self = shift;
    return @{$self->{dists}};
}

=head2 set_dists

  Usage: $phygeboots->set_dists(\@dists);

  Desc: Set an array of Bio::Matrix::PhylipDist object in the PhyGeBoots object

  Ret: None

  Args: An array of Bio::Matrix::PhylipDist objects.

  Side_Effects: Die if no argument is specified, if this argument is not an
                array reference or if at least one of its elements is not a 
                Bio::Matrix::PhylipDist object

  Example: $phygeboots->set_dists(\@dists);

=cut

sub set_dists {
    my $self = shift;
    my $aref = shift ||
	croak("ARG. ERROR: No arg. was supplied to set_dists function");

    my $obj = 'Bio::Matrix';
    unless (ref($aref) eq 'ARRAY') {
	croak("ARG. ERROR: Arg=$aref supplied to set_dists isn't a ARRAYREF");
    }
    else {
	foreach my $elobj (@{$aref}) {
	    unless (ref($elobj) =~ m/$obj/) {
		croak("ARG. ERROR: Element=$elobj for set_dists isnt a $obj");
	    }
	}
    }
    $self->{dists} = $aref;
}


=head2 get_trees

  Usage: my @trees = $self->get_trees();

  Desc: Get an array of Bio::Tree object from PhyGeBoots object

  Ret: An array of Bio::Tree objects

  Args: None

  Side_Effects: None

  Example: my @trees = $self->get_trees();

=cut

sub get_trees {
    my $self = shift;
    return @{$self->{trees}};
}

=head2 set_trees

  Usage: $phygeboots->set_trees(\@trees);

  Desc: Set an array of Bio::Tree object in the PhyGeBoots object

  Ret: None

  Args: An array of Bio::Tree objects.

  Side_Effects: Die if no argument is specified, if this argument is not an
                array reference or if at least one of its elements is not a 
                Bio::Tree object

  Example: $phygeboots->set_trees(\@trees);

=cut

sub set_trees {
    my $self = shift;
    my $aref = shift ||
	croak("ARG. ERROR: No arg. was supplied to set_trees function");

    my $obj = 'Bio::Tree::Tree';
    unless (ref($aref) eq 'ARRAY') {
	croak("ARG. ERROR: Arg=$aref supplied to set_trees isn't a ARRAYREF");
    }
    else {
	foreach my $elobj (@{$aref}) {
	    unless (ref($elobj) eq $obj) {
		croak("ARG. ERROR: Element=$elobj for set_trees isnt a $obj");
	    }
	}
    }
    $self->{trees} = $aref;
}



=head2 get_consensus

  Usage: my $consensus = phygeboots->get_consensus()

  Desc: Get the consensus tree, a Bio::Tree::Tree object from phyGeBoots object

  Ret: A Bio::Tree::Tree object

  Args: None

  Side_Effects: None

  Example: my $consensus = phygeboots->get_consensus()

=cut

sub get_consensus {
    my $self = shift;
    return $self->{consensus};
}

=head2 set_consensus

  Usage: $phygeboots->set_consensus($consensus);

  Desc: Set the consensus tree for a PheGeBoots object

  Ret: None

  Args: A Bio::Tree::Tree object

  Side_Effects: Die if no argument is supplied to the function or if the arg.
                is not a Bio::Tree::Tree object

  Example: $phygeboots->set_consensus($consensus);

=cut

sub set_consensus {
    my $self = shift;
    my $cons = shift ||
	croak("ARG. ERROR: No arg was supplied to set_consensus function");

    my $obj = 'Bio::Tree::Tree';
    unless (ref($cons) eq $obj) {
	croak("ARG. ERROR: Arg=$cons supplied to set_consensus isnt $obj");
    }
    $self->{consensus} = $cons;
}


##########################
## ANALITICAL FUNCTIONS ##
##########################

=head2 run_bootstrap

  Usage: my @aligns = $phygeboots->run_bootstrap($arg_href);

  Desc: Calculates the bootstrapping for a the alignment contained in
        the Bio::Cluster::SequenceFamily object

  Ret: An array of Bio::SimpleAlign objects (bootstrap)

  Args: $args_href, a hash reference with the following options:
        datatype, permute, blocksize, replicates, readweights, readcat & quiet.
        (for more information: 
         http://evolution.genetics.washington.edu/phylip/doc/seqboot.html)
 

  Side_Effects: Died if the arguments are wrong.
                Return undef in there are no Bio::Cluster::SequenceFamily or
                Alignment inside the object
                Set alignments inside PhyGeBoots object
                
  Example: my @aligns = $phygeboots->run_bootstrap(
                                                    { 
                                                      datatype   => 'Sequence',
                                                      replicates => 1000,
                                                    }
                                                  );

=cut

sub run_bootstrap {
    my $self = shift;
    my $args_href = shift ||
	croak("ARG. ERROR: No args. were supplied to run_bootstrap()");
    
    ## Bootstrap tool uses the arguments as a array, so the argument checking
    ## also will change them to this format

    my %perm_args = (
	'datatype' => {
	    'Molecular sequences'               => 1, 
	    'Discrete morphological characters' => 1, 
	    'Restriction sites'                 => 1, 
	    'Gene frequencies'                  => 1,
	    'Sequence'                          => 1, 
	    'Morph'                             => 1,
	    'Rest.'                             => 1, 
            'Gene Freqs'                        => 1,
	}, 
	'permute'     => 'yes|no',
	'blocksize'   => 'int', 
	'replicates'  => 'int', 
	'readweights' => 'yes|no', 
	'readcat'     => 'yes|no', 
	'quiet'       => 'yes|no',
	);

    ## Checking args.

    my @bootstrargs = ();
    unless (ref($args_href) eq 'HASH') {
	croak("ARG. ERROR: Arg. supplied to run_bootstrapping() arent hashref");
    }
    else {
	foreach my $argkey (keys %{$args_href}) {
	    my $val = $args_href->{$argkey};
	    unless (exists $perm_args{$argkey}) {
		my $err = "ARG. ERROR: $argkey is not a permited arg. for ";
		$err .= "run_bootstrapping() function";
		croak($err);
	    }
	    else {
		if (ref($perm_args{$argkey}) eq 'HASH') {
		    my %perm_dt = %{$perm_args{$argkey}};
		    
		    unless (exists $perm_dt{$val}) {
			my $err1 = "ARG. ERROR: $val is not a permited value ";
			$err1 .= "for arg=$argkey with run_bootstrapping()";
			croak($err1);
		    }
		}
		else {
		    if ($perm_args{$argkey} eq 'int') {
			unless ($val =~ m/^\d+$/) {
			    my $err2 = "ARG. ERROR: $val for $argkey is not ";
			    $err2 .= "an integer for run_bootstrapping()";
			    croak($err2);
			}
		    }
		    else {
			unless ($val =~ m/^(yes|no)$/) {
			    my $err3 = "ARG. ERROR: $val for $argkey do not ";
			    $err3 .= "have yes|no val. for run_bootstrapping()";
			    croak($err3);
			}
		    }
		}
		push @bootstrargs, ($argkey, $val);
	    }
	}
    }

    ## Define the variable to catch the align objects

    my @boots_aligns = ();

    ## Now it will get all the run_bootstrapping factory

    my $fact_boots = Bio::Tools::Run::Phylo::Phylip::SeqBoot->new(@bootstrargs);

    ## Get the Bio::Cluster::SequenceFamily object and the alignment from it

    my $seqfam = $self->get_seqfam();
    if (defined $seqfam) {

	my $align = $seqfam->alignment();
	if (defined $align) {
	    
	    @boots_aligns = @{$fact_boots->run($align)};
	}
    }

    ## Finally it will return the array and also it will set the aligns
    ## in the object
    
    $self->set_aligns(\@boots_aligns);
    return @boots_aligns;
}


=head2 run_distances

  Usage: my @dists = $phygeboots->run_distances($arg_href);

  Desc: Calculates the distances for a set of alignment contained in
        PhyGeBoots object as bootstraps

  Ret: An array of Bio::Matrix::Phylip objects

  Args: $args_href, a hash reference with the following options:
        method => a scalar to pass to the Bio::Align::DNAStatistics->distance()
        function (JukesCantor,Uncorrected,F81,Kimura,Tamura or TajimaNei).
        (for more information: 
         http://evolution.genetics.washington.edu/phylip/doc/seqboot.html)
 

  Side_Effects: Died if the arguments are wrong.
                Return undef if there are no aligns inside the PhyGeBoots obj.
                Set alignments inside PhyGeBoots object.
                It will run Kimura by default.
                
  Example: my @dists = $phygeboots->run_distances({method => 'JukesCantor'});

=cut

sub run_distances {
    my $self = shift;
    my $args_href = shift || { method => 'Kimura' }; ## Use Kimura by default.

    my $method = '';
    unless (ref($args_href) eq 'HASH') {
	croak("ARG. ERROR: Arg. supplied to run_distances() arent hashref");
    }
    else {
	$method = $args_href->{method} ||
	    croak("ARG. ERROR: No method arg. was used for run_distances()");

	my %perm_methods = ( 
	    'JukesCantor' => 1,
	    'Uncorrected' => 1,
	    'F81'         => 1,
	    'Kimura'      => 1,
	    'Tamura'      => 1,
	    'TajimaNei'   => 1 );
	
	unless (exists $perm_methods{$method}) {
	    croak("ARG. ERROR: $method isnt permited method for run_distances");
	}
    }

    my @dists = ();

    ## Create the factory to calculate the distances

    my $factory = Bio::Align::DNAStatistics->new();

    my @aligns = $self->get_aligns();
    foreach my $align (@aligns) {

	my $dists = $factory->distance( -align => $align, -method => $method );
	push @dists, $dists;
    }
    
    $self->set_dists(\@dists);
    return @dists;
}


=head2 run_njtrees

  Usage: my @trees = $phygeboots->run_njtrees($arg_href);

  Desc: Calculates the trees for a set of distances contained in
        PhyGeBoots object as bootstraps

  Ret: An array of Bio::Tree::Tree objects

  Args: $args_href, a hash reference with the following options:
        type => NJ|UPGMA, outgroup => int, lowtri => 1|0, uptri => 1|0, 
        subrep => 1|0 or/and jumble => int
        (for more information: 
         http://evolution.genetics.washington.edu/phylip/doc/seqboot.html)
 

  Side_Effects: Died if the arguments are wrong.
                Return undef if there are no dists inside the PhyGeBoots obj.
                Set trees inside PhyGeBoots object.
                
  Example: my @trees = $phygeboots->run_njtrees($arg_href);

=cut

sub run_njtrees {
    my $self = shift;
    my $args_href = shift;

    my %perm_args = ( type     => '(NJ|UPGMA)', 
		      outgroup => '\d+', 
		      lowtri   => '[1|0]', 
		      uptri    => '[1|0]', 
		      subrep   => '[1|0]', 
		      jumble   => '\d+', 
		      quiet    => '[1|0]',
	);

    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ARG. ERROR: Arg. supplied to run_njtrees() arent hashref");
	}

	foreach my $arg (keys %{$args_href}) {
	    unless (exists $perm_args{$arg}) {
		croak("ARG. ERROR: $arg isnt permited arg. for run_njtrees()");
	    }
	    else {
		my $val = $args_href->{$arg};
		my $regexp = '^' . $perm_args{$arg} . '$';
		if ($val !~ m/$regexp/i) {
		    croak("ARG. ERROR: $arg has a non-permited value ($val)")
		}
	    }
	}
    }

    my @args = ();
    foreach my $key (keys %{$args_href}) {
	push @args, ($key, $args_href->{$key});
    }

    ## After check, create the array and the factory

    my @trees = ();
    
    my $factory = Bio::Tools::Run::Phylo::Phylip::Neighbor->new(@args);
    
    ## And now it will run one tree per distance

    my @dists = $self->get_dists();
    foreach my $dist (@dists) {

	my ($tree) = $factory->run($dist);
	push @trees, $tree;
    }
    
    $self->set_trees(\@trees);
    return @trees;
}


=head2 run_mltrees

  Usage: my @trees = $phygeboots->run_mltrees($arg_href);

  Desc: Calculates the trees for a set of distances contained in
        PhyGeBoots object as bootstraps using codeml program

  Ret: An array of Bio::Tree::Tree objects

  Args: $args_href, a hash reference with the following options:
        phyml_arg => {}, $hash ref with phyml arguments
        (for more information: 
         http://www.atgc-montpellier.fr/phyml/usersguide.php)
 

  Side_Effects: Died if the arguments are wrong.
                Return undef if there are no dists inside the PhyGeBoots obj.
                Set trees inside PhyGeBoots object.
                
  Example: 

=cut

sub run_mltrees {
    my $self = shift;
    my $args_href = shift;

    my %perm_args = ( 
	phyml_arg      => 'HASH',
	);

    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ARG. ERROR: Arg. supplied to run_mltrees() arent hashref");
	}

	foreach my $arg (keys %{$args_href}) {
	    unless (exists $perm_args{$arg}) {
		croak("ARG. ERROR: $arg isnt permited arg. for run_mltrees()");
	    }
	    else {
		my $val = $args_href->{$arg};
		my $regexp = '^' . $perm_args{$arg};
		if ($val !~ m/$regexp/i) {
		    croak("ARG. ERROR: $arg has a non-permited value ($val)")
		}
	    }
	}
    }

    ## Default parameters:

    my %phyargs = ();
    if (defined $args_href->{phyml_args}) {
	%phyargs = %{$args_href->{phyml_args}};
    }
    else {
	%phyargs = (
	    -data_type       => 'nt',
	    );
    }

    ## After check, create the array and the factory

    my @trees = ();

    ## And now it will run one tree per alignment
    ## Factory only can deal with one tree per run... so it is necessary
    ## redo de factory as many times as trees to create.

    my @aligns = $self->get_aligns();
    foreach my $align (@aligns) {
	my $factory = Bio::Tools::Run::Phylo::Phyml->new(%phyargs);
	my $tree = $factory->run($align);
	push @trees, $tree;
    }
    
    $self->set_trees(\@trees);
    return @trees;
}

=head2 run_consensus

  Usage: my $consensus = $phygeboots->run_consensus($arg_href);

  Desc: Calculates the consensus for an array of trees from the PhyGeBoots 
        object

  Ret: One Bio::Tree::Tree object

  Args: $args_href, a hash reference with the following options:
        
 

  Side_Effects: Died if the arguments are wrong.
                Return undef if there are no tree inside the PhyGeBoots obj.
                Set trees inside PhyGeBoots object.
                
  Example: 

=cut

sub run_consensus {
    my $self = shift;
    my $args_href = shift;

     my %perm_args = (
	 type     => '\w+',
	 rooted   => '\d+',
	 outgroup => '\d+', 
	 quiet    => '[1|0]',
	 );

    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ARG. ERROR: Arg. supplied to run_consensus() arent hashref");
	}

	foreach my $arg (keys %{$args_href}) {
	    unless (exists $perm_args{$arg}) {
		croak("ARG. ERROR: $arg isnt permited arg.for run_consensus()");
	    }
	    else {
		my $val = $args_href->{$arg};
		my $regexp = '^' . $perm_args{$arg} . '$';
		if ($val !~ m/$regexp/i) {
		    croak("ARG. ERROR: $arg has a non-permited value ($val)")
		}
	    }
	}
    }

    my @args = ();
    foreach my $key (keys %{$args_href}) {
	push @args, ($key, $args_href->{$key});
    }

    ## After check, create the array and the factory

    my $consensus = '';
    
    my $factory = Bio::Tools::Run::Phylo::Phylip::Consense->new(@args);
    
    ## And now it will run one tree per distance

    my @trees = $self->get_trees();

    if (scalar(@trees) > 0) {
	$consensus = $factory->run(\@trees);
    }

    $self->set_consensus($consensus);
    return $consensus;
}



####
1; #
####
