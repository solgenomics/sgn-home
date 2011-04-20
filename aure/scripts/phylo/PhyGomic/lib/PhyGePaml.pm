
package PhyGePaml;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Try::Tiny;
use Math::BigFloat;

use File::Temp qw/ tempfile tempdir/;

use FindBin;
use lib "$FindBin::Bin/../lib";


###############
### PERLDOC ###
###############

=head1 NAME

PhyGePaml.pm
A class to analyze syn/non-syn ratio using PAML

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use PhyGePaml;

  my $phygepaml = PhyGePaml->new({ seqfams     => $hash_ref,
                                   cds         => $hash_ref,
                                   strains     => $hash_ref,
                                   topotypes   => $hash_ref, 
                                   });


  ## Accessors

  my $seqfam_href = $phygepaml->get_seqfams();
  $phygepaml->set_seqfams($seqfam_href);

  my $cds_href = $phygepaml->get_cds();
  $phygepaml->set_cds($cds_href);

  my $strains_href = $phygepaml->get_strains();
  $phygepaml->set_strains($strains_href);

  my $topotypes_href = $phygepaml->get_topotypes();
  $phygepaml->set_topotypes($topotypes_href);

  my $mlmatrix_href = $phygepaml->get_mlmatrix();
  $phygepaml->set_mlmatrix($mlmatrix_href);

  ## Running tools

  $phygepaml->cds_alignment(\%args);
  $phygepaml->run(\%args);
  


=head1 DESCRIPTION

 Object to analyze dN/dS ratio for between different members of the phygecluster
 sequence family.

 It have different methods to create the cds (codon) alignment:

 1- Loading the cds for a consensus sequence for each alignment/sequence family.
    Align it, and trimming the sequence un the alignment.

 2- Calculating the consensus from the alignment and use estscan or getting
    the longest 6 frame for this consensus.


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $phygepaml = PhyGePaml->new($arguments_href);

  Desc: Create a phygepaml object with the specified parameters.

  Ret: a PhyGePaml.pm object

  Args: A hash reference with the following key-value pairs: 
         + seqfams     => a hash ref. with key=id and 
                          value=Bio::Cluster::SequenceFamily obj,
         + strains     => a hash ref. with key=members and value=strains
         + cds         => a hash ref. with key=id and value=consensus_cds (Seq)
         + topotypes   => a hash ref. with key=id and value=Bio::Tree::TopoType 
                          object
        
  Side_Effects: Die if the argument used is not a hash or there are argument
                incompatibility.

  Example: my $phygepaml = PhyGePaml->new($arguments_href);

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %permargs = ( 
	seqfams     => {},
	strains     => {},
	cds         => {},
	topotypes   => {},
	);

    ## Check argument

    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ARGUMENT ERROR: $args_href used for new() is not HASH REF.");
	}
	else {
	    my %args = %{$args_href};
	    foreach my $argkey (keys %args) {
		my $exp = $permargs{$argkey};
		unless (defined $exp) {
		    croak("ARG. ERROR: $argkey isnt a permited arg. for new()");
		}
		else {
		    if (ref($args{$argkey}) ne ref($exp)) {
			croak("ARG. ERROR: $args{$argkey} isnt $exp for new()");
		    }
		}	
	    }
	}
    }

    ## Fill the default values
	
    foreach my $argkey2 (keys %permargs) {
	unless (exists $args_href->{$argkey2}) {
	    $args_href->{$argkey2} = $permargs{$argkey2};
	}
    }

    ## Set vars in the object

    $self->set_seqfams($args_href->{seqfams});
    $self->set_strains($args_href->{strains});
    $self->set_cds($args_href->{cds});
    $self->set_topotypes($args_href->{topotypes});
    
    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get_seqfams/set_seqfams

  Usage: my $seqfams_href = $phygepaml->get_seqfams();
         $phygepaml->set_seqfam($seqfam_href); 

  Desc: Get or set a hash with key=ID and value=Bio::Cluster::SequenceFamily 
        objects

  Ret: Get: Hash reference with key=ID and value=Bio::Cluster::SequenceFamily.
       Set: None

  Args: Get: None
        Set: Hash reference with key=ID and value=Bio::Cluster::SequenceFamily

  Side_Effects: Die if no argument is used or if the object is not a hash ref.
                or if its values arent Bio::Cluster::SequenceFamily objects.

  Example: my %seqfams = %{$phygepaml->get_seqfams()}; 
           $phygepaml->set_seqfam(\%seqfams);

=cut

sub get_seqfams {
    my $self = shift;
    return $self->{seqfams};
}

sub set_seqfams {
    my $self = shift;
    my $seqfam_href = shift;
   
    unless (defined $seqfam_href) {
	croak("ARG. ERROR: No arg. was used for set_seqfams function");
    }
    else {
	if ($seqfam_href =~ m/\w+/) {
	    my $obj = 'Bio::Cluster::SequenceFamily';
	    unless (ref($seqfam_href) eq 'HASH') {
		croak("ARG. ERROR: $seqfam_href set_seqfams() isnt HASHREF");
	    }
	    else {
		my %hash = %{$seqfam_href};
		foreach my $key (keys %hash) {
		    my $val = $hash{$key};
		    if (ref($hash{$key}) ne $obj) {
			croak("ARG. ERROR: $val used set_seqfams isnt $obj");
		    }
		}
	    }
	}
	$self->{seqfams} = $seqfam_href;    
    }
}

=head2 get_strains/set_strains

  Usage: my $strains_href = $phygepaml->get_strains(); 
         $phygepaml->set_strains($strains_href);

  Desc: Get or set a hashref with the members and the strains

  Ret: Get: a hashref with keys=member_id and value=strain
       Set: None

  Args: Get: None
        Set: a hashref with keys=member_id and value=strain

  Side_Effects: Die if no argument is used with set or if this argument isnt a
                hashref.

  Example: my %strains = %{$phygepaml->get_strains()};
           $phygepaml->set_strains($strains_href);

=cut

sub get_strains {
    my $self = shift;
    return $self->{strains};
}

sub set_strains {
    my $self = shift;
    my $strain_href = shift;
   
    unless (defined $strain_href) {
	croak("ARG. ERROR: No arg. was used for set_strains function");
    }
    else {
	unless (ref($strain_href) eq 'HASH') {
	    croak("ARG. ERROR: When arg is not a hash ref. for set_strains");
	}
	$self->{strains} = $strain_href;    
    }
}

=head2 get_cds/set_cds

  Usage: my $cds_href = $phygepaml->get_cds(); 
         $phygepaml->get_cds($cds_href); 

  Desc: Get or set a hashref with: keys=ID and value=Bio::Seq object
        
  Ret: Get: a hashref with keys=ID and value=Bio::Seq object
       Set: none       

  Args: Get: none
        Set: a hashref with keys=ID and value=Bio::Seq object

  Side_Effects: Die if no argument is used with set_cds or if the argument
                used is not a hashref. with Bio::Seq objects as values.

  Example: my %cds = %{$phygepaml->get_cds()};
           $phygepaml->get_cds(\%cds);

=cut

sub get_cds {
    my $self = shift;
    return $self->{cds};
}

sub set_cds {
    my $self = shift;
    my $cds_href = shift;
   
    unless (defined $cds_href) {
	croak("ARG. ERROR: No arg. was used for set_cds function");
    }
    else {
	if ($cds_href =~ m/\w+/) {
	    my $obj = 'Bio::Seq';
	    unless (ref($cds_href) eq 'HASH') {
		croak("ARG. ERROR: $cds_href set_cds() isnt HASHREF");
	    }
	    else {
		my %hash = %{$cds_href};
		foreach my $key (keys %hash) {
		    my $val = $hash{$key};
		    if (ref($hash{$key}) !~ m/Bio::\w*Seq/) {
			croak("ARG. ERROR: $val used set_cds isnt $obj");
		    }
		}
	    }
	}
	$self->{cds} = $cds_href;    
    }
}

=head2 get_topotypes/set_topotypes

  Usage: my $topotypes_href = $phygetopo->get_topotypes();
         $phygetopo->get_topotypes($topotypes_href); 

  Desc: Get or set a hashref with key=ID and value=TopoTreeType object

  Ret: Get: a hashref with keys=ID and value=TopoTreeType
       Set: None

  Args: Get: None
        Set: a hashref with keys=ID and value=TopoTreeType

  Side_Effects: Die if no argument is used, if the argument is not a hashref
                or if its values are not TopoTreeType objects.

  Example: my %topotypes = %{$phygetopo->get_topotypes()};
           $phygetopo->set_topotypes(\%topotypes);

=cut

sub get_topotypes {
    my $self = shift;
    return $self->{topotypes};
}

sub set_topotypes {
    my $self = shift;
    my $topotypes_href = shift;
   
    unless (defined $topotypes_href) {
	croak("ARG. ERROR: No arg. was used for set_topotypes function");
    }
    else {
	unless (ref($topotypes_href) eq 'HASH') {
	    croak("ARG. ERROR: When arg is not a hash ref. for set_topotypes");
	}
	else {
	    my %hash = %{$topotypes_href};
	    my $obj = 'Bio::Tree::TopoType';
	    foreach my $key (keys %hash) {
		unless (ref($hash{$key}) eq $obj) {
		    croak("ARG. ERROR: Values arent $obj for set_topotypes()");
		}
	    }
	}
	$self->{topotypes} = $topotypes_href;    
    }
}


##########################
## ANALYTICAL FUNCTIONS ##
##########################



####
1; #
####
