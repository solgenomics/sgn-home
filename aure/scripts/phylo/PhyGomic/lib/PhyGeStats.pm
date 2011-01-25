
package PhyGeStats;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Math::BigFloat;
use Statistics::R;

use File::Temp qw/ tempfile tempdir/;


###############
### PERLDOC ###
###############

=head1 NAME

PhyGeStats.pm
A class with functions to analyze the results produced by the
PhyGomics modules

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use PhyGeStats;

  my $phygestats = PhyGeStats->new();


  ## Accessors

  my $srh = $phygestast->get_srh();
  $phygestats->set_srh($srh);

  my $phygetopo = $phygestats->get_phygetopo();
  $phygestats->set_phygetopo($phygetopo);

 
  ## Create a new R bridge and set into the phygestats object

  $phygestats->startR();
  $phygestats->stopR();

  ## Pass the Statistical::R methods

  $phygestats->stR($function, $parameter);


=head1 DESCRIPTION

 PhyGeStats is a wrapper for Statistics::R with specific functions to
 analyze the data contained in a phygetopo object. 

 For example:

 my $phygestats = PhyGeStats->new();
 $phygestats->startR();
 
 my $bargraph_filename = $phygestats->bargraph();
 my $treegraph_filename = $phygestats->treegraph();


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $phygestats = PhyGeStats->new($arguments_href);

  Desc: Create a phygestats object with the specified parameters.

  Ret: a PhyGeStats.pm object

  Args: A hash reference with the following key-value pairs: 
         + phygetopo   => a phygetopo object 
         + srh         => a Statistics::R object
        
  Side_Effects: Die if the argument used is not a hash or there are argument
                incompatibility (for example run_trees without run_distances).

  Example: my $phygestats_empty = PhyGeStats->new();
           my $phygestats = PhyGeStats->new(
                                             { 
                                               phygetopo => $phygetopo,
                                               srh       => $r_stats,
                                             }
                                           );

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %permargs = ( 
	phygetopo => 'PhyGeTopo',
	srh       => 'Statistics::R',
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
		    if (ref($args{$argkey}) ne $exp) {
			croak("ARG. ERROR: $args{$argkey} isnt $exp for new()");
		    }
		}	
	    }
	}
    }

    ## As default values ir will create two empty objects

    unless (defined $args_href->{phygetopo}) {
	$args_href->{phygetopo} = PhyGeTopo->new();
    }
    unless (defined $args_href->{srh}) {
	$args_href->{srh} = Statistics::R->new();
    }

    
    ## Set vars in the object

    $self->set_phygetopo($args_href->{phygetopo});
    $self->set_srh($args_href->{srh});
    
    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get_phygetopo

  Usage: my $phygetopo = $phygestats->get_phygetopo(); 

  Desc: Get a PhyGeTopo object contained into the PhyGeStats function

  Ret: A PhyGeTopo object

  Args: None

  Side_Effects: None

  Example:  my $phygetopo = $phygetopo->get_phygestats();

=cut

sub get_phygetopo {
    my $self = shift;
    return $self->{phygetopo};
}

=head2 set_phygetopo

  Usage: $phygestats->set_phygetopo($phygetopo);

  Desc: Set PhyGeTopo object in the phygestats object

  Ret: None

  Args: A PhyGeTopo object

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
                object

  Example: $phygetstats->set_phygetopo($phygetopo);

=cut

sub set_phygetopo {
    my $self = shift;
    my $phygetopo = shift;
   
    unless (defined $phygetopo) {
	croak("ARG. ERROR: No arg. was used for set_phygetopo function");
    }
    else {
	if ($phygetopo =~ m/\w+/) {
	    unless (ref($phygetopo) eq 'PhyGeTopo') {
		croak("ERROR: $phygetopo for set_phygetopo() isnt PhyGeTopo ");
	    }	    
	}
	$self->{phygetopo} = $phygetopo;    
    }
}

=head2 get_srh

  Usage: my $srh = $phygetopo->get_srh(); 

  Desc: Get a Statistics::R object contained into the PhyGeStats function

  Ret: A Statistics::R object

  Args: None

  Side_Effects: None

  Example: my $srh = $phygetopo->get_srh(); 

=cut

sub get_srh {
    my $self = shift;
    return $self->{srh};
}

=head2 set_srh

  Usage: $phygestats->set_srh($srh);

  Desc: Set Statistics::R object in the phygestats object

  Ret: None

  Args: A Statistics::R object

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
                object

  Example: $phygestats->set_srh($srh);

=cut

sub set_srh {
    my $self = shift;
    my $srh = shift;
   
    unless (defined $srh) {
	croak("ARG. ERROR: No arg. was used for set_srh function");
    }
    else {
	if ($srh =~ m/\w+/) {
	    unless (ref($srh) eq 'Statistics::R') {
		croak("ERROR: $srh for set_phygetopo() isnt Statistics::R obj");
	    }	    
	}
	$self->{srh} = $srh;    
    }
}



####
1; #
####
