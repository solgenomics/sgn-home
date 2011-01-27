
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
	phygetopo    => 'ARRAY',
	r_connection => "[1|0|Statistics::R]",
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
		    if (ref($args{$argkey}) =~ m/$exp/) {
			croak("ARG. ERROR: $args{$argkey} isnt $exp for new()");
		    }
		}	
	    }
	}
    }

    ## As default values ir will create two empty objects

    my $phygetopo_aref = $args_href->{phygetopo};
    unless (defined $args_href->{phygetopo}) {
	$phygetopo_aref = [];
    }

    my $srh = $args_href->{r_connection};
    if (defined $srh) {
	if ($srh =~ m/^1$/) {
	    $srh = Statistics::R->new();
	    $srh->start_sharedR();
	}
	elsif ($srh =~ m/^0$/) {
	    $srh = Statistics::R->new();
	}
    }
    else {
	$srh = Statistics::R->new();
	$srh->start_sharedR();
    }
    
    
    ## Set vars in the object

    $self->set_phygetopo($phygetopo_aref);
    $self->set_r_connection($srh);
    
    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get_phygetopo

  Usage: my $phygetopo_aref = $phygestats->get_phygetopo(); 

  Desc: Get a PhyGeTopo objects array ref contained into the PhyGeStats 
        function

  Ret: An array reference of PhyGeTopo objects

  Args: None

  Side_Effects: None

  Example:  my @phygetopo = @{$phygetopo->get_phygestats()};

=cut

sub get_phygetopo {
    my $self = shift;
    return $self->{phygetopo_aref};
}

=head2 set_phygetopo

  Usage: $phygestats->set_phygetopo($phygetopo_aref);

  Desc: Set PhyGeTopo array ref. in the phygestats object

  Ret: None

  Args: A PhyGeTopo object

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
                object

  Example: $phygetstats->set_phygetopo($phygetopo);

=cut

sub set_phygetopo {
    my $self = shift;
    my $phyget_aref = shift;
   
    unless (defined $phyget_aref) {
	croak("ARG. ERROR: No arg. was used for set_phygetopo function");
    }
    else {
	if ($phyget_aref =~ m/\w+/) {
	    unless (ref($phyget_aref) eq 'ARRAY') {
		croak("ERROR: $phyget_aref for set_phygetopo() isnt array ref");
	    }
	    else {
		foreach my $phygetop (@{$phyget_aref}) {
		    unless (ref($phygetop) eq 'PhyGeTopo') {
			my $msg = "ERROR: member $phygetop for set_phygetopo()";
			$msg .= "isnt array ref.";
			croak($msg);
		    }
		}
	    }
	}
	$self->{phygetopo_aref} = $phyget_aref;    
    }
}

=head2 add_phygetopo

  Usage: $phygestats->add_phygetopo($phygetopo_aref);

  Desc: Add PhyGeTopo array ref. in the phygestats object

  Ret: None

  Args: A PhyGeTopo object

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
                object

  Example: $phygetstats->add_phygetopo($phygetopo);

=cut

sub add_phygetopo {
    my $self = shift;
    my $phyget_aref = shift;
   
    unless (defined $phyget_aref) {
	croak("ARG. ERROR: No arg. was used for add_phygetopo function");
    }
    else {
	if ($phyget_aref =~ m/\w+/) {
	    unless (ref($phyget_aref) eq 'ARRAY') {
		croak("ERROR: $phyget_aref for add_phygetopo() isnt array ref");
	    }
	    else {
		foreach my $phygetop (@{$phyget_aref}) {
		    unless (ref($phygetop) eq 'PhyGeTopo') {
			my $msg = "ERROR: member $phygetop for add_phygetopo()";
			$msg .= "isnt array ref.";
			croak($msg);
		    }
		}
	    }
	}
	push @{$self->{phygetopo_aref}}, $phyget_aref;    
    }
}

=head2 delete_phygetopo

  Usage: my $phygetopo_aref = $phygestats->delete_phygetopo();

  Desc: Delete PhyGeTopo array ref. in the phygestats object

  Ret: Deleted PhyGeTopo objects

  Args: None

  Side_Effects: None

  Example: my $phygetopo_aref = $phygestats->delete_phygetopo();

=cut

sub delete_phygetopo {
    my $self = shift;
    
    my @del_phygetopo = ();
    my $phygetopo_aref = $self->get_phygetopo();
    
    my $i = 0;
    foreach my $phygetopo (@{$phygetopo_aref}) {
	push @del_phygetopo, delete @{$phygetopo_aref}[$i];
	$i++;
    }
    return \@del_phygetopo;
}



=head2 get_r_connection

  Usage: my $srh = $phygestats->get_r_connection(); 

  Desc: Get a Statistics::R object contained into the PhyGeStats function

  Ret: A Statistics::R object

  Args: None

  Side_Effects: None

  Example: my $srh = $phygestats->get_r_connection(); 

=cut

sub get_r_connection {
    my $self = shift;
    return $self->{r_connection};
}

=head2 set_r_connection

  Usage: $phygestats->set_r_connection($srh);

  Desc: Set Statistics::R object in the phygestats object

  Ret: None

  Args: A Statistics::R object

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
                object

  Example: $phygestats->set_r_connection($srh);

=cut

sub set_r_connection {
    my $self = shift;
    my $srh = shift;
   
    unless (defined $srh) {
	croak("ARG. ERROR: No arg. was used for set_r_connection function");
    }
    else {
	if ($srh =~ m/\w+/) {
	    unless (ref($srh) eq 'Statistics::R') {
		croak("ERROR: $srh set_r_connection() isnt Statistics::R obj");
	    }	    
	}
	$self->{r_connection} = $srh;    
    }
}


##########################
## ANALYTICAL FUNCTIONS ##
##########################

##
## Topotypes is a hash with the following data:
## key=ID and value=Bio::Tree::TopoType
## Bio::Tree::TopoType has topology and members
## So, it will produce two graphs
## 






####
1; #
####
