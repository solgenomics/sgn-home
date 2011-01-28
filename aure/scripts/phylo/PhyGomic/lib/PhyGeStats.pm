
package PhyGeStats;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Math::BigFloat;
use Statistics::R;

use File::Temp qw/ tempfile tempdir/;
use Cwd;

use Bio::Tree::TopoType;

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
	phygetopo    => 'HASH',
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

    my $phygetopo_href = $args_href->{phygetopo};
    unless (defined $args_href->{phygetopo}) {
	$phygetopo_href = {};
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

    $self->set_phygetopo($phygetopo_href);
    $self->set_r_connection($srh);
    
    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get_phygetopo

  Usage: my $phygetopo_href = $phygestats->get_phygetopo(); 

  Desc: Get a PhyGeTopo objects array ref contained into the PhyGeStats 
        function

  Ret: An hash reference with key=name and value=PhyGeTopo objects

  Args: None

  Side_Effects: None

  Example:  my %phygetopo = %{$phygetopo->get_phygestats()};

=cut

sub get_phygetopo {
    my $self = shift;
    return $self->{phygetopo_href};
}

=head2 set_phygetopo

  Usage: $phygestats->set_phygetopo($phygetopo_aref);

  Desc: Set PhyGeTopo hash ref. in the phygestats object

  Ret: None

  Args: 

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
                object

  Example: $phygetstats->set_phygetopo(\%phygetopos);

=cut

sub set_phygetopo {
    my $self = shift;
    my $phyget_href = shift ||
	croak("ARG. ERROR: No arg. was used for set_phygetopo function");
    
    unless (ref($phyget_href) eq 'HASH') {
	croak("ERROR: $phyget_href for set_phygetopo() isnt hash ref");
    }
    else {
	foreach my $n (keys %{$phyget_href}) {
	    my $phytop = $phyget_href->{$n};
	    unless (ref($phytop) eq 'PhyGeTopo') {
		croak("ERROR: $phytop for set_phygetopo() PhyGeTopo inst obj");
	    }
	}
	$self->{phygetopo_href} = $phyget_href;    
    }
}

=head2 add_phygetopo

  Usage: $phygestats->add_phygetopo($name, $phygetopo);

  Desc: Add PhyGeTopo array ref. in the phygestats object

  Ret: None

  Args: $name, a scalar for PhyGeTopo object
        $phygetopo, a PhyGeTopo object

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
                object

  Example: $phygetstats->add_phygetopo($name, $phygetopo);

=cut

sub add_phygetopo {
    my $self = shift;
    my $name = shift ||
	croak("ARG. ERROR: No arg. was used for add_phygetopo function");
    my $phyget = shift ||
	croak("ARG. ERROR: No phygetopo arg. was used for add_phygetopo()");
   
    unless (ref($phyget) eq 'PhyGeTopo') {
	croak("ERROR: $phyget for add_phygetopo() isnt PhyGeTopo object");
    }
	    
    $self->{phygetopo_href}->{$name} = $phyget;
}

=head2 delete_phygetopo

  Usage: my $phygetopo_href = $phygestats->delete_phygetopo();

  Desc: Delete PhyGeTopo array ref. in the phygestats object

  Ret: Deleted PhyGeTopo objects

  Args: None

  Side_Effects: None

  Example: my $phygetopo_aref = $phygestats->delete_phygetopo();

=cut

sub delete_phygetopo {
    my $self = shift;
    
    my $phygetop_href = $self->get_phygetopo();
    $self->set_phygetopo({});
    return $phygetop_href;
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


=head2 _r_infile_topomemb

  Usage: my $file = $phygestats->_r_infile_topomemb($argument_href); 

  Desc: Print into a file topotype_id and member_count in table format

  Ret: $file, a filename

  Args: optional, a hash reference with the following arguments:
        dirname   => dirname to create the file
        filename  => filename for the file
        tempfile  => [1|0], enable/disable the file creation as temp file
                     (1 by default)
        phygetopo => array ref. with PhyGeTopo objects to be used
                     (all contained in the phygestats object by default)

  Side_Effects: Die if some arguments are wrong

  Example: my $file0 = $phygestats->_r_infile_topomemb();
           my $file1 = $phygestats->_r_infile_topomemb(
                                       { 
                                         phygetopo => $phygetopo1,
                                         dirname   => '/home/user/analysis/',
                                         filename  => 'phygetopo_tab_for_r.txt',
				       }
                                     );


=cut

sub _r_infile_topomemb {
    my $self = shift;
    my $arghref = shift;

    my $file = '';
    my $fh = '';

    ## Check arguments

    my %permarg = (
	dirname   => '\w+',
	filename  => '\w+',
	tempfile  => '[0|1|no|yes]',
	phygetopo => 'PhyGeTopo'
	);


    if (defined $arghref) {
	unless (ref($arghref) eq 'HASH') {
	    croak("ERROR: $arghref isnt a hash ref. for _r_infile_topomemb()");
	}
	else {
	    foreach my $keyarg (keys %{$arghref}) {
		unless (exists $permarg{$keyarg}) {
		    croak("ERROR: $keyarg isnt a permited argument");
		}
		else {
		    my $valarg = $arghref->{$keyarg};
		    unless ($valarg =~ m/$permarg{$keyarg}/i) {
			croak("ERROR: $valarg doesnt have a right value");
		    }
		}
	    }
	}
    }

    

    ## First, create a hash with key=topotype_id and value=aref_member_counts
    ## the first member always will be the headers

    my %topodata = ();
    my %topologies = ();
    my %newick = ();

    my $phyg_href = $arghref->{phygetopo} || $self->get_phygetopo();
    
    unless (ref($phyg_href) eq 'HASH') {
	croak("ERROR: phygetopo value $phyg_href isnt an hash reference");
    }

    foreach my $phygename (sort keys %{$phyg_href}) {
	print STDERR "TEST PHYGENAME: $phygename\n";
	my %topotypes = %{$phyg_href->{$phygename}->get_topotypes()};
	
	foreach my $topo_id (sort keys %topotypes) {
	    print STDERR "TEST TOPOID: $topo_id\n";
	    my $topology = $topotypes{$topo_id}->get_topology();
	    my $newick = $topotypes{$topo_id}->get_topology_as_newick();
	    my @members = @{$topotypes{$topo_id}->get_members()};
	    
	    ## Check if exists that topology
	    my $match;
	    foreach my $ktopo_id (keys %topologies) {
		my $keep_topo = $topologies{$ktopo_id};
		if (Bio::Tree::TopoType::is_same_tree($topology, $keep_topo)) {
		    $match = $ktopo_id;
		}
	    }

	    ## If exists that topology add to the old one, if not, add a new
	    ## one 

	    if (defined $match) {
		print STDERR "MATCH ENABLE\n";
		my $grouptopo_id = $match;
		$topodata{$grouptopo_id}->{$phygename} = scalar(@members);
	    }
	    else {
		print STDERR "MATCH DISABLE\n";
		if (exists $topodata{$topo_id}) {

		    ## It will need a new topology_id
		    $topo_id = $topo_id . '_x';
		}
		$topodata{$topo_id}->{$phygename} = scalar(@members);
		$topologies{$topo_id} = $topology;
		$newick{$topo_id} = $newick;  
		print STDERR "ADDED TOPOLOGY $topology => $topo_id\n";
	    }
	}
    }

    ## Second, create the file and print the hash only if there are any
    ## phygetopo object

    my @phygenames = sort keys %{$phyg_href};
    if (scalar(@phygenames) > 0) {

	## First, create a file

	my $filename = $arghref->{filename};
	my $dirname = $arghref->{dirname} || getcwd();
	my $tempfile = $arghref->{tempfile};
	my $msg = "ERROR: filename or dirname are incompat args. with tempfile";
    
	if (defined $filename) {
	    if (defined $tempfile && $tempfile =~ m/[1|yes]/i) {
		croak($msg);
	    }
	    $file = $dirname . '/' . $filename;
	    open $fh, '>', $file;
	}
	else {	
	    ($fh, $file) = tempfile('r_infile_topom_XXXXXX', TMPDIR => 1);
	}
	
	my $header = "\t" . '"' . join("\"\t\"", @phygenames) . "\"\t\"T\"\n";
	print $fh $header;
    
	foreach my $row_name (sort keys %topodata) {
	    my %row_data = %{$topodata{$row_name}};
	    print $fh '"' . $row_name . '"';
	    foreach my $phname (@phygenames) {
		if (exists $row_data{$phname}) {
		    print $fh "\t$row_data{$phname}";
		}
		else {
		    print $fh "\t0";
		}
	    }
	    print $fh "\t" . '"' . $newick{$row_name} . '"' . "\n";
	}

	## Close only if the filehandle was open
	close($fh);
    }
    
    
    return $file;
}




####
1; #
####
