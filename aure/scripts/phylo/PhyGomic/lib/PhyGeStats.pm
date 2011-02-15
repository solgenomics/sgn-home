
package PhyGeStats;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Math::BigFloat;
use YapRI::Base;
use YapRI::Data::Matrix;

use File::Temp qw/ tempfile tempdir/;
use String::Random qw/ random_regex random_string/;
use Cwd;

use FindBin;
use lib "$FindBin::Bin/../lib";
use Bio::Tree::TopoType qw/ is_same_tree /;

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



=head1 DESCRIPTION

 PhyGeStats is a wrapper for YapRI::Base with specific functions to
 analyze the data contained in a phygetopo object. 

 For example:

 my $phygestats = PhyGeStats->new();

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
         + phygetopo    => a phygetopo object 
         + rbase        => a YapRI::Base object
         + matrix       => a YapRI::Data::Matrix object
        
  Side_Effects: Die if the argument used is not a hash or there are argument
                incompatibility (for example run_trees without run_distances).

  Example: my $phygestats_empty = PhyGeStats->new();
           my $phygestats = PhyGeStats->new(
                                             { 
                                               phygetopo => $phygetopo,
                                               rbase     => $rbase,
                                               matrix    => $matrix,
                                             }
                                           );

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %permargs = ( 
	phygetopo    => 'HASH',
	rbase        => "YapRI::Base",
	matrix       => "YapRI::Data::Matrix",
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

    my $srh = $args_href->{rbase};
    unless (defined $srh) {
	$srh = YapRI::Base->new();
    }

    my $matrix = $args_href->{matrix} || '';
    
    ## Set vars in the object

    $self->set_phygetopo($phygetopo_href);
    $self->set_rbase($srh);
    $self->set_matrix($matrix);
    
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

=head2 get_rbase

  Usage: my $srh = $phygestats->get_rbase(); 

  Desc: Get a YapRI::Base object contained into the PhyGeStats function

  Ret: A YapRI::Base object

  Args: None

  Side_Effects: None

  Example: my $srh = $phygestats->get_rbase(); 

=cut

sub get_rbase {
    my $self = shift;
    return $self->{rbase};
}

=head2 set_rbase

  Usage: $phygestats->set_rbase($srh);

  Desc: Set YapRI::Base object in the phygestats object

  Ret: None

  Args: A YapRI::Base object

  Side_Effects: Die if no argument is used or if the object is not a YapRI::Base
                object

  Example: $phygestats->set_rbase($srh);

=cut

sub set_rbase {
    my $self = shift;
    my $srh = shift;
   
    unless (defined $srh) {
	croak("ARG. ERROR: No arg. was used for set_rbase function");
    }
    else {
	if ($srh =~ m/\w+/) {
	    unless (ref($srh) eq 'YapRI::Base') {
		croak("ERROR: $srh set_rbase() isnt YapRI::Base obj");
	    }	    
	}
	$self->{rbase} = $srh;    
    }
}


=head2 get_matrix

  Usage: my $matrix = $phygestats->get_matrix(); 

  Desc: Get a YapRI::Data::Matrix object contained into the PhyGeStats function

  Ret: A YapRI::Data::Matrix

  Args: None

  Side_Effects: None

  Example: my $matrix = $phygestats->get_matrix();

=cut

sub get_matrix {
    my $self = shift;
    return $self->{matrix};
}

=head2 set_matrix

  Usage: $phygestats->set_matrix($matrix);

  Desc: Set YapRI::Data::Matrix object in the phygestats object

  Ret: None

  Args: A YapRI::Data::Matrix object

  Side_Effects: Die if no argument is used or if it is not a YapRI::Data::Matrix

  Example: $phygestats->set_matrix($matrix);

=cut

sub set_matrix {
    my $self = shift;
    my $matrix = shift;
   
    unless (defined $matrix) {
	croak("ARG. ERROR: No arg. was used for set_matrix function");
    }
    else {
	if ($matrix =~ m/\w+/) {
	    unless (ref($matrix) eq 'YapRI::Data::Matrix') {
		croak("ERROR: $matrix set_matrix() isnt YapRI::Data::Matrix");
	    }	    
	}
	$self->{matrix} = $matrix;    
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


=head2 _compare_phygetopos

  Usage: my %phygetopos = $phygestats->_compare_phygetopos($toponame) 

  Desc: Compare phygetopologies and return a hash with:
          key   = unique_phygetopo_name
          value = hash ref. with key   = phygename
                                 value = topo_id

        For exampple, if it compare two phylotopo:
          'a' with topologies 'topoA', 'topoB' and 'topoC'
          'b' with topologies 'topoA', and 'topoB'

        If 'topoA' for 'a' and 'topoB' for 'b' are the same, it will return:
        %phygetopos = ( 'topo1' => {a => 'topoA', b => 'topoB' },
                        'topo2' => {a => 'topoB', b => undef   },
                        'topo3' => {a => 'topoC', b => undef   },
                        'topo4' => {a => undef,   b => 'topoA' }
        );


  Ret: %phygetopos, a hash (see above)

  Args: $toponame, a base name to use to define the new topologies

  Side_Effects: Die if some arguments are wrong.
                It uses 'topology' by default

  Example: my %phygetopos = $phygestats->_compare_phygetopos('topology') 

=cut

sub _compare_phygetopos {
    my $self = shift;
    my $base = shift ||
	'topology';

    my %phygt_comp = ();
    my %phygt_keep = ();

    ## Get count

    my $maxph = 0;

    my %phygt = %{$self->get_phygetopo()};

    foreach my $phygename (sort keys %phygt) {

	my %topotypes = %{$phygt{$phygename}->get_topotypes()};	
	$maxph += scalar(keys %topotypes);
    }

    my $fcount = length($maxph);

    foreach my $phygename (sort keys %phygt) {
	
	my %topotypes = %{$phygt{$phygename}->get_topotypes()};	
	foreach my $topo_id (sort keys %topotypes) {
	    
	    my $topology = $topotypes{$topo_id}->get_topology();
	    
	    ## Check if exists that topology

	    my $topo_match = 'none';
	    my $comp_match = '';

	    foreach my $comp_id (sort keys %phygt_comp) {
		my %types_comps = %{$phygt_comp{$comp_id}};
		
		## Compare with the topologies stored into %phygt_keep hash

		foreach my $phname (keys %types_comps) {
		    
		    my $ktopo = $phygt_keep{$comp_id}->{$phname}; 
		    if (is_same_tree($topology, $ktopo)) {
			$topo_match = $ktopo;
			$comp_match = $comp_id;
		    }
		}
	    }
	
	    if (ref($topo_match) eq 'Bio::Tree::Tree') {  ## if exists add old
		$phygt_keep{$comp_match} = { $phygename => $topology };
		$phygt_comp{$comp_match} = { $phygename => $topo_id };
	    }
	    else {                                  ## if doesnt exist, to new
		my $n = scalar(keys %phygt_comp) + 1;
		my $new_comp = $base . '_' . sprintf('%0' . $fcount . 's' , $n);
		$phygt_keep{$new_comp} = { $phygename => $topology };
		$phygt_comp{$new_comp} = { $phygename => $topo_id };
	    }
	}
    }
    return %phygt_comp;
}

=head2 _phygt2matrix

  Usage: my $matrix = $self->_phygt2matrix($mtxname, $rownames); 

  Desc: Creates a YapRI::Data::Matrix with the PhyGeTopo data

  Ret: $matrix, a YapRI::Data::Matrix object

  Args: $mtxname, matrix name,
        $rownames, base for the rownames

  Side_Effects: Use 'PhyGeTopo_Comp' as default matrix name.

  Example: my $ymatrix = $self->_phygt2matrix(); 

=cut

sub _phygt2matrix {
    my $self = shift;
    my $mtxname = shift
	|| 'PhyGeTopo_Comp';
    my $rowbase = shift;

    ## Get the phygetopo objects, colnames and coln

    my %phygt = %{$self->get_phygetopo()};
    my @colnames = sort keys %phygt;
    my $coln = scalar(@colnames);

    ## Get the comparison hash, rownames and rown

    my %comp_phygt = $self->_compare_phygetopos($rowbase);
    my @rownames = sort keys %comp_phygt;
    my $rown = scalar(@rownames);


    ## Create a new matrix

    my $mtx = { name     => $mtxname,
		coln     => $coln, 
		rown     => $rown, 
		colnames => \@colnames, 
		rownames => \@rownames,
    };

    my $matrix = YapRI::Data::Matrix->new($mtx);

    ## Get the data (sort always by names). If it is defined for the comparison
    ## hash it will get the members, if isnt defined, it will add 0.

    my @data = ();
    foreach my $rowname (sort keys %comp_phygt) {
	my %datah = %{$comp_phygt{$rowname}};
	foreach my $colname (@colnames) {
	    if (defined $datah{$colname}) {
		my $topo_id = $datah{$colname};
		my %topoty = %{$phygt{$colname}->get_topotypes()};
		my @members = @{$topoty{$topo_id}->get_members()};
		push @data, scalar(@members);
	    }
	    else {
		push @data, 0;
	    }
	}
    }

    ## Add the data to the matrix

    $matrix->set_data(\@data);

    return $matrix;
}

=head2 create_matrix

  Usage: $phystats->create_matrix($mtxname, $rowbasename); 

  Desc: Creates a YapRI::Data::Matrix with the PhyGeTopo data and load it into
        the accessor.

  Ret: None

  Args: $mtxname, matrix name,
        $rowbasename, a basename for rows

  Side_Effects: Die if no phygetopo objects were loaded into the object.
                If there is a matrix into the accessor, it will overwrite it.

  Example: $phystats->create_matrix('topomtx', 'topology');   

=cut

sub create_matrix {
    my $self = shift;
    my $mtxname = shift;
    my $rowbase = shift;

    ## Check that phygetopo contains data.

    my %phygt = %{$self->get_phygetopo()};

    if (scalar( keys %phygt) == 0) {
	croak("ERROR: There isnt any phygetopo data inside phygestat object.");
    }

    ## Run _phygt2matrix and load the result into the object

    my $matrix = $self->_phygt2matrix($mtxname, $rowbase); 

    $self->set_matrix($matrix);
}




####
1; #
####
