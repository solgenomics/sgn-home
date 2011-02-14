
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
         + rbase => a YapRI::BAse object
        
  Side_Effects: Die if the argument used is not a hash or there are argument
                incompatibility (for example run_trees without run_distances).

  Example: my $phygestats_empty = PhyGeStats->new();
           my $phygestats = PhyGeStats->new(
                                             { 
                                               phygetopo => $phygetopo,
                                               rbase => $rbase,
                                             }
                                           );

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %permargs = ( 
	phygetopo    => 'HASH',
	rbase => "YapRI::Base",
	r_dir        => '\w+',
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

    my $r_dir = $args_href->{r_dir} || cwd();

    ## As default values ir will create two empty objects

    my $phygetopo_href = $args_href->{phygetopo};
    unless (defined $args_href->{phygetopo}) {
	$phygetopo_href = {};
    }

    my $srh = $args_href->{rbase};
    unless (defined $srh) {
	$srh = YapRI::Base->new();
    }
    
    ## Set vars in the object

    $self->set_phygetopo($phygetopo_href);
    $self->set_rbase($srh);
    
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

  Side_Effects: Die if no argument is used or if the object is not a PhyGeTopo
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

  Usage: my $matrix = $self->_phygt2matrix(); 

  Desc: Creates a YapRI::Data::Matrix with the PhyGeTopo data

  Ret: $matrix, a YapRI::Data::Matrix object

  Args: None

  Side_Effects: None

  Example: my $ymatrix = $self->_phygt2matrix(); 

=cut

sub _phygt2matrix {



}






=head2 _r_infile_tm

  Usage: my $file = $phygestats->_r_infile_tm($argument_href); 

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

  Example: my $file0 = $phygestats->_r_infile_tm();
           my $file1 = $phygestats->_r_infile_tm(
                                       { 
                                         phygetopo => $phygetopo1,
                                         dirname   => '/home/user/analysis/',
                                         filename  => 'phygetopo_tab_for_r.txt',
				       }
                                     );


=cut

sub _r_infile_tm {
    my $self = shift;
    my $arghref = shift;

    my $file = '';
    my $fh = '';

    ## Check arguments

    my %permarg = (
	dirname   => '\w+',
	filename  => '\w+',
	tempfile  => '[0|1|no|yes]',
	phygetopo => 'HASH'
	);


    if (defined $arghref) {
	unless (ref($arghref) eq 'HASH') {
	    croak("ERROR: $arghref isnt a hash ref. for _r_infile_topomemb()");
	}
	else {
	    foreach my $keyarg (keys %{$arghref}) {
		unless (exists $permarg{$keyarg}) {
		    croak("ERROR: $keyarg isnt a valid argument");
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

	my %topotypes = %{$phyg_href->{$phygename}->get_topotypes()};
	
	foreach my $topo_id (sort keys %topotypes) {

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

		my $grouptopo_id = $match;
		$topodata{$grouptopo_id}->{$phygename} = scalar(@members);
	    }
	    else {

		if (exists $topodata{$topo_id}) {

		    ## It will need a new topology_id, by default it will
		    ## add the phylotopo name if it exists before

		    $topo_id = $topo_id . '_' . $phygename;
		}
		$topodata{$topo_id}->{$phygename} = scalar(@members);
		$topologies{$topo_id} = $topology;
		$newick{$topo_id} = $newick;  
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
	
	my $header = "\t" . '"' . join("\"\t\"", @phygenames);

	## Add a topolovgy column
	$header .= "\"\t\"Topology\"\n";
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

=head2 _r_loadfile

  Usage: my $r_obj_name = $phygestats->_r_loadfile($file, $basename); 

  Desc: Load the topomember file into R as a file

  Ret: $r_obj_name, a R object name that contains the table, and also
       the name of the YapRI::Base block

  Args: $file, filename,
        $basename, for the R object

  Side_Effects: Die if no filename or R object basename is supplied.
                Die if the R connection is not set or is not started

  Example: my $r_obj_name = $phygestats->_r_loadfile($file, 'TopoMembTable'); 


=cut

sub _r_loadfile {
    my $self = shift;
    my $file = shift ||
	croak("ERROR: No file was supplied to _r_loadfile function.");
    my $basename = shift ||
	croak("ERROR: No R basename was supplied to _r_loadfile function.");

    my $srh = $self->get_rbase();

    my $r_obj = $basename . "_" . random_regex('\w\w\w\w\w\w');

    ## Check that exists the connection and that it has been started

    unless (defined $srh) {
	croak("ERROR: rbase is not set for PhyGeStat object");
    }
 
    
    ## Now build the command to load the data

    my $r_cmd = "$r_obj <- read.table(\"$file\", header=TRUE)";

    ## Run the command

    $srh->create_block($r_obj);
    $srh->add_command("$r_cmd", $r_obj);

    ## Return the object

    return $r_obj;
}

=head2 _phto2ymtx

  


=cut





####
1; #
####
