
package PhyGeStats;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Math::BigFloat;
use YapRI::Base qw/ r_var /;
use YapRI::Data::Matrix;
use YapRI::Graph::Simple;

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
  use PhyGeTopo;
  use PhyGeCluster;

  ## Prepare data from PhyGeCluster:
 
  my $phygecluster = PhyGeCluster->new( { blastfile => $blastfile } );
  $phygecluster->load_seqfile( { sequencefile => $seqfile } );
  $phygecluster->load_strainfile( { strainfile => $strainfile } );

  $phygecluster->run_alignments( { program => 'clustalw' } );
  $phygecluster->run_distances( { method => 'Kimura' } );
  $phygecluster->run_mltrees({ dnaml => {}});

  ## Prepare PhyGeTopo as analysis of PhyGeCluster:

  my %seqfams = %{$phygecluster->get_clusters()};
  my %strains = %{$phygecluster->get_strains()};

  my $phygetopo = PhyGeTopo->new({seqfams => \%seqfams, strains => \%strains});
  my %topotypes = $phygetopo->run_topoanalysis();

  ## Create a rbase object:

  my $rbase = YapRI::Base->new();


  ## Create a PhyGeStats:

  my $phygestats = PhyGeStats->new({ rbase     => $rbase, 
                                     phygetopo => { 'ML' => $phygetopo } } );


  ## Create the matrix table with the data

  $phygestats->create_matrix();

  ##  Run analysis and create the result files

  $phygestats->create_composition_table('MyTable.tab');
  $phygestats->create_composition_graph('MyCompositionGraph.bmp');
  $phygestats->create_tree_graph('MyTreeGraph.bmp');


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
        
  Side_Effects: Die if the argument used is not a hash or there are argument
                incompatibility (for example run_trees without run_distances).

  Example: my $phygestats_empty = PhyGeStats->new();
           my $phygestats = PhyGeStats->new(
                                             { 
                                               phygetopo => $phygetopo,
                                               rbase     => $rbase,
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
		    if (ref($args{$argkey}) ne $exp && $args{$argkey} =~ m/./) {
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

		foreach my $phname (sort keys %types_comps) {
		    
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

=head2 _tree_list

  Usage: my %trees = $phygestats->_tree_list($toponame) 

  Desc: Get a hash with key   = global_topology_id
                        value = tree (newick format)


  Ret: %trees, a hash (see above)

  Args: $toponame, a base name to use to define the new topologies

  Side_Effects: Die if some arguments are wrong.
                It uses 'topology' by default

  Example: my %trees = $phygestats->_tree_list('topology') 

=cut

sub _tree_list {
    my $self = shift;
    my $rowbase = shift;

    my %comp_phygt = $self->_compare_phygetopos($rowbase);
    my %phygt = %{$self->get_phygetopo()};

    my %trees = ();

    foreach my $phyname (sort keys %phygt) {

	my %topotypes = %{$phygt{$phyname}->get_topotypes()};
	
	foreach my $topo_id (sort keys %topotypes) {
	    
	    my $topo = $topotypes{$topo_id};

            ## Search the global topology name for this topo_id
	    
	    foreach my $gtopo_cid (sort keys %comp_phygt) {
		my $old_topo_id = $comp_phygt{$gtopo_cid}->{$phyname};
		if (defined $old_topo_id) {
		    if ($old_topo_id eq $topo_id) {
			$trees{$gtopo_cid} = $topo->get_topology_as_newick();
		    }
		}
	    }
	}
    }
    return %trees;
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

    ## Check matrix name

    if ($mtxname =~ m/:/) {
	croak("ERROR: Matrix name can not have non word characters like ':'");
    }

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

  Args: $mtxname, a matrix name,
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

=head2 create_composition_graph

  Usage: $phystats->create_composition_graph($filename, $graph_args); 

  Desc: Creates a file with a barplot graph

  Ret: None

  Args: $filename, name for the graph file (required)
        $graph_args, a hash reference with the arguments for the graph, 
                     according YapRI::Graph::Simple accessors

  Side_Effects: Die if no filename argument is used.
                Die if graph_args isnt a hash ref.
                Die if no rbase was set before run this method.
                Die if no matrix was set before run this method 
                Use default YapRI::Graph::Simple accessors

  Example: $phystats->create_composition_graph('MyFile.bmp');

=cut

sub create_composition_graph {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename was supplied to create_composition_graph.");
    
    my $grhref = shift;

    if (defined $grhref) {
	unless (ref($grhref) eq 'HASH') {
	    croak("ERROR: $grhref supplied to create_comp. isnt a HASHREF.");
	}
    }

    ## Check is exists rbase and matrix.

    my $rbase = $self->get_rbase();
    if (ref($rbase) ne 'YapRI::Base') {
	croak("ERROR: No rbase was set before run create_composition_graph.");
    }
    
    my $rmtx = $self->get_matrix();
    if (ref($rmtx) ne 'YapRI::Data::Matrix') {
	croak("ERROR: No matrix was set before run create_composition_graph.");
    }

    ## It needs the matrix by topotypes (rows), no by methods (cols)

    my $trmtx = $rmtx->transpose();
    my $t = $trmtx->get_name();

    ## Now it will create the default graph parameters.
    
    my $rown = $trmtx->get_rown();  ## To know how many groups are.
    my $coln = $trmtx->get_coln();

    ## Define the object to set the colors in the barblot

    my $barcols = 'barcolors <- terrain.colors(' . $rown . ')'; 

    my %grargs = (
	rbase  => $rbase,
	rdata  => { height => $trmtx },
	grfile => $filename, 
	device => { bmp     => { width  => 800, 
				 height => 600,
		    } 
	},
	grparams => { par => { cex => 1.5 }, 
	},
	sgraph => { barplot => { beside => 'TRUE',
				 main   => 'Topotypes Abundance',
				 xlab   => 'Topotype',
				 ylab   => 'Count',
				 col    => { $barcols => '' }, 
		    }
	},
	gritems => [ 
	    { 
		legend => {
		    x      => "topright",
		    legend => $trmtx->get_rownames(),
		    fill   => { $barcols => '' }, 
		}
	    } 
	],
	);
	
    foreach my $arg (keys %{$grhref}) {
	unless (exists $grargs{$arg}) {
	    croak("ERROR: $arg isnt a valid argument for create_graph");
	}
	else {
	    ## Ignore rbase, rdata and grfile

	    if ($arg !~ m/(rbase|rdata|grfile)/) {
		my %argsh = %{$grargs{$arg}};
		foreach my $karg (keys %argsh) {
		    $grargs{$arg}->{$karg} = $argsh{$karg}; ## remplace them 
		}
	    }
	}
    }

    ## Run the graph commands

    my $rgraph = YapRI::Graph::Simple->new(\%grargs);

    my $block = $rgraph->build_graph();

    $rgraph->run_graph($block);
}

=head2 create_composition_table

  Usage: $phystats->create_composition_table($filename); 

  Desc: Print into a tab delimiter file a table with: 
         + headers  = phygetopo methods (phygenames),
         + rownames = member_id
         + data     = topotype_id after the comparison (see _compare_phygetopos)

  Ret: None

  Args: $filename, a filename to print the table, 
        $rowbase, a basename for rows.

  Side_Effects: Die if no filename argument is used.
                Die if phygetopo accessor is empty.

  Example: $phystats->create_composition_table('MyTable.tab');   

=cut

sub create_composition_table {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename was supplied to create_composition_table.");
    
    my $rowbase = shift;

    ## Check that phygetopo contains data.

    my %phygt = %{$self->get_phygetopo()};

    if (scalar( keys %phygt) == 0) {
	croak("ERROR: There isnt any phygetopo data inside phygestat object.");
    }

    ## Get the composition hash to get the same names that the rest of 
    ## the methods.

    my %comp_phygt = $self->_compare_phygetopos($rowbase);

    ## It will create a hash with:
    ##    keys  = member_id,
    ##    value = hashref. with key = phygename and value = new topology_id

    my %phymembers = ();

    foreach my $phyname (sort keys %phygt) {

	my %topotypes = %{$phygt{$phyname}->get_topotypes()};
	
	foreach my $topo_id (sort keys %topotypes) {
	    
	    my @members = @{$topotypes{$topo_id}->get_members()};
            ## Search the global topology name for this topo_id
	    
	    my $gtopo_id = '';
	    foreach my $gtopo_cid (sort keys %comp_phygt) {
		my $old_topo_id = $comp_phygt{$gtopo_cid}->{$phyname};
		if (defined $old_topo_id) {
		    if ($old_topo_id eq $topo_id) {
			$gtopo_id = $gtopo_cid;
		    }
		}
	    }

	    ## Add to the hash

	    foreach my $member (@members) {
		my $clid = $member->id();
		unless (exists $phymembers{$clid}) {
		    $phymembers{$clid} = { $phyname => $gtopo_id };
		}
		else {
		    $phymembers{$clid}->{$phyname} = $gtopo_id;
		}
	    }
	}
    }

    ## Now all the members should be in the hash, it will print the file

    open my $fh, '>', $filename;

    my @phynames = sort keys %phygt;
    
    my $header = "\t" . join("\t", @phynames);
    print $fh "$header\n";
    foreach my $memb_id (sort keys %phymembers) {
	my $line = $memb_id . "\t";
	my @data = ();
	foreach my $phyname (@phynames) {
	    if (exists $phymembers{$memb_id}->{$phyname}) {
		if ($phymembers{$memb_id}->{$phyname} =~ m/./) {
		    push @data, $phymembers{$memb_id}->{$phyname};
		}
		else {
		    push @data, 'NA';
		}
	    }
	    else {
		push @data, 'NA';
	    }
	}
	$line .= join("\t", @data);
	print $fh "$line\n";
    }
}

=head2 create_tree_graph

  Usage: $phystats->create_tree_graph($filename, $gr_args_href); 

  Desc: Creates a image file with trees

  Ret: None

  Args: $filename, name for the graph file (required),
        $gr_args_href, a hash reference with the following keys:
          - stack   => (horizontal|vertical|matrix), default horizontal.
          - device  => a hash ref. with the R device arguments
          - grparam => a hash ref. with key = par and values = hashref. with 
                       par args. mfrow controls the number of plots inside
                       the device.
          - plot.phylog => a hash ref. with key=plot.phylog and 
                           value=hash ref. with this R function args.

  Side_Effects: Die if no filename argument is used.
                Die if no rbase was set before run this method. 

  Example: $phystats->create_tree_graph('MyFile.bmp');

=cut

sub create_tree_graph {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename was supplied to create_tree_graph.");
    
    my $grhref = shift;

    if (defined $grhref) {
	unless (ref($grhref) eq 'HASH') {
	    croak("ERROR: $grhref supplied create_tree_graph isnt a HASHREF.");
	}
    }
    
    ## Get trees

    my %trees = $self->_tree_list();

    ## Get rbase or die

    my $rbase = $self->get_rbase();
    if (ref($rbase) ne 'YapRI::Base') {
	croak("ERROR: No rbase was set before run create_composition_graph.");
    }

    ## If tree number is 0, die

    my $tree_n = scalar(keys %trees);

    if ($tree_n == 0) {
	croak("ERROR: No trees were are contained at phygestat object.");
    }

    ## 1) Create a block.

    my $block = 'TREEGRAPH_' . random_regex('\w\w\w\w');
    $rbase->create_block($block);

    ## 2) Load the R module to write phylogenetic trees

    $rbase->add_command('suppressPackageStartupMessages(library(ade4))', 
			$block );
    
    ## 4) Get size and par according stack

    my $stack = $grhref->{stack} || "horizontal";
    my $size = { width => 150 * $tree_n, height => 200 };
    my $mfrow = { mfrow => [1, $tree_n] };

    if ($stack eq 'vertical') {
	$size = { width => 200, height => 150 * $tree_n };
	$mfrow = { mfrow => [$tree_n, 1] };
    }
    elsif ($stack eq 'matrix') {
	sqrt($tree_n) =~ m/^(\d)/;
	my $row = $1;
	my $div = $tree_n / $row;
	$div =~ m/^(\d+)/;
	my $col = $1;
	if ($div - $col > 0) {
	    $col++;
	}
	$size = { width => 150 * $col, height => 150 * $row };
	$mfrow = { mfrow => [$row, $col] };
    }

    ## 5) Create the device

    my $defdev_cmd = 'bmp(filename="' . $filename . '", ';
    $defdev_cmd .= 'width = ' . $size->{width} . ', ';
    $defdev_cmd .= 'height = ' . $size->{height} . ')';

    my $dev_cmd = $grhref->{device} || $defdev_cmd;

    $rbase->add_command($dev_cmd, $block);

    ## 6) Add par arguments

    my $defpar = { par => $mfrow };
    $defpar->{par}->{cex} = 1.5;

    my $parhref = $grhref->{grparam} || $defpar;
    
    $rbase->add_command(r_var($parhref), $block);
    
    ## 7) Create the objects

    foreach my $topo_id (sort keys %trees) {
	my $obj_cmd = $topo_id . ' <- "' . $trees{$topo_id} . ';"'; 
	$rbase->add_command($obj_cmd, $block)
    }

    ## 8) Create the plotting commands

    my $def_plot = { 
	'plot.phylog' => { 
	    f          => 0.5, 
	    cnod       => 2, 
	    cleav      => 2, 
	    'clabel.l' => 3,
	    csub       => 3,
	    possub     => "bottomright",
	} 
    };

    my $plot = $grhref->{'plot.phylog'} || $def_plot; 

    foreach my $topo_id (sort keys %trees) {
	
	## It will overwrite each command with x = newick2phylog($topo_id)
	## and sub = $topo_id

	my $obj = 'newick2phylog(' . $topo_id . ')';
	$plot->{'plot.phylog'}->{x} = { $obj => '' };
	$plot->{'plot.phylog'}->{'sub'} = $topo_id;
	
	$rbase->add_command(r_var($plot), $block);
    }

    ## 9) Run commands

    $rbase->run_block($block);
    
}




####
1; #
####
