
package PhyGeTopo;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Try::Tiny;
use Math::BigFloat;


###############
### PERLDOC ###
###############

=head1 NAME

PhyGeTopo.pm
A class to analyze tree topology.

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use PhyGeTopo;

  my $phygetopo = PhyGeTopo->new({ seqfams     => $hash_ref,
                                   strains     => $hash_ref,
                                   annotations => $hash_ref,
                                   topotypes   => $hash_ref, 
                                   });


  ## Accessors

  my $seqfam_href = $phygetopo->get_seqfams();
  $phygetopo->set_seqfams($seqfam_href);

  my $strains_href = $phygetopo->get_strains();
  $phygetopo->set_strains($strains_href);

  my $annotations_href = $phygetopo->get_annotations();
  $phygetopo->set_annotations($annotation_href);

  my $topotypes = $phygetopo->get_topotypes();
  $phygetopo->set_topotypes($topotypes_href);

  ## To run different tools over the phygetopo object.

  my @topotypes = $phygetopo->run_topoanalysis($args_href);
  foreach my $topotype (@topotypes) {
    my $topology_tree = $topotype->get_topology();
    my @tree_members = @{$topotype->get_members()};
    my $member_count = $topotype->member_count();
  }

  ## To print results

  my %files = $phygetopo->out_topologies();


=head1 DESCRIPTION

 Object to analyze the topologies for trees


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $phygetopo = PhyGeTopo->new($arguments_href);

  Desc: Create a phygetopo object with the specified parameters.

  Ret: a PhyGeTopo.pm object

  Args: A hash reference with the following key-value pairs: 
         + seqfams     => a hash ref. with key=id and 
                          value=Bio::Cluster::SequenceFamily obj,
         + strains     => a hash ref. with key=members and value=strains
         + annotations => a hash ref. with key=id and value=annotation
         + topotypes   => a hash ref. with key=id and value=Bio::Tree::TopoType 
                          object
        
  Side_Effects: Die if the argument used is not a hash or there are argument
                incompatibility (for example run_trees without run_distances).

  Example: my $phygetopo = PhyGeTopo->new($arguments_href);

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %permargs = ( 
	seqfams     => {},
	strains     => {},
	annotations => {},
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
    $self->set_annotations($args_href->{annotations});
    $self->set_topotypes($args_href->{topotypes});
    
    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get_seqfams

  Usage: my $seqfams_href = $phygetopo->get_seqfams(); 

  Desc: Get a hash with key=ID and value=Bio::Cluster::SequenceFamily objects

  Ret: A hash reference with key=ID and value=Bio::Cluster::SequenceFamily obj.

  Args: None

  Side_Effects: None

  Example: my %seqfams = %{$phygetopo->get_seqfams()}; 

=cut

sub get_seqfams {
    my $self = shift;
    return $self->{seqfams};
}

=head2 set_seqfams

  Usage: $phygetopo->set_seqfam($seqfam_href);

  Desc: Set Bio::Cluster::SequenceFamily objects in the PhyGeTopo object

  Ret: None

  Args: A hash reference with key=ID and value=Bio::Cluster::SequenceFamily obj.

  Side_Effects: Die if no argument is used or if the object is not a hash ref.
                or if its values arent Bio::Cluster::SequenceFamily objects.

  Example: $phygetopo->set_seqfam(\%seqfams);

=cut

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

=head2 get_strains

  Usage: my $strains_href = $phygetopo->get_strains(); 

  Desc: Get a hashref with the members and the strains

  Ret: a hashref with keys=member_id and value=strain

  Args: None

  Side_Effects: None

  Example: my %strains = %{$phygetopo->get_strains()};

=cut

sub get_strains {
    my $self = shift;
    return $self->{strains};
}

=head2 set_strains

  Usage: $phygetopo->set_strains($strains_href);

  Desc: Set strain information into PhyGeTopo object

  Ret: None

  Args: a hash reference with key=member_id and value=strain

  Side_Effects: Die if no argument is used.

  Example: $phygetopo->set_seqfam(\%strains);

=cut

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

=head2 get_annotations

  Usage: my $annotations_href = $phygetopo->get_annotations(); 

  Desc: Get a hashref with the IDs and the annotations hashrefs.

  Ret: a hashref with keys=ID and value=hashref with key=tag and value=string
       or a hashref with keys=ID and value=string if tag argument is used

  Args: A scalar as $tag (annotation tag like GenBank...)

  Side_Effects: None

  Example: my %annotations = %{$phygetopo->get_annotations()};
           my %genbank_annot = %{$phygetopo->get_annotations('GenBank')};

=cut

sub get_annotations {
    my $self = shift;
    my $tag = shift;

    my $annot_href = {};
    if (defined $tag) {
	my %glob_annot = %{$self->{annotations}};
	foreach my $id (keys %glob_annot) {
	    if (exists $glob_annot{$id}->{$tag}) {
		$annot_href->{$id} = $glob_annot{$id}->{$tag};
	    }
	}
    }
    else {
	$annot_href = $self->{annotations};
    }
    return $annot_href;
}

=head2 set_annotations

  Usage: $phygetopo->set_annotations($annotations_href);

  Desc: Set annotation information into PhyGeTopo object

  Ret: None

  Args: a hash reference with key=member_id and value=hashref with
        key=tag and value=annotation string

  Side_Effects: Die if no argument is used or if the value of the hashref
                is not another hashref.

  Example: $phygetopo->set_annotations(\%annotations);

=cut

sub set_annotations {
    my $self = shift;
    my $annot_href = shift;
   
    unless (defined $annot_href) {
	croak("ARG. ERROR: No arg. was used for set_annotations function");
    }
    else {
	unless (ref($annot_href) eq 'HASH') {
	    croak("ARG. ERROR: When arg. is not a hashref for set_annotations");
	}
	else {
	    my %hash = %{$annot_href};
	    foreach my $key (keys %hash) {
		unless (ref($hash{$key}) eq 'HASH') {
		    croak("ARG. ERROR: When val. isnt hashref set_annotations");
		}
	    }
	}
	$self->{annotations} = $annot_href;    
    }
}

=head2 get_topotypes

  Usage: my $topotypes_href = $phygetopo->get_topotypes(); 

  Desc: Get a hashref with key=ID and value=TopoTreeType object

  Ret: a hashref with keys=ID and value=TopoTreeType

  Args: None

  Side_Effects: None

  Example: my %topotypes = %{$phygetopo->get_topotypes()};

=cut

sub get_topotypes {
    my $self = shift;
    return $self->{topotypes};
}

=head2 set_topotypes

  Usage: $phygetopo->set_topotypes($topotypes_href);

  Desc: Set topotype information into PhyGeTopo object

  Ret: None

  Args: a hash reference with key=ID and value=Bio::Tree::TopoType

  Side_Effects: Die if no argument is used, if the argument is not a hashref
                or if its values are not TopoTreeType objects.

  Example: $phygetopo->set_topotypes(\%topotypes);

=cut

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


=head2 run_topoanalysis

  Usage: %topotypes = $phygetopo->run_topoanalysis($args_href);

  Desc: run a topology analysis over all the trees defined in the
        phygetopo object.

  Ret: An hash with key=Topotype_ID and value=Bio::Tree::TopoTypes object

  Args: a hash reference with keys=argument and value=value, such as:
        branch_cutoffs => hashref. with key=branch_lenght cutoff and 
                          value=new quant. value.
        base_toponame  => a scalar to be used as basename for topologies
                          (topotype by default)

  Side_Effects: Set topotypes in the phygetopo object.

  Example: %topotypes = $phygetopo->run_topoanalysis({ 
                                        branch_cutoffs => { 0.1 => 1 },
                                    });

=cut

sub run_topoanalysis {
    my $self = shift;
    my $args_href = shift;

    my %permargs = ( 
	branch_cutoffs => 'HASH',
	base_toponame  => '\w+',
	);

    ## Check argument

    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ERROR: $args_href used for run_topoanalysis() isnt HASHREF");
	}
	else {
	    my %args = %{$args_href};
	    foreach my $argkey (keys %args) {
		my $exp = $permargs{$argkey};
		unless (defined $exp) {
		    croak("ERROR: $argkey isnt permited for run_topoanalysis");
		}
		else {
		    if (ref($args{$argkey}) !~ m/$exp/) {
			my $err = "ERROR: $args{$argkey} isnt $exp for ";
			$err .= "run_topoanalysis()";
			croak($err);
		    }
		}	
	    }
	}
    }

    ## Define topotypes hash

    my %topotypes = ();
    my %toponewicks = ();
    my %tord = ();  ## To order by members

    ## Now it will get all the trees and strains

    my $strains_href = $self->get_strains();
    my %seqfams = %{$self->get_seqfams()};
    
    foreach my $seqfam_id (keys %seqfams) {
	my $tree = $seqfams{$seqfam_id}->tree();

	if (defined $tree) {

	    ## Create the args and pass brach_cutoffs if it exists

	    my %topoargs = ( tree    => $tree,
			     strains => $strains_href,
 		);
	    if (exists $args_href->{branch_cutoffs}) {
		$topoargs{branch_cutoffs} = $args_href->{branch_cutoffs}
	    }

	    ## Run make_topotype and get the newick line for it

	    my $topotype = Bio::Tree::TopoType::_make_topotype(\%topoargs);
	    
	    if (defined $topotype) {
		my $toponewick = Bio::Tree::TopoType::_tree2newick($topotype);
		
		if (exists $toponewicks{$toponewick}) {
		    push @{$toponewicks{$toponewick}}, $tree;
		    $tord{$toponewick}++;
		}
		else {
		    $toponewicks{$toponewick} = [$tree];
		    $topotypes{$toponewick} = $topotype;
		    $tord{$toponewick} = 1;
		}
	    }
	}
    }

    ## After get all the topotypes it will create one topotype object
    ## for each of them.

    my %topoobjs = ();
    my $i = 0;
    my $base = $args_href->{base_toponame} || "topotype";

    foreach my $topo_nwid (sort {$tord{$b} <=> $tord{$a}} keys %tord) {
	$i++;

	my $topotype_obj = Bio::Tree::TopoType->new(
	    {
		topology => $topotypes{$topo_nwid},
		members  => $toponewicks{$topo_nwid}, 
	    }
	    );
	
	my $name = $base . '_' . $i;
	$topoobjs{$name} = $topotype_obj;
    }

    ## Finally set the object

    $self->set_topotypes(\%topoobjs);
    
    return %topoobjs;			   
}






####
1; #
####
