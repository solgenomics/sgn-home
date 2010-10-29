
package BlastCluster;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Cluster::SequenceFamily;

###############
### PERLDOC ###
###############

=head1 NAME

BlastCluster.pm
a class to cluster sequences based in blast results

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use BlastCluster;

  ## To get each cluster from a blast file

  my $blastcluster = BlastCluster->new({
                                   blastfile     => $filename1,
                                   blastformat   => $format, 
                                   sequencefile  => $filename2,
                                   strainfile    => $filename3,   
                                   clustervalues => { $variable => $condition},
                                   rootname      => $name,
                                   });

  ## Accessors

  $blastcluster->set_clusters(\%clusters);
  my $clusters_href = $blastcluster->get_clusters();

  ## Load sequences into the object

  $blastcluster->load_sequences($seqfile);

  ## Load strains

  $blastcluster->load_strains($strainfile); 


 

  ## Calculate the alignment

  my $blastcluster_aligned = $blastcluster->align_members();



=head1 DESCRIPTION

 Object to parse and cluster sequences based in blast results


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $blastcluster = BlastCluster->new($arguments_href);

  Desc: Create a blastcluster object

  Ret: a BlastCluster.pm object

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, clustervalues and rootname.

  Side_Effects: Die if the argument used is not a hash.
                blastformat = 'blasttable' by default.
                clustervalues = { percent_identity => '> 90', 
                                  length('hit')    => '> 30' }
                rootname = 'cluster'

  Example: my $blastcluster = BlastCluster->new(%arguments);

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %cluster_objs = ();
    my %strains = ();

    if (defined $args_href->{'blastfile'}) {

	## Get the cluster information

	my %clusters = parse_blastfile($args_href);

	## Parse the seqfile if exists

	my %seqmembers = ();
	if (defined $args_href->{'sequencefile'}) {
	    %seqmembers = parse_seqfile($args_href);
	}

	## Parse the strain file if exists

	my %seqstrains = ();
	if (defined $args_href->{'strainfile'}) {
	    %seqstrains = parse_strainfile($args_href);
	}
    
	## Now it will convert the seqid in the cluster hash into
	## Bio::Seq objects, and it will create the Bio::Cluster::SequenceFamily
	
	foreach my $cluster_id (keys %clusters) {
	    
	    my @memberseq = ();
	    my @member_ids = @{$clusters{$cluster_id}};
	    
	    foreach my $member_id (@member_ids) {
		if (exists $seqmembers{$member_id}) {
		    push @memberseq, $seqmembers{$member_id};
		}
		else {
		    my $seq = Bio::Seq->new( -id => $member_id );
		    push @memberseq, $seq;
		}
	    }

	    my $seqcluster = Bio::Cluster::SequenceFamily->new(
		-family_id => $cluster_id,
		-members   => \@memberseq, 
		);
	    $cluster_objs{$cluster_id} = $seqcluster;
	    
	}
    }
    
    ## Finally it will set the clusters in the object
    ## (if none blastfile was supplied, it will be an empty object)

    $self->set_clusters(\%cluster_objs);
    $self->set_strains(\%strains);

    return $self;
}

#################
### ACCESSORS ###
#################

=head2 get_clusters

  Usage: my $clusters_href = $blastcluster->get_clusters();

  Desc: Get a cluster data from BlastCluster object

  Ret: A hash refence with key=clustername and value=Bio::Cluster
       object
 
  Args: None
 
  Side_Effects: None
 
  Example: my $clusters_href = $blastcluster->get_clusters();

=cut

sub get_clusters {
  my $self = shift;

  return $self->{clusters};
}


=head2 set_clusters

  Usage: $blastcluster->set_clusters(\%clusters);

  Desc: Set a cluster data from BlastCluster object

  Ret: None
 
  Args: A hash with key=clustername and value=Bio::Cluster::SequenceFamily
        object
 
  Side_Effects: Die if the argument is not a hash reference
                or if the values for the hash are not 
                Bio::Cluster::SequenceFamily objects
 
  Example: $blastcluster->set_clusters(\%clusters);

=cut

sub set_clusters {
  my $self = shift;
  my $clusters_href = shift;

  if (defined $clusters_href) {
	if (ref($clusters_href) ne 'HASH') {
	    croak("ARGUMENT ERROR: $clusters_href is not a hash reference.\n");
	}
	else {
	    
            ## Check the objects

	    foreach my $cluster_name (keys %{$clusters_href}) {
		my $cluster = $clusters_href->{$cluster_name};
		if (ref($cluster) ne 'Bio::Cluster::SequenceFamily') {
		    croak("ARGUMENT ERROR: value for cluster=$cluster_name IS
                           NOT a Bio::Cluster::SequenceFamily object.\n");
		}
	    }

	    ## If it has survive to the checking it will set the object
	    
	    $self->{clusters} = $clusters_href;
	}
  }
  else {
      croak("ARGUMENT ERROR: No argument was used for set_clusters function\n");
  }
}


=head2 add_cluster

  Usage: $blastcluster->add_cluster($cluster_name, $member_aref);

  Desc: Add a new cluster to the BlastCluster object

  Ret: None
 
  Args: $cluster_name, a scalar
        $member_aref, an array reference with member as Bio::Seq objects
 
  Side_Effects: Die if the arguments are wrong or undef
 
  Example: $blastcluster->add_cluster('cluster_1', ['seq1', 'seq2']);

=cut

sub add_cluster {
  my $self = shift;
  my $cluster_name = shift ||
      croak("ARG. ERROR: No argument was supplied to add_cluster function.\n");
  my $member_aref = shift ||
      croak("ARG. ERROR: None member array ref. was added to add_cluster.\n");

  unless (ref($member_aref) eq 'ARRAY') {
      croak("ARG. ERROR: member array ref. argument is not an ARRAY REF.\n");
  }

  ## Check the members

  foreach my $memberseq (@{$member_aref}) {
      unless (ref($memberseq) =~ m/Bio::Seq/) {
	  croak("ARG. ERROR: member=$memberseq used for add_cluster is not an 
                 Bio::Seq object.\n");
      }
  }

  ## First, create the cluster object

  my $new_cluster = Bio::Cluster::SequenceFamily->new( 
      -family_id => $cluster_name,
      -members   => $member_aref,
      );

  ## Second get the clusters and add a new one
  
  my $clusters_href = $self->get_clusters();
  $clusters_href->{$cluster_name} = $new_cluster;
}


=head2 remove_cluster

  Usage: $blastcluster->remove_cluster($cluster_name);

  Desc: Remove a new cluster to the BlastCluster object

  Ret: None
 
  Args: $cluster_name, a scalar
 
  Side_Effects: Die if the arguments are wrong or undef
 
  Example: $blastcluster->remove_cluster('cluster_1');

=cut

sub remove_cluster {
  my $self = shift;
  my $cluster_name = shift ||
      croak("ARG. ERROR: No argument was supplied to add_cluster function.\n");
  
  my $clusters_href = $self->get_clusters();
  delete($clusters_href->{$cluster_name});
}

=head2 find_cluster

  Usage: my $cluster = $blastcluster->find_cluster($member_id);

  Desc: Find a cluster object by member_id

  Ret: A Bio::Cluster::SequenceFamily object
 
  Args: $member_id, an scalar
 
  Side_Effects: return undef if the member does not exist
 
  Example: my $cluster = $blastcluster->find_cluster('seq_id1');

=cut

sub find_cluster {
  my $self = shift;
  my $member_id = shift;
  
  my $cluster_found;
  my $clusters_href = $self->get_clusters();
  foreach my $cluster_name (keys %{$clusters_href}) {
      my $cluster = $clusters_href->{$cluster_name};
      
      foreach my $member ($cluster->get_members()){
	  if ($member->display_id eq $member_id) {
	      $cluster_found = $cluster;
	  }
      }
  }

  return $cluster_found;
}

=head2 get_strains

  Usage: my $strains_href = $blastcluster->get_strains();

  Desc: Get a strain data from BlastCluster object

  Ret: A hash refence with key=member_id and value=strain
 
  Args: None
 
  Side_Effects: None
 
  Example: my $strains_href = $blastcluster->get_strains();

=cut

sub get_strains {
  my $self = shift;

  return $self->{strains};
}


=head2 set_strains

  Usage: $blastcluster->set_strains(\%strains);

  Desc: Set a strain data from BlastCluster object

  Ret: None
 
  Args: A hash with key=member_id and value=strain
 
  Side_Effects: Die if the argument is not a hash refence
 
  Example: $blastcluster->set_strains(\%strains);

=cut

sub set_strains {
  my $self = shift;
  my $strains_href = shift;

  if (defined $strains_href) {
	if (ref($strains_href) ne 'HASH') {
	    croak("ARGUMENT ERROR: $strains_href is not a hash reference.\n");
	}
	else {	    
	    $self->{strains} = $strains_href;
	}
  }
  else {
      croak("ARGUMENT ERROR: No argument was used for set_strains function\n");
  }
}



#####################
## Other functions ##
#####################

=head2 parse_blastfile

  Usage: my %clusters = parse_blastfile($arguments_href);

  Desc: Parse the blast file using bioperl Bio::SearchIO
        Cluster the sequences based in clustervalues

  Ret: %cluster, a hash with keys=cluster_name and value=array_ref with
       sequences ids.

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, clustervalues and rootname.

  Side_Effects: Die if the argument used is not a hash or if none 
                blastfile is detailed.
                blastformat = 'blasttable' by default.
                clustervalues = { percent_identity => '> 90', 
                                  length('hit')    => '> 30' }
                rootname = 'cluster'

  Example: my %clusters = parse_blastfile($arguments_href);

=cut

sub parse_blastfile {
    my $arg_href = shift;
    
    my %clusters = ();
    my %cluster_members = ();

    if (defined $arg_href) {
	if (ref($arg_href) ne 'HASH') {
	    croak("ARGUMENT ERROR: $arg_href is not a hash reference.\n");
	}
	unless (defined $arg_href->{'blastfile'}) {
	    croak("ARGUMENT ERROR: 'blastfile' isnot defined for $arg_href.\n");
	}
	else {

	    ## 0) Define default values
	    
	    unless (defined $arg_href->{'blastformat'}) {
		$arg_href->{'blastformat'} = 'blasttable';
	    }
	    unless (defined $arg_href->{'clustervalues'}) {
		$arg_href->{'clustervalues'} = { 
		    "percent_identity" => ['>', 90], 
		    "hsp_length"    => ['>', 30] 
		};
	    }
	    unless (defined $arg_href->{'rootname'}) {
		$arg_href->{'rootname'} = 'cluster';
	    }
	    
	    ## 1) Parse the file and create the clusters

	    my %cluster_values = %{$arg_href->{'clustervalues'}};

	    my $searchio = Bio::SearchIO->new( 
		-format => $arg_href->{'blastformat'}, 
		-file   => $arg_href->{'blastfile'}
		);
	    
	    while( my $result = $searchio->next_result() ) {
		
		my $q_name = $result->query_name();

		while( my $hit = $result->next_hit() ) {

		    my $s_name = $hit->name();
   
		    while( my $hsp = $hit->next_hsp() ) {
			
			## define the last cluster n
			
			my $cluster_n = scalar(keys(%clusters));
			my $cluster_name = $arg_href->{'rootname'} . 
			                   '_' . 
					   $cluster_n;

			## define the conditions needed for be in a cluster.
			my $conditions_n = scalar(keys(%cluster_values));
			
			foreach my $cluster_cond (keys %cluster_values) {
			    my $cond = $cluster_values{$cluster_cond}->[0];
			    my $val = $cluster_values{$cluster_cond}->[1];
			    
			    if ($cond =~ m/^<$/) {
				if ($hsp->$cluster_cond < $val) {
				    $conditions_n--;
				}
			    } 
			    elsif ($cond =~ m/^<=$/) {
				if ($hsp->$cluster_cond <= $val) {
				    $conditions_n--;
				}
			    }
			    elsif ($cond =~ m/^==$/) {
				if ($hsp->$cluster_cond == $val) {
				    $conditions_n--;
				}
			    }
			    elsif ($cond =~ m/^>=$/) {
				if ($hsp->$cluster_cond >= $val) {
				    $conditions_n--;
				}
			    }
			    elsif ($cond =~ m/^>$/) {
				if ($hsp->$cluster_cond > $val) {
				    $conditions_n--;
				}
			    }
			    else {
				croak("ARG.ERROR: $cluster_cond hasnt permited
                                       value ('<', '<=', '==', '>=' or '>')");
			    }
			}

			## it will be a cluster if
			if ($conditions_n == 0) {
			    unless (exists $cluster_members{$q_name}) {
				unless (exists $clusters{$cluster_name}) {
				    $clusters{$cluster_name} = [$q_name];
				}
				else {
				    push @{$clusters{$cluster_name}}, $q_name;
				}
				$cluster_members{$q_name} = 1;
			    }
			    unless (exists $cluster_members{$s_name}) {
				unless (exists $clusters{$cluster_name}) {
				    $clusters{$cluster_name} = [$s_name];
				}
				else {
				    push @{$clusters{$cluster_name}}, $s_name;
				}
				$cluster_members{$s_name} = 1;
			    }
			}
			else {
			    ## Create a new cluster name
			    my $new_cluster_name =  $arg_href->{'rootname'} . 
				                    '_' . 
					            ($cluster_n+1);
			    
			    unless (exists $cluster_members{$q_name}) {
				$clusters{$new_cluster_name} = [$q_name];
				$cluster_members{$q_name} = 1;
			    }
			    unless (exists $cluster_members{$s_name}) {
				$clusters{$new_cluster_name} = [$s_name];
				$cluster_members{$s_name} = 1;
			    }
			}
		    }
		}
	    }  
	}
    }
    else {
	croak("ARG. ERROR: No argument was used for parse_blast function\n");
    }
    return %clusters;
}


=head2 parse_seqfile

  Usage: my %seqs = parse_seqfile($arguments_href);

  Desc: Parse the fasta file using bioperl Bio::SeqIO
        and return a hash with keys=seq_id and values=seqobjs

  Ret: %seqs, a hash with keys=seq_id and value= with
       sequences ids.

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, clustervalues and rootname.

  Side_Effects: Die if the argument used is not a hash.

  Example: my %seqs = parse_seqfile($arguments_href);

=cut

sub parse_seqfile {
    my $args_href = shift;
    
    my %seqs = ();
    
    if (defined $args_href) {
	if (ref($args_href) ne 'HASH') {
	    croak("ARG. ERROR: $args_href is not a hash reference.\n");
	}
	unless (defined $args_href->{'sequencefile'}) {
	    croak("ARG. ERROR: 'sequencefile' isnot defined for $args_href.\n");
	}
	else {
	    my $seqio = Bio::SeqIO->new( -format => 'fasta', 
					-file => $args_href->{'sequencefile'} );
	
	    while((my $seqobj = $seqio->next_seq())) {
		my $seqid = $seqobj->display_id();
		$seqs{$seqid} = $seqobj;
	    }
	}
    }
    
    return %seqs;
}

=head2 parse_strainfile

  Usage: my %strains = parse_strainfile($arguments_href);

  Desc: Parse the strain file and return a hash with keys=seq_id 
        and values=strain

  Ret: %strain, a hash with keys=seq_id and value= with
       sequences ids.

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, strainfile, clustervalues and rootname.

  Side_Effects: Die if the argument used is not a hash.
                Die if the if do not exists $args{'strainfile'}

  Example: my %strain = parse_strainfile($arguments_href);

=cut

sub parse_strainfile {
    my $arg_href = shift;
    
    my %strains = ();
    
    if (defined $arg_href) {
	if (ref($arg_href) ne 'HASH') {
	    croak("ARG. ERROR: $arg_href is not a hash reference.\n");
	}
	unless (defined $arg_href->{'strainfile'}) {
	    croak("ARG. ERROR: 'strainfile' isnot defined for $arg_href.\n");
	}
	else {
	    open my $fileio, '<', $arg_href->{'strainfile'};
	
	    while(<$fileio>) {
		chomp($_);
		my @cols = split('\t', $_);
		
		$strains{$cols[0]} = $cols[1];
	    }
	}
    }
    
    return %strains;
}

####
1; #
####
