
package PhyGeCluster;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Math::BigFloat;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Cluster::SequenceFamily;

use Bio::Tools::Run::Alignment::Clustalw;
use Bio::Tools::Run::Alignment::Kalign;
use Bio::Tools::Run::Alignment::MAFFT;
use Bio::Tools::Run::Alignment::Muscle;
use Bio::Tools::Run::Alignment::TCoffee;

###############
### PERLDOC ###
###############

=head1 NAME

PhyGeCluster.pm
a class to cluster sequences based in blast results

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use PhyGeCluster;

  ## To get each cluster from a blast file

  my $phygecluster = PhyGeCluster->new({
                                   blastfile     => $filename1,
                                   blastformat   => $format, 
                                   sequencefile  => $filename2,
                                   strainfile    => $filename3,   
                                   clustervalues => { $variable => $condition},
                                   rootname      => $name,
                                   });

  ## Accessors

  $phygecluster->set_clusters(\%clusters);
  my $clusters_href = $phygecluster->get_clusters();

  ## Load sequences into the object

  $phygecluster->load_seqfile($seqfile);

  ## Load strains

  $phygecluster->load_strainfile($strainfile); 


 

  ## Calculate the alignment

  my $phygecluster_aligned = $phygecluster->align_members();



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

  Usage: my $phygecluster = PhyGeCluster->new($arguments_href);

  Desc: Create a phygecluster object
        Use fastblast_parser argument to use the non-bioperl
        blast parser (faster)

  Ret: a PhyGeCluster.pm object

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, clustervalues and rootname.

  Side_Effects: Die if the argument used is not a hash.
                blastformat = 'blasttable' by default.
                clustervalues = { percent_identity => '> 90', 
                                  length('hit')    => '> 30' }
                rootname = 'cluster'

  Example: my $phygecluster = PhyGeCluster->new(%arguments);

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %cluster_objs = ();
    my %strains = ();

    if (defined $args_href->{'blastfile'}) {

	## Get the cluster information

	my %clusters = ();
	if (defined $args_href->{'fastblast_parser'}) {
	   %clusters = fastparse_blastfile($args_href);
	}
	else {
	    %clusters = parse_blastfile($args_href);
	}

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

  Usage: my $clusters_href = $phygecluster->get_clusters();

  Desc: Get a cluster data from PhyGeCluster object

  Ret: A hash refence with key=clustername and value=Bio::Cluster
       object
 
  Args: None
 
  Side_Effects: None
 
  Example: my $clusters_href = $phygecluster->get_clusters();

=cut

sub get_clusters {
  my $self = shift;

  return $self->{clusters};
}


=head2 set_clusters

  Usage: $phygecluster->set_clusters(\%clusters);

  Desc: Set a cluster data from PhyGeCluster object

  Ret: None
 
  Args: A hash with key=clustername and value=Bio::Cluster::SequenceFamily
        object
 
  Side_Effects: Die if the argument is not a hash reference
                or if the values for the hash are not 
                Bio::Cluster::SequenceFamily objects
 
  Example: $phygecluster->set_clusters(\%clusters);

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

  Usage: $phygecluster->add_cluster($cluster_name, $member_aref);

  Desc: Add a new cluster to the PhyGeCluster object

  Ret: None
 
  Args: $cluster_name, a scalar
        $member_aref, an array reference with member as Bio::Seq objects
 
  Side_Effects: Die if the arguments are wrong or undef
 
  Example: $phygecluster->add_cluster('cluster_1', ['seq1', 'seq2']);

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

  Usage: $phygecluster->remove_cluster($cluster_name);

  Desc: Remove a new cluster to the PhyGeCluster object

  Ret: None
 
  Args: $cluster_name, a scalar
 
  Side_Effects: Die if the arguments are wrong or undef
 
  Example: $phygecluster->remove_cluster('cluster_1');

=cut

sub remove_cluster {
  my $self = shift;
  my $cluster_name = shift ||
      croak("ARG. ERROR: No argument was supplied to add_cluster function.\n");
  
  my $clusters_href = $self->get_clusters();
  delete($clusters_href->{$cluster_name});
}

=head2 find_cluster

  Usage: my $cluster = $phygecluster->find_cluster($member_id);

  Desc: Find a cluster object by member_id

  Ret: A Bio::Cluster::SequenceFamily object
 
  Args: $member_id, an scalar
 
  Side_Effects: return undef if the member does not exist
 
  Example: my $cluster = $phygecluster->find_cluster('seq_id1');

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

  Usage: my $strains_href = $phygecluster->get_strains();

  Desc: Get a strain data from PhyGeCluster object

  Ret: A hash refence with key=member_id and value=strain
 
  Args: None
 
  Side_Effects: None
 
  Example: my $strains_href = $phygecluster->get_strains();

=cut

sub get_strains {
  my $self = shift;

  return $self->{strains};
}


=head2 set_strains

  Usage: $phygecluster->set_strains(\%strains);

  Desc: Set a strain data from PhyGeCluster object

  Ret: None
 
  Args: A hash with key=member_id and value=strain
 
  Side_Effects: Die if the argument is not a hash refence
 
  Example: $phygecluster->set_strains(\%strains);

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


###################################
## OBJECT DATA LOADING FUNCTIONS ##
###################################

=head2 load_seqfile

  Usage: $phygecluster->load_seqfile($arguments_href);

  Desc: Parse the seqquence file using parse_seqfile function and
        load the results into the PhyGeCluster object

  Ret: None

  Args: A hash reference with the following key: sequencefile

  Side_Effects: Die if the argument used is not a hash or if none 
                sequencefile is detailed.
                Ignore the sequences that are not into the 
                Bio::Cluster::SequenceFamily object

  Example: $phygecluster->load_seqfile($arguments_href);

=cut

sub load_seqfile {
    my $self = shift;
    my $args_href = shift;
    
    if (defined $args_href) {
	if (ref($args_href) ne 'HASH') {
	    croak("ARG. ERROR: $args_href is not a hash reference.\n");
	}
	unless (defined $args_href->{'sequencefile'}) {
	    croak("ARG. ERROR: 'sequencefile' is not defined for $args_href\n");
	}
	else {
	    ## First, get the sequences.

	    my %seqs = parse_seqfile->($args_href);

	    ## Get the clusters and the different members

	    my %clusters = %{$self->get_clusters()};
	    foreach my $cluster_id (keys %clusters) {
		my $seqcluster_obj = $clusters{$cluster_id};
		my @members = $seqcluster_obj->get_members();
		
		foreach my $seqmember (@members) {
		    my $seq = $seqs{$seqmember->display_id()};
		    if (defined $seq) {
			$seqmember->seq($seq->seq());
		    }
		}  
	    }
	    
	}
    }
    else {
      	croak("ARG. ERROR: No argument was used for load_seqfile function\n");
    }
}

=head2 load_strainfile

  Usage: $phygecluster->load_strainfile($arguments_href);

  Desc: Parse the strain file using parse_seqfile function and
        load the results into the PhyGeCluster object

  Ret: None

  Args: A hash reference with the following key: strainfile

  Side_Effects: Die if the argument used is not a hash or if none 
                strainfile is detailed.

  Example: $phygecluster->load_strainfile($arguments_href);

=cut

sub load_strainfile {
    my $self = shift;
    my $arg_href = shift;
    
    my %strains = ();
    
    if (defined $arg_href) {
	if (ref($arg_href) ne 'HASH') {
	    croak("ARG. ERROR: $arg_href is not a hash reference.\n");
	}
	unless (defined $arg_href->{'strainfile'}) {
	    croak("ARG. ERROR: 'strainfile' is not defined for $arg_href.\n");
	}
	else {
	    my %strains = parse_strainfile($arg_href);
	    $self->set_strains(\%strains);
	}
    }
    else {
      	croak("ARG. ERROR: No arg. was used for load_strainfile function\n");
    }
}



#######################
## PARSING FUNCTIONS ##
#######################

=head2 parse_blastfile

  Usage: my %clusters = parse_blastfile($arguments_href);

  Desc: Parse the blast file using bioperl Bio::SearchIO
        Cluster the sequences based in clustervalues. 

  Ret: %cluster, a hash with keys=cluster_name and value=array_ref with
       sequences ids.

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, clustervalues and rootname.
        Permited clustervalues: evalue, expect, frac_identical, frac_conserved,
		                gaps, hsp_length,  num_conserved, 
                                num_identical, score, bits, percent_identity

  Side_Effects: Die if the argument used is not a hash, if none 
                blastfile is detailed or if the clustervalues is not
                permitted
                blastformat = 'blasttable' by default.
                clustervalues = { percent_identity => ['>', 90'], 
                                  hsp_length       => ['>', 60'] }
                rootname = 'cluster'
                Print status messages with $arg_href->{'debug'}
                Calculate the percent_identity as frac_identical * 100
                instead to use percent_identity from the blast result
                (it is a Bio::SearchIO

  Example: my %clusters = parse_blastfile($arguments_href);

=cut

sub parse_blastfile {
    my $arg_href = shift;
    
    my %clusters = ();
    my %cluster_members = ();

    ## Define the permited blast variables according the 
    ## HSP methods

    my %blastvar = ( evalue           => 1,
		     expect           => 1,
		     frac_identical   => 1,
		     frac_conserved   => 1, 
		     gaps             => 1,
		     hsp_length       => 1,
		     num_conserved    => 1, 
		     num_identical    => 1,
		     score            => 1,
		     bits             => 1,
		     percent_identity => 1, 
	);


    ## Now parsing...

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
		    "hsp_length"    => ['>', 60] 
		};
	    }
	    else {
		foreach my $clval_key (keys %{$arg_href->{'clustervalues'}}) {
		    unless (exists $blastvar{$clval_key}) {
			my $er = "WRONG BLAST VARIABLE: $clval_key is not a ";
			$er .= "permited var. for parse_blast function.";
			croak($er);
		    }
		}
	    }

	    unless (defined $arg_href->{'rootname'}) {
		$arg_href->{'rootname'} = 'cluster';
	    }
	    
	    ## 1) Parse the file and create the clusters
	    
	    my %clustervals = %{$arg_href->{'clustervalues'}};
	    
	    my $searchio = Bio::SearchIO->new( 
		-format => $arg_href->{'blastformat'}, 
		-file   => $arg_href->{'blastfile'}
		);
	    
	    my $ht = $searchio->result_count();
	    my $hc = 0;

	    ## define the last cluster n
			    
	    my $cluster_n = 0;

	    while( my $result = $searchio->next_result() ) {
		
		my $q_name = $result->query_name();		

		while( my $hit = $result->next_hit() ) {

		    $hc++;
		    my $s_name = $hit->name();
		    if (exists $arg_href->{'report_status'}) {
			print_parsing_status($hc, $ht, 
					   "Percentage of blast file parsed:");
		    }
		    
		    while( my $hsp = $hit->next_hsp() ) {

			## define the last cluster n
			    
			$cluster_n = scalar(keys(%clusters));

			## Define the cluster name according the q_name
		
			my $cluster_name = $cluster_members{$q_name};
			    
			unless ($q_name eq $s_name) {
			    
			    ## define the conditions needed for be in 
			    ## a cluster.
			    my $conditions_n = scalar(keys(%clustervals));

			    foreach my $cluster_cond (keys %clustervals) {
				my $cutc = $clustervals{$cluster_cond};
				
				if (ref($cutc) ne 'ARRAY') {
				    my $er2 = "WRONG clustervalues format ";
				    $er2 .= "($cutc) for parse_blastfile()";
				    croak($er2);
				}
				
				my $cond = $clustervals{$cluster_cond}->[0];
				my $val = $clustervals{$cluster_cond}->[1];
				unless ($val =~ m/^\d+$/) {
				    my $er3 = "WRONG clustervalues val=";
				    $er3 .= "$cluster_cond only int. are ";
				    $er3 .= "permited for parse_blastfile()";
				    croak($er3);
				}
				    
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
				    my $er4 = "WRONG clustervalues val=";
				    $er4 .= "$cluster_cond, it is not a ";
				    $er4 .= "permited condition ('<','<=',";
				    $er4 .= "'==','>=' or '>')\n";
				    croak($er4);
				}
			    }    
			    
			    unless (exists $cluster_members{$s_name}) {
				
				## New cluster condition, it will overwrite
				## the default cluster name 
				
				if ($conditions_n > 0) {
				    $cluster_name = $arg_href->{'rootname'} 
				    . '_' . ($cluster_n+1);
				}
				
				unless (exists $clusters{$cluster_name}) {
				    $clusters{$cluster_name} = [$s_name];
				}
				else {
				    push @{$clusters{$cluster_name}}, 
				    $s_name;
				}
				$cluster_members{$s_name} = $cluster_name;
			    }
			}
			else {  
                       
			    ## Selfblast hit. It always will add that if it
			    ## has not added before as subject (s_name) hit

			    unless (exists $cluster_members{$q_name}) {
			    
				$cluster_name =  $arg_href->{'rootname'} . 
				    '_' . 
				    ($cluster_n+1);
				    
				unless (exists $clusters{$cluster_name}) {
				    $clusters{$cluster_name} = [$q_name];
				}
				else {
				    push @{$clusters{$cluster_name}}, 
				    $q_name;
				}
				$cluster_members{$q_name} = $cluster_name;
			    }
			}
		    }		    
		}
	    }
	}
	if (exists $arg_href->{'report_status'}) {
			print "\n";
	}
    }
    else {
	croak("ARG. ERROR: No argument was used for parse_blast function\n");
    }
    return %clusters;
}


=head2 fastparse_blastfile

  Usage: my %clusters = fastparse_blastfile($arguments_href);

  Desc: Parse the blast file using its own script
        Cluster the sequences based in clustervalues

  Ret: %cluster, a hash with keys=cluster_name and value=array_ref with
       sequences ids.

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, clustervalues and rootname.
        Permited clustervalues: query_id, subject_id, percent_identity, 
                                align_length, mismatches, gaps_openings, 
                                q_start, q_end, s_start, s_end, e_value, 
                                bit_score.

  Side_Effects: Die if the argument used is not a hash or if none 
                blastfile is detailed.
                Die if is using non permited clustervalues
                blastformat = 'blasttable' by default.
                clustervalues = { percent_identity => ['>', 90'], 
                                  align_length     => ['>', 60'] }
                rootname = 'cluster'
                Print status messages with $arg_href->{'debug'}

  Example: my %clusters = parse_blastfile($arguments_href);

=cut

sub fastparse_blastfile {
    my $arg_href = shift;
    
    my %clusters = ();
    my %cluster_members = ();

    ## Define the possible clustervalues

    my %blastvar = (
	            'query_id'         => 0, 
		    'subject_id'       => 1, 
		    'percent_identity' => 2, 
		    'align_length'     => 3,
		    'mismatches'       => 4, 
		    'gaps_openings'    => 5,
		    'q_start'          => 6,
		    'q_end'            => 7, 
		    's_start'          => 8, 
		    's_end'            => 9,
		    'e_value'          => 10,
		    'bit_score'        => 11
	);

    ## Now parsing...

    if (defined $arg_href) {
	
	if (ref($arg_href) ne 'HASH') {
	    croak("ARGUMENT ERROR: $arg_href is not a hash reference.\n");
	}

	unless (defined $arg_href->{'blastfile'}) {
	    croak("ARGUMENT ERROR: 'blastfile' isnot defined for $arg_href.\n");
	}
	else {

	    ## Check fastblastparser compatibilities
	    my $bl_format = $arg_href->{'blastformat'};
	    if (defined $bl_format && $bl_format !~ m/m8|blasttable/) {
		my $er0 = "ARG. ERROR: fastparse_blastfile function only ";
		$er0 .= "can be used with blasttable or m8 formats";
		croak($er0);
	    }
	    
		
	    ## 0) Define default values

	    unless (defined $arg_href->{'clustervalues'}) {
		$arg_href->{'clustervalues'} = { 
		    "percent_identity" => ['>', 90], 
		    "align_length"     => ['>', 60] 
		};
	    }
	    else {
		foreach my $clval_key (keys %{$arg_href->{'clustervalues'}}) {
		    unless (exists $blastvar{$clval_key}) {
			my $er = "WRONG BLAST VARIABLE: $clval_key is not a ";
			$er .= "permited var. for fastparse_blast function.";
			croak($er);
		    }
		}
	    }
	    unless (defined $arg_href->{'rootname'}) {
		$arg_href->{'rootname'} = 'cluster';
	    }
	    
	    ## 1) Parse the file and create the clusters
	    
	    my %clustervals = %{$arg_href->{'clustervalues'}};
	    
	    my $infile = $arg_href->{'blastfile'};
	    
	    my $ht = `cut -f1 $infile | wc -l`;
	    chomp($ht);
	    my $hc = 0;

	    open my $in, '<', $infile;

	    while(<$in>) {
		chomp($_);
		$hc++;

		my %datap = ();
		my @data = split(/\t/, $_);
		foreach my $field (sort {$blastvar{$a} <=> $blastvar{$b}} 
				   keys(%blastvar) ) {

		    $datap{$field} = $data[$blastvar{$field}];
		}

		my $q_name = $datap{'query_id'};
		my $s_name = $datap{'subject_id'};
		
		## define the last cluster n
		
		my $cluster_n = scalar(keys(%clusters));

		## Define the cluster name according the q_name
		
		my $cluster_name = $cluster_members{$q_name};
		
		if (exists $arg_href->{'report_status'}) {
		    print_parsing_status($hc, $ht, 
					 "Percentage of blast file parsed:");
		}
		    
		unless ($q_name eq $s_name) {
		    
		    ## define the conditions needed for be in 
		    ## a cluster.
		    my $conditions_n = scalar(keys(%clustervals));
			
		    foreach my $cluster_cond (keys %clustervals) {
			my $cutc = $clustervals{$cluster_cond};
			if (ref($cutc) ne 'ARRAY') {
			    my $er2 = "WRONG clustervalues format ";
			    $er2 .= "($cutc) for parse_blastfile()";
			    croak($er2);
			}

			my $cond = $clustervals{$cluster_cond}->[0];
			my $val = $clustervals{$cluster_cond}->[1];
			unless ($val =~ m/^\d+$/) {
			    my $er3 = "WRONG clustervalues val=$cluster_cond";
			    $er3 .= "only int. are permited parse_blastfile()";
			    croak($er3);
			}


			if ($cond =~ m/^<$/) {
			    if ($datap{$cluster_cond} < $val) {
				$conditions_n--;
			    }
			} 
			elsif ($cond =~ m/^<=$/) {
			    if ($datap{$cluster_cond} <= $val) {
				$conditions_n--;
			    }
			}
			elsif ($cond =~ m/^==$/) {
			    if ($datap{$cluster_cond} == $val) {
				$conditions_n--;
			    }
			}
			elsif ($cond =~ m/^>=$/) {
			    if ($datap{$cluster_cond} >= $val) {
				$conditions_n--;
			    }
			}
			elsif ($cond =~ m/^>$/) {
			    if ($datap{$cluster_cond} > $val) {
				$conditions_n--;
			    }
			}
			else {
			    my $er4 = "WRONG clustervalues val=";
			    $er4 .= "$cluster_cond, it is not a ";
			    $er4 .= "permited condition ('<','<=',";
			    $er4 .= "'==','>=' or '>')\n";
			    croak($er4);
			}
		    }    
			    
		    unless (exists $cluster_members{$s_name}) {
				
			## New cluster condition, it will overwrite the
			## default cluster name 

			if ($conditions_n > 0) {
			    $cluster_name =  $arg_href->{'rootname'} . 
				'_' . 
				($cluster_n+1);
			}

			unless (exists $clusters{$cluster_name}) {
			    $clusters{$cluster_name} = [$s_name];
			}
			else {
			    push @{$clusters{$cluster_name}}, $s_name;
			}
			$cluster_members{$s_name} = $cluster_name;
		    }
		}
		else {  
		
		    ## Selfblast hit. It always will add that if it has not
		    ## added before as subject (s_name) hit
		
		    unless (exists $cluster_members{$q_name}) {
		    
			$cluster_name =  $arg_href->{'rootname'} . 
			    '_' . 
			    ($cluster_n+1);

			unless (exists $clusters{$cluster_name}) {
			    $clusters{$cluster_name} = [$q_name];
			}
			else {
			    push @{$clusters{$cluster_name}}, $q_name;
			}
			$cluster_members{$q_name} = $cluster_name;
		    }
		}
	    }
	    close $in;
	    if (exists $arg_href->{'report_status'}) {
		print "\n";
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
                Print status messages with $arg_href->{'debug'}

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
	    croak("ARG. ERROR: 'sequencefile' is not defined for $args_href\n");
	}
	else {

	    my $seqfile = $args_href->{'sequencefile'};
	    my $seqio = Bio::SeqIO->new(-format => 'fasta', 
					-file => $seqfile );
	    
	    
	    my $st = `grep -c '>' $seqfile`;
	    chomp($st);
	    my $sc = 0;
	    
	    while((my $seqobj = $seqio->next_seq())) {
		
		$sc++;
		my $seqid = $seqobj->display_id();
		$seqs{$seqid} = $seqobj;
		if (exists $args_href->{'report_status'}) {
		    print_parsing_status($sc, $st, 
					 "Percentage of seq file parsed:");
		}
	    }
	    if (exists $args_href->{'report_status'}) {
		print "\n";
	    }
	}
    }
    else {
      	croak("ARG. ERROR: No argument was used for parse_seqfile function\n");
    }
    
    return %seqs;
}




=head2 parse_strainfile

  Usage: my %strains = parse_strainfile($arguments_href);

  Desc: Parse the strain file and return a hash with keys=seq_id 
        and values=strain

  Ret: %strain, a hash with keys=seq_id and value= strain

  Args: A hash reference with the following keys: blastfile, blastformat,
        sequencefile, strainfile, clustervalues and rootname.

  Side_Effects: Die if the argument used is not a hash.
                Die if the if do not exists $args{'strainfile'}
                Print status messages with $arg_href->{'debug'}

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
	    croak("ARG. ERROR: 'strainfile' is not defined for $arg_href.\n");
	}
	else {
	    open my $fileio, '<', $arg_href->{'strainfile'};
	    my $st = `cut -f1 $arg_href->{'strainfile'} | wc -l`;
	    chomp($st);
	    my $sc = 0;

	    while(<$fileio>) {
		chomp($_);
		$sc++;
		my @cols = split('\t', $_);
		$strains{$cols[0]} = $cols[1];
		if (exists $arg_href->{'report_status'}) {
		    print_parsing_status($sc, $st, 
					 "Percentage of strain file parsed:");
		}
	    }
	    close $fileio;
	    if (exists $arg_href->{'report_status'}) {
		print "\n";
	    }
	}
    }
    else {
      	croak("ARG. ERROR: No arg. was used for parse_strainfile function\n");
    }
    
    return %strains;
}

=head2 print_parsing_status

  Usage: print_parsing_status($progress, $total, $message);

  Desc: Print as STDERR the percentage of the file parsed

  Ret: none

  Args: $progress, a scalar, an integer with the progress
        $total, a scalar, a integer with the total
        $message to print before the percentage

  Side_Effects: Die if 1st and 2nd arguments are undefined or are not integers
                Print as default message: Percentage of the file parsed: ";

  Example: print_parsing_status($progress, $total, $message);

=cut

sub print_parsing_status {
    my $a = shift ||
	croak("ERROR: 1st argument is not defined for print_parsing_status()");
    my $z = shift ||
	croak("ERROR: 2nd argument is not defined for print_parsing_status()");
    my $message = shift ||
	"Percentage of the file parsed:";

    unless ($a =~ m/^\d+$/) {
	croak("ERROR: 1st argument is not an int. for print_parsing_status()");
    }
    unless ($a =~ m/^\d+$/) {
	croak("ERROR: 1st argument is not an int. for print_parsing_status()");
    }
    
    my $perc = $a * 100 / $z;
    my $perc_obj = Math::BigFloat->new($perc);
    my $perc_formated = $perc_obj->bfround(-2);

    print STDERR "\t$message $perc_formated %\r";
}


##########################
## ANALYTICAL FUNCTIONS ##
##########################

=head2 cluster_sizes

  Usage: my %cluster_sizes = phygecluster->cluster_sizes();
         my %cluster_sizes = phygecluster->cluster_sizes($min, $max);

  Desc: return a hash with keys=cluster_id and values=size

  Ret: %cluster_sizes,  with keys=cluster_id and values=size

  Args: $min, $max values to select the sizes ($min <= $size <= $max)

  Side_Effects: If $min value is detailed but $max value is absent, 
                $max = $min.
                If there are not any cluster with the size requeriments,
                it will return undef

  Example: my %cluster_sizes = phygecluster->cluster_sizes();
           my %singlets = phygecluster->cluster_sizes(1,1);
=cut

sub cluster_sizes {
    my $self = shift;
    my $min = shift;
    my $max = shift;

    if ($min && !$max) {
	$max = $min;
    }

    my %cluster_sizes = ();

    my %clusters = %{$self->get_clusters()};

    foreach my $cluster_name (keys %clusters) {
	my $size = $clusters{$cluster_name}->size();

	if (defined $min) {
	    if ($min <= $size && $size <= $max) {
		$cluster_sizes{$cluster_name} = $size;
	    }
	}
	else {
	    $cluster_sizes{$cluster_name} = $size;
	}
    }

    return %cluster_sizes;
}


=head2 run_alignments

  Usage: phygecluster->run_alignments({ program => $prog, 
                                        parameters => $parameters_href});

  Desc: Create the sequence aligns for each cluster and load it into the
        Bio::Cluster::SequenceFamily object

  Ret: none (load the alignments into the Bio::Cluster::SequenceFamily object)

  Args: $args_href, a hash reference with the following permited keys:
        program => $scalar, parameters => $array_ref with program parameters

  Side_Effects: Died if some of the parameters are wrong.
                It does not run alignments over singlets (cluster with 1 member)

  Example:  phygecluster->run_alignments({ program => $prog, 
                                        parameters => $parameters_href});

=cut

sub run_alignments {
    my $self = shift;
    my $args_href = shift ||
	croak("ARG. ERROR: None args. were suplied to run_alignments()");

    ## Check the arguments and create the factory object

    my $factory;
    unless (ref($args_href) eq 'HASH') {
	croak("ARG. ERROR: Arg=$args_href for run_alignments() isn't hashref.");
    }
    else {
	my %args = %{$args_href};
	
	unless (exists $args{program}) {
	    my $er0 = "ARG. ERROR: Args. supplied to run_alignments doesn't ";
	    $er0 .= "contains 'program' specification";
	    croak($er0);
	}
	else {
	    my %perm_prog = (
		clustalw => 'Bio::Tools::Run::Alignment::Clustalw',
		kalign   => 'Bio::Tools::Run::Alignment::Kalign',
		mafft    => 'Bio::Tools::Run::Alignment::MAFFT',
		muscle   => 'Bio::Tools::Run::Alignment::Muscle',
		tcoffee  => 'Bio::Tools::Run::Alignment::TCoffee',
		);

	    ## Check the program

	    my $pprog = join(',', sort(keys %perm_prog));
	    unless (exists $perm_prog{$args{program}}) {
		my $er1 = "ARG. ERROR: program selected for run_alignments()";
		$er1 .= "is not in the list of permited programs ($pprog)";
		croak($er1);
	    }

	    unless (exists $args{parameters}) {
		my $er3 = "ARG. ERROR: Args. supplied to run_alignments ";
		$er3 .= "doesn't contains 'parameters' specification";
		croak($er3);
	    }
	    else {
		unless (ref($args{parameters}) eq 'ARRAY') {
		    my $er4 = "ARG. ERROR: 'parameters' args. supplied to ";
		    $er4 .= "run_alignments()";
		    $er4 .= "doesn't contains 'parameters' specification";
		    croak($er4);
		}
		else {
		    my @parameters = @{$args{parameters}};
		    $factory = $perm_prog{$args{program}}->new();
		}
	    }
	}	    
    }

    if (defined $factory) {
 
	my %clusters = %{$self->get_clusters()};

	foreach my $cluster_id (keys %clusters) {
	 
	    my @seq_members = $clusters{$cluster_id}->get_members();

	    ## Only make sense if the cluster has more than one member

	    if (scalar(@seq_members) > 1) {
	
		my $alignobj = $factory->align(\@seq_members);

		## Also it will set the alinments into the SequenceFamily 
		## object

		$clusters{$cluster_id}->alignment($alignobj);
	    }
	}
    }
}

=head2 run_distance

  Usage: my %distances = $phygecluster->run_distance($method);

  Desc: Calculates a distance matrix for all pairwise distances of
        sequences in an alignment using Bio::Align::DNAStatistics object

  Ret: %distances, a hash with keys=cluster_id and value=Bio::Matrix::PhylipDis
       object

  Args: $method, a scalar to pass to the Bio::Align::DNAStatistics->distance()
        function (JukesCantor,Uncorrected,F81,Kimura,Tamura,F84 or TajimaNei).

  Side_Effects: Died if some of the parameters are wrong.
                Return undef for the clusters that do not have alignment into 
                the Bio::Cluster::SequenceFamily object.
                It will run JukesCantor distance by default.

  Example:  my %distances = $phygecluster->run_distance('Kimura');

=cut

sub run_distance {
    my $self = shift;
    my $method = shift || 'JukesCantor';

    

    
}


=head2 prune_cluster

  Usage: $phygecluster->prune_cluster()

  Desc: Remove sequence members and clusters that have not some conditions

  Ret: none (modify the PhyGeCluster, Bio::Cluster::SequenceFamily and 
       Bio::SimpleAlign objects)

  Args: A hash reference with the following parameters:
        by_simplealign => { functionname => { condition => integer } }
        by_strain_comp => { relation => $strainname_aref }
        
        relations could be: 'random', 'min_distance', 'max_distance', 
                            'min_length' and 'max_length'
        for example: 
           + { random => ['A', 'A', 'B'] }
             will select 2 seq from A strain and 1 from B strain. 
           + { max_length => ['A', 'B', 'A'] }
             will order the sequences by length and will take 2xA + 1xB
           + { min_distance => ['A', 'B'], min_distance => ['A', 'C']}
             will take: 2xA + 1xB + 1xC where A close to B and A' close to C
           + { min_distance => ['A', 'B', 'C']}
             will take: 1xA + 1xB + 1xC with the small distance between them
             calculated as A-B and A-C

  Side_Effects: Died if some of the parameters are wrong.
              

  Example:  $phygecluster->prune_cluster(
               { by_simplealign => { 
                     num_sequences => 5, 
                     length        => 500,
                 },          
               });



=cut

sub prune_cluster {
    my $self = shift;
    my $args_href = shift ||
	croak("ARG. ERROR: None args. were suplied to prune_cluster()");

    ## Check the arguments

    my %perm_args = ( by_simplealign => 'HASH',
		      by_strain_comp => 'HASH',
	);

    foreach my $argkey (keys %{$args_href}) {
	
	unless (exists $perm_args{$argkey}) {
	    croak("ARG. ERROR: $argkey isn't a valid arg. for prune_cluster()");
	}
	else {
	    my $argvar = $args_href->{$argkey};
	    unless (ref($argvar) eq $perm_args{$argkey}) {
		croak("ARG. ERROR: $argvar isn't valid arg. prune_cluster()");
	    }
	}
    }

    ## Now it will get the data and extract the alignments
    ## and run the filters in the order by_simplealign_function,
    ## by_strainpairs_distance and by_strain_count

    ## Store the clusters and members removed into a hash with the reason
    ## for that
    
    my %rmclusters = ();
    my %rmmembers = ();

    my %clusters = %{$self->get_clusters()};
    foreach my $clid (keys %clusters) {
	
	my $seqfam = $clusters{$clid};
	my $align = $seqfam->alignment();

	## First by_simplealign_function, and store the removed clusters
	## to modify the object later.

	if (defined $args_href->{'by_simplealign_function'}) {
	    my %args_simaln = %{$args_href->{'by_simplealign_function'}};
	    
	    foreach my $func (keys %args_simaln) {

		## If the function value is less than the constrains
		## add to the hash

		if ($align->$func < $args_simaln{$func}) {
		    $rmclusters{$clid} = $func;
		}
	    }	    
	}

	## Second detail the constrains
	


	


    }

	
}


####
1; #
####
