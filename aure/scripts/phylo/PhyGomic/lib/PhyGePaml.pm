
package PhyGePaml;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Try::Tiny;
use Math::BigFloat;

use File::Temp qw/ tempfile tempdir/;

use Bio::Tools::Run::StandAloneBlast;
use Bio::Tools::Run::Phylo::PAML::Codeml;

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
	seqfams        => {},
	strains        => {},
	cds            => {},
	topotypes      => {},
	codeml_results => {},
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
    $self->set_codeml_results($args_href->{codeml_results});
    $self->disable_reportstatus();
    
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

=head2 get_codeml_results/set_codeml_results

  Usage: my $codeml_results_href = $phygetopo->get_codeml_results();
         $phygetopo->set_codeml_results($codeml_results_href); 

  Desc: Get or set a hashref with key=ID and value=Bio::Matrix::Generic object
        Bio::Matrix::Generic entries are hashrefs. with the followings keys:
          S, N, kappa, dS, lnL, dN, omega, t

  Ret: Get: a hashref with keys=ID and value=Bio::Matrix::Generic
       Set: None

  Args: Get: None
        Set: a hashref with keys=ID and value=Bio::Matrix::Generic

  Side_Effects: Die if no argument is used, if the argument is not a hashref
                or if its values are not Bio::Matrix::Generic objects.

  Example:  my $codeml_results_href = $phygetopo->get_codeml_results();
            $phygetopo->set_codeml_results($codeml_results_href); 

=cut

sub get_codeml_results {
    my $self = shift;
    return $self->{codeml_results};
}

sub set_codeml_results {
    my $self = shift;
    my $codeml_res_href = shift;
   
    unless (defined $codeml_res_href) {
	croak("ARG. ERROR: No arg. was used for set_codeml_results function");
    }
    else {
	unless (ref($codeml_res_href) eq 'HASH') {
	    croak("ARG. ERROR: When arg isnt a hashref for set_codeml_results");
	}
	else {
	    my %hash = %{$codeml_res_href};
	    my $obj = 'Bio::Matrix::Generic';
	    foreach my $key (keys %hash) {
		unless (ref($hash{$key}) eq $obj) {
		    croak("ARG. ERROR: Values arent $obj set_codeml_results()");
		}
	    }
	}
	$self->{codeml_results} = $codeml_res_href;    
    }
}

################
### SWITCHES ###
################

=head2 enable/disable_reportstatus

  Usage: $phygepaml->enable_reportstatus();
         $phygepaml->disable_reportstatus();

  Desc: Enable/disable the switch report status

  Ret: none
 
  Args: None
 
  Side_Effects: None
 
  Example: $phygepaml->enable_reportstatus();
           $phygepaml->disable_reportstatus();

=cut 

sub enable_reportstatus {
    my $self = shift;
    $self->{reportstatus} = 1;
}

sub disable_reportstatus {
    my $self = shift;
    $self->{reportstatus} = 0;
}

=head2 is_on_reportstatus

  Usage: my $boolean = $phygepaml->is_on_reportstatus();

  Desc: Get 0 or 1 for reportstatus switch

  Ret: a boolean, 0 or 1.
 
  Args: None
 
  Side_Effects: None
 
  Example: if ($phygepaml->is_on_reportstatus() ) {
               print STDERR "REPORT";
           }

=cut 

sub is_on_reportstatus {
    my $self = shift;
    return $self->{reportstatus};
}


=head2 print_run_status

  Usage: print_run_status($progress, $total, $message, $id);

  Desc: Print as STDERR the percentage of the process run

  Ret: none

  Args: $progress, a scalar, an integer with the progress
        $total, a scalar, a integer with the total
        $message to print before the percentage
        $id, to print after the percetage

  Side_Effects: Die if 1st and 2nd arguments are undefined or are not integers
                Print as default message: Percentage of the file parsed: ";

  Example: print_run_status($progress, $total, $message);

=cut

sub print_run_status {
    my $a = shift;

    unless (defined $a) {
        croak("ERROR: 1st argument is not defined for print_parsing_status()");
    }
    my $z = shift;
    
    unless (defined $z) {
        croak("ERROR: 2nd argument is not defined for print_parsing_status()");
    }
        
    my $message = shift ||
        "Percentage of the file parsed:";

    my $id = shift || 'NA';

    unless ($a =~ m/^\d+$/) {
        croak("ERROR: 1st argument is not an int. for print_parsing_status()");
    }
    unless ($z =~ m/^\d+$/) {
        croak("ERROR: 2nd argument is not an int. for print_parsing_status()");
    }
    
    my $perc = $a * 100 / $z;
    my $perc_obj = Math::BigFloat->new($perc);
    my $perc_formated = $perc_obj->bfround(-2);

    print STDERR "\t$message $perc_formated %    \t(processing:$id)       \r";
}



##########################
## ANALYTICAL FUNCTIONS ##
##########################

=head2 translate_cds

  Usage: my %protein_objs = $phygepaml->translate_cds();

  Desc: Translate the cds to protein with the ORF +1
        
  Ret: %protein_objs, a hash ref. with key=id and values=Bio::Seq object. 

  Args: None

  Side_Effects: None

  Example: my %protein_objs = $phygepaml->translate_cds();

=cut

sub translate_cds {
    my $self = shift;
    
    my %proteins = ();
    my %cds = %{$self->get_cds()};

    foreach my $cds_id (keys %cds) {
	my $prot = $cds{$cds_id}->translate();
	$proteins{$cds_id} = $prot;
    }
    return %proteins;
}

=head2 align_member_cds

  Usage: my %seqfam_cds_objs = $phygepaml->align_member_cds();

  Desc: Using cds as reference, align it with the consensus for the family
        get the coordinates and modify the members of the family to get its cds.
        
  Ret: %seqfam_cds_objs, a hash ref. with key=id and 
                                          values=Bio::Cluster::SequenceFamily

  Args: None

  Side_Effects: None

  Example: my %seqfam_cds_objs = $phygepaml->align_member_cds();

=cut

sub align_member_cds {
    my $self = shift;
    my $params = shift;
    
    my @params;
    if (defined $params) {
	if (ref($params) ne 'HASH') {
	    croak("ERROR: Parameters used for align_member_cds isnt a HASHREF");
	}
	else {
	    foreach my $p (keys %{$params}) {
		push @params, ($p, $params->{$p});
	    }
	}
    }
    else {
	@params = ('-program', 'blastx', '-W', 2);
    }

    my %new_seqfam = ();

    my %cds = %{$self->get_cds()};
    my %seqfam = %{$self->get_seqfams()};

    ## Create a blast factory to run bl2seq

    my $blast_fc = Bio::Tools::Run::StandAloneBlast->new(@params);

    ## First get the consensus for each of the sequences and compare with
    ## the cds.

    my $new_align;
    
    foreach my $id (sort keys %seqfam) {
	
	my $aln = $seqfam{$id}->alignment();
	if (defined $aln && defined $cds{$id}) {
	    
	    my $cons = $aln->consensus_string();

	    ## If the cds ends in a stop codon remove it.

	    my $pepseq = $cds{$id}->translate()->seq();
	    if ($pepseq =~ m/\*$/i) {
		$pepseq =~ s/\*$//;
	    }

	    my $cds_obj = Bio::Seq->new( -id       => $id, 
					 -seq      => $pepseq, 
					 -alphabet => 'protein');
	    my $cons_obj = Bio::Seq->new( -id       => $id . '_cons', 
					  -seq      => $cons,
					  -alphabet => 'dna',
		);

	    ## It can give problems if the seq lenght is to small

	    if ($cds_obj->length() > 10 && $cons_obj->length() > 30) {

		my $blast_rep = $blast_fc->bl2seq($cons_obj, $cds_obj);

		## Take the first result

		my @blast_results = $blast_rep->next_result;
		my @blast_hits = $blast_results[0]->next_hit;
		if (scalar(@blast_hits) > 0) {
		    my @blast_hsps = $blast_hits[0]->next_hsp;
		    my $que_st = $blast_hsps[0]->start('query');
		    my $que_en = $blast_hsps[0]->end('query');
		    
		    $new_align = $aln->slice($que_st, $que_en);
		}
	    }
	}

	## Second create a new family with the new alignment if it is defined

	if (defined $new_align) {
	    my $new_fam = Bio::Cluster::SequenceFamily->new(
		-family_id => $id,
		-members   => [$seqfam{$id}->get_members]
		);
	    $new_fam->alignment($new_align);
	    $new_seqfam{$id} = $new_fam;
	}
    }

    return %new_seqfam;
}


=head2 set_seqfam_cds

  Usage: $phygepaml->set_seqfam_cds($bl2seq_parameters);

  Desc: Replace the seqfams alignments with the alignments of the cds of the
        members.
        
  Ret: None

  Args: $bl2seq_parameters, a hashref. with the parameters to use for the
                            alignment between the original consensus of the
                            seqfam and the cds.

  Side_Effects: Die if the $bl2seq_parameters isnt a hashref.

  Example: $phygepaml->set_seqfam_cds();

=cut

sub set_seqfam_cds {
    my $self = shift;
    my $params = shift;

    my %seqfam_cds_objs = $self->align_member_cds($params);
    $self->set_seqfams(\%seqfam_cds_objs);
}

=head2 predict_cds

  Usage: $phygepaml->predict_cds($parameters);

  Desc: Predict the cds of an alignment using the consensus sequence.
        Some methods are available: 
         -longest6frame
         -estscan
         -proteinalign
        
  Ret: None

  Args: $parameters, a hashref. with key=method/arguments and 
                                     value=method_string/argument_href

  Side_Effects: Die if some arguments are wrong

  Example: $phygepaml->predict_cds($parameters);

=cut

sub predict_cds {
    my $self = shift;
    my $params = shift ||
	croak("ERROR: No parameters were specified for predict_cds()");

    ## Check variables

    if (ref($params) ne 'HASH') {
	croak("ERROR: Parameters specified for predict_cds is not a hashref.");
    }
    my %params = %{$params};
    
    my %perm_params = (
	method    => '(longest6frame|estscan|proteinalign)',
	arguments => {},
	);

    foreach my $pa (keys %params) {
	unless (exists $perm_params{$pa}) {
	    croak("ERROR: $pa isnt a permited parameter for predict_cds");
	}
	else {
	    if (ref($perm_params{$pa})) {
		if (ref($params{$pa}) ne ref($perm_params{$pa})) {
		    croak("ERROR: $pa parameters has a non permited value.");
		}
	    }
	    else {
	    	if ($params{$pa} !~ m/^$perm_params{$pa}$/) {
		    croak("ERROR: $pa parameters has a non permited value.");
		}
	    }
	}
    }
    unless (exists $params{arguments}) {
	$params{arguments} = {};
    }

    ## Now it will predict the cds for each contigs
    
    my %cds = ();

    my %seqfams = %{$self->get_seqfams()};
    foreach my $id (sort keys %seqfams) {
	my $aln = $seqfams{$id}->alignment();
	my $consensus = $aln->consensus_iupac();
	my $conseq = Bio::Seq->new( -id       => $id, 
				    -seq      => $consensus, 
				    -alphabet => 'dna' );

	my ($cds_seq, $frame, $start, $end);

	if ($params{method} eq 'estscan') {  
	   ## $cds_seq = _cds_by_estscan($conseq, $params{arguments});
	}
	elsif ($params{method} eq 'proteinalign') {   
	    $cds_seq = _cds_by_proteinalign($conseq, 
					    $params{arguments}->{blastdb});
	}
	else {  ## By default
	    ($cds_seq, $frame, $start, $end) = _cds_by_longest6frame(
		$conseq,
		$params{arguments}->{force_firstmet},
		); 
	}

	if (defined $cds_seq) {
	    $cds{$id} = $cds_seq;
	}
    }
    
    $self->set_cds(\%cds);
}

=head2 _cds_by_longest6frame

  Usage: my ($cds, $frame, $start, $end) = _cds_by_longest6frame($seq);

  Desc: Calculate the 6 ORF for a sequences and get the longest sequence.
        If more than one ORF have the longest sequence, it will take the
        first (see translation function at Bio::PrimarySeqI for more details)
        If force_firstmet is enabled it will get the longest6frame where appears
        a first metionine.        

  Ret: $cds, a Bio::Seq object with the predicted cds,
       $frame, frame for the predicted cds
       $start, start coord. for the predicted cds refering to the original seq.
       $end, end coord. for predicted cds refering to the original seq.
    
  Args: $seq, a Bio::Seq object
        $force_firstmet, a boolean.        

  Side_Effects: Die if some arguments are wrong

  Example: my $cds = _cds_by_longest6frame($seq);

=cut

sub _cds_by_longest6frame {
    my $seq = shift ||
	croak("ERROR: No seq. object was supplied to _cds_by_longest6frame");
    my $firstmet = shift || 0;

    if ($firstmet !~ m/^(1|0)$/) {
	croak("ERROR: Wrong value for force_firstmet. It only can be boolean");
    }

    if (ref($seq) !~ m/Bio::\w*Seq/) {
	croak("ERROR: Argument supplied to _cds_by_longest6frame isnt seq obj");
    }

    my $longest_orf;

    ## Create the 6 ORFs (orf1 = +1, orf2 = +2, orf3 = +3, orf4 = -1, 
    ## orf5 = -2, orf6 = -3)

    my %frames = (
	orf1 => 1,
	orf2 => 2,
	orf3 => 3,
	orf4 => -1,
	orf5 => -2,
	orf6 => -3,
	);

    my %orfs = (
	orf1 => $seq->translate(),
	orf2 => $seq->translate( -frame => 1),
	orf3 => $seq->translate( -frame => 2),
	orf4 => $seq->revcom->translate(),
	orf5 => $seq->revcom->translate( -frame => 1),
	orf6 => $seq->revcom->translate( -frame => 2),
	);

    ## Evaluate fragments
    
    my %fragms = ();
    my %fragms_seq = ();
    
    foreach my $orf (sort keys %orfs) {

	my @frags = split(/\*/, $orfs{$orf}->seq());

	## That will remove the stop codon from the sequence.

	my $st = 1;
	my $en = 0;

	my $n = 0;
	foreach my $fragm (@frags) {

	    $n++;
	    ## Find the first metionine

	    my $met;
	    if ($fragm =~ m/M/) {
		$met = index($fragm, 'M');
		$st += ($met + 1) * 3 - 3;
		
		$fragm = substr($fragm, $met);
	    }
	    
	    my $len = length($fragm);
	    $en = $st + $len * 3 - 1;

	    if ($n < scalar(@frags)) { ## If it is not the last add the stop
		$en += 3;
	    }
	    
	    my $frag_id = $orf . ':' . $st . '-' . $en;
	    if ($len > 0) {
		$fragms{$frag_id} = $len;
		$fragms_seq{$frag_id} = $fragm;
	    }

	    $st = $en + 3;
	}
    }

    ## Order fragments

    my @frags_ordered = ();
    
    my @frags = sort { $fragms{$b} <=> $fragms{$a} } keys %fragms;
    if ($firstmet == 1) {
	
	foreach my $frag_id (@frags) {
	    if ($fragms_seq{$frag_id} =~ m/^M/) {
		push @frags_ordered, $frag_id;
	    }
	}
    }
    else {
	@frags_ordered = @frags;
    }

    ## Get the longest fragment and return the data

    my $longest = $frags_ordered[0];
    
    if (defined $longest && $longest =~ m/(orf\d+):(\d+)-(\d+)/) {
	my $orf_sel = $1;
	my $start = $2;
	my $end = $3;
	my $cds;
    
	if ($frames{$orf_sel} > 0) {
	    $cds = $seq->trunc($start, $end);
	}
	else {
	    $cds = $seq->revcom()->trunc($start, $end);
	}
    
	return ($cds, $frames{$orf_sel}, $start, $end);
    }
    else {
	return ();
    }
}

=head2 _cds_by_proteinalign

  Usage: my $cds = _cds_by_proteinalign($seq, $protein_blast_database);

  Desc: Calculate the CDS of a sequence based in the match with a know protein.
        Add N to the cds sequence when some gaps are finded in the match.

  Ret: $cds, a Bio::Seq object with the predicted cds,
    
  Args: $seq, a Bio::Seq object,
        $protein_blast_database, a protein blast database name.

  Side_Effects: Die if some arguments are wrong

  Example: my $cds = _cds_by_proteinalign($seq);

=cut

sub _cds_by_proteinalign {
    my $seq = shift ||
	croak("ERROR: No seq. object was supplied to _cds_by_proteinalign");
    my $protein_blastdb = shift ||
	croak("ERROR: No protein blastdb was supplied to _cds_by_proteinalign");

    if (ref($seq) !~ m/Bio::\w*Seq/) {
	croak("ERROR: Argument supplied to _cds_by_proteinalign isnt seq obj");
    }

    my $cds;
    
    return $cds;
}


##################
### PAML TOOLS ###
##################

=head2 run_codeml

  Usage: $phygepaml->run_codeml($parameters_href);

  Desc: Run codeml program (from PAML) over the aligments of the phygepaml
        object.

  Ret: None.
    
  Args: $parameters_href, a hashref with key=parameter, value=value to run 
        codeml (see Bio::Tools::Run::Phylo::PAML::Baseml for more details)

  Side_Effects: Die if some arguments are wrong

  Example: $phygepaml->run_codeml($parameters_href);

=cut

sub run_codeml {
    my $self = shift;
    my $params = shift;

    if (defined $params && ref($params) ne 'HASH') {
	croak("ERROR: $params supplied to run_codeml isnt a hashref.");
    }

    my %seqfams = %{$self->get_seqfams()};
    my %paml = ();
    my $t = scalar(keys %seqfams);
    my $a = 0;


    foreach my $seqfam_id (sort keys %seqfams) {
	
	my $aln = $seqfams{$seqfam_id}->alignment();
	my $tree = $seqfams{$seqfam_id}->tree();

	$a++;
	if ($self->is_on_reportstatus) {
	    print_run_status($a, $t, "\t\tPerc. of codeml run:", $seqfam_id );
	}

	## PAML give problems with the stop codon, so it will check first
	## if there are any problems with that

	my $stopcodon = 0;
	foreach my $seqmem ($aln->each_seq()) {
	    
	    my $pepseq = $seqmem->translate()->seq();
	    if ($pepseq =~ m/\*/) {		
		$stopcodon = 1;
	    }	    
	}
	
	if (defined $aln && $stopcodon == 0) {

	    my $codeml = Bio::Tools::Run::Phylo::PAML::Codeml->new();

	    ## Pass the parameters

	    foreach my $par (keys %{$params}) {
		$codeml->set_parameter($par, $params->{$par});
	    }

	    ## To have the right parser verbose must to be 0, 
	    ## So this line will overwrite that

	    $codeml->set_parameter('verbose', 0);
	    $codeml->alignment($aln);

	    my ($rc, $parser) = $codeml->run();

	    if (defined $rc && $rc == 1) { ## Success
		my $result = $parser->next_result;

		if (defined $result) {

		    my $MLmatrix = $result->get_MLmatrix();

		    ## The sequences in the alignment are in the same order than
		    ## in the input/output PAML files

		    my @members = $aln->each_seq();
		    my @seqnames = ();

		    foreach my $seq (@members) {
			my $seqid = $seq->display_id();
			push @seqnames, $seqid;
		    }

		    ## The matrix is symmetric, but it is incomplete
		
		    my $i = 0;
		    foreach my $iname (@seqnames) {
		    
			my $j = 0;
			foreach my $jname (@seqnames) {
		    
			    unless (defined $MLmatrix->[$i]->[$j]) {
				$MLmatrix->[$i]->[$j] = $MLmatrix->[$j]->[$i];
			    }
			    $j++;
			}
		    $i++;
		    }
       
		    my $mtx = Bio::Matrix::Generic->new(
			-values   => $MLmatrix,
			-rownames => \@seqnames,
			-colnames => \@seqnames,
			-matrix_name => $seqfam_id,
			);

		    $paml{$seqfam_id} = $mtx;
		}
		else {
		    my $error = $codeml->error_string();
		    warn("\ncodeml fails:\n$error\n");
		}
	    }
	    else {
		my $error = $codeml->error_string();
		warn("\ncodeml fails:\n$error\n");
	    }
	}
    }
    
    ## Finally it will set the results into the object

    $self->set_codeml_results(\%paml);
}

=head2 out_codeml

  Usage: $phygepaml->out_codeml($filename, $header);

  Desc: Produce an output file with the following columns:
        -f1:  cluster_id
        -f2:  topology_id
        -f3:  pair_str1
        -f4:  pair_str2
        -f5:  pair_id1
        -f6:  pair_id2
        -f7:  N
        -f8:  S
        -f9:  dN
        -f10: dS
        -f11: omega
        -f12: kappa
        -f13: lnL
        -f14: t

  Ret: None.
    
  Args: $filename, a filename to print the data
        $header, aboolean to print/no print a header

  Side_Effects: Die no filename is supplied

  Example: $phygepaml->out_codeml($filename);

=cut

sub out_codeml {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename was supplied to out_codeml function.");

    my $header = shift || 1;
    if ($header !~ m/0/) { ## Overwrite different things than a boolean
	$header = 1;
    }


    my %seqfams = %{$self->get_seqfams()};
    my %strains = %{$self->get_strains()};
    my %topologies = %{$self->get_topotypes()};
    my %codeml = %{$self->get_codeml_results()};

    ## convert topologies

    my %fams2topo = ();
    foreach my $topology_id (sort keys %topologies) {
	my @members = @{$topologies{$topology_id}->get_members()};
	foreach my $members (@members) {
	    $fams2topo{$members->id()} = $topology_id;
 	}
    }

    ## open the file to print the data

    open my $ofh, '>', $filename;
    
    if ($header == 1) {
	print $ofh "CL_ID\tTOPOLOGY_ID\tSTR1\tSTR2\tSEQ1\tSEQ2\tN\tS\tdN\tdS\t";
	print $ofh "OMEGA\tKAPPA\tLnL\tT\n";
    }

    foreach my $seqfam_id (sort keys %seqfams) {

	my $topo_id = $fams2topo{$seqfam_id} || 'NA';	
	
	## Get the members
    
	my $mtx = $codeml{$seqfam_id};

	if (defined $mtx) {
	
	    my @names = $mtx->column_names();
	    
	    ## It will order the data by strains

	    my %bystr = ();
	    foreach my $name (@names) {
		$bystr{$name} = $strains{$name};
	    }
	    my @ord_names = sort {$bystr{$a} cmp $bystr{$b}} keys %bystr;

	    my %outdata = ();  ## To remove the redundancy
	    foreach my $oname1 (@ord_names) {
		foreach my $oname2 (@ord_names) {
		    
		    if ($oname1 ne $oname2) {

			my $id = join(',', sort ($oname1, $oname2));
			unless (exists $outdata{$id}) {
			    print $ofh "$seqfam_id\t$topo_id\t";
			    print $ofh "$strains{$oname1}\t$strains{$oname2}\t";
			    print $ofh "$oname1\t$oname2\t";
			    
			    my %entry = %{$mtx->entry($oname1, $oname2)};
			    print $ofh "$entry{N}\t$entry{S}\t";
			    print $ofh "$entry{dN}\t$entry{dS}\t";
			    print $ofh "$entry{omega}\t$entry{kappa}\t";
			    print $ofh "$entry{lnL}\t$entry{t}\n";
			    $outdata{$id} = 1;
			}		    
		    }
		}
	    }

	}
	else {
	    print $ofh "$seqfam_id\t$topo_id\t";
	    print $ofh "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\n";
	}
    }
}


####
1; #
####
