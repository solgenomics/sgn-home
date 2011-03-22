
package Strain::Composition;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;



###############
### PERLDOC ###
###############

=head1 NAME

Strain::Composition.pm
An object to count and check strain compositions.

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use Strain::Composition;

  my $comp = Strain::Composition->new({ 
                                        strains     => \%strains,
                                        composition => \%composition,
                                        members     => \@members,
  });


  ## Accessors

  my $strains_href = $comp->get_strains();
  $comp->set_strains($strains_href);

  my $composition_href = $comp->get_composition();
  $comp->set_composition($composition_href);

  my $members_href = $comp->get_members();
  $comp->set_members($members_href);

  ## Tools

  unless ($comp->is_complete()) {
     my $sucess = $comp->add_member($member_id);
  }

  $comp->reset_members();


=head1 DESCRIPTION

 Object to manage the strain composition.


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $comp = Strain::Composition->new($arguments_href);

  Desc: Create a Strain::Composition object with the specified parameters.

  Ret: a Strain::Composition.pm object

  Args: A hash reference with the following key-value pairs: 
         strains     => a hashref. with key=member_id, value=strains
         composition => a hashref. with key=strain, value=amount
         members     => a hashref. with key=strain, value=arrayref. with
                        member_ids.
        
  Side_Effects: Die if the argument used is not a hashref. or if 
                the arguments used dont have the right format

  Example: my $composition = Strain::Composition->new({ strains => \%strains });

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class );                         
    
    my %permargs = ( 
	strains     => {},
	composition => {},
	members     => {},
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

    $self->set_strains($args_href->{strains});
    $self->set_composition($args_href->{composition});
    $self->set_members($args_href->{members});
    
    return $self;
}


#################
### ACCESSORS ###
#################

=head2 get/set_strains

  Usage: my $strains_href = $comp->get_strains(); 
         $comp->set_strains($strains_href);

  Desc: Get/Set a hashref with the members and the strains

  Ret: Get: a hashref with keys=member_id and value=strain
       Set: None

  Args: Get: None
        Set: A hash reference with key=member_id and value=strain

  Side_Effects: Get: None
                Set: Die if no argument is used.

  Example: my %strains = %{$comp->get_strains()};
           $comp->set_strains(\%strains);

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


=head2 get/set_composition

  Usage: my $composition_href = $comp->get_composition(); 
         $comp->set_composition($composition_href);

  Desc: Get/Set a hashref with the strains and composition

  Ret: Get: a hashref with keys=strain and value=strain member count
       Set: None

  Args: Get: None
        Set: A hash reference with key=strain and value=strain member count

  Side_Effects: Get: None
                Set: Die if no argument is used.

  Example: my %composition = %{$comp->get_composition()};
           $comp->set_composition(\%composition);

=cut

sub get_composition {
    my $self = shift;
    return $self->{composition};
}

sub set_composition {
    my $self = shift;
    my $comp_href = shift;
   
    unless (defined $comp_href) {
	croak("ARG. ERROR: No arg. was used for set_composition function");
    }
    else {
	unless (ref($comp_href) eq 'HASH') {
	    croak("ARG. ERROR: When arg is not a hashref. for set_composition");
	}
	else {
	    foreach my $str (keys %{$comp_href}) {
		if ($comp_href->{$str} !~ m/^\d+$/) {
		    croak("ARG. ERROR: value for $str isnt an integer.");
		}
	    }
	
	}
	$self->{composition} = $comp_href;    
    }
}


=head2 get/set_members

  Usage: my $members_href = $comp->get_members(); 
         $comp->set_members($members_href);

  Desc: Get/Set a hashref with the members

  Ret: Get: an hashref. with keys=strains and value=arrayref. with member_ids
       Set: None

  Args: Get: None
        Set:  an hashref. with keys=strains and value=arrayref. with member_ids

  Side_Effects: Get: None
                Set: Die if no argument is used.
                     Die if the strain used don't exists at composition.
                     Die if the hash values are not arrayrefs.
                     Die if more than the number of strain at composition
                     weer used in the array.
                     Die if one of the members doesnt exists into the strains
                     accessor.

  Example: my %members = %{$comp->get_members()};
           $comp->set_members(\%members);

=cut

sub get_members {
    my $self = shift;
    return $self->{members};
}

sub set_members {
    my $self = shift;
    my $members_href = shift;
   
    unless (defined $members_href) {
	croak("ARG. ERROR: No arg. was used for set_members function");
    }
    else {
	unless (ref($members_href) eq 'HASH') {
	    croak("ARG. ERROR: When arg is not a hashref. for set_members");
	}

	

	my %strains = %{$self->get_strains()};
	my %comp = %{$self->get_composition()};

	foreach my $str (sort keys %{$members_href}) {

	    ## Check that the strain used exists at composition accessor

	    unless (exists $comp{$str}) {
		croak("ERROR: $str used for set_members doesnt exist for comp");
	    }
	    else {
		
		## Check if the value is an arrayref.

		if (ref($members_href->{$str}) ne 'ARRAY') {
		    croak("ERROR: Value for $str isnt an arrayref.");
		}
		else {

		    ## Check if the members exists in the strain accessor.

		    my @members = @{$members_href->{$str}};

		    if (scalar(@members) > $comp{$str}) {
			my $err1 = "ERROR: More than $comp{$str} members where";
			$err1 .= " used for strain $str";
			croak($err1);
		    }

		    foreach my $member (@members) {
			unless (exists $strains{$member}) {
			    croak("ERROR: $member doesnt exists at strains");
			}
			else {
			    if ($strains{$member} ne $str) {
				my $err2 = "ERROR: Strain for member $member ";
				$err2 .= "isnt strain than the strain accessor";
				croak($err2);
			    }
			}
		    }
		}
	    }
	}

	$self->{members} = $members_href;    
    }
}

=head2 add_member

  Usage: my $str = $comp->add_member($member_id); 

  Desc: Add a new member to the member accessor

  Ret: $strain, strain for the member_id

  Args: $member_id, a member id

  Side_Effects: Die if no argument is used.
                Die if if the member doesnt exists into the strain accessor
                If the number of the members for a strain is equal to 
                composition it will not add the member and it will return undef.

  Example: my $str = $comp->add_member($member_id); 

=cut

sub add_member {
    my $self = shift;
    my $member_id = shift;

    my %strains = %{$self->get_strains()};
    unless (exists $strains{$member_id}) {
	croak("ERROR: $member_id doesnt exist into the strain accessor.");
    }
    else {
	unless ( $self->is_complete($strains{$member_id}) ) {
	    
	    my $members_aref = $self->get_members();
	    
	    if (exists $members_aref->{$strains{$member_id}}) {
		push @{$members_aref->{$strains{$member_id}}}, $member_id;
	    }
	    else {
		$members_aref->{$strains{$member_id}} = [$member_id];
	    }
	    
	    ## Reset with the same hashref. to re-check arguments
	    $self->set_members($members_aref);
	
	    return $strains{$member_id};
	}
	else {
	    return undef;
	}
    }
}

=head2 delete_members

  Usage: my %members = $comp->delete_members();
         my @members = $comp->delete_members($strain); 

  Desc: Delete all the members for all the strains, or just all the members for
        a concrete strain if a strain is used.

  Ret: %members, a hash with key=strains, values=arrayref. of members.
       @members, an array with members ids for a concrete strain.

  Args: $strain, a strain

  Side_Effects: Die if the strain used doesnt exists into composition accessor.
                Return empty if there are no members.

  Example: my %members = $comp->delete_members();
           my @members = $comp->delete_members($strain); 

=cut

sub delete_members {
    my $self = shift;
    my $strain = shift;

    my %composition = %{$self->get_composition()};
    my $members_href = $self->get_members();

    if (defined $strain) {
	
	my @del_members = ();
	unless (exists $composition{$strain}) {
	    croak("ERROR: Strain $strain doesnt exist into composition acc.");
	}
	
	## delete only if exists

	if (exists $members_href->{$strain}) {
	    my $members_aref = delete($members_href->{$strain});
	    @del_members = @{$members_aref};
	}

	return @del_members;
    }
    else {
	
	## Set to empty and return members hash

	$self->set_members({});
	return %{$members_href};
    }
}





#############
### TOOLS ###
#############

=head2 is_complete

  Usage: my $boolean = $comp->is_complete(); 
         my $boolean = $comp->is_complete($strain); 

  Desc: Check if exists the same number of members than the count number for
        composition accessor. If a strain is used, it will check only that 
        strain.

  Ret: $boolen, 0 FALSE and 1 TRUE

  Args: $strain, a strain name [optional]

  Side_Effects: Die if the strain used doesnt exists into composition accessor.

  Example: if ($comp->is_complete()) {
               print "Is complete.";
           }
 
=cut

sub is_complete {
    my $self = shift;
    my $strain = shift;

    ## declare boolean (by default it will be true)
    
    my $boolean = 1;
    
    my %composition = %{$self->get_composition()};
    my %memb = %{$self->get_members()};

    ## Check strain

    if (defined $strain) {
	
	unless (exists $composition{$strain}) {
	    croak("ERROR: Strain $strain doesnt exist into composition acc.");
	}
	
	if (exists $memb{$strain}) {
	    if (scalar(@{$memb{$strain}}) < $composition{$strain}) {
		$boolean = 0;
	    }
	}
	else {
	    $boolean = 0;
	}
    }
    else {
	
	foreach my $str (keys %composition) {
	    if (exists $memb{$str}) {
		if (scalar(@{$memb{$str}}) < $composition{$str}) {
		    $boolean = 0;
		}
	    }
	    else {
		$boolean = 0;
	    }
	} 
    }
    
    return $boolean;
}

=head2 count_members

  Usage: my %count_members = $comp->count_members();

  Desc: Count how many members there are by strain

  Ret: %count_members, a hash with key=strain, value=count

  Args: None

  Side_Effects: None

  Example: my %count_members = $comp->count_members();
 
=cut

sub count_members {
    my $self = shift;

    my %count_members = ();
    my %members = %{$self->get_members()};

    foreach my $strain (keys %members) {
	$count_members{$strain} = scalar(@{$members{$strain}});
    }
    return %count_members;
}

=head2 total_members

  Usage: my $total_members = $comp->count_members();

  Desc: Count how many members there are in total

  Ret: $total_members, total number of members in the object

  Args: None

  Side_Effects: None

  Example: my $total_members = $comp->total_members();
 
=cut

sub total_members {
    my $self = shift;
    
    my $total_members = 0;
    my %count_members = $self->count_members();
   
    foreach my $strain (keys %count_members) {
	$total_members += $count_members{$strain};
    }
   
    return $total_members;
}


####
1; #
####
