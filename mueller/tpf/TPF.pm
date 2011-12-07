
package TPF;

use Moose;

use File::Slurp;

has 'tpf' => (is => 'rw',
	      isa => 'ArrayRef'
    );

has 'file' => (is => 'rw',
	       isa => 'Str'
    );

sub parse { 
    my $self = shift;
    my $file = shift || $self->file;
    my @tpf = read_file($file);

    my @parsed;
    foreach (@tpf) { 
	chomp;
	push @parsed, [ split /\t/ ];
    }
    $self->tpf(\@parsed);
}

=head2 remove

 Usage:        $tpf->remove($name1, $name2);
 Desc:         removes the elements between $name1 and $name2
               ($name1 and $name2 remain in the file)
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub remove { 
    my $self = shift;
    my $start_id = shift;
    my $end_id = shift;

    $self->replace($start_id, $end_id);
}

=head2 insert

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub insert { 
    my $self = shift;
    my $start_id = shift;
    my @group = @_;

    my $start_index = $self->index($start_id);
    $self->replace_coords($start_index, $start_index+1, @group);
}

=head2 replace

 Usage:        $tpf->replace($start_id, $end_id, @lines_to_insert)
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub replace { 
    my $self = shift;
    my $start_id = shift;
    my $end_id = shift;
    my @new_tpf_lines = @_;

    my $start_index = $self->index($start_id);
    my $end_index = $self->index($end_id);

    if (!defined($start_index) || !defined($end_index)) { 
	warn "$start_id or $end_id not found. Sorry.";
	warn "The following lines could not be integrated: ".(join ", ", map { $_->[0] } @new_tpf_lines)."\n";
    }
    
    return $self->replace_coords($start_index, $end_index, @new_tpf_lines);
    
    

}

sub replace_coords { 
    my $self = shift;
    my $start_index = shift;
    my $end_index = shift;
    my @new_tpf_lines = @_;

    if ($start_index > $end_index) { 
	warn "start_index must be smaller than end_index";
        return 1;
    }

    my $tpf = $self->tpf();

    my $length = $end_index -$start_index -1;

    if ($length < 0 ) { 
	warn "Length is negative! not proceeding.\n";
	$length = 0; 

    }

    splice @$tpf, $start_index+1, $length, @new_tpf_lines;

    return 0;
}

=head2 index

 Usage:        my $index = $tpf->index($element_id)
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub index { 
    my $self = shift;
    my $element_id = shift;
    
    if (!defined($element_id)) { return undef; }
    my $tpf = $self->tpf();
    for(my $i=0; $i<@$tpf; $i++) {
	if (defined($tpf->[$i]) && defined($tpf->[$i]->[0]) &&  ($tpf->[$i]->[0] eq $element_id)) { 
	    return $i;
	}
    }
    return undef;
}


=head2 interval

 Usage:        my @tpf = $tpf->interval($start_id, $end_id);
 Desc:         returns the section of the tpf file between the given ids
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub interval { 
    my $self = shift;
    my $start_id = shift;
    my $end_id = shift;
    
    my $start_index = $self->index($start_id);
    my $end_index   = $self->index($end_id);

    if ((!defined($start_index) || !defined($end_index)) || 
	($start_index > $end_index)) { 
	warn "start_index=$start_index ($start_id), end_index=$end_index ($end_id) does not define useable interval.";
	return undef;
    }
    
    my $tpf = $self->tpf();

    my @interval;
    foreach my $n ($start_index .. $end_index) { 
	push @interval, $tpf->[$n];
    }
    return @interval;
}

=head2 interval_ids

 Usage:
 Desc:         returns a list of identifiers in the interval
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub interval_ids { 
    my $self = shift;
    my $start_id = shift;
    my $end_id = shift;
    my @interval = $self->interval($start_id, $end_id);
    return map { $_->[0] } @interval;
}

=head2 size

 Usage:        my $s = $tpf->size()
 Desc:         returns the current size, in lines, of the currently stored
               tpf file
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub size { 
    my $self = shift;
    return scalar(@{$self->tpf()});
}

=head2 render

 Usage:
 Desc:
 Ret:
 Args:
 Side Effects:
 Example:

=cut

sub render { 
    my $self = shift;
    foreach my $t (@{$self->tpf()}) { 
	print join "\t", @$t;
	print "\n";
    }
}

1;
