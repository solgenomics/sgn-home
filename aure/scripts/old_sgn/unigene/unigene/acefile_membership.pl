#!/usr/bin/perl -w
use strict;
use CXGN::BioTools::FastaParser;

if (!$ARGV[0] || $ARGV[0] eq "help") {
    print STDERR <<EOF;

    Script to parse up the output of an SGN unigene assembly. Reads the 
    output directory for .ace files, CAP3 singlets, and preclustering
    singlets.

    Usage: <output directory name>

EOF
    exit(0);
}
{
my ($output_directory) = @ARGV;

open FINDPIPE, "find $output_directory/ -name \"*.ace\" |" 
  or die "Failed to open pipe for find command ($!)";

foreach ( <FINDPIPE> ) {
    my ($contig_name, $contig_no);
    if (! m/cluster-([0-9]+).+\.cap\.ace/) {
	die "File \"$_\" unrecognized\n";
    }
    my $cluster_no = $1;

    open ACEFILE, "<$_"
	or die "Can't open ace file \"$_\" ($!)\n";

    my $line = <ACEFILE>;
    die "Not ACE format $_\n" if ($line !~ m/^AS/);
    my (undef, $n_contigs) = split/\s/,$line;
    if ($n_contigs == 0) {
      close ACEFILE;
      next;
    }

    while(<ACEFILE>) {	  
	if (m/^CO/) {
	    (undef, $contig_name) = split;
	    ($contig_no) = $contig_name =~ m/Contig([0-9]+)/;
	    last;
	}
    }

    my %components = ();

  again:

    foreach my $k (keys (%components)) {
      delete ($components{$k});
    }

    while(<ACEFILE>) {
	my ($readname, $ort, $start_pos, $padded_length, $first_qbase, $last_qbase);
	if (m/^AF/) {
	    (undef, $readname, $ort, $start_pos) = split;
	    if (defined($components{$readname})) {
		warn "Warning: component $readname already read from file\n";
	    }
	    if ($ort eq "C") {
		$ort = "-";
	    } else {
		$ort = "+";
	    }
	    $components{$readname} = { ort => $ort, start_pos => $start_pos };
	}

	if (m/^RD/) {
	    (undef, $readname, $padded_length) = split;
	    if (! defined($components{$readname})) {
		warn "Found RD line for component $readname, but there was no AF line\n";
	    }
	    $components{$readname}->{padded_length} = $padded_length;
	    while(<ACEFILE>) {
		if (m/^QA/) {
		    (undef, undef, undef, $first_qbase, $last_qbase) = split;
		    $components{$readname}->{start_trim} = $first_qbase - 1;
		    $components{$readname}->{start_pos} += $first_qbase - 1;
		    $components{$readname}->{end_trim} =
			$padded_length - $last_qbase;
		    $components{$readname}->{end_pos} = 
		      $components{$readname}->{start_pos} +
			($last_qbase - $first_qbase);
		    last;
		}
	    }
	}

	if (m/^CO/) {
	    my $line = $_;
	    dump_output($cluster_no, $contig_no, \%components);  
	    (undef, $contig_name) = split/\s+/,$line;
	    ($contig_no) = $contig_name =~ m/Contig([0-9]+)/;
	    goto again;
	}
    }

#End of File condition, process last contig
    close ACEFILE;
    dump_output($cluster_no, $contig_no, \%components);
}

close FINDPIPE;

# Read unigene.seq file for singlets.
my ($seq_id);
my $ufh =
  CXGN::BioTools::FastaParser->load_file("unigene.seq");

foreach $seq_id ( $ufh->get_seqids() ) {
  if ($seq_id =~ m/^([X0-9]+)-S-(.+)|(SGN-E\d+)/) {
      my $cluster_no = $1 ? $1 : 0;
      my $est_id = $2 ? $2 : $3;
      
    my $seq_length = length($ufh->get_sequence($seq_id));
    if ($cluster_no eq "X") {
      print "Singlet:\t",join("\t",-1,-1,1,$seq_length,$seq_length),"\n";
      print join("\t",$est_id,"+",1,$seq_length,0,0),"\n";
    } else {
      print "Singlet:\t",join("\t",$cluster_no,-1,1,$seq_length,$seq_length),"\n";
      print join("\t",$est_id,"+",1,$seq_length,0,0),"\n";
    }
  }
}

}

sub dump_output {
    my ($cluster_no, $contig_no, $component_ref) = @_;

    my @reads = sort { $component_ref->{$a}->{start_pos} <=> $component_ref->{$b}->{start_pos} } keys %{$component_ref};

    my ($max_e, $min_s, $max_l) = (0,0,0);
    foreach ( @reads ) {
	my ($s, $e, $l);
	$s = $component_ref->{$_}->{start_pos}-$component_ref->{$_}->{start_trim};
	$e = $component_ref->{$_}->{end_pos}+$component_ref->{$_}->{end_trim};
	$l = $component_ref->{$_}->{end_pos};

	$min_s = $s if ($s < $min_s);
	$max_e = $e if ($e > $max_e);
	$max_l = $l if ($l > $max_l);
    }

    my $n_components = @reads;

    # untrimmed length is untrimmed end - untrimmed start 
    # Format: cluster# contig# n_components untrimmed_length trimmed_length
    print "Contig:\t",join("\t",$cluster_no,$contig_no,$n_components,$max_e-$min_s,
	       $max_l),"\n";
    foreach ( @reads ) {
	print join("\t",$_,$component_ref->{$_}->{ort},$component_ref->{$_}->{start_pos},$component_ref->{$_}->{end_pos}, $component_ref->{$_}->{start_trim},$component_ref->{$_}->{end_trim}),"\n";
    }
}





# QA padded start padded end
# Length used = padded end - padded start
# CO Aligned END = ALign Start + Length used
# Trim end = padded_length - padded_end
