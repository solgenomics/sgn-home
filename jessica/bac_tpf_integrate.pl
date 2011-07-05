#!usr/bin/perl

use strict;
use warnings;


#similify bacs, only use the "best" ones, close to total length and good correctness

#produces a tpf file with bac's incorprated
	    
my $bac_file = shift;
my $agpchr = shift;
#my $tpfchr = shift;
my $tpfchr = "$agpchr.tpf";
#be sure the agp and tpf are named the SAME AND in the SAME location, to fit formatting above.

my %chr  = ( #chromosome ID translation.
    "CM001064.1" => "CHR1",
    "CM001065.1" => "CHR2",
    "CM001066.1" => "CHR3",
    "CM001067.1" => "CHR4",
    "CM001068.1" => "CHR5",
    "CM001069.1" => "CHR6",
    "CM001070.1" => "CHR7",
    "CM001071.1" => "CHR8",
    "CM001072.1" => "CHR9",
    "CM001073.1" => "CHR10",
    "CM001074.1" => "CHR11",
    "CM001075.1" => "CHR12"
    );

my $start_contig = 'start_contig';
my $end_contig = 'end_contig';
my $id_contig = 'id_contig';
my $chr_contig = 'chr_contig';


my %position_convert;
my @start_contig;
my @end_contig;
my @id_contig;
my @chr_contig;
my %chr_defign;
open(my $F,"<",$agpchr) || die "can't open $agpchr";
while(<$F>){
    chomp;
    if(/^$/ || /^\#/){next;} #skip blank and comment lines
    my ($obj, $start_obj, $end_obj, $id, $type, $comp_id, $length_agp, $orientation_agp) = split/\t/;
 
    my $chr_current = $chr{$obj} || '';
    if($start_obj && $comp_id && $comp_id =~m/^AEKE/){
	$comp_id=~s/(.*?)\.\d+/$1/; #remove version number
	$position_convert{$comp_id}=$start_obj;

	$chr_defign{$comp_id}=$chr_current;
    }
    if($start_obj && $end_obj && $comp_id && $chr_current){
    push (@start_contig, $start_obj || '');
    push (@end_contig, $end_obj || '');
    push (@id_contig, $comp_id || '');
    push (@chr_contig, $chr_current || '');
}	
    else{
	print STDERR "start_obj=$start_obj && end_obj$end_obj= && comp_id=$comp_id && chr_current=$chr_current\n";
    }
}

open($F,"<",$tpfchr) || die "can't open $tpfchr";
my %tpf_output_inital;
my %scaffold_convert;

while(<$F>){
    chomp;
    if(/^$/ || /^\#/){next;} #skip blank and comment lines
    my ($component, $data, $scaffold, $orientation_tpf, $blank) = split/\t/;

    if ($component && $component =~m/^AEKE/){  #exist, not gap (gaps do not have associated scaffold information)
	$scaffold_convert{$scaffold} = $component;
    }

    if($component =~m/^AEKE/){ #component as key (contig id)

$tpf_output_inital{$component} = [$component, $data, $scaffold, $orientation_tpf];
    }

    if($component eq 'GAP'){#scaffold as key (gap id)
	$tpf_output_inital{'s'.$scaffold} = [$component, $data, $scaffold, $orientation_tpf];
    }
}

my %bac_scaffold;
#my %rejects;
my $bac_begin;
my $bac_end;
my $bac_chr;
my $bac_position;
open(my $line, "<", $bac_file) || die "can't open $bac_file.";
while(<$line>){
    print STDERR ".";
    chomp;
    if(/^$/ || /^\#/){next;} #skip blank and comment lines
    my ($q_id, $s_id, $per_id, $len, $mismatch, $gaps, $q_start, $q_end, $s_start, $s_end, $e_val, $bit) = split /\t/;
    if (exists($scaffold_convert{$s_id}) && 
	exists($position_convert{$scaffold_convert{$s_id}})){

	$bac_position =($position_convert{$scaffold_convert{$s_id}} + $q_start);
	$bac_chr =  $chr_defign{$scaffold_convert{$s_id}};
	$bac_begin = $s_start + $q_start - 1;
	$bac_end = $s_end + $q_start - 1;

	my $i = 0;
	my $contig_count = @chr_contig;

	while($i < $contig_count){
#	      if($bac_chr eq $chr_contig[$i] && $bac_begin <= $end_contig[$i]){
		  #check bac has crossover
	    if($bac_begin >= $start_contig[$i] &&
	       $bac_begin <= $end_contig[$i]){#around bac_start
		if($id_contig[$i] =~m/^[0-9]/){#gap (contig start with 'A' and BACs start with 'C')
		    #alter to remove gap where bac overlaps
		    $end_contig[$i] = $bac_begin;
		}	
		#insert bac after obj in array
		$bac_scaffold{$q_id} = $s_id;
		splice(@start_contig,$i,0,$bac_begin || '');
		splice(@end_contig,$i,0,$bac_end || '');
		splice(@chr_contig,$i,0,$bac_chr || '');
		splice(@id_contig,$i,0,$q_id || '');
		# $i++ #to correctly calculate next contig/gap
	    }

	    elsif($bac_begin < $start_contig[$i] && 
		  $bac_end > $end_contig[$i]){#within bac
		#place obj in list (named after array)

#		$rejects{$q_id} = [ 
#		    $start_contig => $start_contig[$i+1],
#		    $end_contig => $end_contig[$i+1],
#		    $chr_contig => $chr_contig[$i+1],
#		    $id_contig => $id_contig[$i+1]
#		    ];
		
		#remove contig/gap from each array
		splice(@start_contig,$i+1,1);
		splice(@end_contig,$i+1,1);
		splice(@chr_contig,$i+1,1);
		splice(@id_contig,$i+1,1);
	    }
	    
	    elsif($bac_begin <= $start_contig[$i] && 
		  $bac_end <= $end_contig[$i]){
		if(!($id_contig[$i] =~m/^[0-9]/)){
		    #alter to remove gap where bac overlaps
		    $start_contig[$i] = $bac_end;
		}
		$i= $contig_count; #move on(next bac), ends while i
	    }
#####!!!!!!
##if not found check the rejects!!!Do we need to? does it check current bacs?
#####!!!!!!

	      }
	      if($bac_end <= $start_contig[$i]){#past bac
		  $i= $contig_count; #move on(next bac), ends while i
	      }
	      $i++; #moves to next contig/gap
	      
	}
    }
}

my $bac_orient;  #NEED TO SET/USE!!!!
my $id_size = @id_contig;
my $current_line = 0;
print STDERR "-";
while($id_contig[$current_line] && $current_line <= $id_size){#contig or gap
    print STDERR ":";
    if($tpf_output_inital{$id_contig[$current_line]}){#there is tpf for the id
#Use of uninitialized value in string eq at bac_integration/bac_chr_position.pl line 187, <$line> line 618.  BELLOW
	if($tpf_output_inital{$id_contig[$current_line]}[0] eq 'GAP'){#gap
	    my $contig_length = $end_contig[$current_line] - $start_contig[$current_line] + 1;
	    print join("\t",('GAP', #true for all gaps
		       $tpf_output_inital{'s'.$id_contig[$current_line]}[1], #data (gap type)
		       $contig_length,
#		       ($end_contig[$current_line] - $start_contig[$current_line] + 1), #calculate length, incase change
		       $tpf_output_inital{'s'.$id_contig[$current_line]}[3], #type
		       '' #true for all gaps
		       ));
	    print "\n";
	}
	else{#contig
#Use of uninitialized value in join or string at bac_integration/bac_chr_position.pl line 202, <$line> line 618.      BELLOW
	    print join("\t",($id_contig[$current_line], 
		      $tpf_output_inital{$id_contig[$current_line]}[1],#data
 		      $tpf_output_inital{$id_contig[$current_line]}[2],#scaffold
		      $tpf_output_inital{$id_contig[$current_line]}[3]#$orientation_tpf
		       ));
	    print "\n";
#remember: %tpf_output_inital{id}=@($component,$data,$scaffold,$orientation_tpf)
	}
    }

    elsif($id_contig[$current_line] =~m/^C\d\dHBa/){#bac's start with C##HBa..., and are not in the previous tpf file
	my $previous_line = $current_line - 1;
	my $next_line = $current_line + 1;
	if($tpf_output_inital{$id_contig[$previous_line]}[0] eq 'GAP'){#bac began in a gap
	    if($tpf_output_inital{$id_contig[$previous_line - 1]}[3]=$tpf_output_inital{$id_contig[$next_line]}[3]){#if orientation of last contig = orientation of next contig: set orientation
		$bac_orient = $tpf_output_inital{$id_contig[$next_line]}[3];
		print STDERR "BAC good, and oriented";
	    }
	    elsif($tpf_output_inital{$id_contig[$previous_line - 1]}[3]!=$tpf_output_inital{$id_contig[$next_line]}[3]){#if previous contig orientation diffrent from following contig orientation, orientation of BAC not known
		print STDERR "ERROR!!!!!    THUS NO BAC ORIENTATION!!!!\n$tpf_output_inital{$id_contig[$previous_line - 1]}[3]!=$tpf_output_inital{$id_contig[$next_line]}[3])\n\n\n\n";
	    }	
	}
	if($tpf_output_inital{$id_contig[$next_line]}[0] eq 'GAP'){#bac began in a gap
	    if($tpf_output_inital{$id_contig[$previous_line]}[3]=$tpf_output_inital{$id_contig[$next_line+1]}[3]){#if orientation of last contig = orientation of next contig: set orientation
		$bac_orient = $tpf_output_inital{$id_contig[$previous_line]}[3];
		print STDERR "BAC good, and oriented";
	    }
	    elsif($tpf_output_inital{$id_contig[$previous_line]}[3]!=$tpf_output_inital{$id_contig[$next_line+1]}[3]){#if previous contig orientation diffrent from following contig orientation, orientation of BAC not known
		print STDERR "ERROR!!!!!    THUS NO BAC ORIENTATION!!!!\n$tpf_output_inital{$id_contig[$previous_line]}[3]!=$tpf_output_inital{$id_contig[$next_line+1]}[3])\n\n\n\n";
	    }
	    print join("\t",($id_contig[$current_line], '?', $bac_scaffold{$id_contig[$current_line]}, $bac_orient || ''));
#!!!!!!!!!!!
	    print "\n";
	}

	else{#was not following a gap, thus overlaped with scaffold and can print normaly, the orientation will be filled automaticaly with conversion of file.
	    print join("\t",($id_contig[$current_line], '?', $bac_scaffold{$id_contig[$current_line]}, ''));
	    print "\n";
##!!!!!!
#do we put a ? like in scaffolds? or length like in gaps?
	}
    }
    $current_line++;
}
	    #BLAST BACs -> scaffolds
 	    #-> tomato scaffold ID, start & end coord.
 	    # |
 	    # \/
	    #use combined *.chr.agp to find
	    #scaffold ID.
	    #retrieve chr id + chromosome coordinate
	    #of match (chr-scaffold start + match start -1)
	    #and end
	    #go through chr*.comp.agp and find contig
	    #matching to start and end points.
	    #replace intervening gaps and contigs.
	    #keep hash of replaced contigs in a
	    #hash with BACid as key.
	    #replace tomato
	    #BACid with genbank BACids.
