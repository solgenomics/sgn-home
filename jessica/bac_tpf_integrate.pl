#!usr/bin/perl

use strict;
use warnings;

#similify bacs before running, only use the "best" ones.
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

my @start_contig;
my @end_contig;
my @id_contig;
my @chr_contig;
my @turnout;
my %chr_defign;
#upload agp fie into several arrays to allow for easy access and editing
open(my $F,"<",$agpchr) || die "can't open $agpchr";

while(<$F>){
    chomp;
    if(/^$/ || /^\#/){next;} #skip blank and comment lines
    my ($obj, $start_obj, $end_obj, $id, $type, $comp_id, $length_agp, $orientation_agp) = split/\t/;
    
    my $chr_current = $chr{$obj} || '';
    if($start_obj && $comp_id && $comp_id =~m/^AEKE/){#contig
	$comp_id=~s/(.*?)\.\d+/$1/; #remove version number
	$chr_defign{$comp_id}=$chr_current;
    }
    
    if($start_obj && $end_obj && $comp_id && $chr_current){
	push (@start_contig, $start_obj);
	push (@end_contig, $end_obj);
	push (@id_contig, $comp_id);
	push (@chr_contig, $chr_current);
	push (@turnout, '');
    }
    else{
	print STDERR "One or more of the following does not exist:\nstart_obj=$start_obj && end_obj$end_obj= && comp_id=$comp_id && chr_current=$chr_current\n";
    }
}

#input of tpf file into a hash of lists
my @tpf_output_initial;
my %scaffold_convert;
open($F,"<",$tpfchr) || die "can't open $tpfchr";
while(<$F>){
    chomp;
    if(/^$/ || /^\#/){next;} #skip blank and comment lines
    my ($component, $data, $scaffold, $orientation_tpf, $blank) = split/\t/;
    
    if ($component && $component =~m/^AEKE/){  #exist, not gap (gaps do not have associated scaffold information)
	$scaffold_convert{$scaffold} = $component;
    }
    push (@tpf_output_initial, [$component, $data, $scaffold, $orientation_tpf]);
}

my %bac_scaffold;
my $bac_begin;
my $bac_end;
my $bac_chr;
my $no_version;
#run through bac(blast) file to find where they fit in tpf/agp and replace/remove information accordingly
open(my $line, "<", $bac_file) || die "can't open $bac_file.";
while(<$line>){
    print STDERR ".";
    chomp;
    if(/^$/ || /^\#/){next;} #skip blank and comment lines
    my ($q_id, $s_id, $per_id, $len, $mismatch, $gaps, $q_start, $q_end, $s_start, $s_end, $e_val, $bit) = split /\t/;
    if (exists($scaffold_convert{$s_id})
	){
	$bac_chr =  $chr_defign{$scaffold_convert{$s_id}};
	$bac_begin = $s_start + $q_start - 1;
	$bac_end = $s_end + $q_start - 1;
	
	my $i = 0;
	my $contig_count = @chr_contig;
	
	while($i < $contig_count){
	    if($id_contig[$i] eq $q_id){#current contig/gap is  current bac    
		#move to next contig/gap
	    }
	    elsif($bac_chr eq $chr_contig[$i] && $bac_begin <= $end_contig[$i]){#close on same chromosome
		#check bac has crossover
		if($bac_end <= $start_contig[$i]){#past bac
		    $i= $contig_count; #move on(next bac), ends while i
		}

		elsif($bac_begin >= $start_contig[$i]){#around bac_start
		    if($id_contig[$i] =~m/^[0-9]/){#gap (contig start with 'A' and BACs start with 'C')
			#alter to remove gap where bac overlaps
			$end_contig[$i] = $bac_begin-1;#gap set to end at  bac start
		    }	
		    #insert bac after obj in array
		    $bac_scaffold{$q_id} = $s_id;
		    splice(@start_contig,$i+1,0,$bac_begin || '');
		    splice(@end_contig,$i+1,0,$bac_end || '');
		    splice(@chr_contig,$i+1,0,$bac_chr || '');
		    splice(@id_contig,$i+1,0,$q_id || '');
		    splice(@tpf_output_initial,$i+1,0,['error','error','error','error']);
		    splice(@turnout,$i+1,0,'');
		    if($id_contig[$i] !~m/^[0-9]/){#not gap
			$no_version = $id_contig[$i];
			$no_version =~s/(.*?)\.\d+/$1/; #remove version number
			$turnout[$i+1] = 'CONTAINED_TURNOUT'."\t".$no_version;
			if($id_contig[$i+2] =~m/^C\d\dHBa/){#current bac followed by bac
			    if($start_contig[$i+2] > $bac_begin){
				$no_version = $q_id;
				$no_version =~s/(.*?)\.\d+/$1/; #remove version number
				$turnout[$i+2] = 'CONTAINED_TURNOUT'."\t".$no_version;
			    }
			    elsif($start_contig[$i+2] < $bac_begin){
				#SWITCH PLACES, AND SET CURRENT TO IN OLD
				#replace both locations with old bac info
				splice(@start_contig,$i+1,1,$start_contig[$i+2]);
				splice(@end_contig,$i+1,1,$end_contig[$i+2]);
				splice(@chr_contig,$i+1,1,$chr_contig[$i+2]);
				splice(@id_contig,$i+1,1,$id_contig[$i+2]);

				#change second bac to current information
				splice(@start_contig,$i+2,1,$bac_begin || '');
				splice(@end_contig,$i+2,1,$bac_end || '');
				splice(@chr_contig,$i+2,1,$bac_chr || '');
				splice(@id_contig,$i+2,1,$q_id || '');
				$no_version = $id_contig[$i];
				$no_version =~s/(.*?)\.\d+/$1/; #remove version number
				$turnout[$i+1] = 'CONTAINED_TURNOUT'."\t".$no_version;
				$no_version = $id_contig[$i+1];
				$no_version =~s/(.*?)\.\d+/$1/; #remove version number
				$turnout[$i+2] = 'CONTAINED_TURNOUT'."\t".$no_version;
			    }
			}
		    }
		    $i++ #to correctly calculate next contig/gap
		}

		elsif($bac_begin <= $start_contig[$i] && 
		      $bac_end >= $end_contig[$i]){#within bac
		    #remove contig/gap from each array 
		    splice(@start_contig,$i,1);
		    splice(@end_contig,$i,1);
		    splice(@chr_contig,$i,1);
		    splice(@id_contig,$i,1);
		    splice(@tpf_output_initial,$i,1);
		    splice(@turnout,$i,1);
		    $i--; #re-set i due to removed space
		}

		elsif($bac_begin <= $start_contig[$i] && 
		      $bac_end <= $end_contig[$i]){#around bac end
		    if(!($id_contig[$i] =~m/^[0-9]/)){#gap
			#alter to remove gap where bac overlaps
			$start_contig[$i] = $bac_end; #set gap start at bac end
		    }
		    else{
			$no_version = $id_contig[$i-1];
			$no_version =~s/(.*?)\.\d+/$1/; #remove version number
			$turnout[$i] = 'CONTAINED_TURNOUT'."\t".$no_version;
		    }
		    $i= $contig_count; #move on(next bac), ends while i
		}
	    }
	    $i++; #moves to next contig/gap
	}
    }
}

my $bac_orient;
my $id_size = @id_contig;
my $current_line = 0;
while($id_contig[$current_line] && $current_line <= $id_size){#contig or gap
    if($tpf_output_initial[$current_line][0] eq 'GAP'){#gap
	my $contig_length = $end_contig[$current_line] - $start_contig[$current_line] + 1;
	print join("\t",('GAP', #true for all gaps
			 $tpf_output_initial[$current_line][1], #data (gap type)
			 $contig_length,
			 $tpf_output_initial[$current_line][3], #type
			 '', #true for all gaps
			 ));
	print "\n";
    }
    elsif($id_contig[$current_line] =~m/^AEKE/){#contig
	if($turnout[$current_line] eq ''){
	    print join("\t",($id_contig[$current_line], 
			     $tpf_output_initial[$current_line][1],#data
			     $tpf_output_initial[$current_line][2],#scaffold
			     $tpf_output_initial[$current_line][3],#$orientation
		       ));
	}
	else{
	    print join("\t",($id_contig[$current_line], 
			     $tpf_output_initial[$current_line][1],#data
			     $tpf_output_initial[$current_line][2],#scaffold
			     $turnout[$current_line]#turnout in place of orientation
		       ));
	}
	print "\n";
#remember: %tpf_output_initial{id}=@($component,$data,$scaffold,$orientation_tpf)
    }

    elsif($id_contig[$current_line] =~m/^C\d\dHBa/){#bac's start with C##HBa..., and are not in the previous tpf file
	my $previous_line = $current_line - 1;
	my $next_line = $current_line + 1;
	if($id_contig[$previous_line] !~m/^AEKE/){#bac began in a gap or another bac
	    if($id_contig[$previous_line - 1]!~m/^AEKE/){#Bac is proceded by bac that starts at a gap
		$next_line++;
	    }
	    if($tpf_output_initial[$previous_line - 1][3]=$tpf_output_initial[$next_line][3]){#if orientation of last contig = orientation of next contig: set orientation
		$bac_orient = $tpf_output_initial[$next_line][3];
	    }
	    elsif($tpf_output_initial[$previous_line - 1][3] ne $tpf_output_initial[$next_line][3]){#if previous contig orientation diffrent from following contig orientation, orientation of BAC not known
		print STDERR "ERROR!!!!!     NO BAC ORIENTATION!!!!\nbefore=$tpf_output_initial[$previous_line - 1][3]!=$tpf_output_initial[$next_line][3]=after)\n\n\n\n";
	    }	
	}
	if($id_contig[$next_line] !~m/^AEKE/){#bac end at non-contig (either gap or another BAC
	    if($id_contig[$next_line + 1]!~m/^AEKE/){#Bac is following bac, followed by gap
		$next_line++;
	    }
	    if($tpf_output_initial[$previous_line][3]=$tpf_output_initial[$next_line+1][3]){#if orientation of last contig = orientation of next contig: set orientation
		$bac_orient = $tpf_output_initial[$previous_line][3];
	    }
	    elsif($tpf_output_initial[$previous_line][3] ne $tpf_output_initial[$next_line+1][3]){#if previous contig orientation diffrent from following contig orientation, orientation of BAC not known
		print STDERR "ERROR!!!!!    NO BAC ORIENTATION!!!!\nbefore=$tpf_output_initial[$previous_line][3]!=$tpf_output_initial[$next_line+1][3]=afer)\n\n\n\n";
	    }
	    if($turnout[$current_line] eq ''){
		print join("\t",($id_contig[$current_line], '?', $bac_scaffold{$id_contig[$current_line]}, $bac_orient || ''));
	    }
	    else{
		print join("\t",($id_contig[$current_line], '?', $bac_scaffold{$id_contig[$current_line]}, $turnout[$current_line]));
	    }
	    print "\n";
	}

	else{#was not following a gap, thus overlaped with scaffold and can print normaly, the orientation will be filled automaticaly with conversion of file.
	    if($turnout[$current_line] eq ''){
		print join("\t",($id_contig[$current_line], '?', $bac_scaffold{$id_contig[$current_line]}, ''));
	    }
	    else{
		print join("\t",($id_contig[$current_line], '?', $bac_scaffold{$id_contig[$current_line]}, $turnout[$current_line]));
	    }		
	    print "\n";
	}
    }
    $current_line++;
}
