#! usr/bin/perl

use strict;
use warnings;
use autodie;

my $bac_file = shift;
my $dir = shift;
#be sure the agp and tpf are named the SAME AND in the location given.



my %bac2scaf;
open (my $F, "<", $bac_file);
print STDERR "Loading $bac_file information...\n";
while(<$F>){#load BAC file
    
    chomp;
    if(/^$/ || /^\#/){next;} #skip blank and comment lines
    my ($bac_id, $s_id, $per_id, $len, $mismatch, $gaps, 
	$bac_start, $bac_end, $s_start, $s_end, $e_val, $bit) = split /\t/;
    my @bac_info = ($s_id, $bac_start, $len, $s_start, $s_end);
    $bac2scaf{$bac_id} = \@bac_info;
    
}#end BAC load

my $current_orientation;

foreach my $agp (glob $dir."*.comp.agp") {
    
    my $tpf = $agp.".tpf";
    
    print STDERR "opening files $agp and $tpf...\n";

    $agp =~ s/(.*?)\.agp$/$1/;
    open (my $OUTTPF, ">", $agp.".bac.agp.tpf");
    open (my $OUTAGP, ">", $agp.".bac.agp");
    open ($F, "<", $agp.".agp");

    my %contig2chr;
    my %start2chr;
    my $contig_input;
    while(<$F>){#look at each agp line
	chomp;

	if(/^$/ || /^\#/){ #print header of agp
	    print $OUTAGP $_."\n";
	    next;
	}

	my @data = split/\t/;	
	$contig_input = $data[5];
	$contig_input =~ s/(.*?)\.\d+/$1/; #remove version number
	
	$data[9]= '';
	$contig2chr{$contig_input} = \@data;
	$start2chr{$data[1]} = \@data;


    }#end agp loop
    
    my %contig2scaf;
    my %scaf2start;
    my %scaf2chr;
    my $current_scaf_id;
    my $current_contig_id;
    open ($F, "<", $tpf);
    
    while(<$F>){ #each tpf line
	
	if(/^\#\#[A-Z]/){
	    print $OUTTPF $_;
	    next;
	}
	
	if(/^$/ || /^\#/){next;} #skip blank and comment lines
	
	chomp; 
	my ($contig_id, $data_type, $scaffold, $orientation_tpf, 
	    $blank) = split/\t/;
	
	if($contig_id !~ m/^GAP$/){
	    $contig2scaf{$contig_id} = $scaffold;

	    if(!$scaf2start{$scaffold}){
	    $scaf2start{$scaffold} = $contig2chr{$contig_id}->[1];
	    }

	    if(!$scaf2chr{$scaffold}) {
		$scaf2chr{$scaffold} = $contig2chr{$contig_id}->[0];
	    }
	}
    }

    my $scaf_start = '';
    my $scaf_end = '';


    foreach my $bac_id ( keys %bac2scaf){ #each bac line data
	
	if($bac2scaf{$bac_id}->[0] && $scaf2chr{$bac2scaf{$bac_id}->[0]}){ 
	    #if scaffold exists

	    $contig2scaf{$bac_id} = $bac2scaf{$bac_id}->[0];
	    
	    if ($bac2scaf{$bac_id}->[4] < $bac2scaf{$bac_id}->[3]){
		#start after end
		
		#reverse orientation
		$scaf_start = $bac2scaf{$bac_id}->[4] + #end in scaffold
		    $scaf2start{$bac2scaf{$bac_id}->[0]} - 1;
		    #start of scaffold

		$scaf_end = $bac2scaf{$bac_id}->[3] + #start in scaffold
		    $scaf2start{$bac2scaf{$bac_id}->[0]} - 1;
		    #start of scaffold
		$current_orientation = 'REV';
		
	    }
	    
	    else{ #correct orientation
		
		$scaf_end = $bac2scaf{$bac_id}->[4] + #end in scaffold
		    $scaf2start{$bac2scaf{$bac_id}->[0]} - 1;
		    #start of scaffold

		$scaf_start = 	$bac2scaf{$bac_id}->[3] + #start in scaffold
		    $scaf2start{$bac2scaf{$bac_id}->[0]} - 1;
		    #start of scaffold	

		$current_orientation = '';
		
	    }

	    my @current_bac_data = (
		$scaf2chr{$bac2scaf{$bac_id}->[0]}, #chr_id
		$scaf_start || '', #s_start
		$scaf_end || '', #s_end
		'', #number to be filled in later
		'O', #type for other, aka: BAC
		$bac_id,
		$bac2scaf{$bac_id}->[1], #start in scaffold
		$bac2scaf{$bac_id}->[2] + $bac2scaf{$bac_id}->[1] - 1,
		 #end in scaffold (length + o_start (normaly 1) -1
		$current_orientation || '', #orientation
		'' #contained?
		);
	    
	    $contig2chr{$bac_id} = \@current_bac_data;
	    $start2chr{$current_bac_data[1]} = \@current_bac_data;

	}
    }#end addition of bacs to contig2chr
    
    
    print $OUTTPF "\n##=== Begining of TPF Data ===\n";
    
    
    for my $key_test ( keys %start2chr ) {#each hashkey of hash start2chr
	if($key_test !~ m/^\d+$/){#if key not all numbers
	    delete $start2chr{$key_test}; #remove key and values from hash
	}
    }

    my @ordered_starts = sort {$a <=> $b} keys %start2chr;
        my $current_id = 0;
    my $i;
    my $next_end;
    my $previous_end;
    my $current_length = 0;
    my $current_position_move = 0;
    my $count = 0;


    while ($ordered_starts[$current_id]){

	if($start2chr{$ordered_starts[$current_id]}->[5] =~ m/^C/){ #is bac
	    
 
	    my $current_bac_end =
		$start2chr{$ordered_starts[$current_id]}->[2];
	    
	    if ($ordered_starts[$current_id + 1]){

		my $next_start = $ordered_starts[$current_id + 1];

		$next_end =
		    $start2chr{$ordered_starts[$current_id + 1]}->[2];

		$previous_end = 0;
		if ($ordered_starts[$current_id - 1]){

		    $previous_end =
			$start2chr{$ordered_starts[$current_id - 1]}->[2];
		    #gapLength;
		}
		
		#look at previous data, if gap, will be accounted for in output
		if($start2chr{
		    $ordered_starts[$current_id - 1]}->[6] =~ m/^\d+$/){
		    #prevous is not gap
		    
		    if($previous_end > $ordered_starts[$current_id]){
			#around bac start

			$count = 1;
			while($start2chr{
			    $ordered_starts[$current_id - $count]}->[9] =~
			      m/^CONTAINED/i &&
			      $start2chr{$ordered_starts[
					     $current_id - $count]}->[9] !~
			      m/^CONTAINED TURNOUT/i){
			    #find furthest/largest contained
			    
			    $count++;
			}
		
			if($start2chr{
			    $ordered_starts[$current_id - 
					    $count]}->[2] > $current_bac_end){
			    #overlaping contained
			    
			    $start2chr{$ordered_starts[$current_id]}->[9] = 
				'CONTAINED'."\t".
				$start2chr{
				    $ordered_starts[$current_id - $count]}->[5]
			}
			
			else{			
			    #make current back contained turnout of previous
			    $start2chr{$ordered_starts[$current_id]}->[9] = 
				'CONTAINED TURNOUT'."\t".
				$start2chr{
				    $ordered_starts[$current_id - 1]}->[5];
			}
		    }#end around bac start
		    
		
		    if($previous_end > $current_bac_end){#bac in previous contig


			$start2chr{$ordered_starts[$current_id]}->[9] = 
			    'CONTAINED'."\t".$start2chr{
				$ordered_starts[$current_id - 1]}->[5];
		    }
		} #end prev not gap

		$i = 1;
		while($i > 0 && $ordered_starts[$current_id + $i] &&
		      $start2chr{
			  $ordered_starts[$current_id + $i]}->[1]){
			
			
		    #reset start
		    my $next_start =$ordered_starts[$current_id + $i];
			
		    $next_end = $start2chr{
			$ordered_starts[$current_id + $i]}->[2];

			
		    if ($ordered_starts[$current_id] < $next_start && 
			$next_end < $current_bac_end){#in bac
			    
			if ($start2chr{
			    $ordered_starts[$current_id + $i]}->[4] eq 'O'){
			    #bac
				
			    $start2chr{
				$ordered_starts[$current_id + $i]}->[9] =
				    'CONTAINED'."\t".$start2chr{
					$ordered_starts[$current_id]}->[5];
			} #end if bac
			
			else{#remove
				
			    if($start2chr{
				$ordered_starts[$current_id + $i]}->[4] eq
			       'U'){ #fill unknown gap!  MAYBE....
				#Could be error here and following 
				#due to estimation  not accurate,
				#tests needed to varify.

				print STDERR "FILLED UNKNOWN GAP! ";

				if($start2chr{
				    $ordered_starts[$current_id+$i]}->[7] =~ 
				      m/^yes$/i){#unknown clone gap (linkage)

				    print STDERR "IN: ".
					$contig2scaf{$start2chr{
					    $ordered_starts[$current_id]}->[5]};
				}

				elsif($start2chr{
				    $ordered_starts[$current_id+$i]}->[7] =~ 
				      m/^no$/i){
				    #unknown contig gap (no linkage)
				    
				    $count=$i+1;
				    while (!$contig2scaf{
					$start2chr{
					    $ordered_starts[
						$current_id+$count]}->[5]}){
					
					$count++;
				    }
				    
				    print STDERR "JOINING: ".
					$contig2scaf{$start2chr{
					    $ordered_starts[$current_id]}->[5]}.
					" and ".
					$contig2scaf{$start2chr{
					    $ordered_starts[
						$current_id + $count]}->[5]};

				}

				print STDERR " from $agp...Where BAC_id=".
				    $start2chr{
					$ordered_starts[$current_id]}->[5].
				    "\n";
			    }

			    splice(@ordered_starts,$current_id + $i,1);

			    $i--;#cancel out addition at end of while


			} #end remove, not bac

			if($start2chr{
			    $ordered_starts[$current_id - 1]}->[9] =~
			   m/^CONTAINED/i &&
			   $start2chr{$ordered_starts[
					  $current_id]}->[9] !~
			   m/^CONTAINED TURNOUT/i){#if prev contained

			    $count = 1;
			    while($start2chr{
				$ordered_starts[$current_id - $count]}->[9] =~
				  m/^CONTAINED/i &&
				  $start2chr{$ordered_starts[
						 $current_id - $count]}->[9] !~
				  m/^CONTAINED TURNOUT/i){
			    #find furthest/largest contained
				
				$count++;
			    }
			    
			    if($start2chr{
				$ordered_starts[
				    $current_id - 
				    $count]}->[2] > $current_bac_end){
				#overlaping contained
				
				$start2chr{$ordered_starts[$current_id]}->[9] = 
				'CONTAINED'."\t".
				$start2chr{
				    $ordered_starts[$current_id - $count]}->[5]
			    }
			} #end previous contained
		    }#end in bac
		    
		    elsif($next_end > $current_bac_end && 
			  $next_start < $current_bac_end){#end of bac

			if($start2chr{
			    $ordered_starts[$current_id]}->[9] eq 
			   'CONTAINED'){

			    $count = 1;
			    while($start2chr{
				$ordered_starts[$current_id - $count]}->[9] =~
				  m/^CONTAINED/i &&
				  $start2chr{
				      $ordered_starts[
					  $current_id - $count]}->[9] !~
				  m/^CONTAINED TURNOUT/i){
				#find furthest/largest contained
				
				$count++;
			    }
			    
			    if($start2chr{$ordered_starts[
					      $current_id - 
					      $count]}->[2] > $ordered_starts[
				   $current_id]){ #overlaping contained
				
				$start2chr{$ordered_starts[$current_id]}->[9] = 
				    'CONTAINED'."\t".
				    $start2chr{$ordered_starts[
						   $current_id - $count]}->[5]
			    }

			    else{
				$start2chr{
				    $ordered_starts[$current_id + $i]}->[9] =
					'CONTAINED TURNOUT'."\t".
					$start2chr{
					    $ordered_starts[$current_id]}->[5];
				#will result in turnout of most recent BAC.
			    }
			}#end if turnout of contained
			
			elsif($start2chr{$ordered_starts[
				$current_id + $i]}->[6] !~ m/^\d+$/ && 
			      $start2chr{$ordered_starts[
				$current_id]}->[2] > $ordered_starts[
				  $current_id+2]){
			    #previous gap that containes overlapping data,

			    #thus gap needs to be removed:
			    splice(@ordered_starts,$current_id + 1,1);
			    $i--; #correct addition at end of while
			}

			else{
			    $start2chr{
				$ordered_starts[$current_id + $i]}->[9] =
				    'CONTAINED TURNOUT'."\t".
				    $start2chr{
					$ordered_starts[$current_id]}->[5];
			    #will result in turnout of most recent BAC
			}
			
		    } #end elsif at end of bac
		    
		    $i++;
			
			
		} #end of while
		    
		
	    }#end while next exists

	    if(!$start2chr{$ordered_starts[$current_id]}->[9] && 
	       $start2chr{$ordered_starts[
			      $current_id-1]}->[9] =~ m/^CONTAINED/i && 
	       $start2chr{$ordered_starts[
			      $current_id - 1]}->[9] !~ m/^CONTAINED TURNOUT/){
		#no contained info, prev contained

			$count = 1;
			while($start2chr{
			    $ordered_starts[$current_id - $count]}->[9] =~
			      m/^CONTAINED/i &&
			      $start2chr{$ordered_starts[
					     $current_id - $count]}->[9] !~
			      m/^CONTAINED TURNOUT/i){
			    #find furthest/largest contained
			    
			    $count++;
			}
		
			if($start2chr{
			    $ordered_starts[$current_id - 
					    $count]}->[2] > $current_bac_end){
			    #overlaping contained
			    
			    $start2chr{$ordered_starts[$current_id]}->[9] = 
				'CONTAINED'."\t".
				$start2chr{
				    $ordered_starts[$current_id - $count]}->[5]
			}
	    }#end if previous contained
	    
	    my $after = $current_id + 1;
	    my $before = $current_id - 1;
	    

	    while($ordered_starts[$after] &&
		  $start2chr{$ordered_starts[$after]}->[4] ne 'W'){#gap

		$after++;
	    }
	    
	    while($ordered_starts[$before] &&
	       $start2chr{$ordered_starts[$before]}->[4] ne 'W'){#gap
		$before--;
	    }
	    
	    
	    #print bac
	    if ($ordered_starts[$after] && $ordered_starts[$before] &&
		$start2chr{$ordered_starts[$after]}->[8] &&
		$start2chr{$ordered_starts[$before]}->[8] &&
		$start2chr{$ordered_starts[$after]}->[8] eq 
		$start2chr{$ordered_starts[$before]}->[8]){
		#before and after orientation same
		if ($start2chr{$ordered_starts[$current_id]}->[8] eq 'REV'){

		    if($start2chr{$ordered_starts[$before]}->[8] eq '+'){

			$start2chr{$ordered_starts[$current_id]}->[8] = '-';
		    }
		    
		    if($start2chr{$ordered_starts[$before]}->[8] eq '-'){

			$start2chr{$ordered_starts[$current_id]}->[8] = '+';
		    }

		} #end reverse

		else{
		    $start2chr{$ordered_starts[$current_id]}->[8] = 
			$start2chr{$ordered_starts[$after]}->[8];
		}
		
	    } #end same orientation around
	    
	    elsif(!$ordered_starts[$after] ||
		  $start2chr{$ordered_starts[$current_id + 1]}->[4] eq 'U'){
		
		if ($start2chr{$ordered_starts[$current_id]}->[8] eq 'REV'){

		    if($start2chr{$ordered_starts[$before]}->[8] eq '+'){

			$start2chr{$ordered_starts[$current_id]}->[8] = '-';
		    }
		    
		    if($start2chr{$ordered_starts[$before]}->[8] eq '-'){

			$start2chr{$ordered_starts[$current_id]}->[8] = '+';
		    }
		}

		else{
		    $start2chr{$ordered_starts[$current_id]}->[8] = 
			$start2chr{$ordered_starts[$before]}->[8];
		}
		
	    }
	    
	    elsif(!$ordered_starts[$before] ||
		  $start2chr{$ordered_starts[$current_id - 1]}->[4] eq 'U'){

		if ($start2chr{$ordered_starts[$current_id]}->[8] eq 'REV'){

		    if($start2chr{$ordered_starts[$after]}->[8] eq '+'){

			$start2chr{$ordered_starts[$current_id]}->[8] = '-';
		    }
		    
		    if($start2chr{$ordered_starts[$after]}->[8] eq '-'){

			$start2chr{$ordered_starts[$current_id]}->[8] = '+';
		    }
		}
		else{
		    $start2chr{$ordered_starts[$current_id]}->[8] = 
			$start2chr{$ordered_starts[$after]}->[8];
		}
	    }
	    
	    else{
		$start2chr{$ordered_starts[$current_id]}->[8] = 'ERROR';
	    }



	    if($start2chr{$ordered_starts[$current_id]}->[8] eq '-'){
		$current_orientation = 'MINUS';
	    }
	    
	    elsif($start2chr{$ordered_starts[$current_id]}->[8] eq '+'){
		$current_orientation = 'PLUS';
	    }

	    else{ #ERROR
		$current_orientation = 'ERROR';
	    }

	    $current_contig_id = $start2chr{$ordered_starts[$current_id]}->[5];
	    $current_contig_id =~ s/(.*?)\.\d+/$1/; #remove version number

	    print $OUTTPF join(
		"\t", 
		$current_contig_id, #bac_id
		'?',
		$contig2scaf{$start2chr{$ordered_starts[$current_id]}->[5]},
		#s_id
		); #end of print join


	    $current_length = $start2chr{$ordered_starts[$current_id]}->[2] -
		$ordered_starts[$current_id];


	    if($start2chr{$ordered_starts[$current_id]}->[9]){
		print $OUTTPF "\t".
		    $start2chr{$ordered_starts[$current_id]}->[9].
		    "\n";
	    } #end of if has contained info
	    
	    else{
		print $OUTTPF "\t".
		    $current_orientation. #orient
		    "\n";
	    }#end else, print orientation


	    
	    print $OUTAGP join(
		"\t",
		$start2chr{$ordered_starts[$current_id]}->[0], #chr_id
		($ordered_starts[$current_id] +
		     $current_position_move), #start
		($start2chr{$ordered_starts[$current_id]}->[2] +
		     $current_position_move), #end
		$current_id + 1, #count
		$start2chr{$ordered_starts[$current_id]}->[4],#type
		$start2chr{$ordered_starts[$current_id]}->[5], #id
		$start2chr{$ordered_starts[$current_id]}->[6],#start_obj
		$start2chr{$ordered_starts[$current_id]}->[7]#end_obj
		);

	    if ($start2chr{$ordered_starts[$current_id]}->[9]){#contained

		print $OUTAGP "\t".
		    $start2chr{$ordered_starts[$current_id]}->[9].
		    "\n";
	    }

	    else{#orientation
		print $OUTAGP "\t".
		    $start2chr{$ordered_starts[$current_id]}->[8].
		    "\n";
	    }



	}# end if(bac)
	
	
	else{ #not bac



	    if($start2chr{$ordered_starts[$current_id]}->[6] =~ m/^\d+$/){
		#not gap

		
		$current_length = $start2chr{
		    $ordered_starts[$current_id]}->[2] -
			$ordered_starts[$current_id];


		print $OUTAGP join(
		    "\t",
		    $start2chr{$ordered_starts[$current_id]}->[0],
		    ($ordered_starts[$current_id] +
		     $current_position_move),
		    ($start2chr{$ordered_starts[$current_id]}->[2] +
		     $current_position_move) || '',
		    $current_id + 1,
		    $start2chr{$ordered_starts[$current_id]}->[4],
		    $start2chr{$ordered_starts[$current_id]}->[5],
		    $start2chr{$ordered_starts[$current_id]}->[6],
		    $start2chr{$ordered_starts[$current_id]}->[7]
		    );
		
		if ($start2chr{$ordered_starts[$current_id]}->[9]){#contained
		    
		    print $OUTAGP "\t".
			$start2chr{$ordered_starts[$current_id]}->[9].
			"\n";
		}
		
		else{#orientation
		    print $OUTAGP "\t".
			$start2chr{$ordered_starts[$current_id]}->[8].
			"\n";
		}
		
		
		if($start2chr{$ordered_starts[$current_id]}->[8] eq '-'){
		    $current_orientation = 'MINUS';
		}
		
		if($start2chr{$ordered_starts[$current_id]}->[8] eq '+'){
		    $current_orientation = 'PLUS';
		}


		$current_contig_id = $start2chr{
		    $ordered_starts[$current_id]}->[5];

		$current_contig_id =~ s/(.*?)\.\d+/$1/; #remove version number

		print $OUTTPF join("\t",
			 $current_contig_id, #id
			 '?',
			 $contig2scaf{$current_contig_id},#s_id
		    );
		
		if($start2chr{$ordered_starts[$current_id]}->[9]){
		    
		    print $OUTTPF "\t".
			$start2chr{$ordered_starts[$current_id]}->[9].
			"\n";
		}
		
		elsif($start2chr{$ordered_starts[$current_id]}->[8]){
		    
		    print $OUTTPF "\t".
			$current_orientation.
			"\n";
		}

		else{
		    print $OUTTPF "\t".''."\n";
		}
		
	    }#end if not gap

	    elsif(#0 < $current_id < @ordered_starts && #before end
		  $start2chr{
		      $ordered_starts[$current_id + 1]}->[6] =~ m/^\d+$/ &&
		  $start2chr{
		      $ordered_starts[$current_id - 1]}->[6] =~ m/^\d+$/){
		#gap, surounded by not gap


		
		print $OUTAGP join(
		    "\t",
		    $start2chr{$ordered_starts[$current_id]}->[0],
		    ($ordered_starts[$current_id] +
		     $current_position_move),
		    ($start2chr{$ordered_starts[$current_id]}->[2] +
		     $current_position_move) || '',
		    $current_id + 1,
		    $start2chr{$ordered_starts[$current_id]}->[4]
		    );
		
		$current_length = #next start:
		    $ordered_starts[$current_id + 1] -
		    #prev end:
		    $start2chr{$ordered_starts[$current_id - 1]}->[2];  
	

		if(0 > $current_length){
		    print"ERROR at start=$current_id!!Length=$current_length\n";
		}
		
		if($start2chr{$ordered_starts[$current_id]}->[4] eq 'U'){
		    #unknown gap
		    print $OUTAGP "\t100\t";
		    $current_position_move += (100 - $current_length)
		}
		else{#known gap
		    print $OUTAGP "\t".$current_length."\t";
		}
		
		print $OUTAGP join(
		    "\t",
		    $start2chr{$ordered_starts[$current_id]}->[6],
		    $start2chr{$ordered_starts[$current_id]}->[7],
		    ''
		    ), "\n";
		
		
		if($start2chr{$ordered_starts[$current_id]}->[4] eq 'U'){
		    
		    my $gap_type = "TYPE-2";
		    if ($start2chr{
			$ordered_starts[$current_id]}->[7] =~ m/no/i){
			
			#unknown contig gap
			$gap_type = "TYPE-3";
		    }
		    
		    print $OUTTPF
			join("\t", 
			     'GAP',
			     $gap_type || '',
			     '', #length of gap, not known
			     '', #type of gap, blank for not known
			     ''
			),"\n"; #end print statment

		}#end if unknown gap
		

		else{ #known gap
		    
		    print $OUTTPF
		    join("\t", 
			 'GAP',
			 "TYPE-2",
			 $current_length, #length of gap
			 "PAIRED ENDS", #type of gap (PAIRED ENDS?)
			 ''
		    ),"\n"; #end print statment

		}
		
	    }  #end if gap print 
	    
	    else{
		print $OUTTPF "***************ERROR*********************\n";
		print $OUTAGP "***************ERROR*********************\n";

	    }
	} #end else of if(bac)


	$current_id++;

    } #end while.    
    
    print $OUTTPF "##=== End of TPF Data ===\n";


}  #end file


