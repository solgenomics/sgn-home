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
    my @bac_info = ($s_id, $bac_start, $len);
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

	if($data[6] =~ m/^\d+/){#contig
	    my $contig_start = $data[1] + $data[6] -1;
	$start2chr{$contig_start} = \@data;
	}

	else{#gap
	    $start2chr{$data[1]} = \@data;
	}

    }#end agp loop
    
    my %contig2scaf;
    my %scaf2contig;
    my %scaf2chr;
    my $current_scaf_id;
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
	    $scaf2contig{$scaffold} = $contig_id;
	    if(!$scaf2chr{$scaffold}) {
		$scaf2chr{$scaffold} = $contig2chr{$contig_id}->[1];
	    }
	}
    }

    my $scaf_start = '';
    my $scaf_end = '';


    foreach my $bac_id ( keys %bac2scaf){ #each bac line data
	
	if($bac2scaf{$bac_id}->[0] && $scaf2chr{$bac2scaf{$bac_id}->[0]}){ 
	    #if scaffold exists

	    $contig2scaf{$bac_id} = $bac2scaf{$bac_id}->[0];
	    
	    if ($contig2chr{$scaf2contig{$bac2scaf{$bac_id}->[0]}}->[1] < 
		$contig2chr{$scaf2contig{$bac2scaf{$bac_id}->[0]}}->[2]){
		#start after end
		
		#correct orientation
		$scaf_start = 
		    $contig2chr{$scaf2contig{$bac2scaf{$bac_id}->[0]}}->[1];
		
		$scaf_end = 
		    $contig2chr{$scaf2contig{$bac2scaf{$bac_id}->[0]}}->[2];
		
		$current_orientation = '';
		
	    }
	    
	    else{ #reverse orientation
		
		$scaf_end = 
		    $contig2chr{$scaf2contig{$bac2scaf{$bac_id}->[0]}}->[1];
		
		$scaf_start = 
		    $contig2chr{$scaf2contig{$bac2scaf{$bac_id}->[0]}}->[2];
		
	       $current_orientation = 'REV';
		
	    }
	    
	    my @current_bac_data = (
		$start2chr{$scaf2chr{$bac2scaf{$bac_id}->[0]}}->[0], #chr_id
		$scaf_start || '', #s_start
		$scaf_end || '', #s_end
		'', #number to be filled in later
		'O', #type for other, aka: BAC
		$bac_id,
		$bac2scaf{$bac_id}->[1], #start in scaffold
		$bac2scaf{$bac_id}->[2], #end in scaffold
		$current_orientation || '', #orientation
		'' #contained?
		);
	    
	    $contig2chr{$bac_id} = \@current_bac_data;
	    my $bac_begin = $current_bac_data[1] + $current_bac_data[6] -1;
	    $start2chr{$bac_begin} = \@current_bac_data;

	}
    }#end addition of bacs to contig2chr
    
    
    print $OUTTPF "\n##=== Begining of TPF Data ===\n";
    
    
    my @ordered_starts = sort {$a <=> $b} keys %start2chr;
        my $current_id = 0;
    my $i;
    my $next_end;
    my $previous_end;
    my $current_length = 0;
    my $current_position = 0;


    while ($ordered_starts[$current_id]){

	$current_position ++;
	if($start2chr{$ordered_starts[$current_id]}->[5] =~ m/^C/){ #is bac
	    
	    
	    my $current_bac_start = 
		$start2chr{$ordered_starts[$current_id]}->[1];
	    
	    my $current_bac_end = $start2chr{
		$ordered_starts[$current_id]}->[1] +
		    $bac2scaf{$start2chr{
			$ordered_starts[$current_id]}->[5]}->[2];
	    
	    if ($ordered_starts[$current_id + 1]){

		my $next_start = $start2chr{
		    $ordered_starts[$current_id + 1]}->[1];

		if($start2chr{
		    $ordered_starts[$current_id + 1]}->[6] !~m/^\d+$/){
		    #next is gap
		    
		    $next_end = $next_start + 
			($start2chr{$ordered_starts[$current_id + 1]}->[5]);
		    #gap length
		}
		
		else{ #next is not gap
		    
		    $next_start += $start2chr{
			$ordered_starts[$current_id + 1]}->[6] - 1;
		    
		    $next_end = $next_start + 
		    $start2chr{$ordered_starts[$current_id + 1]}->[7] -
		    $start2chr{$ordered_starts[$current_id + 1]}->[6];
		    
		}
		
		if($start2chr{
		    $ordered_starts[$current_id - 1]}->[6] =~m/^\d+$/){
		    #not gap
		    
		    $previous_end = $start2chr{
			$ordered_starts[$current_id - 1]} +#start
			$start2chr{$ordered_starts[$current_id - 1]}->[7] -
			$start2chr{$ordered_starts[$current_id - 1]}->[6];
		    #cotig length
		    
		}
		
		else{#gap
		    
		    $previous_end=$start2chr{
			$ordered_starts[$current_id - 1]}->[1] +
		    #start scaf
			$start2chr{$ordered_starts[$current_id - 1]}->[5];
		    #gapLength;
		}
		
		#look at previous data, if gap, will be accounted for in output
		if($start2chr{
		    $ordered_starts[$current_id - 1]}->[6] =~ m/^\d+$/){
		    #prevous is not gap

		    if($start2chr{$ordered_starts[$current_id-1]}->[2] > 
		       $start2chr{$ordered_starts[$current_id]}->[1] ){
			#around bac start
			
			#make current back contained turnout of previous
			$start2chr{$ordered_starts[$current_id]}->[9] = 
			    'CONTAINED TURNOUT'."\t".
			    $start2chr{$ordered_starts[$current_id - 1]}->[5];
		    }#end around bac start
		    
		    
		    if($previous_end > $current_bac_end){#bac in previous contig
			$start2chr{$ordered_starts[$current_id]}->[9] = 
			'CONTAINED'."\t".$start2chr{
			    $ordered_starts[$current_id - 1]}->[5];
		    }

		    $i = 1;
		    while($i > 0 && $ordered_starts[$current_id + $i] &&
			  $start2chr{
			      $ordered_starts[$current_id + $i]}->[1]){
			
			
			#reset start
			my $next_start = 
			    $start2chr{$ordered_starts[$current_id+$i]}->[1];
			
			$next_end = '';#reset end

			if($start2chr{$ordered_starts[$current_id + $i]}->[6] =~
			   m/^\d+$/){#not gap
			    $next_end = $next_start + 
				$start2chr{
				    $ordered_starts[$current_id + $i]}->[7] -
					$start2chr{$ordered_starts[
						       $current_id + $i]}->[6];
			    
			}
			elsif($start2chr{
			    $ordered_starts[$current_id + $i]}->[5] =~
			      m/^\d+$/){ #gap

			    $next_end = $next_start + 
				$start2chr{$ordered_starts[
					       $current_id + $i]}->[5];
			    #gap length
			}
			
			if ($current_bac_start < $next_start && 
			    $next_end < $current_bac_end){#in bac
			    
			    if ($start2chr{$ordered_starts[
					       $current_id + $i]}->[4] eq 'O'){
				               #bac
				
				$start2chr{
				    $ordered_starts[$current_id + $i]}->[9] =
					'CONTAINED'."\t".$start2chr{
					    $ordered_starts[$current_id]}->[5];
			    } #end if bac
			    
			    else{#remove
				
				if($start2chr{$ordered_starts[
						  $current_id + $i]}->[4] eq
				   'U'){ #fill unknown gap!  MAYBE....
				    my $k=$i+1;
				    while (!$contig2scaf{
					$start2chr{$ordered_starts[
						       $current_id+$k]}->[5]}){
					
					$k++;
					   }

				    print "JOINED SCAFFOLDS(".
					$contig2scaf{$start2chr{
					    $ordered_starts[$current_id]}->[5]}.
				        " and ".
					$contig2scaf{$start2chr{
					    $ordered_starts[
						$current_id + $k]}->[5]}.
					")!!!in $agp...\twhere BAC_id=".
					$start2chr{
					    $ordered_starts[$current_id]}->[5].
					"\n";

				}

				splice(@ordered_starts,$current_id + $i,1);

				$i--;#cancel out addition at end of while


			    } #end remove, not bac
			    
			}#end in bac
			
			elsif($next_end > $current_bac_end && 
			      $next_start < $current_bac_end){#end of bac

			    for(my $j=0; $j<$i; $j++){#find turnout_id

			    $start2chr{$ordered_starts[$current_id + $i]}->[9] =
				'CONTAINED TURNOUT'."\t".
				$start2chr{
				    $ordered_starts[$current_id + $j]}->[5];

			    } #end for turnout


			} #end elsif at end of bac
			
			$i++;
			
			
		    } #end of while
		    
		    
		    
		} #end if not gap

		
	    }#end while next exists
	    
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

	    print $OUTTPF join(
		"\t", 
		$start2chr{$ordered_starts[$current_id]}->[5], #bac_id
		'?',
		$bac2scaf{$start2chr{$ordered_starts[$current_id]}->[5]}->[0],
		#s_id
		); #end of print join


	    $current_length = $start2chr{$ordered_starts[$current_id]}->[2] -
		$start2chr{$ordered_starts[$current_id]}->[1];


	    print $OUTAGP join(
		"\t",
		$start2chr{$ordered_starts[$current_id]}->[0],
		$current_position,
		($current_position + $current_length) || '',
		$current_id,
		$start2chr{$ordered_starts[$current_id]}->[4],
		$start2chr{$ordered_starts[$current_id]}->[5],
		$start2chr{$ordered_starts[$current_id]}->[6],
		$start2chr{$ordered_starts[$current_id]}->[7],
		$start2chr{$ordered_starts[$current_id]}->[8]
		), "\n";
	    
	    $current_position += $current_length;

	    if($start2chr{$ordered_starts[$current_id]}->[9]){
		print $OUTTPF "\t".
		    $start2chr{$ordered_starts[$current_id]}->[9].
		    "\n";
	    } #end of if has contained info
	    
	    else{
		print $OUTTPF "\t".
		    $current_orientation. #orient
		    "\n";
	    } #end else, print orientation


	}# end if(bac)
	
	
	else{ #not bac



	    if($start2chr{$ordered_starts[$current_id]}->[6] =~ m/^\d+$/){
		#not gap

		
		$current_length = $start2chr{
		    $ordered_starts[$current_id]}->[2] -
			$start2chr{$ordered_starts[$current_id]}->[1];


		print $OUTAGP join(
		    "\t",
		    $start2chr{$ordered_starts[$current_id]}->[0],
		    $current_position,
		    ($current_position + $current_length) || '',
		    $current_id,
		    $start2chr{$ordered_starts[$current_id]}->[4],
		    $start2chr{$ordered_starts[$current_id]}->[5],
		    $start2chr{$ordered_starts[$current_id]}->[6],
		    $start2chr{$ordered_starts[$current_id]}->[7],
		    $start2chr{$ordered_starts[$current_id]}->[8]
		    ), "\n";
		
		$current_position += $current_length;


		if($start2chr{$ordered_starts[$current_id]}->[8] eq '-'){
		    $current_orientation = 'MINUS';
		}
		
		if($start2chr{$ordered_starts[$current_id]}->[8] eq '+'){
		    $current_orientation = 'PLUS';
		}


		my $current_contig_id = (
		    $start2chr{$ordered_starts[$current_id]}->[5]);

		$current_contig_id =~ s/(.*?)\.\d+/$1/; #remove version number

		print $OUTTPF join("\t",
			 $start2chr{$ordered_starts[$current_id]}->[5], #id
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
		    print $OUTTPF "\n";
		}
		
	    }
	    elsif(
		$start2chr{$ordered_starts[$current_id + 1]}->[6] =~ m/^\d+$/ &&
		$start2chr{$ordered_starts[$current_id - 1]}->[6] =~ m/^\d+$/){
		#gap

		if($start2chr{$ordered_starts[$current_id]}->[4] eq 'U'){
		    #unknown gap
		    $current_length = 100;
		}

		else{
		   $current_length = $start2chr{
		       $ordered_starts[$current_id]}->[2] -
			   $start2chr{$ordered_starts[$current_id]}->[1];
		}
		
		
		print $OUTAGP join(
		    "\t",
		    $start2chr{$ordered_starts[$current_id]}->[0],
		    $current_position,
		    ($current_position + $current_length) || '',
		    $current_id,
		    $start2chr{$ordered_starts[$current_id]}->[4],
		    $current_length || '',
		    $start2chr{$ordered_starts[$current_id]}->[6],
		    $start2chr{$ordered_starts[$current_id]}->[7],
		    ''
		    ), "\n";
		
		$current_position += $current_length;
		
		
		my $gap_length = #end:
		    ($start2chr{$ordered_starts[$current_id + 1]}->[1] + 
		        #next s_start
		     $start2chr{$ordered_starts[$current_id + 1]}->[6] - 
		        #next start on scaf
		     1) - #correct for counting
		     #start:
		     #prev s_start:
		     ($start2chr{$ordered_starts[$current_id - 1]}->[1] + 
		      #prev end:
		      $start2chr{$ordered_starts[$current_id - 1]}->[7] - 
		      #prev start on staf:
		      $start2chr{$ordered_starts[$current_id - 1]}->[6]);  
		#end gap_length calculation


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
			 $gap_length, #length of gap
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


