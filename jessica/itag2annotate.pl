#!usr/bin/perl

use warnings;
use strict;

my $itag_file = shift;
#my $nucliotide_file = shift;
my @details;
my %hash_notes;
my $seq_current = 'start';
my %current_data;
my $current_parent;
my %testing;
#my %print_details = ();

open (my $F, "<", $itag_file);
open (my $OUT, ">", $itag_file.".out");
foreach(<$F>){
    if(/^$/ || m/^\#\#\#/ || (/^\#/ && !/^\#\#/)){next;}     
   #skip empty and comment lines, but not the details (/^\#\#/)

    if(/^\#\#/){
	push(@details, "$'");
	next;
    }
#    print "\nDETAILS:\n".join("\t",@details)."\n\n";
    chomp;
    my ($seqname, $sorce, $feature, $start, $end, $score, 
	$strand, $frame, $attrubutes) = split/\t/; #split input data

	my @notes = split(/;/,$attrubutes);
	foreach my $current(@notes){

	    my ($key,$value) = split(/=/,$current);
	    $hash_notes{$key}=$value;
	    
	}

	my ($ID,$chr,$gene_num,$element_num) = split('.', $hash_notes{'ID'});


	if(!$current_parent || $current_parent ne $hash_notes{'Parent'}){
	    #new parent

	    if($seqname ne $seq_current){#new seqence/feature
	
		print $OUT ">Feature ".$seqname."\n"; #start new feature
		$seq_current = $seqname;
		print STDERR ".";
	    }

#	    print "NEW PARRENT\n";
	    if($current_parent){
		my ($parent_feature,$parent_id) = split (/:/, $current_parent);
		my $count = 0;
#		print"$current_parent\t$parent_feature\t$parent_id\n";
#		print "Keys:\t".%current_data."\n\n";
#		print join ("\t", keys %current_data);
#		print "\n\n$current_data{$parent_feature}\n";
		if($current_data{$parent_feature}->[$count]->[1] && 
		   $current_data{$parent_feature}->[$count]->[1] =~ m/^\d+$/){
		    #if three is data saved in the current parent feature,
		    #and a end position exists all numbers.

		    print"current: $current_data{$parent_feature}->[0]->[0] , $current_data{$parent_feature}->[0]->[1]\n\n";

#		    print $OUT join ("\t",
#				     @{$current_data{$parent_feature}->[0]}).
#					 "\n";#print first element of type
#		    delete($current_data{$parent_feature}->[0]);
		    my $count = 1;
		    print STDERR "WORKING Error in printing Parent feature???\ncurrent_data{parent_feature}->[count]->[0]=current_data{parent_feature}->[$count]->[0]=current_data{$parent_feature}->[$count]->[0]=$current_data{$parent_feature}->[$count]->[0] ::   $current_data{$parent_feature} ->[$count]->[0]:$current_data{$parent_feature}->[$count] ->[0]:$current_data{$parent_feature}->[$count]->[0]:\n\n";
#		    print STDERR "ERROR in printing Parent feature???\ncurrent_data{parent_feature}->[count]->[1]=current_data{parent_feature}->[$count]->[1]=current_data{$parent_feature}->[$count]->[1]=$current_data{$parent_feature}->[$count]->[1]\n\n";
#		    print STDERR "TEST:::  current_data{$parent_feature}->[$count] = ".$current_data{$parent_feature}->[$count];
		    while($current_data{$parent_feature}->[$count]->[1]){
			print "rest of parent feature...\n\n";
			print $OUT 
			    $current_data{$parent_feature}->[$count]->[0].
			    "\t".
			    $current_data{$parent_feature}->[$count]->[1].
			    "\n";

#			delete($current_data{$parent_feature}->[$count]);
#			@{$current_data{$parent_feature}->[$count]} = ();
			$count++;

		    }#end while in parent_feature

		    foreach my $current_feature(%current_data){
#			print "rest of data...\n";
			if($current_feature ne $parent_feature && 
			   $current_data{$current_feature}){

			    my $count=1;
			    my $total = @{$current_data{$current_feature}};
#			    print STDERR "total= $total\n";
			    print $OUT join (
				"\t",
				@{$current_data{$current_feature}->[0]}).
				"\t".$ID.
				"\n";#print first element of type

			    while($count<$total){
				
				print $OUT $current_data{
				    $current_feature}->[$count]->[0].
					"\t".
					$current_data{
					    $current_feature}->[$count]->[1].
						"\n";
				@{$current_data{$current_feature}->[$count]} = ();
				$count++;
			    }#end while in current_feature

			}#end not parent, already printed
		    }#end all data in hash
		    %current_data = ();

#		    foreach my $current_details_feature(%print_details{}){
#			foreach my $current_details(
#			    %print_detials{$current_details_feature}{}){
#
#			    print OUT "\t\t\t$current_details_feature\t";
#			    print OUT "$print_details{$current_details_feature}{current_details}\n";
#			}
#		    }



		} #end parent_feature exists
		else{ #first/only element in Feature
#		    print "Thinks first element..current_data{$parent_feature}.\n";
		    print $OUT join (
			"\t",
			$start,
			$end,
		    $feature,
		    $hash_notes{'Name'}
			)."\n";
#		    print "first/only element in Feature\n";
		    %current_data = ();
		}

#	    foreach my $printing($current_data{$parent_feature}){
#		print "$printing\n";
	#	print $OUT join ("\t", @{$printing})."\n";
		
#	    }
	    }
	    else{#no previous parent
		
		print $OUT join (
		    "\t",
		    $start,
		    $end,
		    $feature,
		    $hash_notes{'Name'}
		    )."\n";
#		print "no previous parent\n";
		 %current_data = ();
	    }
	  #  print "\n$current_parent\t";
	    if($hash_notes{'Parent'}){
#		print STDERR "Set Parrent::: current_parent = $hash_notes{'Parent'}\n\n";
		$current_parent = $hash_notes{'Parent'};
	    }
	    else{
#		print STDERR "No parent?????????? \t $start,\t$end,\t$feature,\t$hash_notes{'Name'}\t\n\n";
		$current_parent = '';
	    }
#	    %hash_notes = ();
#	    @details = ();
#	    %current_data = ();
#	    print "\n$current_parent\n\n";

	} #new or diffrent parent


    else{ #same parent, add data

#	    print "$current_parent eq $hash_notes{'Parent'}\n";


#	    print"pushing\t";
	    my @to_push = (
		$start,
		$end
		);

	    #%current_data{@feature[all to be printed data]}

	    if($current_data{$feature}->[0]->[1]){
#		print "old feature: $feature\n";
		push (@{$current_data{$feature}},[@to_push]);
	    }
	    else{
#		print "new feature: $feature\n";
		$current_data{$feature}->[0] = [];
		push (@{$current_data{$feature}},[@to_push]);
		#@to_push does not work here
		#as it will give array containing the array,
		#or it will give length not the array itself.
	    }


	}


#   foreach my $type(keys %hash_notes){
#	if(!$testing{$type}){ #find all unique types in notes
#	    print $type."\t";
#	    $testing{$type}="....................";
#	    print $testing{$type}."\n";
    ##############################################################
    #ID,#################Alias,#################length,###########
    #from_BOGAS,#########Name,##################Note,#############
    #Parent,#############Ontology_term,#########interpro2go_term,#
    #nb_exon,############eugene_evidence_code,##sifter_term#######
    ##############################################################
#	}
#    } 

#    * Column 1: Start location of feature
#    * Column 2: Stop location of feature
#    * Column 3: Feature key
#    * Column 4: Qualifier key
#    * Column 5: Qualifier value


}
#print join ("\n", @details)."\n\n";



