
package PhyGeConfig;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Config::General;

## To export some functions

use Exporter qw( import );

our @EXPORT_OK = qw( write_conf_file read_conf_file );


###############
### PERLDOC ###
###############

=head1 NAME

PhyGeConfig.pm
A module to write/read/check the configuration file for phygomics.pl

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

  use PhyGeConf qw( write_conf_file read_conf_file);

  ## To print an empty config file with one path.

  write_conf_file($filename);

  ## To print an example with 3 paths

  write_conf_file($filename, { examples => 1, path => 3 });

  ## To read a configuration file

  my %conf = read_conf_file($filename);


=head1 DESCRIPTION

 A module to write/read/check the configuration file for phygomics.pl


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 




=head2 write_conf_file

  Usage: write_conf_file($filename, $args_href);

  Desc: Create an empty phygomics configuration file

  Ret: None

  Args: $filename, name of the configuration file (phygomics.conf by default)
        $args_href, a hash ref. with the following pairs:
         example => $boolean (to print with examples)
         paths   => $int (number of path to print in the conf. file)        

  Side_Effects: None

  Example:  write_conf_file($filename, $args_href);

=cut

sub write_conf_file {
    my $filename = shift || 'phygomics.conf';
    my $args_href = shift;

    my %args = ();
    if (defined $args_href && length($args_href) > 0) {
	if (ref($args_href) ne 'HASH') {
	    croak("ERROR: $args_href supplied to write_conf_file isnt HASHREF");
	}
	else {
	    %args = %{$args_href};
	}
    }

    my %perm = (
	example => '(0|1)',
	paths   => '\d+',
	);

    foreach my $arg (keys %args) {
	unless (exists $perm{$arg}) {
	    croak("ERROR: $arg isnt a valid argument for write_conf_file.");
	}
	else {
	    if ($args{$arg} !~ m/^$perm{$arg}$/) {
		croak("ERROR: $arg doesnt have valid val. for write_conf_file");
	    }
	}
    }

    ## Get the fields
    
    my %fields = _confields();

    ## Create the conf. file and write it.

    open my $cfh, '>', $filename;
    
    my $date = `date`;
    chomp($date);

    ## Get <general>

    print $cfh qq/

#####################################
##        CONFIGURATION FILE       ##
##      FOR PhyGomics Pipeline     ##
## $date    ##
#####################################
## Configuration file tips:        ##
##  + use '#' for comment          ## 
##  + optional variables can be    ##
##    deleted or just leave empty  ##
##  + (<) for a variable means     ##
##    it will use as lower cutoff  ##
##  + (>) for a variable means     ##
##    it will use as upper cutoff  ##
#####################################

#####################################
## GENERAL CONFIGURATION VARIABLES ##
#####################################
<general>
\t<input_filenames>
\t\tsource\t\t\t=\t## Mandatory, example: 'MyAssembly.ace'
\t\tsource_filetype\t\t=\t## Mandatory, example: 'ace'
\t\tmemberseq\t\t=\t## Mandatory, example: 'MyMembers.fasta'
\t\tmemberstrain\t\t=\t## Mandatory, example: 'MyStrains.tab'
\t<\/input_filenames>

\t<cluster_values>\t\t\t## INCOMPATIBLE w. <fastcluster_values>
\t\tevalue\t\t\t=\t## Optional, (<) example: '1e-10'
\t\texpect\t\t\t=\t## Optional, (<) example: '1e-10'
\t\tfrac_identical\t\t=\t## Optional, (>) example: '0.80'
\t\tgaps\t\t\t=\t## Optional, (<) example: '5'
\t\thsp_length\t\t=\t## Optional, (>) example: '100'
\t\tnum_conserved\t\t=\t## Optional, (>) example: '75'
\t\tnum_identical\t\t=\t## Optional, (>) example: '75'
\t\tscore\t\t\t=\t## Optional, (>) example: '120'
\t\tbits\t\t\t=\t## Optional, (>) example: '200'
\t\tpercent_identity\t=\t## Optional, (>) example: '80'
\t\tmax_cluster_members\t=\t## Optional, (<) example: '4'
\t<\/cluster_values>
	
\t<fastcluster_values>\t\t\t## INCOMPATIBLE W. <cluster_values>
\t\tquery_id\t\t=\t## Optional, example: 'Str_'
\t\tsubject_id\t\t=\t## Optional, example: 'Str_'
\t\tpercent_identity\t=\t## Optional, (>) example: '80'
\t\talign_length\t\t=\t## Optional, (>) example: '100'
\t\tmismatches\t\t=\t## Optional, (<) example: '5'
\t\tgaps_openings\t\t=\t## Optional, (<) example: '10'
\t\tq_start\t\t\t=\t## Optional, (>) example: '10'
\t\tq_end\t\t\t=\t## Optional, (<) example: '1000'
\t\ts_start\t\t\t=\t## Optional, (>) example: '50'
\t\ts_end\t\t\t=\t## Optional, (<) example: '5000'
\t\te_value\t\t\t=\t## Optional, (<) example: '1e-10'
\t\tbit_score\t\t=\t## Optional, (>) example: '250'
\t\tmax_cluster_members\t=\t## Optional, (<) example: '6'
\t<\/fastcluster_values>
<\/general>

/;

print $cfh qq/

##################################
## PATH CONFIGURATION VARIABLES ##
##################################	

/;

    my $path_n = $args{paths} || 1;

    my $n = 0;
    while($n < $path_n) {
	$n++;

	print $cfh qq/

<path $n>	    
\t<path_data>\t\t\t\t## Mandatory section
\t\tpath_name\t\t=\t## Mandatory, example: 'ML analysis'
\t<\/path_data>

\t<homologous_search>\t\t\t## Optional section
\t\tblast_arguments\t\t=\t## Optional, example: '-a 4 -e 1e-20'
\t\tblast_database\t\t=\t## Optional, example: 'MyRefStrain.fasta'
\t\thomologous_strain\t=\t## Optional, example: 'Str1'
\t\tblastmatch_filter\t=\t## Optional, example: 'hsp_length > 100'
\t<\/homologous_search>

\t<alignment>\t\t\t\t## Mandatory section
\t\tprogram\t\t\t=\t## Mandatory, example: 'clustalw'
\t\t<parameters>\t\t\t## Optional check different alignment
\t\t\t\t\t\t## programs to get the possible 
\t\t\t\t\t\t## parameters (quiet)
\t\t\t\t\t\t## example:
\t\t\t\t\t\t##   quiet   = 1
\t\t\t\t\t\t##   matrix  = BLOSUM
\t\t\t\t\t\t##   ktuple  = 1
\t\t\t\t\t\t##   pairgap = 1
\t\t<\/parameters>\t\t\t##   maxdiv  = 99
\t<\/alignment>

\t<prune_alignments>\t\t\t## Optional section
\t\tscore\t\t\t=\t## Optional, (>) example: '120'
\t\tlength\t\t\t=\t## Optional, (>) example: '100'
\t\tnum_residues\t\t=\t## Optional, (>) example: '100'
\t\tnum_members\t\t=\t## Optional, (<) example: '6'
\t\tpercentage_identity\t=\t## Optional, (>) example: '80'
\t<\/prune_alignments>

\t\<prune_overlaps>\t\t\t## Optional section
\t\tcomposition\t\t=\t## Optional, example: 'Str1=1,Str2=1'
\t\tmethod\t\t\t=\t## Optional, example: 'ovlscore'
\t\tevalseed\t\t=\t## Optional, example: '3'
\t\tfilter\t\t\t=\t## Optional, example: 'identity > 80'
\t\ttrim\t\t\t=\t## Optional, example: '1'
\t\tremovegaps\t\t=\t## Optional, example: '0'
\t<\/prune_overlaps>

\t<distance>\t\t\t\t## Mandatory section (for NJ\/UPGMA)
\t\tmethod\t\t\t=\t## Optional, example: 'Kimura'
\t\tquiet\t\t\t=\t## Optional, example: '1'
\t<\/distance>

\t<prune_distances>\t\t\t## Optional section
\t\tcomposition\t\t=\t## Optional, example: 'Str1=1,Str2=1'
\t\tmin_distance\t\t=\t## Optional, example: 'Str1=Str2'
\t\tmax_distance\t\t=\t## Optional, example: 'Str1=Str2'
\t<\/prune_distances>

\t<tree>\t\t\t\t\t## Mandatory section
\t\tmethod\t\t\t=\t## Mandatory, example: 'ML'
\t\tprogram\t\t\t=\t## Optional, (only ML) example: 'dnaml'
\t\tquiet\t\t\t=\t## Optional, example: '1'
\t\toutgroup\t\t=\t## Optional, example: 'str3'
\t\t<nj_parameters>\t\t\t## Optional parameters to pass to phylip
\t\t\t\t\t\t## example:
\t\t\t\t\t\t##   lowtri  = 1
\t\t\t\t\t\t##   jumble  = 10
\t\t<\/nj_parameters>
\t\t<ml_parameters>\t\t\t## Optional parameters to pass to phylip or dnaml
\t\t\t\t\t\t## example:
\t\t\t\t\t\t##   trans_ratio  = 1.5
\t\t<\/ml_parameters>
\t<\/tree>

\t<reroot_tree>\t\t\t\t## Optional section
\t\tmidpoint\t\t=\t## Optional, example: '1'
\t\tstrainref\t\t=\t## Optional, example: 'Str3'
\t\tlongest\t\t\t=\t## Optional, example: '0'
\t<\/reroot_tree>

\t<bootstrapping>\t\t\t\t## Optional section
\t\treplicates\t\t=\t## Optional, example: '100'
\t\tquiet\t\t\t=\t## Optional, example: '1'
\t\toutgroup\t\t=\t## Optional, example: 'str3'
\t\tfilter\t\t\t=\t## Optional, (>) example: '75'
\t\tmidpoint\t\t=\t## Optional, example: '0'
\t\tnormalized\t\t=\t## Optional, example: '1'
\t<\/bootstrapping>
 
\t<topoanalysis>\t\t\t\t## Optional section
\t\tbranch_cutoff\t\t=\t## Optional, example: '0.01-1'
\t<\/topoanalysis>

\t<allele_identification>\t\t\t## Optional section
\t\ttarget\t\t\t=\t## Optional, example: 'str1'
\t\tparents\t\t\t=\t## Optional, example: 'str2=str1-2,str3=str1-3'
\t\tfilter_length\t\t=\t## Optional, (>) example: '100'
\t\tfilter_identity\t\t=\t## Optional, (>) example: '75'
\t<\/allele_identification>

\t<codeml_analysis>\t\t\t## Optional section
\t\tsource_cds\t\t=\t## Optional, example: 'file'
\t\tfile_cds\t\t=\t## Optional, example: 'MyCDS.fasta'
\t\t<codeml_parameters>\t\t## Optional parameters to pass to codeml
\t\t\t\t\t\t## example:
\t\t\t\t\t\t##   model  = 3
\t\t<\/codeml_parameters>
\t<\/codeml_analysis>

<\/path $n>

/;
    }
    
    close($cfh);
}



=head2 _confields

  Usage: my %fields = _confields();

  Desc: Function with the field data

  Ret: An array with the field data

  Args: None   

  Side_Effects: None

  Example: my %fields = _confields();

=cut

sub _confields {

    ## Define default vars.
    
    my %int_opt = ( regexp => '\d+', request => 'optional',  default => undef );
    my %int_man = ( regexp => '\d+', request => 'mandatory', default => undef );
    my %wrd_opt = ( regexp => '\w+', request => 'optional',  default => undef );
    my %wrd_man = ( regexp => '\w+', request => 'mandatory', default => undef );
    my %fil_man = ( regexp => '.+',  request => 'mandatory', default => undef );
    my %fre_opt = ( regexp => '.+',  request => 'optional',  default => undef );
    my %fre_man = ( regexp => '.+',  request => 'mandatory', default => undef );
    my %bol_opt = ( regexp  => '(1|0)',  
		    request => 'optional', 
		    default => undef );
    

    ## Create the field hash

    my %fields = (
	general => {
	    input_filenames    => { 
		source              => \%fil_man,
		source_filetype     => { regexp  => '(ace|blast)', 
					 request => 'mandatory', 
					 default => undef },
		memberseq           => \%fil_man,
		memberstrain        => \%fil_man,
	    },
	    cluster_values     => {
		evalue              => { regexp  => '.+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef },
		expect              => { regexp  => '.+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef },
		frac_identical      => { regexp  => '\d\.?\d*', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		gaps                => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef },
		hsp_length          => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => 100 },
		num_conserved       => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		num_identical       => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		score               => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		bits                => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		percent_identity    => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => 75 },
		max_cluster_members => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => 20 },
	    },
	    fastcluster_values => {
		query_id            => \%fre_opt, 
		subject_id          => \%fre_opt, 
		percent_identity    => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => 75 }, 
		align_length        => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => 100 },
		mismatches          => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef }, 
		gaps_openings       => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef }, 
		q_start             => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef }, 
		q_end               => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef }, 
		s_start             => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },  
		s_end               => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef }, 
		e_value             => { regexp  => '.+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef }, 
		bit_score           => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef }, 
		max_cluster_members => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => 20 },
	    },
	},
	path => {
	    path_data => {
		path_name  => \%fre_man,
	    },
	    homologous_search => {
		blast_arguments   => \%fre_opt,
		blast_database    => \%fre_opt,
		homologous_strain => \%fre_opt,
		blastmatch_filter => \%fre_opt,
	    },
	    alignment => {
		program    => { regexp  => '\w+', 
				request => 'mandatory', 
				default => 'clustalw' },
		parameters => { regexp  => undef,
				request => 'optional' },
	    },
	    prune_alignments => {
		score               => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		length              => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		num_residues        => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
		num_members         => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '<',
					 default => undef },
		percentage_identity => { regexp  => '\d+', 
					 request => 'optional',
					 cutoff  => '>',
					 default => undef },
	    },
	    prune_overlaps => {
		composition => \%fre_opt,
		method      => \%fre_opt,
		evalseed    => \%int_opt,
		filter      => \%fre_opt,
		trim        => \%bol_opt,
		removegaps  => \%bol_opt,
	    },
	    distance => {
		method => { regexp  => '\w+', 
			    request => 'optional', 
			    default => 'Kimura' },
		quiet => \%bol_opt,
	    },
	    prune_distances => {
		composition  => \%fre_opt,
		min_distance => \%fre_opt,
		max_distance => \%fre_opt,
	    },
	    tree => {
		method          => { regexp  => '(NJ|UPGMA|ML)', 
				     request => 'mandatory', 
				     default => 'ML' },
		program         => \%wrd_opt,
                outgroup        => \%fre_opt,
		quiet           => \%bol_opt,
		nj_parameters   => { regexp  => undef,
				     request => 'optional' },
		ml_parameters   => { regexp  => undef,
				     request => 'optional' },
	    },
	    reroot_tree => {
		midpoint  => \%bol_opt,
		strainref => \%fre_opt,
		longest   => \%bol_opt,
	    },
	    bootstrapping => {
		replicates => \%int_opt,
		quiet      => \%bol_opt,
		outgroup   => \%fre_opt,
		midpoint   => \%bol_opt,
		normalized => \%bol_opt,
		filter     => \%int_opt,
	    },
	    topoanalysis => {
		branch_cutoff => \%fre_opt,
	    },
	    allele_identification => {
		target          => \%fre_opt,
                parents         => \%fre_opt,
		filter_length   => \%fre_opt,
                filter_identity => \%fre_opt,
	    },
	    codeml_analysis => {
		source_cds        => \%wrd_opt,
                file_cds          => { regexp => '.+',  
				       request => 'optional', 
				       default => undef },
		codeml_parameters => { regexp  => undef,
				       request => 'optional' },
	    },
	},
	);
    return %fields;
}


=head2 read_conf_file

  Usage: my %conf = read_conf_file($filename); 
         my %conf = read_conf_file($conf_options_href);

  Desc: Read the configuration file and return the configuration variables

  Ret: %conf, a hash with the configuration variables as a hash of hashes
       %conf = ( first_level => { second_level => { third_level => $value } })
       (see Config::General for more information)

  Args: $filename, name of the configuration file (phygomics.conf by default)
        $conf_href, a hash ref. with the Config::General options 
        (for example -LowerCaseName or -ConfigFile)       

  Side_Effects: Die if no argument is used.
                Check the arguments and die if a non permited parameter
                is used.

  Example:  my %conf = read_conf_file($filename);

=cut

sub read_conf_file {
    my $arg = shift ||
	croak("ERROR: No argument was supplied to read_conf_file() function.");

    my %conf = Config::General->new($arg)->getall();
    my %fields = _confields();

    ## Check the variables

    my $err = "ERROR reading conf. file: ";

    foreach my $lv1 (sort keys %conf) {
    
	unless (exists $fields{$lv1}) {
	    croak($err . " $lv1 isnt a permitted level 1 config.");
	}

	if ($lv1 ne 'path') {
	
	    foreach my $lv2 (sort keys %{$conf{$lv1}}) {
	    
		unless (exists $fields{$lv1}->{$lv2}) {
		    $err .= "$lv1>$lv2";
		    croak($err . " isnt a permitted level 2 config.");
		}

		foreach my $lv3 (sort keys %{$conf{$lv1}->{$lv2}}) {
	    
		    unless (exists $fields{$lv1}->{$lv2}->{$lv3}) {
			$err .= "$lv1>$lv2>$lv3";
			croak($err . " isnt a permitted level3 config.");
		    }

		    my $regexp = $fields{$lv1}->{$lv2}->{$lv3}->{regexp};
		    my $request = $fields{$lv1}->{$lv2}->{$lv3}->{request};
		    my $cutoff = $fields{$lv1}->{$lv2}->{$lv3}->{cutoff};

		    if ($request eq 'mandatory') {
			unless (defined $conf{$lv1}->{$lv2}->{$lv3}) {
			    $err .= "$lv1>$lv2>$lv3";
			    $err .= " mandatory config. parameter was not";
			    croak($err . " supplied.");
			}
		    }

		    if (defined $regexp) {
		    
			if ($conf{$lv1}->{$lv2}->{$lv3} !~ m/^$regexp$/) {
			    $err .= "$lv1>$lv2>$lv3";
			    $err .= " doesnt have a permitted value ($regexp)";
			    croak($err);
			}
			else {  ## Modify the argument to add the cutoff
			    
			    if (defined $cutoff) {
				my $val = $conf{$lv1}->{$lv2}->{$lv3};
				$conf{$lv1}->{$lv2}->{$lv3} = [$cutoff, $val];
			    }			
			}
		    }
		}
	    }
	}
	else {
	    
	    foreach my $p (sort keys %{$conf{$lv1}}) {
	    
		foreach my $lv2 (sort keys %{$conf{$lv1}->{$p}}) {

		    unless (exists $fields{$lv1}->{$lv2}) {
			$err .= "path $p>$lv2";
			croak($err . " isnt a permitted level 2 config.");
		    }

		    foreach my $lv3 (sort keys %{$conf{$lv1}->{$p}->{$lv2}}) {
			
			unless (exists $fields{$lv1}->{$lv2}->{$lv3}) {
			    $err .= "path $p>$lv2>$lv3";
			    croak($err . " isnt a permitted level3 config.");
			}

			my $regexp = $fields{$lv1}->{$lv2}->{$lv3}->{regexp};
			my $request = $fields{$lv1}->{$lv2}->{$lv3}->{request};
			my $cutoff = $fields{$lv1}->{$lv2}->{$lv3}->{cutoff};

			if ($request eq 'mandatory') {
			    unless (defined $conf{$lv1}->{$p}->{$lv2}->{$lv3}) {
				$err .= "path $p>$lv2>$lv3";
				$err .= " mandatory config. parameter is undef";
				croak($err);
			    }
			}
			if (defined $regexp) {
			    
			    my $val = $conf{$lv1}->{$p}->{$lv2}->{$lv3};
			    if ($val !~ m/^$regexp$/) {
				$err .= "path $p>$lv2>$lv3 (value=$val)";
				$err .=" doesnt have permitted value ($regexp)";
				croak($err);
			    }
			    else {  ## Modify the argument to add the cutoff
			    
				if (defined $cutoff) {
				    my $cutval = [$cutoff,$val];
				    $conf{$lv1}->{$p}->{$lv2}->{$lv3} = $cutval;
				}			
			    }
			}
		    }
		}
	    }
	}
    }

    ## Check if the mandatory variables are present

    foreach my $fv1 (sort keys %fields) {
	
	if ($fv1 ne 'path') {
	    foreach my $fv2 (sort keys %{$fields{$fv1}}) {
		foreach my $fv3 (sort keys %{$fields{$fv1}->{$fv2}}) {
		    my $req = $fields{$fv1}->{$fv2}->{$fv3}->{request};
		    if ($req eq 'mandatory') {
			unless (defined $conf{$fv1}->{$fv2}->{$fv3}) {
			    $err .= "$fv1>$fv2>$fv3 mandatory config parameter";
			    croak($err . " was not supplied.");
			}
		    }
		}
	    }
	}
	else {
	    foreach my $fv2 (sort keys %{$fields{$fv1}}) {
		foreach my $fv3 (sort keys %{$fields{$fv1}->{$fv2}}) {
		    
		    my $req = $fields{$fv1}->{$fv2}->{$fv3}->{request};
		    if ($req eq 'mandatory') {
			my @paths = keys %{$conf{path}};

			foreach my $p (sort @paths) {
			    unless (defined $conf{$fv1}->{$p}->{$fv2}->{$fv3}) {
				$err .= "$fv1 $p>$fv2>$fv3 mandatory config";
				croak($err . " parameter was not supplied.");
			    }
			}
		    }
		}
	    }
	}
    }
    
    return %conf;
}






####
1; #
####
