#!/usr/bin/perl

=head1 NAME

  phygeconfig.t
  A piece of code to test PhyGeConfig.pm

=cut

=head1 SYNOPSIS

 perl phygeconfig.t
 prove phygeconfig.t

=head1 DESCRIPTION

 Test PhyGeConfig.pm, a module to read/write configuration files.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 178;
use Test::Exception;

use File::Temp qw/ tempfile tempdir /;
use File::Spec;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1

BEGIN {
    use_ok('PhyGeConfig', qw( write_conf_file read_conf_file ));
}

## Prepare the data

my $tempdir = tempdir( CLEANUP => 1 );
my $fileconf = File::Spec->catfile('phygeconf.test.conf'); 

write_conf_file($fileconf, { paths => 1});

## TEST 2 to 5

is(-f $fileconf, 1,
   "testing write_conf_file, checking file")
    or diag("Looks like this has failed");

throws_ok { write_conf_file(undef, []) } qr/ERROR: ARRAY/,
    "TESTING DIE ERROR: when argument used for write_conf_file isnt HASHREF";

throws_ok { write_conf_file(undef, { fake => 1 }) } qr/ERROR: fake isnt/,
    "TESTING DIE ERROR: when argument used for write_conf_file isnt valid arg.";

throws_ok { write_conf_file(undef, { paths => 'fake'}) } qr/ERROR: paths does/,
    "TESTING DIE ERROR: when argument used for write_conf_file isnt valid val.";

## Get the possible fields

my %fields1st = PhyGeConfig::_confields();
my %fields2nd = ();
my %fields3rd = ();

foreach my $fd1 (keys %fields1st) {
    foreach my $fd2 (keys %{$fields1st{$fd1}}) {
	
	$fields2nd{$fd2} = 1;
	foreach my $fd3 (keys %{$fields1st{$fd1}->{$fd2}}) {
	    
	    $fields3rd{$fd3} = 1;
	}
    }
}

my %ffile1st = ();
my %ffile2nd = ();
my %ffile3rd = ();

open my $fh, '<', $fileconf;
while(<$fh>) {
    chomp($_);
    if ($_ =~ m/^<(\w+)\s*\d*>/) {
	$ffile1st{$1} = 1;
    }
    elsif ($_ =~ m/^\t<(\w.+?)>/) {
	$ffile2nd{$1} = 1;
    }
    elsif ($_ =~ m/^\t\t(\w+_?\w*)\s+=.+/) {
	$ffile3rd{$1} = 1;
    }
    elsif ($_ =~ m/^\t\t<(\w+)>/) {
	$ffile3rd{$1} = 1;
    }
}
close($fh);

## TEST 6 to 81

foreach my $arg1 (sort keys %fields1st) {
    is ($ffile1st{$arg1}, 1, 
	"testing write_conf_file, checking variables of first level ($arg1)")
	or diag("Looks like this has failed");
}

foreach my $arg2 (sort keys %fields2nd) {
    is ($ffile2nd{$arg2}, 1, 
	"testing write_conf_file, checking variables of second level ($arg2)")
	or diag("Looks like this has failed");
}

foreach my $arg3 (sort keys %fields3rd) {
    is ($ffile3rd{$arg3}, 1, 
	"testing write_conf_file, checking variables of third level ($arg3)")
	or diag("Looks like this has failed");
}


## TEST 82 to 174

my $fileconf2 = File::Spec->catfile($tempdir, 'phygeconf.example.conf'); 

open my $fh2, '<', $fileconf;
open my $fh3, '>', $fileconf2;

while(<$fh2>) {
    chomp($_);
    $_ =~ s/\t\t\t\t\t\t##\s\s+/\t\t\t/;
    $_ =~ s/## Mandatory,//;
    $_ =~ s/## Optional,//;
    $_ =~ s/ \(.+?\) //;
    $_ =~ s/example: //;
    $_ =~ s/'//g;
    $_ =~ s/\s+##.+$//;
    $_ =~ s/=\s+/= /;
    print $fh3 "$_\n";
}
close($fh2);
close($fh3);

my %config = read_conf_file($fileconf2);

## This the data.

my %exp_data = ( 
    general => {
        input_filenames => {
	    source                  => 'MyAssembly.ace',
	    source_filetype         => 'ace',
	    memberseq               => 'MyMembers.fasta',
	    memberstrain            => 'MyStrains.tab',
        },
        cluster_values => {
                evalue                  => ['<', 1e-10],
                expect                  => ['<', 1e-10],
                frac_identical          => ['>', '0.80'],
                gaps                    => ['<', 5],
                hsp_length              => ['>', 100],
                num_conserved           => ['>', 75],
                num_identical           => ['>', 75],
                score                   => ['>', 120],
                bits                    => ['>', 200],
                percent_identity        => ['>', 80],
                max_cluster_members     => ['<', 4],
        },       
        fastcluster_values => {
                query_id                => 'Str_',
                subject_id              => 'Str_',
                percent_identity        => ['>', 80],
                align_length            => ['>', 100],
                mismatches              => ['<', 5],
                gaps_openings           => ['<', 10],
                q_start                 => ['>', 10],
                q_end                   => ['<', 1000],
                s_start                 => ['>', 50],
                s_end                   => ['<', 5000],
                e_value                 => ['<', 1e-10],
                bit_score               => ['>', 250],
                max_cluster_members     => ['<', 6],
	},
    },
    path => { 
	1 => {            
	    path_data => {
		path_name                   => 'ML analysis',
	    },
	    homologous_search => {
                blast_arguments         => '-a 4 -e 1e-20',
                blast_database          => 'MyRefStrain.fasta',
                homologous_strain       => 'Str1',
                blastmatch_filter       => 'hsp_length > 100',
	    },
	    alignment => {
                program                 => 'clustalw',
		parameters => {
		    quiet                   => 1,
		    matrix                  => 'BLOSUM',
		    ktuple                  => 1,
		    pairgap                 => 1,
                },
	    },
	    prune_alignments => {
                score                   => ['>', 120],
                length                  => ['>', 100],
                num_residues            => ['>', 100],
                num_members             => ['<', 6],
                percentage_identity     => ['>', 80],
	    },
	    prune_overlaps => {
                composition             => 'Str1=1,Str2=1',
                method                  => 'ovlscore',
                evalseed                => 3,
                filter                  => 'identity > 80',
                trim                    => 1,
                removegaps              => 0,
	    },
	    distance => {
                method                  => 'Kimura',
                quiet                   => 1,
	    },
	    prune_distances => {
                composition             => 'Str1=1,Str2=1',
                min_distance            => 'Str1=Str2',
                max_distance            => 'Str1=Str2',
	    },
	    tree => {
                method                  => 'ML',
                program                 => 'dnaml',
                quiet                   => 1,
                outgroup                => 'str3',
		nj_parameters           => {
		    lowtri              => 1,
		    jumble              => 10,
		},
		ml_parameters           => {
		    trans_ratio         => 1.5
		},
	    },
	    reroot_tree => {
                midpoint                => 1,
                strainref               => 'Str3',
                longest                 => 0,
	    },
	    bootstrapping => {
		replicates              => 100,
                quiet                   => 1,
                outgroup                => 'str3',
                filter                  => 75,
		midpoint                => 0,
		normalized              => 1,
	    }, 
	    topoanalysis => {
		branch_cutoff               => '0.01-1',
	    },
	    allele_identification => {
		target                  => 'str1',
		parents                 => 'str2=str1-2,str3=str1-3',
		filter_length           => 100,
		filter_identity         => 75,
	    },
	},
    },
    );

## Now it will check that every conf. has these values

foreach my $par1lv (sort keys %exp_data) {
    
    is(ref($config{$par1lv}), 'HASH', 
	"testing read_conf_file, checking 1st level arg ($par1lv)")
	or diag("Looks like this has failed");

    foreach my $par2lv (sort keys %{$exp_data{$par1lv}}) {
	
	is(ref($config{$par1lv}->{$par2lv}), 'HASH', 
	"testing read_conf_file, checking 2st level arg ($par2lv)")
	or diag("Looks like this has failed");
 
	foreach my $par3lv (sort keys %{$exp_data{$par1lv}->{$par2lv}}) {
	
	    if ($par1lv eq 'general') {
		
		my $expval3 = $exp_data{$par1lv}->{$par2lv}->{$par3lv};
		my $obtval3 = $config{$par1lv}->{$par2lv}->{$par3lv};

		if (ref($expval3) eq 'ARRAY' && ref($obtval3) eq 'ARRAY' ) {

		    is(join('', @{$obtval3}), join('', @{$expval3}), 
		       "testing read_conf_file, checking 3st ($par3lv)")
			or diag("Looks like this has failed");
		}
		else {
		    is($obtval3, $expval3, 
		       "testing read_conf_file, checking 3st level ($par3lv)")
			or diag("Looks like this has failed");
		}
	    }
	    else {

		is(ref($config{$par1lv}->{$par2lv}->{$par3lv}), 'HASH', 
		   "testing read_conf_file, checking 3st level value ($par3lv)")
		    or diag("Looks like this has failed");
	    
	    
		my %expdata4lv = %{$exp_data{$par1lv}->{$par2lv}->{$par3lv}};
		foreach my $par4lv (sort keys %expdata4lv) {
		    
		    if ($par3lv eq 'alignment' && $par4lv eq 'parameters' ||
			$par3lv eq 'tree' && $par4lv eq 'nj_parameters'   ||
			$par3lv eq 'tree' && $par4lv eq 'ml_parameters') {
		
			my %expdata5lv = %{$expdata4lv{$par4lv}};
			foreach my $par5lv (sort keys %expdata5lv) {
			
			    my $expval5 = $expdata5lv{$par5lv};
			    my $obtval5 = $config{$par1lv}->{$par2lv}
			                                  ->{$par3lv}
			                                  ->{$par4lv}
			                                  ->{$par5lv};

			
			    is($obtval5, $expval5, 
			       "testing read_conf_file, checking 5st ($par5lv)")
				or diag("Looks like this has failed");
			}
		    }
		    else {
		
			my $expval4 = $expdata4lv{$par4lv};
			my $obtval4 = $config{$par1lv}->{$par2lv}
                                                      ->{$par3lv}
                                                      ->{$par4lv};

			if (ref($expval4) eq 'ARRAY' &&
			    ref($obtval4) eq 'ARRAY' ) {

			    is(join('', @{$obtval4}), join('', @{$expval4}), 
			       "testing read_conf_file, checking 4st ($par4lv)")
				or diag("Looks like this has failed");
			}
			else {
			    is($obtval4, $expval4, 
			       "testing read_conf_file, checking 4st ($par4lv)")
				or diag("Looks like this has failed");
			}
		    }
		}
	    }
	}
    }
}

## TEST 175

throws_ok { read_conf_file() } qr/ERROR: No argument/,
    "TESTING DIE ERROR: when no argument was used for read_conf_file";

## test some regexp, TEST 176 to 178

my $fileconf3 = File::Spec->catfile($tempdir,, 'phygeconf.failcheck1.conf'); 

open my $fh4, '<', $fileconf2;
open my $fh5, '>', $fileconf3;

while(<$fh4>) {
    chomp($_);
    $_ =~ s/source\s+=\sMyAssembly\.ace//;
    print $fh5 "$_\n";
}
close($fh4);
close($fh5);

throws_ok { read_conf_file($fileconf3) } qr/general>input_filenames>source ma/,
    "TESTING DIE ERROR: when mandatory arg. didnt was used for read_conf_file";

my $fileconf4 = File::Spec->catfile($tempdir,, 'phygeconf.failcheck2.conf'); 

open my $fh6, '<', $fileconf2;
open my $fh7, '>', $fileconf4;

while(<$fh6>) {
    chomp($_);
    $_ =~ s/evalseed\s+=\s\d+/evalseed\t\t=\tfake/;
    print $fh7 "$_\n";
}
close($fh6);
close($fh7);

throws_ok { read_conf_file($fileconf4) } qr/path 1>prune_overlaps>evalseed/,
    "TESTING DIE ERROR: when wrong value (int) is used for read_conf_file";

my $fileconf5 = File::Spec->catfile($tempdir,, 'phygeconf.failcheck3.conf'); 

open my $fh8, '<', $fileconf2;
open my $fh9, '>', $fileconf5;

while(<$fh8>) {
    chomp($_);
    $_ =~ s/trim\s+=\s\d+/trim\t\t=\t8/;
    print $fh9 "$_\n";
}
close($fh8);
close($fh9);

throws_ok { read_conf_file($fileconf5) } qr/path 1>prune_overlaps>trim/,
    "TESTING DIE ERROR: when wrong value (boolean) is used for read_conf_file";

####
1; #
####
