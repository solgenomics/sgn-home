
=head1

load_pedigrees_in_database.pl

=head1 SYNOPSIS

    $this_script.pl -H [dbhost] -D [dbname] [-t]

=head1 COMMAND-LINE OPTIONS

 -H  host name
 -D  database name
 -i infile
 -u sgn user name
 -t  Test run . Rolling back at the end.


=head2 DESCRIPTION



=head2 AUTHORS

Jeremy Edwards (jde22@cornell.edu)
Naama Menda (nm249@cornell.edu)

June 2012

=cut


#!/usr/bin/perl
use strict;
use Getopt::Std;
use CXGN::Tools::File::Spreadsheet;
use CXGN::People::Person;

use Bio::Chado::Schema;
use CXGN::DB::InsertDBH;
use Carp qw /croak/ ;
##
##


our ($opt_H, $opt_D, $opt_i, $opt_t, $opt_u);

getopts('H:i:tD:u:');

my $dbhost = $opt_H;
my $dbname = $opt_D;
my $file = $opt_i;



my $dbh = CXGN::DB::InsertDBH->new( { dbhost=>$dbhost,
				      dbname=>$dbname,
				      dbargs => {AutoCommit => 0,
						 RaiseError => 1}
				    }
    );
my $schema= Bio::Chado::Schema->connect(  sub { $dbh->get_actual_dbh() } );
$dbh->do('SET search_path TO public');

#getting the last database ids for resetting at the end in case of rolling back
###############
my $last_nd_experiment_id = $schema->resultset('NaturalDiversity::NdExperiment')->get_column('nd_experiment_id')->max;
my $last_cvterm_id = $schema->resultset('Cv::Cvterm')->get_column('cvterm_id')->max;
my $last_nd_experiment_project_id = $schema->resultset('NaturalDiversity::NdExperimentProject')->get_column('nd_experiment_project_id')->max;
my $last_nd_experiment_stock_id = $schema->resultset('NaturalDiversity::NdExperimentStock')->get_column('nd_experiment_stock_id')->max;
#my $last_nd_experiment_phenotype_id = $schema->resultset('NaturalDiversity::NdExperimentPhenotype')->get_column('nd_experiment_phenotype_id')->max;
#my $last_phenotype_id = $schema->resultset('Phenotype::Phenotype')->get_column('phenotype_id')->max;
my $last_stock_id = $schema->resultset('Stock::Stock')->get_column('stock_id')->max;
my $last_stock_relationship_id = $schema->resultset('Stock::StockRelationship')->get_column('stock_relationship_id')->max;
my $last_project_id = $schema->resultset('Project::Project')->get_column('project_id')->max;
my $last_nd_geolocation_id = $schema->resultset('NaturalDiversity::NdGeolocation')->get_column('nd_geolocation_id')->max;
my $last_geoprop_id = $schema->resultset('NaturalDiversity::NdGeolocationprop')->get_column('nd_geolocationprop_id')->max;
my $last_projectprop_id = $schema->resultset('Project::Projectprop')->get_column('projectprop_id')->max;


my %seq  = (
    'nd_experiment_nd_experiment_id_seq' => $last_nd_experiment_id,
    'cvterm_cvterm_id_seq' => $last_cvterm_id,
    'nd_experiment_project_nd_experiment_project_id_seq' => $last_nd_experiment_project_id,
    'nd_experiment_stock_nd_experiment_stock_id_seq' => $last_nd_experiment_stock_id,
#    'nd_experiment_phenotype_nd_experiment_phenotype_id_seq' => $last_nd_experiment_phenotype_id,
#    'phenotype_phenotype_id_seq' => $last_phenotype_id,
    'stock_stock_id_seq'         => $last_stock_id,
    'stock_relationship_stock_relationship_id_seq'  => $last_stock_relationship_id,
    'project_project_id_seq'     => $last_project_id,
    'nd_geolocation_nd_geolocation_id_seq'          => $last_nd_geolocation_id,
    'nd_geolocationprop_nd_geolocationprop_id_seq'  => $last_geoprop_id,
    'projectprop_projectprop_id_seq'                => $last_projectprop_id,
    );

# find the cvterm for a genotyping experiment
my $pedigree_cvterm = $schema->resultset('Cv::Cvterm')->create_with(
    { name   => 'pedigree experiment',
      cv     => 'experiment type',
      db     => 'null',
      dbxref => 'pedigree experiment',
    });

#my $username = $opt_u || 'kulakow' ; #'cassavabase';
my $username = 'kulakow' ; #'cassavabase';
my $sp_person_id= CXGN::People::Person->get_person_by_username($dbh, $username);

die "User $username for cassavabase must be pre-loaded in the database! \n" if !$sp_person_id ;




my $accession_cvterm = $schema->resultset("Cv::Cvterm")->create_with(
    { name   => 'accession',
      cv     => 'stock type',
      db     => 'null',
      dbxref => 'accession',
    });





my $female_parent = $schema->resultset("Cv::Cvterm")->create_with(
    { name   => 'female_parent',
      cv     => 'stock relationship',
      db     => 'null',
      dbxref => 'female_parent',
    });

my $male_parent = $schema->resultset("Cv::Cvterm")->create_with(
    { name   => 'male_parent',
      cv     => 'stock relationship',
      db     => 'null',
      dbxref => 'male_parent',
    });



########################

#new spreadsheet, skip first column
my $spreadsheet=CXGN::Tools::File::Spreadsheet->new($file, 0);


my $organism = $schema->resultset("Organism::Organism")->find_or_create(
    {
	genus   => 'Manihot',
	species => 'Manihot esculenta',
    } );
my $organism_id = $organism->organism_id();

my @rows = $spreadsheet->row_labels();
my @columns = $spreadsheet->column_labels();

my $project_name = "IITA_cassava_breeding_program";
my $experiment_name = "pedigree_load_experiment";
my $location = "SGN";


eval {
    foreach my $num (@rows ) {
	my $stock_name = $spreadsheet->value_at($num , "accession_name");
	my $cross = $spreadsheet->value_at($num , "pedigree");
	my $cross_character = '/';
	my @parents = split($cross_character,$cross);

	if (!defined($parents[1])) {
	    print STDERR "Incorrectly formatted cross for: ".$stock_name."\t".$cross."\n";
	    next;
	}
	
	print STDERR $stock_name."\t".$parents[0]."\t".$parents[1]."\n";

	my $accession_stock = $schema->resultset("Stock::Stock")->find(
	    { name => $stock_name,
	    } );

	print "accession stock name: ".$accession_stock->name()."\n";

	my $female_parent_stock = $schema->resultset("Stock::Stock")->find(
            { name       => $parents[0],
            } );

	print STDERR $female_parent_stock->name()."\n".$female_parent->cvterm_id()."\n".$female_parent_stock->stock_id()."\n";
	

	$accession_stock->find_or_create_related('stock_relationship_objects', {
		type_id => $female_parent->cvterm_id(),
		object_id => $accession_stock->stock_id(),
		subject_id => $female_parent_stock->stock_id(),
	 					  } );



	my $male_parent_name = $parents[1];
	if ($parents[1] eq '.') {
	    $male_parent_name = $parents[0];
	}


	#print STDERR $male_parent_stock->name()."\n".$male_parent->cvterm_id()."\n".$male_parent_stock->stock_id()."\n"; 

	unless ($parents[1] eq '?') {

	    my $male_parent_stock = $schema->resultset("Stock::Stock")->find(
		{ uniquename       => $male_parent_name,
		} );
	    
	    
	    $accession_stock->find_or_create_related('stock_relationship_objects', {
		type_id => $male_parent->cvterm_id(),
		object_id => $accession_stock->stock_id(),
		subject_id => $male_parent_stock->stock_id(),
						  } );

	}
    }
};

if ($@) { print "An error occured! Rolling backl!\n\n $@ \n\n "; }
else {
    print "Transaction succeeded! Commiting ! \n\n";
    $dbh->commit();
}
