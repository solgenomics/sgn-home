
=head1

verify_accessions_in_pedigrees_exist_in_database.pl

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

my $visible_to_role_cvterm = $schema->resultset("Cv::Cvterm")->create_with(
    { name   => 'visible_to_role',
      cv => 'local',
      db => 'null',
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
my $experiment_name = "loading_crosses";
my $location = "unspecified";


eval {
    foreach my $num (@rows ) {
	my $stock_name = $spreadsheet->value_at($num , "accession_name");
	#my $stock_name = "tester";
        my $geolocation = $schema->resultset("NaturalDiversity::NdGeolocation")->find_or_create(
           {
                description => $location,
           } ) ;

	my $project = $schema->resultset("Project::Project")->find_or_create(
            {
                name => $project_name,
                description => "testing pedigree loading",
            } ) ;
	my $accession_stock = $schema->resultset("Stock::Stock")->find_or_create(
            { organism_id => $organism_id,
              name       => $stock_name,
              uniquename => $stock_name,
              type_id     => $accession_cvterm->cvterm_id,
            } );

	my $accession_stock_prop = $schema->resultset("Stock::Stockprop")->find_or_create(
	       { type_id =>$visible_to_role_cvterm->cvterm_id(),
		 value => "IITA",
		 stock_id => $accession_stock->stock_id()
		 });

	# $accession_stock->find_or_create_related('stock_relationship_objects', {
	# 	    type_id => $female_parent_of->cvterm_id(),
	# 	    subject_id => $female_accession_stock->stock_id(),
	# 					  } );
	#add the owner for this stock
	#check first if it exists
	my $owner_insert = "INSERT INTO phenome.stock_owner (sp_person_id, stock_id) VALUES (?,?)";
	my $sth = $dbh->prepare($owner_insert);
	my $check_query = "SELECT sp_person_id FROM phenome.stock_owner WHERE ( sp_person_id = ? AND stock_id = ? )";
	# my $person_ids = $dbh->selectcol_arrayref($check_query, undef, ($sp_person_id, $plot_stock->stock_id) );
	my $person_ids = $dbh->selectcol_arrayref($check_query, undef, ( $sp_person_id, $accession_stock->stock_id) );
	if (!@$person_ids) {
	$sth->execute($sp_person_id, $accession_stock->stock_id);
	}
	my $experiment = $schema->resultset('NaturalDiversity::NdExperiment')->create(
            {
                nd_geolocation_id => $geolocation->nd_geolocation_id(),
                type_id => $accession_cvterm->cvterm_id(),
            } );
	#link to the project
	$experiment->find_or_create_related('nd_experiment_projects', {
	    project_id => $project->project_id()
                                            } );
	#link the experiment to the stock
	$experiment->find_or_create_related('nd_experiment_stocks' , {
	    stock_id => $accession_stock->stock_id(),
	    type_id  =>  $accession_cvterm->cvterm_id(),
                                           });
    }
};

if ($@) { print "An error occured! Rolling backl!\n\n $@ \n\n "; }
else {
    print "Transaction succeeded! Commiting ! \n\n";
    $dbh->commit();
}
