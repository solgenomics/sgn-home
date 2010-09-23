#!/usr/bin/env perl


##########################
## 1) POD               ##
##########################

=head1 NAME

 ChadoStock_v1_7_dbversion.pl

=head1 SYNOPSIS

  example_with_dbversion.pl [options]

  Options:

    -D <dbname> (mandatory)
      dbname to load into

    -H <dbhost> (mandatory)
      dbhost to load into

    -p <script_executor_user> (mandatory)
      username to run the script

    -F force to run this script and don't stop it by 
       missing previous db_patches

  Note: If the first time that you run this script, obviously
        you have not any previous dbversion row in the md_dbversion
        table, so you need to force the execution of this script 
        using -F

=head1 DESCRIPTION

 Dbpatch to create the chado stock module with dbversion.

=head1 AUTHOR

Aureliano Bombarely<ab782@cornell.edu>
Naama Menda<nm249@cornell.edu>

=head1 COPYRIGHT & LICENSE

Copyright 2009 Boyce Thompson Institute for Plant Research

This program is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.

=cut

use strict;
use warnings;

use Pod::Usage;
use Getopt::Std;
use CXGN::DB::InsertDBH;
use CXGN::Metadata::Dbversion;   ### Module to interact with the metadata.md_dbversion table


## Declaration of the parameters used to run the script

our ($opt_H, $opt_D, $opt_p, $opt_F, $opt_h);
getopts("H:D:p:Fh");


##########################
## 2) ARGUMENT CHECKING ##
##########################

## If is used -h <help> or none parameters is detailed print pod

if (!$opt_H && !$opt_D && !$opt_p && !$opt_F && !$opt_h) {
    print STDOUT "There are n\'t any tags. Print help\n\n";
    pod2usage(1);
} 
elsif ($opt_h) {
    pod2usage(1);
} 

## Specify the mandatory parameters

if (!$opt_H || !$opt_D) {
    print STDOUT "\nMANDATORY PARAMETER ERROR: -D <db_name> or/and -H <db_host> parameters has not been specified for $patch_name.\n";
    pod2usage(1);

} 

if (!$opt_p) {
    print STDOUT "\nMANDATORY PARAMETER ERROR: -p <script_executor_user> parameter has not been specified for $patch_name.\n";
    pod2usage(1);
}

##########################
## 3) DB CONNECTION     ##
##########################

## Create the $schema object for the db_version object
## This should be replace for CXGN::DB::DBICFactory as soon as it can use CXGN::DB::InsertDBH

my $dbh =  CXGN::DB::InsertDBH->new(
                                     { 
                                         dbname => $opt_D, 
                                         dbhost => $opt_H 
                                     }
                                   )->get_actual_dbh();

print STDOUT "\nCreating the Metadata Schema object.\n";

my $metadata_schema = CXGN::Metadata::Schema->connect(   
                                                       sub { $dbh },
                                                      { on_connect_do => ['SET search_path TO metadata;'] },
                                                      );


##########################
## 4) DBVERSION DATA    ##
##########################

## Declaration of the name of the script and the description

my $patch_name = 'ChadoStock_v1_7_dbversion.pl';
my $patch_descr = 'Dbpatch to create chado stock module v1.7 in a postgreSQL database';

print STDOUT "\n+--------------------------------------------------------------------------------------------------+\n";
print STDOUT "Executing the patch:\n   $patch_name.\n\nDescription:\n  $patch_descr.\n\nExecuted by:\n  $opt_p.";
print STDOUT "\n+--------------------------------------------------------------------------------------------------+\n\n";

## And the requeriments if you want not use all
##
my @previous_requested_patches = (); ## ADD HERE

print STDOUT "\nChecking if this db_patch was executed before or if have been executed the previous db_patches.\n";

### Now it will check if you have runned this patch or the previous patches

my $dbversion = CXGN::Metadata::Dbversion->new($metadata_schema)
                                         ->complete_checking( { 
                                                                 patch_name  => $patch_name,
                                                                 patch_descr => $patch_descr, 
                                                                 prepatch_req => \@previous_requested_patches,
                                                                 force => $opt_F 
                                                              } 
                                                             );


##########################
## 5) METADBDATA        ##
##########################

### CREATE AN METADATA OBJECT and a new metadata_id in the database for this data

my $metadata = CXGN::Metadata::Metadbdata->new($metadata_schema, $opt_p);

### Get a new metadata_id (if you are using store function you only need to supply $metadbdata object)

my $metadata_id = $metadata->store()
                           ->get_metadata_id();


##########################
## 6) DB CHANGES        ##
##########################

### Now you can insert the data using different options:
##
##  1- By sql queryes using $dbh->do(<<EOSQL); and detailing in the tag the queries
##
##  2- Using objects with the store function
##
##  3- Using DBIx::Class first level objects
##

## In this case we will use the SQL tag

print STDOUT "\nExecuting the SQL commands.\n";

print STDOUT "\nExecuting the SQL commands.\n";

$dbh->do(<<EOSQL);

-- $Id: stock.sql,v 1.7 2007-03-23 15:18:03 scottcain Exp $
-- ==========================================
-- Chado stock module
--
-- DEPENDENCIES
-- ============
-- :import cvterm from cv
-- :import pub from pub
-- :import dbxref from general
-- :import organism from organism
-- :import genotype from genetic
-- :import contact from contact

-- ================================================
-- TABLE: stock
-- ================================================

create table stock (
       stock_id serial not null,
       primary key (stock_id),
       dbxref_id int,
       foreign key (dbxref_id) references dbxref (dbxref_id) on delete set null INITIALLY DEFERRED,
       organism_id int,
       foreign key (organism_id) references organism (organism_id) on delete cascade INITIALLY DEFERRED,
       name varchar(255),
       uniquename text not null,
       description text,
       type_id int not null,
       foreign key (type_id) references cvterm (cvterm_id) on delete cascade INITIALLY DEFERRED,
       is_obsolete boolean not null default 'false',
       constraint stock_c1 unique (organism_id,uniquename,type_id)
);
create index stock_name_ind1 on stock (name);
create index stock_idx1 on stock (dbxref_id);
create index stock_idx2 on stock (organism_id);
create index stock_idx3 on stock (type_id);
create index stock_idx4 on stock (uniquename);

COMMENT ON TABLE stock IS 'Any stock can be globally identified by the
combination of organism, uniquename and stock type. A stock is the physical entities, either living or preserved, held by collections. Stocks belong to a collection; they have IDs, type, organism, description and may have a genotype.';
COMMENT ON COLUMN stock.dbxref_id IS 'The dbxref_id is an optional primary stable identifier for this stock. Secondary indentifiers and external dbxrefs go in table: stock_dbxref.';
COMMENT ON COLUMN stock.organism_id IS 'The organism_id is the organism to which the stock belongs. This column should only be left blank if the organism cannot be determined.';
COMMENT ON COLUMN stock.type_id IS 'The type_id foreign key links to a controlled vocabulary of stock types. The would include living stock, genomic DNA, preserved specimen. Secondary cvterms for stocks would go in stock_cvterm.';
COMMENT ON COLUMN stock.description IS 'The description is the genetic description provided in the stock list.';
COMMENT ON COLUMN stock.name IS 'The name is a human-readable local name for a stock.';


-- ================================================
-- TABLE: stock_pub
-- ================================================

create table stock_pub (
       stock_pub_id serial not null,
       primary key (stock_pub_id),
       stock_id int not null,
       foreign key (stock_id) references stock (stock_id)  on delete cascade INITIALLY DEFERRED,
       pub_id int not null,
       foreign key (pub_id) references pub (pub_id) on delete cascade INITIALLY DEFERRED,
       constraint stock_pub_c1 unique (stock_id,pub_id)
);
create index stock_pub_idx1 on stock_pub (stock_id);
create index stock_pub_idx2 on stock_pub (pub_id);

COMMENT ON TABLE stock_pub IS 'Provenance. Linking table between stocks and, for example, a stocklist computer file.';


-- ================================================
-- TABLE: stockprop
-- ================================================

create table stockprop (
       stockprop_id serial not null,
       primary key (stockprop_id),
       stock_id int not null,
       foreign key (stock_id) references stock (stock_id) on delete cascade INITIALLY DEFERRED,
       type_id int not null,
       foreign key (type_id) references cvterm (cvterm_id) on delete cascade INITIALLY DEFERRED,
       value text null,
       rank int not null default 0,
       constraint stockprop_c1 unique (stock_id,type_id,rank)
);
create index stockprop_idx1 on stockprop (stock_id);
create index stockprop_idx2 on stockprop (type_id);

COMMENT ON TABLE stockprop IS 'A stock can have any number of
slot-value property tags attached to it. This is an alternative to
hardcoding a list of columns in the relational schema, and is
completely extensible. There is a unique constraint, stockprop_c1, for
the combination of stock_id, rank, and type_id. Multivalued property-value pairs must be differentiated by rank.';


-- ================================================
-- TABLE: stockprop_pub
-- ================================================

create table stockprop_pub (
     stockprop_pub_id serial not null,
     primary key (stockprop_pub_id),
     stockprop_id int not null,
     foreign key (stockprop_id) references stockprop (stockprop_id) on delete cascade INITIALLY DEFERRED,
     pub_id int not null,
     foreign key (pub_id) references pub (pub_id) on delete cascade INITIALLY DEFERRED,
     constraint stockprop_pub_c1 unique (stockprop_id,pub_id)
);
create index stockprop_pub_idx1 on stockprop_pub (stockprop_id);
create index stockprop_pub_idx2 on stockprop_pub (pub_id); 

COMMENT ON TABLE stockprop_pub IS 'Provenance. Any stockprop assignment can optionally be supported by a publication.';


-- ================================================
-- TABLE: stock_relationship
-- ================================================

create table stock_relationship (
       stock_relationship_id serial not null,
       primary key (stock_relationship_id),
       subject_id int not null,
       foreign key (subject_id) references stock (stock_id) on delete cascade INITIALLY DEFERRED,
       object_id int not null,
       foreign key (object_id) references stock (stock_id) on delete cascade INITIALLY DEFERRED,
       type_id int not null,
       foreign key (type_id) references cvterm (cvterm_id) on delete cascade INITIALLY DEFERRED,
       value text null,
       rank int not null default 0,
       constraint stock_relationship_c1 unique (subject_id,object_id,type_id,rank)
);
create index stock_relationship_idx1 on stock_relationship (subject_id);
create index stock_relationship_idx2 on stock_relationship (object_id);
create index stock_relationship_idx3 on stock_relationship (type_id);

COMMENT ON COLUMN stock_relationship.subject_id IS 'stock_relationship.subject_id is the subject of the subj-predicate-obj sentence. This is typically the substock.';
COMMENT ON COLUMN stock_relationship.object_id IS 'stock_relationship.object_id is the object of the subj-predicate-obj sentence. This is typically the container stock.';
COMMENT ON COLUMN stock_relationship.type_id IS 'stock_relationship.type_id is relationship type between subject and object. This is a cvterm, typically from the OBO relationship ontology, although other relationship types are allowed.';
COMMENT ON COLUMN stock_relationship.rank IS 'stock_relationship.rank is the ordering of subject stocks with respect to the object stock may be important where rank is used to order these; starts from zero.';
COMMENT ON COLUMN stock_relationship.value IS 'stock_relationship.value is for additional notes or comments.';



-- ================================================
-- TABLE: stock_relationship_cvterm
-- ================================================

CREATE TABLE stock_relationship_cvterm (
	stock_relationship_cvterm_id SERIAL NOT NULL,
	PRIMARY KEY (stock_relationship_cvterm_id),
	stock_relatiohship_id integer NOT NULL,
	--FOREIGN KEY (stock_relationship_id) references stock_relationship (stock_relationship_id) ON DELETE CASCADE INITIALLY DEFERRED,
	cvterm_id integer NOT NULL,
	FOREIGN KEY (cvterm_id) REFERENCES cvterm (cvterm_id) ON DELETE RESTRICT,
	pub_id integer,
	FOREIGN KEY (pub_id) REFERENCES pub (pub_id) ON DELETE RESTRICT
);
COMMENT ON TABLE stock_relationship_cvterm is 'For germplasm maintenance and pedigree data, stock_relationship. type_id will record cvterms such as "is a female parent of", "a parent for mutation", "is a group_id of", "is a source_id of", etc The cvterms for higher categories such as "generative", "derivative" or "maintenance" can be stored in table stock_relationship_cvterm';


-- ================================================
-- TABLE: stock_relationship_pub
-- ================================================

create table stock_relationship_pub (
      stock_relationship_pub_id serial not null,
      primary key (stock_relationship_pub_id),
      stock_relationship_id integer not null,
      foreign key (stock_relationship_id) references stock_relationship (stock_relationship_id) on delete cascade INITIALLY DEFERRED,
      pub_id int not null,
      foreign key (pub_id) references pub (pub_id) on delete cascade INITIALLY DEFERRED,
      constraint stock_relationship_pub_c1 unique (stock_relationship_id,pub_id)
);
create index stock_relationship_pub_idx1 on stock_relationship_pub (stock_relationship_id);
create index stock_relationship_pub_idx2 on stock_relationship_pub (pub_id);

COMMENT ON TABLE stock_relationship_pub IS 'Provenance. Attach optional evidence to a stock_relationship in the form of a publication.';


-- ================================================
-- TABLE: stock_dbxref
-- ================================================

create table stock_dbxref (
     stock_dbxref_id serial not null,
     primary key (stock_dbxref_id),
     stock_id int not null,
     foreign key (stock_id) references stock (stock_id) on delete cascade INITIALLY DEFERRED,
     dbxref_id int not null,
     foreign key (dbxref_id) references dbxref (dbxref_id) on delete cascade INITIALLY DEFERRED,
     is_current boolean not null default 'true',
     constraint stock_dbxref_c1 unique (stock_id,dbxref_id)
);
create index stock_dbxref_idx1 on stock_dbxref (stock_id);
create index stock_dbxref_idx2 on stock_dbxref (dbxref_id);

COMMENT ON TABLE stock_dbxref IS 'stock_dbxref links a stock to dbxrefs. This is for secondary identifiers; primary identifiers should use stock.dbxref_id.';
COMMENT ON COLUMN stock_dbxref.is_current IS 'The is_current boolean indicates whether the linked dbxref is the current -official- dbxref for the linked stock.';


-- ================================================
-- TABLE: stock_cvterm
-- ================================================

create table stock_cvterm (
     stock_cvterm_id serial not null,
     primary key (stock_cvterm_id),
     stock_id int not null,
     foreign key (stock_id) references stock (stock_id) on delete cascade INITIALLY DEFERRED,
     cvterm_id int not null,
     foreign key (cvterm_id) references cvterm (cvterm_id) on delete cascade INITIALLY DEFERRED,
     pub_id int not null,
     foreign key (pub_id) references pub (pub_id) on delete cascade INITIALLY DEFERRED,
     constraint stock_cvterm_c1 unique (stock_id,cvterm_id,pub_id)
);
create index stock_cvterm_idx1 on stock_cvterm (stock_id);
create index stock_cvterm_idx2 on stock_cvterm (cvterm_id);
create index stock_cvterm_idx3 on stock_cvterm (pub_id);

COMMENT ON TABLE stock_cvterm IS 'stock_cvterm links a stock to cvterms. This is for secondary cvterms; primary cvterms should use stock.type_id.';


-- ================================================
-- TABLE: stock_genotype
-- ================================================

create table stock_genotype (
       stock_genotype_id serial not null,
       primary key (stock_genotype_id),
       stock_id int not null,
       foreign key (stock_id) references stock (stock_id) on delete cascade,
       genotype_id int not null,
       foreign key (genotype_id) references genotype (genotype_id) on delete cascade,
       constraint stock_genotype_c1 unique (stock_id, genotype_id)
);
create index stock_genotype_idx1 on stock_genotype (stock_id);
create index stock_genotype_idx2 on stock_genotype (genotype_id);

COMMENT ON TABLE stock_genotype IS 'Simple table linking a stock to
a genotype. Features with genotypes can be linked to stocks thru feature_genotype -> genotype -> stock_genotype -> stock.';


-- ================================================
-- TABLE: stockcollection
-- ================================================

create table stockcollection (
	stockcollection_id serial not null, 
        primary key (stockcollection_id),
	type_id int not null,
        foreign key (type_id) references cvterm (cvterm_id) on delete cascade,
        contact_id int null,
        foreign key (contact_id) references contact (contact_id) on delete set null INITIALLY DEFERRED,
	name varchar(255),
	uniquename text not null,
	constraint stockcollection_c1 unique (uniquename,type_id)
);
create index stockcollection_name_ind1 on stockcollection (name);
create index stockcollection_idx1 on stockcollection (contact_id);
create index stockcollection_idx2 on stockcollection (type_id);
create index stockcollection_idx3 on stockcollection (uniquename);

COMMENT ON TABLE stockcollection IS 'The lab or stock center distributing the stocks in their collection.';
COMMENT ON COLUMN stockcollection.uniquename IS 'uniqename is the value of the collection cv.';
COMMENT ON COLUMN stockcollection.type_id IS 'type_id is the collection type cv.';
COMMENT ON COLUMN stockcollection.name IS 'name is the collection.';
COMMENT ON COLUMN stockcollection.contact_id IS 'contact_id links to the contact information for the collection.';


-- ================================================
-- TABLE: stockcollectionprop
-- ================================================

create table stockcollectionprop (
    stockcollectionprop_id serial not null,
    primary key (stockcollectionprop_id),
    stockcollection_id int not null,
    foreign key (stockcollection_id) references stockcollection (stockcollection_id) on delete cascade INITIALLY DEFERRED,
    type_id int not null,
    foreign key (type_id) references cvterm (cvterm_id),
    value text null,
    rank int not null default 0,
    constraint stockcollectionprop_c1 unique (stockcollection_id,type_id,rank)
);
create index stockcollectionprop_idx1 on stockcollectionprop (stockcollection_id);
create index stockcollectionprop_idx2 on stockcollectionprop (type_id);

COMMENT ON TABLE stockcollectionprop IS 'The table stockcollectionprop
contains the value of the stock collection such as website/email URLs;
the value of the stock collection order URLs.';
COMMENT ON COLUMN stockcollectionprop.type_id IS 'The cv for the type_id is "stockcollection property type".';


-- ================================================
-- TABLE: stockcollection_stock
-- ================================================

create table stockcollection_stock (
    stockcollection_stock_id serial not null,
    primary key (stockcollection_stock_id),
    stockcollection_id int not null,
    foreign key (stockcollection_id) references stockcollection (stockcollection_id) on delete cascade INITIALLY DEFERRED,
    stock_id int not null,
    foreign key (stock_id) references stock (stock_id) on delete cascade INITIALLY DEFERRED,
    constraint stockcollection_stock_c1 unique (stockcollection_id,stock_id)
);
create index stockcollection_stock_idx1 on stockcollection_stock (stockcollection_id);
create index stockcollection_stock_idx2 on stockcollection_stock (stock_id);

COMMENT ON TABLE stockcollection_stock IS 'stockcollection_stock links
a stock collection to the stocks which are contained in the collection.';

EOSQL

## Now it will add this new patch information to the md_version table.  It did the dbversion object before and
## set the patch_name and the patch_description, so it only need to store it.



##########################
## 7) STORE AND COMMIT  ##
##########################   

$dbversion->store($metadata);

$dbh->commit;

__END__

