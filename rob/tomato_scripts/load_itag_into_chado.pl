#!/usr/bin/env perl
use strict;
use warnings;
use autodie ':all';

use File::Temp;

use IPC::Cmd 'can_run';
use IPC::Open2;

use Path::Class;

use Bio::Chado::Loader::FASTA;
use Bio::Chado::Loader::FixGMODBulkGFF3Polypeptides;

sub tee(@) {
    print join ' ', @_;
    print "\n";
    return @_;
}

$ENV{GMOD_ROOT} = make_fake_gmod_root();

my $gmod_bulkloader = 'gmod_bulk_load_gff3.pl';
die "$gmod_bulkloader must be in your PATH"
    unless can_run( $gmod_bulkloader );

my ( $dbhost, $dbname, $dbuser, $dbpass ) = @ARGV;

# load the assembly and gene models
load_gff3(
    glob('ITAG*_assembly.gff3'),
    glob('ITAG*_gene_models.gff3'),
   );


# fix the polypeptide features
Bio::Chado::Loader::FixGMODBulkGFF3Polypeptides->new(
    'db_dsn'  => "dbi:Pg:dbhost=$dbhost;dbname=$dbname",
    'db_user' => $dbuser,
    'db_pass' => $dbpass,
    'mrna_name_like' => 'Solyc%',
    );


# load the genomic sequences
my $fasta_l = Bio::Chado::Loader::FASTA->new(
    'organism_name' => 'Solanum _ Solanum lycopersicum',
    'type_name' => 'ultracontig',
    'analysis_name' => 'SL2.40_assembly',
    'source'  => 'SL2.40_assembly',
    'db_dsn'  => "dbi:Pg:dbhost=$dbhost;dbname=$dbname",
    'db_user' => $dbuser,
    'db_pass' => $dbpass,
    );
$fasta_l->run(
    glob('ITAG*_genomic.fasta'),
    );

    # glob(q|ITAG*_proteins.fasta|),
    # glob(q|ITAG*_genomic.fasta|),
    # glob(q|ITAG*_cdna.fasta|),
    # glob(q|ITAG*_cds.fasta|),



sub load_gff3 {
    my @files = @_;

    my $loader_stdin;
    my $loader_stdout;
    my $loader_pid = open2(
        #IO::String->new( \$loader_out ),
        \*STDOUT,
        $loader_stdin,
        $gmod_bulkloader,
        '--analysis',
        '--random_tmp_dir',
        '--dbxref',
        '--noexon',
        '--organism' => 'tomato',
        '--dbname'   => $dbname,
        '--dbuser'   => $dbuser,
        '--dbpass'   => $dbpass,
        '--save_tmpfiles',
        );

    for my $f ( @files ) {
        open my $fh, '<', $f;
        while( <$fh> ) {
            next unless /^SL2.40ch01/;
            $loader_stdin->print( $_ );
        }
    }
    close $loader_stdin;
    waitpid $loader_pid, 0;

    print $loader_stdout;
}

sub make_fake_gmod_root {
    my $tempdir = File::Temp->newdir;

    my $gmod    = file( $tempdir, 'conf', 'gmod.conf' );
    my $default = file( $tempdir, 'conf', 'default.conf' );
    $gmod->dir->mkpath;

    $gmod->openw->print(<<'');
CONF=/home/rob/dev/gmod/gmod_install_root/conf
TMP=/home/rob/dev/gmod/gmod_install_root/tmp

    $default->openw->print(<<'');
DBUSER=postgres
DBHOST=nightshade
DEFAULT=y
DBDRIVER=Pg
DBNAME=nonexistent_database!
DBPORT=5432
DBORGANISM=tomato
DBPASS=lalala
SQLFILE=/home/rob/dev/gmod/gmod_install_root/src/chado/modules/complete.sql

    return $tempdir;
}
