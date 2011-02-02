#!/usr/bin/perl

=head1 NAME

  base.t
  A piece of code to test the YapRI::Base module used 
  for PhylGomic pipeline

=cut

=head1 SYNOPSIS

 perl base.t
 prove base.t

=head1 DESCRIPTION

 Test YapRI::Base module used by PhylGomic pipeline.

=cut

=head1 AUTHORS

 Aureliano Bombarely Gomez
 (ab782@cornell.edu)

=cut

use strict;
use warnings;
use autodie;

use Data::Dumper;
use Test::More tests => 70;
use Test::Exception;
use Test::Warn;

use File::stat;
use Image::Size;
use Cwd;

use FindBin;
use lib "$FindBin::Bin/../lib";

## TEST 1

BEGIN {
    use_ok('YapRI::Base');
}

## Add the object created to an array to clean them at the end of the script

my @rih_objs = ();

## Create an empty object and test the possible die functions. TEST 2 to 6

my $rih0 = YapRI::Base->new();
push @rih_objs, $rih0;

is(ref($rih0), 'YapRI::Base', 
   "Test new function for an empty object; Checking object ref.")
    or diag("Looks like this has failed");

## By default it will create an empty temp dir


throws_ok { YapRI::Base->new(['fake']) } qr/ARGUMENT ERROR: Arg./, 
    'TESTING DIE ERROR when arg. supplied new() function is not hash ref.';

throws_ok { YapRI::Base->new({ fake => {} }) } qr/ARGUMENT ERROR: fake/, 
    'TESTING DIE ERROR for new() when arg. is not a permited arg.';

throws_ok { YapRI::Base->new({ cmddir => undef }) } qr/ARGUMENT ERROR: value/, 
    'TESTING DIE ERROR for new() when arg. has undef value';

throws_ok { YapRI::Base->new({ cmdfiles => []}) } qr/ARGUMENT ERROR: ARRAY/, 
    'TESTING DIE ERROR for new() when arg. doesnt have permited value';


###############
## ACCESSORS ##
###############

## Testing accessors for cmddir, TEST 7 to 11

my $currdir = getcwd;
my $testdir = $currdir . '/test';
mkdir($testdir);

$rih0->set_cmddir($testdir);

is($rih0->get_cmddir(), $testdir, 
    "testing get/set_cmddir, checking test dirname")
    or diag("Looks like this has failed");

throws_ok { $rih0->set_cmddir() } qr/ERROR: cmddir argument/, 
    'TESTING DIE ERROR when no arg. was supplied to set_cmddir()';

throws_ok { $rih0->set_cmddir('fake') } qr/ERROR: dir arg./, 
    'TESTING DIE ERROR when dir. arg. used doesnt exist in the system';

$rih0->delete_cmddir(); 

is($rih0->get_cmddir(), undef, 
    "testing delete_cmddir, checking if cmddir has been deleted from object")
    or diag("Looks like this has failed");

is(-f $testdir, undef,
    "testing delete_cmddir, checking if the dir has been deleted from system")
    or diag("Looks like this has failed");

## Once the cmddir has been deleted, it can be reset with set_default_cmddir.
## TEST 12 and 13

$rih0->set_default_cmddir();

is($rih0->get_cmddir() =~ m/RiPerldir_/, 1, 
   "testing set_default_cmddir, checking the default dirname")
   or diag("Looks like this has failed");

is(-d $rih0->get_cmddir(), 1, 
   "testing set_default_cmddir, checking that the default dir has been created")
   or diag("Looks like this has failed");

## Now it will create a file
## testing cmdfiles accessors, TEST 14 to 18


my $testfile0 = $rih0->get_cmddir() . '/testfile_for_ribase0.txt';
open my $testfh0, '>', $testfile0;

$rih0->set_cmdfiles({ $testfile0 => $testfh0 });
my $cmdfile0 = $rih0->get_cmdfiles($testfile0);

is($cmdfile0, $testfh0, 
    "testing get/set_cmdfiles, checking filehandle identity")
    or diag("Looks like this has failed");

throws_ok { $rih0->set_cmdfiles() } qr/ERROR: No cmdfile argument/, 
    'TESTING DIE ERROR when no arg. was supplied to set_cmdfiles()';

throws_ok { $rih0->set_cmdfiles('fake') } qr/ERROR: cmdfiles arg./, 
    'TESTING DIE ERROR when arg. used for set_cmdfiles isnt a HASHREF';

throws_ok { $rih0->set_cmdfiles({ test => 'fake'}) } qr/: cmdfiles key/, 
    'TESTING DIE ERROR when value arg. used for set_cmdfiles isnt a GLOB';

throws_ok { $rih0->set_cmdfiles({ $testfile0 => 'fake'}) } qr/: cmdfiles val/, 
    'TESTING DIE ERROR when value arg. used for set_cmdfiles isnt a GLOB';


## Testing add_cmdfile function, TEST 19 to 23

my $testfile1 = $rih0->get_cmddir() . '/testfile_for_ribase1.txt';
open my $testfh1, '>', $testfile1;

$rih0->add_cmdfile($testfile1, $testfh1);
my %cmdfiles0 = %{$rih0->get_cmdfiles()};

is(scalar(keys %cmdfiles0), 2,
    "testing add_cmdfile, checking number of files.")
    or diag("Looks like this has failed");

throws_ok { $rih0->add_cmdfile() } qr/ERROR: No filename argument/, 
    'TESTING DIE ERROR when no arg. was supplied to add_cmdfile()';

throws_ok { $rih0->add_cmdfile('fake') } qr/ERROR: cmdfile fake/, 
    'TESTING DIE ERROR when arg. file used for add_cmdfile doesnt exist';

throws_ok { $rih0->add_cmdfile($testfile1) } qr/ERROR: No filehandle arg./, 
    'TESTING DIE ERROR when no fh arg. was used for add_cmdfile';

throws_ok { $rih0->add_cmdfile($testfile1, 'fake') } qr/ERROR: cmdfiles value/, 
    'TESTING DIE ERROR when value arg. used for add_cmdfile isnt a GLOB';

## Testing delete_cmdfile function, TEST 24 to 27

is( -e $testfile1, 1,
    "testing delete_cmdfile, checking that the file exists previous deletion")
    or diag("Looks like this has failed");

$rih0->delete_cmdfile($testfile1);
my %cmdfiles1 = %{$rih0->get_cmdfiles()};

is(scalar(keys %cmdfiles1), 1,
    "testing delete_cmdfile, checking number of files.")
    or diag("Looks like this has failed");

is( -e $testfile1, undef,
    "testing delete_cmdfile, checking that the file has been deleted")
    or diag("Looks like this has failed");

throws_ok { $rih0->delete_cmdfile() } qr/ERROR: No filename argument/, 
    'TESTING DIE ERROR when no arg. was supplied to delete_cmdfile()';

## Test add/get_default_cmdfile functions, TEST 28 to 30

$rih0->delete_cmdfile($testfile0);

$rih0->add_default_cmdfile();
my ($deffile, $deffh) = $rih0->get_default_cmdfile();

is($deffile =~ m/RiPerlcmd_/, 1,
    "testing add/get_default_cmdfile, testing default filename")
    or diag("Looks like this has failed");

is(-e $deffile, 1, 
   "testing add/get_default_cmdfile, testing that a file has been created")
    or diag("Looks like this has failed");

## Test the cleanup, TEST 31 and 32

my $cmddir0 = $rih0->get_cmddir();
$rih0->cleanup();
is(-d $cmddir0, undef, 
    "Testing cleanup function, checking removing of the  def. cmddir")
    or diag("Looks like this has failed");

throws_ok { $rih0->add_default_cmdfile() } qr/ERROR: Default/, 
    'TESTING DIE ERROR when cmddir isnt set for add_default_cmdfile()';


### TESTING add_commands, TEST 32 to 40

my $rih1 = YapRI::Base->new({ use_defaults => 1 });
push @rih_objs, $rih1;

my @r_commands = (
    'x <- c(2)',
    'y <- c(3)',
    'x * y',
);

foreach my $r_cmd (@r_commands) {
    $rih1->add_command($r_cmd);
}

my ($def_cmdfile, $def_cmdfh)  = $rih1->get_default_cmdfile();
close($def_cmdfh);

open my $newfh, '<', $def_cmdfile;

my $l = 0;
while (<$newfh>) {
    chomp($_);
    is($_, $r_commands[$l], 
	"testing add_command, checking command lines in default file")
	or diag("Looks like this has failed");
    $l++;
}

throws_ok { $rih1->add_command() } qr/ERROR: No command/, 
    'TESTING DIE ERROR when no command is added add_command()';

throws_ok { $rih1->add_command('x <- c(9)', 'fake') } qr/ERROR: No filehandle/, 
    'TESTING DIE ERROR when filename added to add_command() doesnt exists';

my @g_commands = $rih1->get_commands();
my $n = 0;
foreach my $g_cmd (@g_commands) {
    is($g_cmd, $r_commands[$n], 
	"testing get_commands, checking command lines in default file")
	or diag("Looks like this has failed");
    $n++;
}

throws_ok { $rih1->get_commands('fake') } qr/ERROR: No filehandle/, 
    'TESTING DIE ERROR when filename added to get_commands() doesnt exists';

## Test if i can add more commands after read the file, TEST 41 to 44

my $new_r_cmd = 'y + x';
push @r_commands, $new_r_cmd;

$rih1->add_command($new_r_cmd);
my @ag_commands = $rih1->get_commands();
my $m = 0;

is(scalar(@ag_commands), 4,
   "testing add/get_commands after read the file, checking number of commands")
    or diag("Looks like this has failed");

foreach my $ag_cmd (@ag_commands) {
    is($ag_cmd, $r_commands[$m], 
	"testing add/get_commands after read the file, checking command lines")
	or diag("Looks like this has failed");
    $m++;
}

#####################################
## TEST Accessors for resultfiles  ##
#####################################

## 1) Create a test file

my $testfile2 = $rih1->get_cmddir() . '/testfile_for_ribase2.txt';
open my $testfh2, '+>', $testfile2;
close($testfh2);

## Get/Set resultfiles function. TEST 46 to 50

$rih1->set_resultfiles({ $def_cmdfile => $testfile2 });
my $get_resultfile = $rih1->get_resultfiles($def_cmdfile);

is($get_resultfile, $testfile2,
    "testing get/set_resultfile, checking filename")
    or diag("Looks like this has failed");

throws_ok { $rih1->set_resultfiles() } qr/ERROR: No resultfile arg./, 
    'TESTING DIE ERROR when no arg. was supplied to set_resultfiles()';

throws_ok { $rih1->set_resultfiles('fake') } qr/ERROR: resultfiles used/, 
    'TESTING DIE ERROR when arg. used for set_resultfiles isnt a HASHREF';

throws_ok { $rih1->set_resultfiles({ test => 'fake'}) } qr/: cmdfiles key/, 
    'TESTING DIE ERROR when key arg. used for set_resultfiles doesnt exist';

throws_ok { $rih1->set_resultfiles({ $def_cmdfile => 'fake'}) } qr/resultfile/, 
    'TESTING DIE ERROR when value arg. used for set_resultfiles doesnt exist';

## Delete resultfiles, TEST 51 to 52

$rih1->delete_resultfile($def_cmdfile);
my $get_resultfile2 = $rih1->get_resultfiles($def_cmdfile);

is($get_resultfile2, undef, 
   "testing delete_resultfile, checking resultfile in the object")
    or diag("Looks like this has failed");

is(-f $testfile2, undef, 
    "testing delete_resultfile, checking filename deletion")
    or diag("Looks like this has failed");

throws_ok { $rih1->delete_resultfile() } qr/ERROR: No filename arg./, 
    'TESTING DIE ERROR when no arg. was supplied to delete_resultfiles()';

## Test add resultfiles, TEST 54 to 58

my $testfile3 = $rih1->get_cmddir() . '/testfile_for_ribase3.txt';
open my $testfh3, '+>', $testfile3;
close($testfh3);

$rih1->add_resultfile($def_cmdfile, $testfile3);
my $get_resultfile3 = $rih1->get_resultfiles($def_cmdfile);

is($get_resultfile3, $testfile3,
    "testing add_resultfile, checking filename")
    or diag("Looks like this has failed");

throws_ok { $rih1->add_resultfile() } qr/ERROR: No filename arg./, 
    'TESTING DIE ERROR when no arg. was supplied to add_resultfile()';

throws_ok { $rih1->add_resultfile('fake') } qr/ERROR: cmdfile fake/, 
    'TESTING DIE ERROR when cmdfile used for add_resultfile doesnt exist';

throws_ok { $rih1->add_resultfile($def_cmdfile) } qr/ERROR: No resultfile/, 
    'TESTING DIE ERROR when no resultfile was used for add_resultfile';

throws_ok { $rih1->add_resultfile({ $def_cmdfile, 'fake'}) } qr/resultfile/, 
    'TESTING DIE ERROR when resultfile used for add_resultfile doesnt exist';


## Test accessors for r_opts_pass

$rih1->set_r_opts_pass('--verbose');
my $r_opts_pass = $rih1->get_r_opts_pass();

is($r_opts_pass, '--verbose', 
    "testing get/set_r_opts_pass, checking r_opts_pass variable")
    or diag("Looks like this has failed");

warning_like { $rih1->set_r_opts_pass('--slave --vanilla --file=test') } 
qr/WARNING: --file/i, 
    "TESTING WARNING when --file= is used for set_r_opts_pass";
    

##########################
## TEST RUNNING COMMAND ##
##########################

## Lets create a new object to test something more complex

my $rih2 = YapRI::Base->new({ use_defaults => 1 });
push @rih_objs, $rih2;

## Add the commands to enable a graph device and check that it exists

my $grfile1 = $rih2->get_cmddir() . "/TestMyGraph.bmp";
$rih2->add_command('bmp(filename="' . $grfile1 . '", width=600, height=800)');
$rih2->add_command('dev.list()');
$rih2->add_command('plot(c(1, 5, 10), type = "l")');
$rih2->add_command('dev.off()');

## Get the command file, and run it

my ($def_cmdfile2, $fh_cmdfh2) = $rih2->get_default_cmdfile();
$rih2->run_command({ debug => 0, cmdfile => $def_cmdfile2});

## Get the file

my $get_result_file2 = $rih2->get_resultfiles($def_cmdfile2);

## So, it will check different things, TEST 61
## 1) Does the output (result file) have the right data ?
##    It should contains: 
##    bmp            ## For bmp enable
##      2
##    null device    ## For bmp disable
##              1

my $filecontent_check = 0;
open my $check_fh1, '<', $get_result_file2;
while(<$check_fh1>) {
    if ($_ =~ m/bmp|null device|\s+1|\s+2/) {
	$filecontent_check++; 
    }
}

is($filecontent_check, 4, 
    "testing run_command, checking result file content")
    or diag("Looks like this has failed");

## Now it will check that the image file was created
## with the right size

## Put the image in the Image object, TEST 62 and 63

my ($img_x, $img_y) = Image::Size::imgsize($grfile1);

is($img_x, 600, 
    "testing run_command, checking image size (width)")
    or diag("Looks like this has failed");

is($img_y, 800, 
    "testing run_command, checking image size (heigth)")
    or diag("Looks like this has failed");

## Check die for run_command

throws_ok  { $rih2->run_command('fake') } qr/ERROR: Arg. used/, 
    'TESTING DIE ERROR when arg. used for run_command isnt a HASHREF';

throws_ok  { $rih2->run_command({ fake => 1}) } qr/ERROR: Key=fake/, 
    'TESTING DIE ERROR when key arg. used for run_command isnt valid';

throws_ok  { $rih2->run_command({ debug => 'please'}) } qr/ERROR: Value=pl/, 
    'TESTING DIE ERROR when value arg. used for run_command isnt valid';

throws_ok  { $rih2->run_command({ cmdfile => 'fake'}) } qr/ERROR: cmdfile=/, 
    'TESTING DIE ERROR when cmdfile used for run_command doesnt exist';

$rih2->cleanup();

throws_ok  { $rih2->run_command() } qr/ERROR: cmddir isnt set/, 
    'TESTING DIE ERROR when cmddir doesnt exist for run_command';

$rih2->set_default_cmddir(); 

throws_ok  { $rih2->run_command() } qr/ERROR: No default cmdfile/, 
    'TESTING DIE ERROR when no default cmdfile exists for run_command doesnt';

$rih2->add_default_cmdfile(); 

$rih2->add_command('bmp(filename="' . $grfile1 . '", width=600, height=800)');
$rih2->add_command('dev.list()');

##Add a non-specified file will make the command fail

$rih2->set_r_opts_pass('--file=');

throws_ok  { $rih2->run_command() } qr/SYSTEM FAILS running R/, 
    'TESTING DIE ERROR when system fail running run_command function';

$rih2->set_r_opts_pass('--slave --vanilla');



##############################################################
## Finally it will clean the files produced during the test ##
##############################################################

foreach my $clean_rih (@rih_objs) {
    $clean_rih->cleanup()
}
  
####
1; #
####
