
package YapRI::Base;

use strict;
use warnings;
use autodie;

use Carp qw| croak cluck |;
use Math::BigFloat;
use File::Temp qw/ tempfile tempdir /;
use File::Path qw/make_path remove_tree/;
use File::stat;


###############
### PERLDOC ###
###############

=head1 NAME

YapRI::Base.pm
A wrapper to interact with R/

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

 


=head1 DESCRIPTION

 Another yet perl wrapper to interact with R


=head1 AUTHOR

Aureliano Bombarely <ab782@cornell.edu>


=head1 CLASS METHODS

The following class methods are implemented:

=cut 



############################
### GENERAL CONSTRUCTORS ###
############################

=head2 constructor new

  Usage: my $rih = YapRI::Base->new($arguments_href);

  Desc: Create a new R interfase object.

  Ret: a YapRI::Base object

  Args: A hash reference with the following parameters:
        tempdir => $tempdir to store the files
        
        
  Side_Effects: Die if the argument used is not a hash or its values arent 
                right.

  Example: 

=cut

sub new {
    my $class = shift;
    my $args_href = shift;

    my $self = bless( {}, $class ); 

    my %permargs = (
	cmddir       => '\w+',
	cmdfiles     => {},
	r_opts_pass  => '-{1,2}\w+',
	use_defaults => '1|yes|',
	);

    ## Check variables.

    my %args = ();
    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ARGUMENT ERROR: Arg. supplied to new() isnt HASHREF");
	}
	else {
	    %args = %{$args_href}
	}
    }

    foreach my $arg (keys %args) {
	unless (exists $permargs{$arg}) {
	    croak("ARGUMENT ERROR: $arg isnt permited arg for new() function");
	}
	else {
	    unless (defined $args{$arg}) {
		croak("ARGUMENT ERROR: value for $arg isnt defined for new()");
	    }
	    else {
		if (ref($permargs{$arg})) {
		    unless (ref($permargs{$arg}) eq ref($args{$arg})) {
			croak("ARGUMENT ERROR: $args{$arg} isnt permited val.");
		    }
		}
		else {
		    if ($args{$arg} !~ m/$permargs{$arg}/) {
			croak("ARGUMENT ERROR: $args{$arg} isnt permited val.");
		    }
		}
	    }
	}
    }


    ## Set the dir to put all the commands

    my $cmddir = $args{cmddir} || '';  ## Empty var by default
    $self->set_cmddir($cmddir);
    if (exists $args{use_defaults}) {
	$self->set_default_cmddir();
    }

    my $cmdfiles_href = $args{cmdfiles} || {}; ## Empty hashref by default  
    $self->set_cmdfiles($cmdfiles_href);
    if (exists $args{use_defaults}) {
	$self->add_default_cmdfile();
    }

    my $resfiles_href = $args{resultfiles} || {}; ## Empty hashref by default  
    $self->set_resultfiles($resfiles_href);

    my $r_optspass = $args{r_opts_pass} || '';  ## Empty scalar by default 
    if (exists $args{use_defaults}) {
	$self->set_default_r_opts_pass();
    }
    else {
	$self->set_r_opts_pass($r_optspass);
    }

    return $self;
}

=head2 cleanup

  Usage: my $deleted_data_href = $rih->cleanup();

  Desc: Close all the filehandles and remove the cmddir with all the files

  Ret: A hash reference with key=datatype and value=datadeleted
       keys = cmddir and cmdfiles; 

  Args: None        
        
  Side_Effects: None

  Example: $rih->cleanup();

=cut

sub cleanup {
    my $self = shift;

    ## First close all the filehandles and remove the file hash from the object

    my %cmdfiles = %{$self->get_cmdfiles()};
    foreach my $file (keys %cmdfiles) {
	my $fh = $cmdfiles{$file};
	if (fileno($fh)) {
	    close($fh);
	}
    }
    
    $self->set_cmdfiles({});

    ## Remove the cmddir with all the files and the cmddir from the object

    my $cmddir = $self->delete_cmddir();

    return ( { cmddir => $cmddir, cmdfiles => \%cmdfiles} );
}



#################
### ACCESSORS ###
#################

=head2 get_cmddir

  Usage: my $cmddir = $rih->get_cmddir(); 

  Desc: Get the command dir used by the r interfase object

  Ret: $cmddir, a scalar

  Args: None

  Side_Effects: None

  Example: my $cmddir = $rih->get_cmddir();   

=cut

sub get_cmddir {
    my $self = shift;
    return $self->{cmddir};
}

=head2 set_cmddir

  Usage: $rih->set_cmddir($cmddir); 

  Desc: Set the command dir used by the r interfase object

  Ret: None

  Args: $cmddir, a scalar

  Side_Effects: Die if no argument is used.
                Die if the cmddir doesnt exists

  Example: $rih->set_cmddir($cmddir); 

=cut

sub set_cmddir {
    my $self = shift;
    my $cmddir = shift;
 
    unless (defined $cmddir) {
	croak("ERROR: cmddir argument used for set_cmddir function is undef.");
    }
    else {
	if ($cmddir =~ m/\w+/) {  ## If there are something check if exists
	    unless (defined(-f $cmddir)) {
		croak("ERROR: dir arg. used for set_cmddir() doesnt exists");
	    }
	}
    }
    
    $self->{cmddir} = $cmddir;
}

=head2 set_default_cmddir

  Usage: $rih->set_default_cmddir(); 

  Desc: Set the command dir used by the r interfase object with a default value
        such as RiPerldir_XXXXXXXX

  Ret: None

  Args: None

  Side_Effects: Create the perl_ri_XXXXXXXX folder in the tmp dir

  Example: $rih->set_default_cmddir(); 

=cut

sub set_default_cmddir {
    my $self = shift;

    my $cmddir = tempdir('RiPerldir_XXXXXXXX', TMPDIR => 1);

    $self->{cmddir} = $cmddir;
}


=head2 delete_cmddir

  Usage: my $cmddir = $rih->delete_cmddir(); 

  Desc: Delete the command dir used by the r interfase object

  Ret: $cmddir, deleted cmddir

  Args: None

  Side_Effects: Die if no argument is used.

  Example: $rih->delete_cmddir(); 

=cut

sub delete_cmddir {
    my $self = shift;

    my $cmddir = $self->get_cmddir();
    remove_tree($cmddir);
    return (delete($self->{cmddir}));
}

=head2 get_cmdfiles

  Usage: my $cmdfiles_href = $rih->get_cmdfiles(); 

  Desc: Get the command files used by the r interfase object

  Ret: $filesdir_href, a hash reference with key=filename and value=fh

  Args: $filename [optional]

  Side_Effects: None

  Example: my %cmdfiles = %{$rih->get_cmdfiles()};
           my $cmdfh = $rih->get_cmdfiles($filename);

=cut

sub get_cmdfiles {
    my $self = shift;
    my $filename = shift;

    if (defined $filename) {
	return $self->{cmdfiles}->{$filename};
    }
    else { 
	return $self->{cmdfiles};
    }
}

=head2 get_default_cmdfile

  Usage: my ($filename, $fh) = $rih->get_default_cmdfile(); 

  Desc: Get the command default file used by the r interfase object

  Ret: $filename, default filename
       $fh, default fh

  Args: None

  Side_Effects: None

  Example: my ($filename, $fh) = $rih->get_default_cmdfile(); 

=cut

sub get_default_cmdfile {
    my $self = shift;
    
    my %cmdfiles = %{$self->get_cmdfiles()};
    foreach my $file (keys %cmdfiles) {
	if ($file =~ m/RiPerlcmd/) {
	    return ($file, $cmdfiles{$file});
	}
    }
}


=head2 set_cmdfiles

  Usage: $rih->set_cmdfiles($cmdfiles_href); 

  Desc: Set the command files used by the r interfase object

  Ret: None

  Args: $cmdfile_hashref, an hash reference with key=filename and value=fh

  Side_Effects: Die if no argument is used.
                Die if the argument is not a hash reference
                Die if the value is not a filehandle

  Example: $rih->set_cmdfiles($cmdfiles_href); 

=cut

sub set_cmdfiles {
    my $self = shift;
    my $cmdfiles_href = shift ||
	croak("ERROR: No cmdfile argument was used for set_cmdfiles function");

    unless(ref($cmdfiles_href) eq 'HASH') {
	croak("ERROR: cmdfiles arg. used for set_cmdfiles isnt a hashref.");
    }
    else {
	foreach my $filename (keys %{$cmdfiles_href}) {

	    unless (-f $filename) {
		croak("ERROR: cmdfiles key=$filename doesnt exist");
	    }

	    my $fh = $cmdfiles_href->{$filename};
	    unless (ref($fh) eq 'GLOB') {
		croak("ERROR: cmdfiles value=$fh isnt a FILEHANDLE.");
	    }
	}
    }
    $self->{cmdfiles} = $cmdfiles_href;
}


=head2 add_cmdfile

  Usage: $rih->add_cmdfile($filename, $fh); 

  Desc: Add a new the command file for the r interfase object

  Ret: None

  Args: $filename, a scalar
        $fh, a filehandle

  Side_Effects: Die if no argument is used.
                Die if the filehandle is not a filehandle

  Example: $rih->add_cmdfile($filename, $fh); 

=cut

sub add_cmdfile {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename argument was used for add_cmdfile function");

    unless (-f $filename) {
	 croak("ERROR: cmdfile $filename doesnt exist");
    }

    my $fh = shift ||
	croak("ERROR: No filehandle arg. was used for add_cmdfile function");

    unless (ref($fh) eq 'GLOB') {
	croak("ERROR: cmdfiles value=$fh isnt a FILEHANDLE.");
    }
    
    $self->{cmdfiles}->{$filename} = $fh;
}

=head2 add_default_cmdfile

  Usage: $rih->add_default_cmdfile(); 

  Desc: Add a default the command file for the r interfase object

  Ret: None

  Args: None

  Side_Effects: Create a RiPerlcmd_XXXXXXXX file in the cmddir folder and
                it will open a new filehandle.
                If exists a default cmdfile in this object, it will not create 
                a new filename/fh pair

  Example: $rih->add_default_cmdfile();

=cut

sub add_default_cmdfile {
    my $self = shift;

    ## First, check there are any default file

    my ($filename, $fh) = $self->get_default_cmdfile();

    ## Second create a default file

    if (defined $filename) {
	my $cmddir = $self->get_cmddir();
	unless (defined $cmddir) {
	    croak("ERROR: Default cmdfile cant be created if cmddir isnt set");
	}
	($fh, $filename) = tempfile("RiPerlcmd_XXXXXXXX", DIR => $cmddir);
	$self->add_cmdfile($filename, $fh);
    }
}

=head2 delete_cmdfile

  Usage: $rih->delete_cmdfile($filename); 

  Desc: Delete a new the command file for the r interfase object, close the
        fh and delete the file.

  Ret: None

  Args: $filename, a scalar

  Side_Effects: Die if no argument is used.

  Example: $rih->delete_cmdfile($filename);

=cut

sub delete_cmdfile {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename argument was used for delete_cmdfile");

    ## 1) delete from the object

    my $fh = delete($self->{cmdfiles}->{$filename});

    ## 2) close the fh if it is opened

    if (fileno($fh)) {
	close($fh);
    }

    ## 3) delete the file from the system
    remove_tree($filename);    
}

=head2 get_resultfiles

  Usage: my $resultfiles_href = $rih->get_resultfiles(); 

  Desc: Get the result files used by the r interfase object

  Ret: $result_href, a hash reference with key=filename and value=resultfile

  Args: $filename [optional]

  Side_Effects: None

  Example: my %resfiles = %{$rih->get_resultfiles()};
           my $resultfile = $rih->get_resultfiles($filename);

=cut

sub get_resultfiles {
    my $self = shift;
    my $filename = shift;

    if (defined $filename) {
	return $self->{resultfiles}->{$filename};
    }
    else { 
	return $self->{resultfiles};
    }
}


=head2 set_resultfiles

  Usage: $rih->set_resultfiles($resfiles_href); 

  Desc: Set the result files used by the r interfase object

  Ret: None

  Args: $resfile_hashref, an hash reference with key=filename and 
        value=resultfile

  Side_Effects: Die if no argument is used.
                Die if the argument is not a hash reference
                Die if the value is a resultfile that doesnt exists

  Example: $rih->set_cmdfiles($cmdfiles_href); 

=cut

sub set_resultfiles {
    my $self = shift;
    my $resfiles_href = shift ||
	croak("ERROR: No resultfile arg. was used for set_resultfiles()");

    unless(ref($resfiles_href) eq 'HASH') {
	croak("ERROR: resultfiles used for set_resultfiles isnt a hashref.");
    }
    else {
	foreach my $filename (keys %{$resfiles_href}) {

	    unless (-f $filename) {
		croak("ERROR: cmdfiles key=$filename doesnt exist");
	    }

	    my $resultfile = $resfiles_href->{$filename};
	    unless (-f $resultfile) {
		croak("ERROR: resultfile $resultfile doesnt exist");
	    }
	}
    }
    $self->{resultfiles} = $resfiles_href;
}


=head2 add_resultfile

  Usage: $rih->add_resultfile($filename, $outfile); 

  Desc: Add a new resultfile associated to a cmdfile

  Ret: None

  Args: $filename, a scalar with a filepath
        $outfile, a scalar with a filepath

  Side_Effects: Die if no argument is used.
                Die if doesnt exists the outfile

  Example: $rih->add_cmdfile($filename, $fh); 

=cut

sub add_resultfile {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename arg. was used for add_resultfile function");

    unless (-f $filename) {
	 croak("ERROR: cmdfile $filename doesnt exist for add_resultfile");
    }

    my $outfile = shift ||
	croak("ERROR: No resultfile arg. was used for add_resultfile function");

    unless (-f $outfile) {
	 croak("ERROR: resultfile $outfile doesnt exist for add_resultfile");
    }
    
    $self->{resultfiles}->{$filename} = $outfile;
}


=head2 delete_resultfile

  Usage: $rih->result_cmdfile($filename); 

  Desc: Delete a result file associated with a cmdfile

  Ret: None

  Args: $filename, a scalar

  Side_Effects: Die if no argument is used.

  Example: $rih->delete_resultfile($filename);

=cut

sub delete_resultfile {
    my $self = shift;
    my $filename = shift ||
	croak("ERROR: No filename argument was used for delete_resultfile");

    ## 1) delete from the object

    my $outfile = delete($self->{resultfiles}->{$filename});

    ## 3) delete the file from the system
    remove_tree($outfile);    
}

=head2 get_r_opts_pass

  Usage: my $r_opts_pass = $rih->get_r_opts_pass(); 

  Desc: Get the r_opts_pass variable (options used with the R command)
        when run_command function is used

  Ret: $r_opts_pass, a string

  Args: None

  Side_Effects: None

  Example: my $r_opts_pass = $rih->get_r_opts_pass(); 
           if ($r_opts_pass !~ m/vanilla/) {
              $r_opts_pass .= ' --vanilla';
           }

=cut

sub get_r_opts_pass {
    my $self = shift;
    return $self->{r_opts_pass};
}

=head2 set_r_opts_pass

  Usage: $rih->set_r_opts_pass($r_opts_pass); 

  Desc: Set the r_opts_pass variable (options used with the R command)
        when run_command function is used. Use R -help for more info.
        The most common options used:
        --save                Do save workspace at the end of the session
        --no-save             Don't save it
        --no-environ          Don't read the site and user environment files
        --no-site-file        Don't read the site-wide Rprofile
        --no-init-file        Don't read the user R profile
        --restore             Do restore previously saved objects at startup
        --no-restore-data     Don't restore previously saved objects
        --no-restore-history  Don't restore the R history file
        --no-restore          Don't restore anything
        --vanilla             Combine --no-save, --no-restore, --no-site-file,
                              --no-init-file and --no-environ
        -q, --quiet           Don't print startup message
        --silent              Same as --quiet
        --slave               Make R run as quietly as possible
        --interactive         Force an interactive session
        --verbose             Print more information about progress

        The only opt that can not be set using set_r_opts_pass is --file, 
        it is defined by the commands stored as cmdfiles.

  Ret: None

  Args: $r_opts_pass, a string

  Side_Effects: Remove '--file=' from the r_opts_pass string

  Example: $rih->set_r_opts_pass('--verbose');

=cut

sub set_r_opts_pass {
    my $self = shift;
    my $r_opts_pass = shift;
 
    if ($r_opts_pass =~ m/(--file=.+)\s*/) {  ## If it exists, remove it
	carp("WARNING: --file opt. will be ignore for set_r_opts_pass()");
	$r_opts_pass =~ s/--file=.+\s*/ /g;
    }
    
    $self->{r_opts_pass} = $r_opts_pass;
}

=head2 set_default_r_opts_pass

  Usage: $rih->set_default_r_opts_pass(); 

  Desc: Set the default r_opts_pass for YapRI::Base (R --slave --vanilla)

  Ret: None

  Args: None

  Side_Effects: None

  Example: $rih->set_default_r_opts_pass(); 

=cut

sub set_default_r_opts_pass {
    my $self = shift;

    my $def_r_opts_pass = '--slave --vanilla';

    $self->{r_opts_pass} = $def_r_opts_pass;
}



##############################
## OUTPUT PARSING FUNCTIONS ##
##############################

=head parse_result 

  Usage: my %outparse = $rih->parse_result($resultfile); 

  Desc: Parse the output and return an hash with:
        key   => r_command line from cmdfile
        value => hashref with key=>line_number and value=value
                 or undef if the r_command did give an output

  Ret: A result hash (see description)

  Args: $resultfile, the resultfilename 

  Side_Effects: Die if resultfile doesnt exist or undef value is used

  Example: my %outparse = $rih->parse_result($resultfile); 
           my $result = $outparse{'y * x'};

=cut

sub parse_result {
    my $self = shift;
    my $resultfile = shift ||
	croak("ERROR: No result file was supplied to parse_result");

    unless (defined -f $resultfile) {
	croak("ERROR: $resultfile supplied to parse_result doesnt exist");
    }

    open my $rfh, '<', $resultfile;

    my %results = ();

    while(<$rfh>) {
	chomp($_);       
    }
}







#################
## CMD OPTIONS ##
#################

=head2 add_command

  Usage: $rih->add_command($r_command, $filename); 

  Desc: Add a R command line to a cmdfile.
        If a filename is specified the command line will be added to a 
        concrete $filename in the base object

  Ret: None

  Args: $r_command, a string with a R command
        $filename, a cmdfile to add the command

  Side_Effects: Add the command to the default if no filename is specified

  Example: $rih->add_command('x <- c(10, 9, 8, 5)')

=cut

sub add_command {
    my $self = shift;
    my $command = shift ||
	croak("ERROR: No command line was added to add_command function");
    my $cmdfile = shift;

    my ($fh, $filename);
    if (defined $cmdfile) {
	$fh = $self->get_cmdfiles($cmdfile);
    }
    else {
	($filename, $fh) = $self->get_default_cmdfile();
	unless (defined $filename) {
	    $self->add_default_cmdfile();
	    ($filename, $fh) = $self->get_default_cmdfile();
	}
    }

    unless(defined $fh) {
	croak("ERROR: No filehandle is associated with the filename used");
    }
    else {
	if (eof($fh) == 1) { ## If the filehandle is close, reopen
	    open $fh, '+>>', $filename;
	}
	print $fh "$command\n";
    }
}

=head2 get_commands

  Usage: my @commands = $rih->get_commands($filename); 

  Desc: Read the $filename and rerun

  Ret: None

  Args: $filename, a cmdfile to add the command

  Side_Effects: None

  Example: my @commands = $rih->get_commands($filename); 

=cut

sub get_commands {
    my $self = shift;
    my $cmdfile = shift;

    my @commands = ();
    my ($filename, $fh);
    unless (defined $cmdfile) {
	($filename, $fh) =  $self->get_default_cmdfile();
    }
    else {
	$fh = $self->get_cmdfiles($cmdfile);
    }

    unless(defined $fh) {
	croak("ERROR: No filehandle is associated with the filename used");
    }
    else {
	unless (eof($fh)) { ## If the filehandle is open
	    seek($fh, 0, 1);
	}
	else {  ## If the filehandle is closed
	    open $fh, '+<', $filename;
	}
	while(<$fh>) {
	    chomp($_);
	    push @commands, $_;
	}
    }
    return @commands;
}

=head2 run_command

  Usage: $rih->run_command($args_href); 

  Desc: Run as command line the R command file

  Ret: None

  Args: $args_href with the following keys/values pair

  Side_Effects: None

  Example: $rih->run_command($args_href); 

=cut

sub run_command {
    my $self = shift;
    my $args_href = shift;

    my %permargs = (
	cmdfile => '\w+',
	debug   => '1|0|yes|no'
	);
    
    my $base_cmd = 'R ';

    ## Add the running opts
    my $r_opts_pass = $self->get_r_opts_pass();

    $base_cmd .= $r_opts_pass;

    ## Check args

    my %args = ();
    if (defined $args_href) {
	unless (ref($args_href) eq 'HASH') {
	    croak("ERROR: Arg. used for run_commands isnt a HASHREF.");
	}
	else {
	    %args = %{$args_href};
	    foreach my $key (keys %args) {
		unless (exists $permargs{$key}) {
		    croak("ERROR: Key=$key isnt permited arg. for run_command");
		}
		else {
		    if ($args{$key} !~ m/$permargs{$key}/i) {
			my $err = "ERROR: Value=$args{$key} isnt a permited ";
			$err .= "value for key=$key at run_command function";
			croak($err);
		    }
		}
	    }	    
	}
    }   
    
    my $err = 'Aborting run_command.';

    ## Check cmddir

    my $cmddir = $self->get_cmddir();
    unless (defined $cmddir) {
	croak("ERROR: cmddir isnt set. Result files cannot be created. $err");
    }

    ## Get the cmdfile
    
    my $cmdfile = $args{cmdfile};
    unless (defined $cmdfile) {
	my ($filename, $fh) = $self->get_default_cmdfile();
	unless (defined $filename && $filename =~ m/\w+/) {
	    croak("ERROR: No default cmdfile was found. $err")
	}
	else {
	    $cmdfile = $filename;
	}
    }
    else {
	unless (-s $cmdfile) {
	    croak("ERROR: cmdfile=$cmdfile doesnt exist. $err");
	}
    }

    $base_cmd .= " --file=$cmdfile";

    ## Create a tempfile to store the results
    
    my (undef, $resultfile) = tempfile( "RiPerlresult_XXXXXXXX", 
					DIR => $cmddir,
					OPEN => 0,
	);
	
    $base_cmd .= " > $resultfile";

    ## finally it will run the command

    if (defined $args{debug} && $args{debug} =~ m/1|yes/i) {
	print STDERR "RUNNING COMMAND:\n$base_cmd\n";
    }

    my $run = system($base_cmd);
    
       
    if ($run == 0) {   ## It means success
	$self->add_resultfile($cmdfile, $resultfile);	
    }
    else {
	croak("\nSYSTEM FAILS running R:\n$run\n\n");
    }
}




####
1; #
####
