#!/usr/bin/perl
# Load unigene build data.  Rewritten from Koni's unigene-db-load.pl.

# This program loads assembled unigene data into the relevant tables
# in SGN/CXGN.  Due to awful, awful legacy database structure, this is
# a complex process, crossing through the ``group'' table several
# times.

# The specific steps this program takes are as follows:

# (1) ensure that the organism group exists,
# (2) creates a new unigene build record in unigene_build,
# (3) create a new ``input group'' in groups and group_linkage containing
#     references to each EST id that went into the build.
# (4) for each unigene, insert a row into unigene with the unigene_build_id,
#     consensus_id, cluster number, contig number and number of members.
# (5) for each unigene member, insert a row into unigene_member with
#     the appropriate unigene_id.

use strict;
no strict "refs";
use Getopt::Std;

use CXGN::DB::InsertDBH;
use CXGN::BioTools::FastaParser;

#DBI->trace(4, '/tmp/rmk-trace.txt');
# These are a list of constants that were sprinkled throughout the old
# program as literals.  Unfortunately, Koni did not believe in naming
# things sensibly; you're just supposed to know that groups of ESTs
# meant to be used as inputs to a unigene assembly is a group of type
# 2.
use constant {
  unigene_build_method => 2,
  input_group_type => 2,
  input_group_member_type => 7,
  # These two defaults are workarounds for a flaw in the insert
  # rules on sgn.unigene.
  default_database_type => "SGN",
  default_sequence_name => 0
};


 our ($opt_a, $opt_e, $opt_g, $opt_h, $opt_i, $opt_o, $opt_u, $opt_U, $opt_Q, $opt_D, $opt_H);
  getopts("a:d:e:g:hi:o:u:U:Q:D:H:");


  my $dbh = CXGN::DB::InsertDBH::connect({dbname=>$opt_D ,dbhost=>$opt_H ,
#  my $dbh = CXGN::DB::InsertDBH::connect({dbname=>"sandbox",dbhost=>"scopolamine",
					  dbargs=>{ RaiseError => 1,
						    PrintError => 1 }})
    or die "Failed to connect to database ($DBI::errstr)";


# We do all of this in a block so that our lexical variables don't
# leak into subsequent subroutine definitions.
main: {
  my $queries = find_queries();

# unless (@ARGV) {
#   help ();
# }
  print STDERR <<EOF;

We need a database connection to validate command-line arguments.
Don't worry, if your arguments are bogus, we'll catch them.

EOF
#

print STDERR "Setting up crummy database api...\n";
  # Prepare all our queries and define functions to execute them.
  setup_database_api ($dbh, $queries);

  # Getopts, ensure files exist and are openable, etc.
print STDERR "Validating the world...\n";
  my ($organism_group_id, $est_infh, $unigene_infh, $membership_infh, $id_outfh, $unigene_build_id, $input_group_id, $dbname, $dbhost) =
    validate_the_world();
#  exit();


  


  # If we're here, we should be okay to run.  All the actual work
  # goes on right in this eval.
  eval {
    print STDERR "starting transaction\n";
    # Create a new unigene_build_id, unless a unigene_build_id is already
    # supplied
    unless ($unigene_build_id) {
      $unigene_build_id = load_new_build ($organism_group_id);
    }

    # Create a new input data group and fill group_linkage with crap,
    # unless an input group is supplied.
    unless ($input_group_id) {
      load_new_input_group ($organism_group_id, $unigene_build_id, $est_infh);
    }

    # Create new unigenes for each contig and singleton, and unigene
    # members for each EST in any unigene.
    load_acefile_membership_data ($unigene_build_id, $unigene_infh, $membership_infh, $id_outfh);

    commit_prompt ($dbh);
  };
  if ($@) {
    $dbh->rollback;
    die ($@);
  }
  
  #close files
  #
  close($membership_infh);
  close($unigene_infh);


#  exit (0)
}
# Done and done.

## Subroutines.
# 
sub load_new_build {
  my ($organism_group_id) = @_;
  unless ($organism_group_id) {
    die ("no organism_group_id supplied to new_build");
  }
  my ($unigene_build_id) = insert_new_build ($organism_group_id);
  print "Inserted a new unigene build (unigene_build_id: $unigene_build_id).\n";
  return ($unigene_build_id);
}

sub load_new_input_group {
  my ($organism_group_id, $unigene_build_id, $est_infh) = @_;
  unless (@_ > 0) {
    die ("no organism_group_id supplied to new_input_group");
  }
  unless ($organism_group_id) {
    die ("undefined organism_group_id supplied to new_input_group");
  }
  unless (@_ > 1) {
    die ("no unigene_build_id supplied to new_input_group");
  }
  unless ($unigene_build_id) {
    die ("undefined unigene_build_id supplied to new_input_group");
  }
  unless (@_ > 2) {
    die ("no est_infh supplied to new_input_group");
  }
  unless ($est_infh) {
    die ("undefined est_infh supplied to new_input_group");
  }
  my $input_group_id;
  {
    my $rows;
    ($input_group_id, $rows) = insert_new_input_group ($organism_group_id);
    unless ($rows) {
      die ("insert_new_input_group ($organism_group_id) inserted $rows rows, not 1.");
    }

    print "Inserted a new input group (group_id: $input_group_id).\n";
  }

  # Add the input group id to the unigene build row for this build.
  my $rows = update_input_group_for_unigene_build ($input_group_id, $unigene_build_id);
  unless ($rows) {
    die ("update_input_group_for_unigene_build($input_group_id, $unigene_build_id) updated zero rows");
  }
  print "Added input group to unigene build.\n";

  # Insert a row into group_linkage for each EST in the input.
  my $count = 0;
  foreach my $est_id ($est_infh->get_seqids()) {
    $est_id = extract_est_id ($est_id);
    my ($input_group_id) = insert_input_group_member ($input_group_id, $est_id);
    $count++;
  }
  print "Inserted $count ESTs.\n";
}

sub load_acefile_membership_data {
  my ($unigene_build_id, $unigene_infh, $membership_infh, $id_outfh) = @_;
  # Now process the acefile membership file.
  my %before;
  my %counted;
  my %after;
  my @counts = ("consensus", "unigene", "member");
  my $res = select_maxima();
  @before{@counts} = $res->fetchrow_array;
  @counted{@counts} = (0, 0, 0);
  print "Starting unigene load.\n";
  with_acefile_membership
    ($membership_infh,
     # What to do for each unigene.
     sub {
       my ($cluster, $contig, $nmemb) = @_;
       my ($consensus_id, $unigene_id, $seq, @qual);
       # If it's not a singleton, get the quality scores and
       # insert a record into unigene_consensi.
       if ($contig) {
	 my $seq_id = "$cluster-$contig";
	 $seq = $unigene_infh->get_sequence ($seq_id, \@qual);
	 unless ($seq) {
	   die ("Can't find sequence for cluster $cluster, contig $contig");
	 }
	 ($consensus_id) = insert_consensus ($seq, join(" ", @qual));
	 $counted{consensus}++;
       }
       # Singletons and contigs both get inserted into unigene.
       my ($unigene_id, $rows) = insert_unigene ($unigene_build_id, $consensus_id, $cluster, $contig, $nmemb);
       unless ($rows == 1) {
	 die ("inserted $rows rows in insert_unigene");
       }
       $counted{unigene}++;
       return (($unigene_id));
     },
     # What to do for each unigene member.
     sub {
       my ($est_id, $dir, $start, $end, $start_trim, $end_trim, $unigene_id) = @_;
       $est_id = extract_est_id ($est_id);
       my (undef, $rows)= 
	 insert_unigene_member ($unigene_id, $est_id, $start, $end,
				$start + $start_trim, $end - $end_trim,
				$dir);
       unless ($rows == 1) {
	 die ("inserted $rows rows in insert_unigene_member");
       }
       printf ($id_outfh "%s\t%s\n", $unigene_id, $est_id);
       $counted{member}++;
     }
    );
  $res = select_maxima();
  @after{@counts} = $res->fetchrow_array;
  foreach my $count (@counts) {
    my $delta = $after{$count} - $before{$count};
    unless ($delta == $counted{$count}) {
      die("Number of $count ids added to database $delta doesn't match number of insert queries $counted{$count}");
    }
  }
  close($id_outfh);
  print <<EOF;
Loaded $counted{consensus} proper unigenes.
Loaded $counted{unigene} unigenes (contigs and singletons).
Loaded $counted{member} unigene members.
EOF
}

# The intent here is to start factoring things like this out so that
# programs that have to agree about it can use operators like these to
# ensure the agreement.
## sub find_est_seq_file {
##   return (sprintf ("%s/basecall/qc-passed.seq", $_[0]));
## }
## sub find_unigene_seq_file {
##   return (sprintf ("%s/unigene/unigene.seq", $_[0]));
## }
sub extract_est_id {
  my ($est_id) = shift;
  if ($est_id =~ m/^\d+$/) {
    return ($est_id);
  }
  if ($est_id =~ m/^(SGN-E|E|sgn\|E)(\d+)$/) {
    return ($2);
  }
  die ("can't extract est_id from \"$est_id\"");
}

sub commit_prompt {
  my ($dbh, $prompt_message, $after_message) = @_;
  unless ($prompt_message) {
    $prompt_message = "Commit?\n(yes|no, default no)> ";
  }
  print $prompt_message;
  if (<STDIN> =~ m/^y(es)/i) {
    print "Committing...";
    $dbh->commit;
    print "okay.\n";
    if ($after_message) {
      print $after_message;
    }
  } else {
    print "Rolling back...";
    $dbh->rollback;
    print "done.\n";
  }
};

sub help {
  print STDERR <<EOF;
$0: load unigene sequence, membership, input EST set, etc.
Flags:
-H database hostname    for example localhost or db.sgn.cornell.edu
-D database name        sandbox or cxgn etc
-a PARSED_ACE_FILE	Parsed ACE membership file (mandatory)
-e EST_FILE		EST sequence FASTA file (mandatory)
-g INPUT_GROUP		Input group id for ESTs (optional, see below)
-h			What you see is what you get
-i FILENAME             Output file for mapping unigene ids to
  			cluster/contig numbers (mandatory)
-o ORGANISM_GROUP	Organism group id for species in build
-u UNIGENE_BUILD	Unigene build id (optional, see below)
-U UNIGENE_FILE         Unigene sequence FASTA file (mandatory).
-Q UNIGENE_QUALITY_FILE Unigene quality FASTA file (mandatory!!!!).

Options -g and -u cause this script to skip over the loading of the
input group and unigene_build, so that multiple parsed ACE membership
files can be loaded for one build (this was required once because of
flaws in the ACE file parsing).
EOF
exit (1);
}
sub validate_the_world {
 

  if ($opt_h) {
    help ();
  }

  unless ($opt_a) {
    die ("Required argument -a <acefile membership file> was not supplied.");
  }
  unless ($opt_e) {
    die ("Required argument -e <EST sequence file> was not supplied.");
  }
  unless ($opt_o) {
    die ("Required argument -o <organism group> was not supplied.");
  }
  unless ($opt_i) {
    die ("Required argument -i <unigene id map file> was not supplied.");
  }
  unless ($opt_U) {
    die ("Required argument -U <unigene sequence file> was not supplied.");
  }
  unless ($opt_Q) { 
    die ("Required argument -Q <unigene quality score file> was not supplied.");
  }

  # Make sure the membership file opens, parses.
  open(my $membership_infh, $opt_a)
    or die ("$opt_a does not name a file that exists.");
  # Do a pass over the acefile membership, to see that it looks okay.
  with_acefile_membership($membership_infh);
  seek ($membership_infh, 0, 0)
    or die ("failed to lseek() to position 0 in $opt_a: $!");

  # Make sure the organism group exists.
  my $sth = select_organism_name ($opt_o);
  unless ($sth->rows) {
    die ("$opt_o is not an organism group in group_linkage.");
  }

  # Open and parse the EST input file.
  my $est_file = $opt_e;
  my $est_infh = CXGN::BioTools::FastaParser->load_file($est_file)
    or die ("$est_file does not name an EST sequence file.");

  # Open and parse the assembly output file.
  my $unigene_file = $opt_U;
  my $unigene_qual_file = $opt_Q;
  my $unigene_infh = CXGN::BioTools::FastaParser->load_file($unigene_file, $unigene_qual_file)
      or die ("$est_file does not name an EST sequence file.");

  my $id_file = $opt_i;
  if (-f $id_file && -s $id_file) {
    die ("$id_file already exists and is not empty (I won't open it for writing)");
  }
  open (my $id_outfh, ">$id_file")
    or die ("failed to open $id_file: $!");

  # If the optional unigene build id is supplied, make sure there's
  # such a build id for this organism group.
  if ($opt_u) {
    my $sth = select_unigene_build ($opt_u, $opt_o);
    unless ($sth->rows > 0) {
      die ("there's no unigene build with id $opt_u for organism group $opt_o");
    }
  }

  # If the optional organism group id is supplied, check that there's
  # such a group.
  if ($opt_g) {
    my $sth = select_organism_group ($opt_g);
    unless ($sth->rows > 0) {
      die ("there's no group with id $opt_g (supposed to be an organism group)");
    }
  }

  if (!$opt_D || !$opt_H) { 
      die "need -D and -H options";
  }

  return ($opt_o, $est_infh, $unigene_infh, $membership_infh, $id_outfh, $opt_u, $opt_g, $opt_D, $opt_H);
}

# Given a database handle and a hash table of query names/SQL statements,
# do the following:
# (1) prepare the SQL statement,
# (2) make the query name name a function that when called executes the
#     prepared statement, and return:
#     * a statement handle if the statement is a SELECT,
#     * an insert id if the statement is an INSERT,
#     * the number of rows affected by an UPDATE or a DELETE.
sub setup_database_api {
  my ($dbh, $q) = @_;
  # Prepare all the q.
  map { my $x = $_;

	# Ensure that we can process what we care about in  the query
	# string /before/ the user tries to run it.
	unless (($q->{$x} =~ m/^select/i) || ($q->{$x} =~ m/^insert/i) ||
		($q->{$x} =~ m/^insert\s+into\s+((\w)|"([^\"]+)")/i) ||
		($q->{$x} =~ m/^(update|delete)/i)) {
	  die ("bug: the cheesy SQL processing in setup_database_api is too lame to handle $q->{$x}");
	}
	my $sth = $dbh->prepare ($q->{$_})
	  or die ($dbh->errstr);
	my $sub = sub {
	   print STDERR "executing $q->{$x}\nArgs: ".join (", ", @_)."\n";
	  # execute() returns rows affected for non-select statements
	  my $rows = $sth->execute(@_);
	  if ($q->{$x} =~ m/^select/i) {
	    return ($sth);
	  } elsif ($q->{$x} =~ m/^insert/i) {
	    if ($q->{$x} =~ m/^insert\s+into\s+((\w+)|"([^\"]+)")/i) {
	      my $ret = $dbh->last_insert_id("$1", "sgn");
	      return (($ret, $rows));
	    } else {
	      die ("bug: lame attempt to parse INSERT query failed.  sorry.");
	    }
	  } elsif ($q->{$x} =~ m/^(update|delete)/i) {
	    return ($rows);
	  }
	};
	*$x = $sub;
      }
    keys(%$q);
}

# This is essentially control structure, implemented as a higher-order
# function.  When called with just a filehandle argument, validates
# the acefile membership file.  When called with various subroutine
# arguments, calls the first argument once per unigene (contigs and
# singletons alike), the second argument once per unigene member
# (contigs and singletons alike, again).
sub with_acefile_membership {
  my ($infh, $unigene_sub, $member_sub) = @_;
  unless ($infh) {
    die ("with_acefile_membership called with no filehandle.");
  }
  unless ($unigene_sub) {
    $unigene_sub = sub { return ((1, 2, 3)); }
  }
  unless ($member_sub) {
    $member_sub = sub {
      unless (@_ >= 6) {
	die ("bad acefile record (too few fields) at $.");
      }
      if (grep { !defined($_); } @_) {
	die ("bad acefile record (undefined fields) at $.");
      }
    }
  }
  # The acefile membership is group/record/field structured: groups
  # are indicated by header lines, which contain the number of records
  # in the group.  Records are lines, and consist of
  # whitespace-separated fields.  There are two kinds of group:
  # contigs and singletons.
  while (<$infh>) {
    chomp;
    split;
    my ($cluster, $contig, $nmemb, $is_singleton);
    if (m/^Contig:/) {
      (undef, $cluster, $contig, $nmemb) = @_;
    } elsif (m/^Singlet:/) {
      (undef, $cluster) = @_;
      $nmemb = 1;
    } else {
      die ("bogus line in header position in acefile membership file\n$_");
    }
    # Do whatever needs to be done once per unigene.
    my @rest = $unigene_sub->($cluster, $contig, $nmemb);
    # Do whatever needs to be done once per unigene member.
    for (my $i = 0; $i < $nmemb; $i++) {
      $_ = <$infh>;
      chomp;
      split;
      $member_sub->(@_, @rest);
    }
  }
}

## SQL queries.  Each key/value pair in this hash will be turned
## into a subroutine in the current package by setup_database_api.
## The subroutine names will be the keys, and the subroutines will
## execute the SQL statements with the subroutine arguments in place
## of the placeholders.  SQL INSERT statements return an insert id,
## SQL SELECT statements return a statement handle, and UPDATE and 
## DELETE statements return the number of rows affected.
sub find_queries {
  my %q =
    (
     # Get the organism name for the organism group.
     select_organism_name =>
     "SELECT o.organism_name
    FROM organism o
    JOIN group_linkage l ON (l.member_id = o.organism_id)
   WHERE l.group_id = ?",

     # To check whether a putative group_id exists.
     select_organism_group => 
     "SELECT group_id FROM groups WHERE group_id = ?",

     # Insert a new unigene build.  The build_nr column probably ought
     # to be a serial number, but isn't.
     insert_new_build =>
     sprintf ("INSERT INTO unigene_build (organism_group_id, build_nr, method_id)
                  SELECT ?, (SELECT 1+MAX(build_nr) FROM unigene_build), %d",
	      unigene_build_method),

     # Insert a new group.  This is really disgusting: we store the unigene build id
     # and the organism group id as /strings/ in a comment in the groups table.
     # This really, really needs to be revised.
     insert_new_input_group =>
     sprintf ("INSERT INTO groups (type, comment)
                  SELECT %d,
                         'Input sequence data group for unigene build (unigene id '
                         ||MAX(build_nr)|| '), organism group ' || ?
                    FROM unigene_build",
	      unigene_build_method),

     # For checking whether a number is a unigene_build_id.
     select_unigene_build => 
     "SELECT unigene_build_id FROM unigene_build
       WHERE unigene_build_id = ? AND organism_group_id = ?",
     # Reccord the input data group id in the unigene build table.
     update_input_group_for_unigene_build =>
     "UPDATE unigene_build
     SET source_data_group_id = ?
   WHERE unigene_build_id = ?",

     # Add est ids to the group linkage table (This is how ests are associated
     # with unigene /builds/.  Unigene membership is represented slightly
     # more sanely.)
     insert_input_group_member => 
     sprintf ("INSERT INTO group_linkage (group_id ,member_id, member_type)
           VALUES (?, ?, %d)", input_group_member_type),

     # Store a sequence and quality score.  The id of this row is used in 
     # the insert_unigene table.
     insert_consensus =>
     "INSERT INTO unigene_consensi (seq, qscores) VALUES (?, ?)",

     # Insert a new unigene.  consensi_id comes from the insert_consensus query,
     # cluster_no, contig_no, and nr_members comes from the parsed acefile membership
     # table, and the database name and sequence name are silly columns used to
     # denote where the data came from (this was introduced when CGN data was added
     # to SGN's database).
     insert_unigene => 
     sprintf ("INSERT INTO unigene (unigene_build_id, consensi_id,
                                cluster_no, contig_no, nr_members,
                                database_name, sequence_name)
                       SELECT ?, ?, ?, ?, ?, '%s', %d",
	      default_database_type,
	      default_sequence_name),

     # Insert a unigene member.  The unigene id comes from the insert_unigene query,
     # est_id, start, stop, qstart, qend, dir come from non-header lines in
     # the parsed acefile membership.
     insert_unigene_member =>
     "INSERT INTO unigene_member (unigene_id, est_id, start, stop, qstart, qend, dir)
              VALUES (?, ?, ?, ?, ?, ?, ?)",
     select_maxima =>
     "SELECT (select last_value from sgn.unigene_consensi_consensi_id_seq),
           
             (SELECT  last_value from sgn.unigene_unigene_id_seq),
             (SELECT last_value FROM sgn.unigene_member_unigene_member_id_seq)"
    );
  return (\%q);
}
