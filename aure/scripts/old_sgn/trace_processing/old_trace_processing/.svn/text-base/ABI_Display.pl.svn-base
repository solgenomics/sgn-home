#!/usr/bin/perl -w
use strict;
use Class::Struct;
use CXGN::Chromatogram;
use CXGN::Chromatogram::Draw;
# Koni 2003 August 04
#
# This program originally written by Jay Gadgil
#
# - Modified to make code a bit more readable & robust
# - switched diagnositic prints to STDERR
#
# Virtually all comments are mine -- some included to explain what was not
# obvious to me at first glance, other to indicate what I changed for
# correctness or readibility.
#
# Much code was superfluous -- incorrect statements with subsequent statements
# that corrected the errors. I removed the incorrect statements and made any
# necessary adjustments. I also removed commented out statements that appeared
# to be incorrect and left over from debugging.

# This enclosure with main() is not necessary. I'm not sure if it will make
# things easier in mod_perl since normally this design would force all 
# variables to be passed to subroutines -- however, I found a few instances
# where variables were not being passed, but were being used, and strict is on
# and this script has been used on CGN for some time, before I touched it,
# without any difficulties.
&main (\@ARGV);

sub main () {
    my @ARG = @{$_[0]};

    my $ABIStart_Position=0; # Maybe at 0 or 128
    my $buffer;

    my ($numbyt, $numwrd);
    my (@Bases, @BasePos);
    my ($tp, $nb);

    my ($ABIFile, $OutFile, $Verbose, $TempDir, $PHDfile, $Segment_Height, 
	$Segment_Width, $phred_path, $hilight, $use_phred) = "";

    # By default, use the ABI base callings. If a PHD file is specified or
    # a path to phred, then phred basecallings are turned on automatically
    # -P forces phred basecallings.
    $use_phred = 0;

    for(my $i=0;$i<@ARG;$i++) {
      if ($ARG[$i] eq "-a") {
	$ABIFile = $ARG[++$i];
      } elsif ($ARG[$i] eq "-o") {
	$OutFile = $ARG[++$i];
      } elsif ($ARG[$i] eq "-v") {
	$Verbose = 1;
      } elsif ($ARG[$i] eq "-t") {
	$TempDir = $ARG[++$i];
	$TempDir .= "/" if ($TempDir !~ m/\/$/);
      } elsif ($ARG[$i] eq "-p") {
	$PHDfile = $ARG[++$i];
	$use_phred = 1;
      } elsif ($ARG[$i] eq "-h") {
	$Segment_Height = $ARG[++$i];
      } elsif ($ARG[$i] eq "-w") {
	$Segment_Width = $ARG[++$i];
      } elsif ($ARG[$i] eq "-phred_path") {
	$phred_path = $ARG[++$i];
	$use_phred = 1;
      } elsif ($ARG[$i] eq "-hilight") {
	$hilight = $ARG[++$i];
      } elsif ($ARG[$i] eq "-P") {
	$use_phred = 1;
      } else {
	# Added this
	print STDERR "$0: Unknown option $ARG[$i]\n";
      }
    }

#set defaults

    $OutFile ||= "/tmp/traceview";
    $TempDir ||= "/tmp/";
    $phred_path ||= "/usr/local/bin/phred";
    $Segment_Height ||= 120;
    $Segment_Width ||= 800;

    if ($ABIFile eq "") {
      print STDERR "You must specify an ABI chromatogram file with -a\n";
      exit -1;
    }

    if ( ! -f "$ABIFile") {
      print STDERR "ABI chromatogram file \"$ABIFile\" does not exist\n";
      exit -1;
    }

    CXGN::Chromatogram::is_abi_file($ABIFile) or die ("$ABIFile is not an ABI file.\n");
    open ABIFile, $ABIFile or die "Can't open file '$ABIFile' ($!)";
    binmode ABIFile, ":raw";

    # Added backslashes to this declaration because the combination of dollar
    # sign and single quote was messing up xemacs CPerl mode (how typical)
    struct( DataIndex => {
			  occur => "\$",
			  sernum => "\$",
			  offset => "\$",
			  sizwrd => "\$",
			  numbyt => "\$",
			  numwrd => "\$",
			  label => "\$"
			 });


# Data Indices -- these are effectively used as constants to index the
#                 vectors of data stored in DataIndex structs (defined above)
    my $TRACE1 = 0;
    my $TRACE2 = 1;
    my $TRACE3 = 2;
    my $TRACE4 = 3;
    my $BASMAP = 4;
    my $BASES =  5;
    my $BASPOS = 6;
    my $SIGSTR = 7;
    my $AVGSPC = 8;
    my $PRIPOS = 9;
    my $MCHNAM = 10;
    my $DYEPRI = 11;
    my $SMPNAM = 12;
    my $THMPRT = 13;
    my $LANENM = 14;
    my $GELNAM = 15;
    my $COMMNT = 16;

    my @labels = ("DATA", "DATA", "DATA", "DATA", 
		  "FWO_", "PBAS", "PLOC", "S/N%",
		  "SPAC", "PPOS", "MCHN", "PDMF",
		  "SMPL", "THUM", "LANE", "GELN",
		  "CMNT");

    # ?
    my @blockSerNum = (9,10,11,12,
		       1, 1, 1, 1,
		       1, 1, 1, 1,
		       1, 1, 1, 1,
		       1);
    my @dataIndex;

    my ($block_size, $n_blocks, $start_offset);

    print STDERR "---------------------\n" if ($Verbose);

    # Load meta information which will dictate the interpretation of the rest
    # of the file
    seek ABIFile, 16, 0;
    read ABIFile, $block_size, 2;
    seek ABIFile, 18, 0;
    read ABIFile, $n_blocks, 4;
    seek ABIFile, 26, 0;
    read ABIFile, $start_offset, 4;

    # Intergers stored in the ABI file are in big-endian format. These unpacks
    # ensure the integer is extracted from the bytes read, into the machines
    # native format. On Intel processors, this reverses the byte order.
    $block_size = unpack("n", $block_size);
    $n_blocks = unpack("N", $n_blocks);
    $start_offset = unpack("N", $start_offset);

    if ($Verbose) {
      print STDERR "Block size = $block_size\n";
      print STDERR "# of Blocks: $n_blocks\n";
      print STDERR "First block begins at block $start_offset\n";
    }

    # There seemed to be an error here, $i<16 was used, but there are 
    # 17 elements in @labels declared above
    for (my $i=0;$i<@labels;$i++) {
      $dataIndex[$i] = new DataIndex;
      $dataIndex[$i]->occur(-1);
      $dataIndex[$i]->offset(0);
      $dataIndex[$i]->numbyt(0);
      $dataIndex[$i]->numwrd(0);
      $dataIndex[$i]->label($labels[$i]);
      $dataIndex[$i]->sernum($blockSerNum[$i]);
    }


    my $curOffset = $start_offset;

    # Its not clear why we bother loading anything but trace data blocks
    # The others, with the possible exception the the BASES block, are not
    # used.
    for (my $i=0; $i<$n_blocks; $i++) {
	&readDataIndex($curOffset, \@dataIndex, $Verbose);
	$curOffset = $curOffset + $block_size;
    }


    if ($Verbose) {
      print STDERR "Number of Trace points => " .
	$dataIndex[$TRACE1]->numwrd . "\n";
      print STDERR "Number of bases => " .
	$dataIndex[$BASES]->numwrd . "\n";
    }

    # OK - so this is strange. I am guessing this indicates how to map
    # the index numbers of the traces to A, C, G, or T. How will that work
    # if phred is used? For now, I am leaving this as it is -- it doesn't
    # seem to be used anyway.
    if ($dataIndex[$BASMAP]->occur != -1) {
      if($Verbose){
	print STDERR "Number of bytes in Basemap: " . 
	  $dataIndex[$BASMAP]->numbyt . "\n";
      }
      if ($dataIndex[$BASMAP]->numbyt <= 4) {
	my ($bM, @baseMap);
	$bM = $dataIndex[$BASMAP]->offset;
	$baseMap[0] = chr($bM >> 24 & 0xff);
	$baseMap[1] = chr($bM >> 16 & 0xff);
	$baseMap[2] = chr($bM >> 8 & 0xff);
	$baseMap[3] = chr($bM & 0xff);
	
	if ($Verbose) {
	  print STDERR "baseMap: " .  $dataIndex[$BASMAP]->offset . "\n";
	  print STDERR "baseMap[0]: " . $baseMap[0] . "\n";
	  print STDERR "baseMap[1]: " . $baseMap[1] . "\n";
	  print STDERR "baseMap[2]: " . $baseMap[2] . "\n";
	  print STDERR "baseMap[3]: " . $baseMap[3] . "\n";
	}
      } else {
	# Apparently (makes sense if what I think this is (see above) is
	# true, this is never expected to be more than 4 bytes. I switched
	# this to a die()
	die("Basemap is gt 4 bytes\n");
      }
    }

    if ($use_phred) {
      if (!$PHDfile) {
	$PHDfile = runPhred($ABIFile, $TempDir, $phred_path);
      }
      PhredBaseCalls($PHDfile, \@Bases, \@BasePos);
    } else {
      # Use the ABI base calls
      if ($dataIndex[$BASES]->occur != -1) {
	if ($Verbose) {
	  print STDERR "----Bases-------\n";
	  print STDERR "Number of bytes in Bases: " . 
	    $dataIndex[$BASES]->numbyt . "\n";
	}	
	if ($dataIndex[$BASES]->numbyt <= 4) {
	  $buffer = $dataIndex[$BASES]->offset;
	} else {
	  $buffer = &readStr($dataIndex[$BASES]->offset, 
			     $dataIndex[$BASES]->numbyt);
	}
	$numbyt = $dataIndex[$BASES]->numbyt;
	@Bases = unpack("C$numbyt", $buffer);	
      }

      if ($dataIndex[$BASPOS]->occur != -1) {
	$numbyt = $dataIndex[$BASPOS]->numbyt;
	$numwrd = $dataIndex[$BASPOS]->numwrd;
	if ($Verbose) {
	  print "----Base Positions----\n";
	  print "Number of bytes in Base Positions: " . $numwrd . "\n";
	}
	if ($numbyt <= 4) {
	    $buffer = $dataIndex[$BASPOS]->offset;
	} else {
	    $buffer = &readStr($dataIndex[$BASPOS]->offset, $numbyt);
	}
	@BasePos = unpack("n*", $buffer);
      }
    }

    # Switched this around so that references are returned directly by the
    # readTrace function. Previously 4 seperate full copy arrays were being
    # returned.
    my @Traces;
    $Traces[0] = &readTrace($dataIndex[$TRACE1], $Verbose);
    $Traces[1] = &readTrace($dataIndex[$TRACE2], $Verbose);
    $Traces[2] = &readTrace($dataIndex[$TRACE3], $Verbose);
    $Traces[3] = &readTrace($dataIndex[$TRACE4], $Verbose);


    $tp = $dataIndex[$TRACE1]->numwrd;

    # This was an error -- it was using $dataIndex[$BASES]->numwrd, but this
    # is wrong if phred base calling is used.
    $nb = @Bases;

    &CXGN::Chromatogram::Draw::DrawTrace($tp,$nb, \@Traces, \@BasePos, \@Bases, $OutFile, 
			   $Segment_Height, $Segment_Width, $Verbose);

    close ABIFile;
    print STDERR "$0: Finished\n" if($Verbose);
}

########################################

sub readStr() {
    my ($offset, $numbytes) = @_;
    my $tmp_read;
    seek ABIFile, $offset, 0;
    read ABIFile, $tmp_read, $numbytes;
    return $tmp_read;
}

# readBE = Read Big Endian
sub readBE() {
    my ($offset, $numbytes) = @_;
    my $tmp_read;
    my $output;
    seek ABIFile, $offset, 0;
    read ABIFile, $tmp_read, $numbytes;

    # Originally, this was coded with the 4 byte case as the else statement
    # allowing it to be potentially mistakenly used when neither 2 nor 4 bytes
    # were specified, but some other number. This is clearly not intended
    # to happen -- perhaps the else case should call die() rather than this
    # fail-safe setup here.
    if ($numbytes == 2) {
	$output = unpack("n", $tmp_read);
    } elsif ($numbytes == 4) {
	$output = unpack("N", $tmp_read);
    } else {
      $output = $tmp_read;
    }

    return $output;
}


sub readDataIndex() {
  my ($offset, $dataindex_ref, $Verbose) = @_;

  my ($Label, $SerNum, $DatType, $SizWrd, $NumBytes, $NumWrd, $Data);

    # This results in wasted calls to seek. Nothing will happen in that seek
    # but it wastes the interpreters' time figuring out that nothing need be
    # done.
    $Label = &readStr($offset, 4);
    $SerNum = &readBE($offset+4, 4);
    $DatType = &readBE($offset+8, 2);
    $SizWrd = &readBE($offset+10, 2);
    $NumWrd = &readBE($offset+12, 4);
    $NumBytes = &readBE($offset+16,4);
    $Data = &readBE($offset+20,4);

    my $found = 0;
    my $foundindex=0;

  # This is a little weird... 17 times we read dataIndex blocks from the file
  # and for each of those 17 times we'll search through an array of 17 to find
  # the right spot to store the one under our nose...
  for (my $i=0;$i<@{$dataindex_ref};$i++) {
    if (($dataindex_ref->[$i]->label eq $Label) 
	and ($dataindex_ref->[$i]->sernum == $SerNum)) {
      $found = 1;
      $foundindex = $i;
      last;
    }
  }

  if ($found) {
    if ($Verbose) {
      print STDERR <<EOF;

Offset:   $offset - FoundIndex: $foundindex
Label:    $Label
SerNum:   $SerNum
DatType:  $DatType
SizWrd:   $SizWrd
NumWrd:   $NumWrd
NumBytes: $NumBytes
Data:     $Data
----------------------------------------------

EOF
    }

    $dataindex_ref->[$foundindex]->occur(1);
    $dataindex_ref->[$foundindex]->sizwrd($SizWrd);
    $dataindex_ref->[$foundindex]->numwrd($NumWrd);
    $dataindex_ref->[$foundindex]->numbyt($NumBytes);
    $dataindex_ref->[$foundindex]->offset($Data);
  } else {
    # Silently ignore unknown fields, unless verbosity is requested.
    print STDERR "Unknown dataIndex block \"$Label\" serial number $SerNum\n"
      if $Verbose;
  }
}



sub readTrace() {
  my ($dataindex_obj, $Verbose) = @_;

  my $offset = $dataindex_obj->offset;
  my $numbyt = $dataindex_obj->numbyt;
  my $numwrd = $dataindex_obj->numwrd;
  my $occur = $dataindex_obj->occur;

  my $buffer;
  my @Trace;

  if ($occur != -1) {
    if ($numbyt <= 4) {
      # This occurs all over the place in this code -- if 4 bytes or less
      # are needed to encode the data referenced, the reference itself (the
      # offset is used to store the data. Its not clear however how the padding
      # works. If there is less than 4 trace point observations however, we
      # might as well just piss into the wind....
      $buffer = $offset;
    } else {
      seek(ABIFile, $offset, 0) 
	or die "Can't seek in file ABIFILE: $!";
      my $j = read (ABIFile, $buffer, $numbyt) 
	or die "Can't read from file ABIFILE";

      print "$j bytes read   $! \t $^E \n" if $Verbose;
    }

    @Trace = unpack("n$numwrd", $buffer);
  }

  # This was returning the entire array which is a deep copy operation. I have
  # switched this to return a reference instead since the calling code simply
  # takes a reference to the returned copy and stuffs that reference into
  # another array.
  return \@Trace;
}



###
#
# Run Phred if PHD file not given
#
###

sub runPhred() {
  my($ABIFile, $tmp_dir, $phredpath) = @_;

  if ( ! -d "$tmp_dir" ) {
    die "Temporary directory $tmp_dir does not exist";
  }

  if ($phredpath eq "") {
    $phredpath = "phred";
  } else {
    if ($phredpath !~ m/phred$/) {
      $phredpath .= "phred";
    }
  }

  # Fixed this to not set this variable if it is already set.
  if (!defined($ENV{PHRED_PARAMETER_FILE})) {
    $ENV{PHRED_PARAMETER_FILE} = "/usr/local/etc/phredpar.dat";
  }

  # Since we will unlink immediately, it is not neceessary for this to
  # have a more random name -- no two processes on the system have the same
  # process id. If the directory was shared between different systems, then
  # this would not work.
  my $FileList = $tmp_dir . "FileList-$$";

  my (undef, $ABI_basename) = $ABIFile =~ m/(\S+\/|^)([^\/]+)$/;
  my $PHD_File = $tmp_dir . "${ABI_basename}.phd.1";

  my $phredoptions = "-if $FileList -sd $tmp_dir -qd $tmp_dir -pd $tmp_dir -process_nomatch";

  open (FILELIST, ">$FileList") 
    or die "Couldn't open temporary file \"$FileList\" ($!)";
  print FILELIST "$ABIFile\n";
  close FILELIST;

  # Note that normal output (stdout) from phred is supressed while anything
  # it might write to stderr should come through. This is desireable as we
  # care not about what goes as expected.
  system("$phredpath $phredoptions > /dev/null")
    and die "Calling phred failed ($! - $?) using $phredpath $phredoptions";

  unlink("$FileList");
  return $PHD_File;
}



###
#
# Get Phred Base Calls
#
###

sub PhredBaseCalls() {
  my ($PHD_File, $base_ref, $basepos_ref) = @_;

  open PHD_FILE, "$PHD_File" 
    or die "Can't open PHD File \"$PHD_File\" ($!)\n";

  my $count=0;

  while(<PHD_FILE>) {
    chomp;
    last if $_ eq "BEGIN_DNA";
  }
  while(<PHD_FILE>) {
    chomp;
    last if $_ eq "END_DNA";
    my ($base, $quality, $position) = split;
    $base_ref->[$count] = uc $base;
    $basepos_ref->[$count] = $position;
    $count++;
  }
  close PHD_FILE;

}
