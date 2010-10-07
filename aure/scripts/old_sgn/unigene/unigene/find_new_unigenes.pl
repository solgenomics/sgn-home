#!/usr/bin/perl
use Data::Dumper;
use strict;

my ($prebuild, $newbuild) = @ARGV;
my %oldestuni=();
my %newestuni=();
my %olduniest=();
my %oldcount=();
my %newcount=();

if(!$prebuild || !$newbuild) {
  die "usage: foo <prebuild> <newbuild>\n";
}

open PRE, "<$prebuild" or die "Cannot open input file $prebuild\n";
open NEW, "<$newbuild" or die "Cannot open input file $newbuild\n";

while(<PRE>) {
  chomp;
  my ($ug_id, $est_id, undef, undef)=split "\t";
  if($oldestuni{$est_id}) {
    if($oldestuni{$est_id} != $ug_id) {
      die "Clone $est_id is in unigene $oldestuni{$est_id} and $ug_id\n";
    }
    if($oldestuni{$est_id} eq $ug_id) {
      die "Death comes beacause $est_id makes $oldestuni{$est_id} eq $ug_id\n";
    }
  }

  my $est_found=0;
  foreach my $bar (@{$olduniest{$ug_id}}) {
    if ($bar == $est_id) {
      print STDERR "CRITICAL ERROR!!!!!\n";
      $est_found = 1;
    }
  }

  if(!$est_found) {
    push @{$olduniest{$ug_id}}, $est_id;
  }

  $oldestuni{$est_id} = $ug_id;
  if($oldcount{$ug_id}) {
    $oldcount{$ug_id} ++;
  } else {
    $oldcount{$ug_id} =1;
  }

#  print STDERR "$ug_id, $est_id\n";
}
close PRE;

print STDERR "Loaded old unigene -> est mapping\n";

while(<NEW>) {
  chomp;
  my ($ug_id, $est_id, undef, undef)=split "\t";
  if($newestuni{$est_id}) {
    if($newestuni{$est_id} != $ug_id) {
      die "Est $est_id is in unigene $newestuni{$est_id} and $ug_id\n";
    }
    if($newestuni{$est_id} eq $ug_id) {
      die "Death comes beacause $est_id makes $newestuni{$est_id} eq $ug_id\n";
    }
  }
  $newestuni{$est_id} = $ug_id;
  if($newcount{$ug_id}) {
    $newcount{$ug_id} ++;
  } else {
    $newcount{$ug_id} = 1;
  }
#  print STDERR "$ug_id, $est_id\n";
}
close NEW;

print STDERR "Loaded new unigene -> est mapping\n";


foreach my $olduni (sort keys(%olduniest)) {
#  print STDERR "olduni\n";
  my %newunigene = ();
  my $newuni;

  foreach my $candidate (@{$olduniest{$olduni}}) {
    if($newuni=$newestuni{$candidate}) {
      if($newunigene{$newuni}) {
	$newunigene{$newuni}++;
      } else {
	$newunigene{$newuni}=1;
      }

# && $newunigene != $newestuni{$candidate}) {
#	die "Old unigene $foo names ests used in more than one new unigene\n";
#      }
    } else {
      print STDERR "New unigene set abandons est $candidate used in old unigene $olduni\n";
    }
  }
  my $count=0;
  foreach my $newuni (sort keys %newunigene) {
    $count++;
    print "$olduni\t$newuni";
    print "\t$oldcount{$olduni}\t$newcount{$newuni}\t$newunigene{$newuni}\n";
  }
  if($count>1) { 
    print STDERR "Old unigene $olduni split into $count new unigenes\n"
  }

}
print STDERR "\n";
