#!/usr/bin/perl -w
#
# $Id: joinby.pl,v 1.8 2001/12/29 02:06:29 gcavet Exp $
#
# joins text files based on common values in key columns

use strict;

my( $option,
    $indelimiter,
    $outdelimiter,
    @fields,
    %hash,
    $key,
    );

$indelimiter = "\t";
$outdelimiter = "\t";

if ($#ARGV < 0) { helpfulexit(); }

# LOAD UP DATA FROM INPUT FILES
open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";


$key = "";
while (<IN>) {
    
    # BREAK INPUT LINE INTO FIELDS
    chomp;
    @fields = split(/$indelimiter/, $_);
    # CHECK FOR EMPTY TRAILING FIELD
    while (substr($_, -(length($indelimiter))) eq $indelimiter) {
      push(@fields, '');
      $_ = substr($_, 0, (length($_) - 1));
    }

    # PULL OUT KEY COLUMN
    $key = join( "|",$fields[0],$fields[1]);
    if (defined($hash{$key})) {
	next;
    }
    $hash{$key} = 1;
   
    print join($outdelimiter, @fields), "\n";
  
}

exit;

sub helpfulexit {

  print <<EOS;

$0: keep the best HSP hits for each pair of genes

usage: $0 FILE

FILE HAS TO BE THE BLAST OUTPUT FILE
EOS

}
