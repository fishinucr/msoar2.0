#!/usr/bin/perl -w
#
# $Id: joinby.pl,v 1.8 2001/12/29 02:06:29 gcavet Exp $
#
# joins text files based on common values in key columns

use strict;

my( $option,
    $indelimiter,
    $outdelimiter,
    $assignable,
    %hash,
    $key,
    );

$indelimiter = "\t";
$outdelimiter = "\t";

if ($#ARGV < 0) { helpfulexit(); }

# LOAD UP DATA FROM INPUT FILES
open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";

while (<IN>) {
    
    # BREAK INPUT LINE INTO FIELDS
    chomp;
    tr/a-z/A-Z/;
    $key= $_;

    # PULL OUT KEY COLUMN
    $hash{$key} = 1;
}

open(IN, '<'.$ARGV[1]) ||
    die "$0: error, can't open input file ", $ARGV[1], "\n";

$assignable = 0;

while (<IN>) {

    # BREAK INPUT LINE INTO FIELDS
    chomp;
    tr/a-z/A-Z/;	
    $key= $_;
    if (defined($hash{$key})){
	print "$key","\n";
	$assignable ++;
    }	
}

print "assignable=",$assignable,"\n";
exit;

sub helpfulexit {

  print <<EOS;

usage: $0 FILE1	FILE2


EOS

}
