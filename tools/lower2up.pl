#!/usr/bin/perl -w
#
# change all the lower case to up case

use strict;

open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";
shift(@ARGV);

 
while (<IN>) {

    chomp;
    tr/a-z/A-Z/;
    print $_,"\n";
}		    

exit;

