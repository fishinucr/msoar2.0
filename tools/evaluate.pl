#!/usr/bin/perl -w
#
# check ortholog assignment result and calculate the number of
# true positives and unkonwns 

use strict;

my( $TP,
    $UK,	
    $total,
    $indelimiter,
    @fields,
    );

$indelimiter = "\t";


  open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";
  shift(@ARGV);

$TP = 0;
$UK = 0;
$total = 0;
 
while (<IN>) {

    # BREAK INPUT LINE INTO FIELDS
    chomp;
    tr/a-z/A-Z/;
    @fields = split(/$indelimiter/, $_);
    $total ++;	    	
    if( $fields[1] eq $fields[3] ){
	$TP ++; 
    }
    elsif( $_ =~ /LOC/ ||  $_ =~ /RIK/ ||$fields[0] =~ /^AB\d\d\d\d\d\d/ || $fields[1] =~ /^AB\d\d\d\d\d\d/ ||$fields[0] =~ /^BC\d\d\d\d\d\d/ || $fields[1] =~ /^BC\d\d\d\d\d\d/ || $fields[1] eq "-" ||$fields[0] =~ /^AF\d\d\d\d\d\d/ || $fields[1] =~ /^AF\d\d\d\d\d\d/ ||$fields[0] =~ /^AK\d\d\d\d\d\d/ || $fields[1] =~ /^AK\d\d\d\d\d\d/ ||$fields[0] =~ /^AL\d\d\d\d\d\d/ || $fields[1] =~ /^AL\d\d\d\d\d\d/ ||$fields[0] =~ /^AI\d\d\d\d\d\d/ || $fields[1] =~ /^AI\d\d\d\d\d\d/||$fields[0] =~ /^AU\d\d\d\d\d\d/ || $fields[1] =~ /^AU\d\d\d\d\d\d/ ||$fields[0] =~ /^AW\d\d\d\d\d\d/ || $fields[1] =~ /^AW\d\d\d\d\d\d/ ||$fields[0] =~ /^AY\d\d\d\d\d\d/ || $fields[1] =~ /^AY\d\d\d\d\d\d/ || $fields[1] =~ /^DQ\d\d\d\d\d\d/ ||$fields[0] =~ /^DQ\d\d\d\d\d\d/ ){
	$UK ++;    
    }	    
	    
}

print "TP=", $TP,"\n","UK=", $UK,"\n", "FP=", $total - $TP - $UK, "\n","Total=",$total,"\n";

exit;

