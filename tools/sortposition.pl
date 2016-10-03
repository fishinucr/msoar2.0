#!/usr/bin/perl -w
#
# $Id: joinby.pl,v 1.8 2001/12/29 02:06:29 gcavet Exp $
#
# sort gene position by its start coordinates

use strict;

my( $option,
    $f,
    $indelimiter,
    $outdelimiter,
    @fields,
    %hash,
    %pos,
    $key,
    $lineno,
    @tmp,
    $id,
    @sorted_keys,
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

    # PULL OUT KEY COLUMN
    $key = join(":",$fields[0],$fields[1]);
    $pos{$key} = $fields[6];
    $hash{$key} = join($outdelimiter, @fields);
}
@sorted_keys = sort { $pos{$a} <=> $pos{$b} } keys %pos;
$lineno = 0;
$id ="";

foreach  $key( @sorted_keys ){
    @tmp = split(/:/, $key);
 
    if( $tmp[0] ne $id ){
	$lineno = $lineno + 1;
	$id = $tmp[0];
    }
    print $lineno,"\t",$hash{$key},"\n";
}

exit;



sub helpfulexit {

  print <<EOS;

$0 sort gene by its start coordinate
Usage: $0 FILE


EOS
     
  exit;
}
