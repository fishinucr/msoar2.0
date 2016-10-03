#!/usr/bin/perl -w
#
# $Id: joinby.pl,v 1.8 2001/12/29 02:06:29 gcavet Exp $
#
# sort gene position by its chromosome id and start coordinates

use strict;

my( $option,
    $f,
    $indelimiter,
    $outdelimiter,
    @fields,
    %hash,
    %pos,
    %chrom,
    $chromid,
    $start,
    $lineno,
    @tmp,
    $id,
    @sorted1,
    @sorted2,
    $aa,
    $bb,
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

    # PULL OUT KEY COLUMN
   $key = $fields[0];
   if (!defined($pos{$key})){
	$pos{$key} = $fields[4];
   }
    if( !defined($chrom{$key})){
	$chrom{$key} = $fields[2];
    }
}

# SORT START COORDINATE
@sorted1 = sort { $pos{$a} <=> $pos{$b} } keys %pos;

# SORT CHROMOSOME ID
@sorted2 = sort { 
		$aa = substr $chrom{$a}, 3;
		$bb = substr $chrom{$b}, 3;
		if($aa eq 'X') { $aa = 50; }
		elsif($aa eq 'Y') { $aa = 51; }
		if($bb eq 'X') { $bb = 50; }
        	elsif($bb eq 'Y') { $bb = 51; }
		return $aa <=> $bb; 	
	        } @sorted1;

$lineno = 0;
$id ="";

foreach  $key( @sorted2 ){
    $lineno = $lineno + 1;
    print $lineno,"\t",$key,"\n";
}

exit;



sub helpfulexit {

  print <<EOS;

$0 sort gene by its start coordinate
Usage: $0 FILE


EOS
     
  exit;
}
