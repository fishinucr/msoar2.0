#!/usr/bin/perl -w
#
# $Id: joinby.pl,v 1.8 2001/12/29 02:06:29 gcavet Exp $
#
# joins text files based on common values in key columns

use strict;

my(
    $indelimiter,
    $outdelimiter,
    @fields,
    $key,
    $prekey,	
    %hash,
    );

$indelimiter = "\t";
$outdelimiter = "\t";
$prekey = "initial";

if ($#ARGV < 0) { helpfulexit(); }

# LOAD UP DATA FROM INPUT FILES
open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";
shift(@ARGV);
    
 
while (<IN>) {

    # BREAK INPUT LINE INTO FIELDS
    chomp;
    @fields = split(/$indelimiter/, $_);

    # PULL OUT KEY COLUMN
    $key = join("|",@fields[0,1]);
    if( !defined( $hash{$key}  )){
	$hash{$key} = join($outdelimiter, @fields);	
    }
  
}

for $key ( keys %hash ) {
	print $hash{$key},"\n";
}
						
close IN;

exit;

sub helpfulexit {

  print <<EOS;

$0: Only keep the best HSPs from blast -m 8 output file

usage: $0 [-options] FILE

EOS

}
