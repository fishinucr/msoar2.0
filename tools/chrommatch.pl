#!/usr/bin/perl -w
# truncate certain field of the input file

use strict;

my( $option,
    $indelimiter,
    $key,
    %hits,
    $key2,
    $chomd1,
    $chomd2,
    $num,
    @fields,
    $a1,
    $a2,
    );

$indelimiter = "\t";

if ($#ARGV < 0) { helpfulexit(); }


open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";
shift(@ARGV);


while (<IN>) {
    # BREAK INPUT LINE INTO FIELDS
    chomp;
    @fields = split(/$indelimiter/, $_);
    $chomd1 = substr( $fields[0], 3,length($fields[0])-1 );
    $chomd2 = substr( $fields[1], 3,length($fields[1])-1 );
    if( $chomd1 eq "X" ){
	    $chomd1 = 23;	    
       }
    if( $chomd1 eq "Y" ){
           $chomd1 = 24;
     }
    if( $chomd2 eq "X" ){
          $chomd2 = 20;
     }
    if( $chomd2 eq "Y" ){
           $chomd2 = 21;
     }
    if ((!defined($hits{$chomd1}{$chomd2}))) {
	       $hits{$chomd1}{$chomd2} = 1;
     }
     else{
	     $hits{$chomd1}{$chomd2}++;
     }
												 
			 
}

											       
for $key (sort {$a<=>$b} keys %hits){
	for $key2 (sort{$a <=> $b} keys %{$hits{$key}} ){
		$a1 = $key;
		$a2 = $key2;
		$num = $hits{$key}{$key2};
		if( $key == 23 ){
			$a1 = "X";
		}
		elsif( $key == 24 ){
			$a1 = "Y";
		}
		
		if( $key2 == 20 ){
			$a2 = "X";
		}
		elsif( $key2 == 21 ){
			$a2 = "Y";
		}
		
		if( $num > 250 ){
			 print $a1." vs. ".$a2,"\t",$num,"\n";
	 	}
		else{
			print " ","\t",$num,"\n";		
		}
	 }
	 
}

exit;

sub helpfulexit {

  print <<EOS;

$0: truncate certain field of the input file

usage: $0 -f|b POSITION [-Options] FILE COLUMN

Options:

-i DELIMITER - set input column delimiter
-d DELIMITER - set output column delimiter
-h           - print this help information

EOS
exit;
}
