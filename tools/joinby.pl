#!/usr/bin/perl -w
#
# joins text files based on common values in key columns
#
# Zheng Fu
# Algorithms and Computational Biology Lab
# University of California, Riverside


use strict;

my( $option,
    $indelimiter,
    $outdelimiter,
    $firstflag,     # whether we're loading the first file
    $outerflag,
    $column,
    $filearrayref,
    $linearrayref,
    @strings,
    @listrefs,
    $listref,
    @newlistrefs,
    $filei,
    @fields,
    %hash,
    $key,
    );

$indelimiter = "\t";
$outdelimiter = "\t";
$filei = -1;
$outerflag = 0;

while (($#ARGV >=0) && ($ARGV[0] =~ /^-(\w+)/)) {

  $option = $1;
  shift(@ARGV);

  if ($option =~ /o/) { $outerflag = 1; }
  if ($option =~ /i/) { $indelimiter = shift(@ARGV); }
  if ($option =~ /d/) { $outdelimiter = shift(@ARGV); }
  if ($option =~ /h/) { helpfulexit(); }
}

if ($#ARGV < 0) { helpfulexit(); }

# LOAD UP DATA FROM INPUT FILES
while ($#ARGV >= 0) {

  if ($#ARGV == 0) { 
    die "$0: error, arguments must include pairs composed of filename and column number\n";
  }

  open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";
  shift(@ARGV);

  $column = shift(@ARGV);
  unless ($column =~ /^[1-9]\d*/) {
    die "$0: error, column indices must be positive integers (the offending value was $column)\n";
  }

  $column -= 1; # CORRECT FOR ZERO-BASED ARRAY INDEXING
    
  $filei++; # COUNT WHICH FILE WE'RE ON
 
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
    $key = splice(@fields, $column, 1);

    # ADD NEW DATA TO HASH ENTRY IF FIRST FILE OR FIND MATCH
    if (($filei == 0) || (defined($hash{$key}))) {

      # ADD DATA TO ARRAY OF ARRAYS IN ELEMENT FOR THIS FILE
      push(@{$hash{$key}->[$filei]}, [ @fields ]);
    }
  }

  # ADD 'NULLS' IF OUTER JOIN IS DESIRED
  if ($outerflag) { 

    # BUILD LIST OF APPROPRIATE NUMBER OF EMPTY STRINGS
    foreach (@fields) { $_ = ''; }

    while (($key, $filearrayref) = each(%hash)) {

      if (!defined($filearrayref->[$filei])) {
	
	push(@{$filearrayref->[$filei]}, [ @fields ]);
      }
    }
  }
}

# NOW ALL FILES HAVE BEEN READ IN AND CAN DO OUTPUT
while (($key, $filearrayref) = each(%hash)) {

  # BUILD ARRAY OF STRINGS TO BE OUTPUT
  @listrefs = ( [ "$key" ] );
  @newlistrefs = ();
  
  for (0 .. $filei) {

    # IF THERE ARE ANY ENTRIES FOR THIS KEY FROM THIS FILE
    if (defined($filearrayref->[$_])) {

      foreach $listref (@listrefs) {
      
	foreach $linearrayref (@{$filearrayref->[$_]}) {

	  push(@newlistrefs, [ @{$listref}, @{$linearrayref} ] );
	}
      }
    }
    @listrefs = @newlistrefs;
    @newlistrefs = ();
  }

  # OUTPUT ALL LISTS FOR THIS KEY
  foreach $listref (@listrefs) {
    print join($outdelimiter, @{$listref}), "\n";
  }
}

exit;

sub helpfulexit {

  print <<EOS;

$0: joins delimited text files on specified key columns

usage: $0 [-o] FILE1 COL1 [FILE2 COL2 ... [FILEn COLn] ]

Options:

-o           - outer join with respect to the first file
-i DELIMITER - set input column delimiter
-d DELIMITER - set output column delimiter
-h           - print this help information

EOS

}
