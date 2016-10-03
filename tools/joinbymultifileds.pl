#!/usr/bin/perl -w
#
# $Id: joinby.pl,v 1.8 2001/12/29 02:06:29 gcavet Exp $
#
# joins text files based on common values in key columns

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
    $field,
    @indices,
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
  if ($option =~ /f/) { # want numbered columns
    if (length($')) { $option = $'; }
    else            { $option = shift(@ARGV); }

    while ($option =~ /[\d-]+/) {

      $field = $&;
      $option = $';

      if ($field =~ /^\d+$/) {
	
	push(@indices, ($field-1)); 
      
      } elsif ($field =~ /^(\d+)-(\d+)$/) {

	push(@indices, (($1-1) .. ($2-1))); 

      } else {

	die "$0: error, can't parse field specification $field\n";
      }
    }
  }
}

if ($#ARGV < 0) { helpfulexit(); }

# LOAD UP DATA FROM INPUT FILES
while ($#ARGV >= 0) {

  open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";
  shift(@ARGV);
    
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
    $key = join("|",@fields[@indices]);
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
   shift( @{$listref} ); 
    print join($outdelimiter, @{$listref}), "\n";
  }
}

exit;

sub helpfulexit {

  print <<EOS;

$0: joins delimited text files on specified key columns

usage: $0 [-options] -f COLS FILE1 [FILE2 ... [FILEn] ]

Options:

-o           - outer join with respect to the first file
-i DELIMITER - set input column delimiter
-d DELIMITER - set output column delimiter
-h           - print this help information
-f           - join by multiple columns

EOS

}
