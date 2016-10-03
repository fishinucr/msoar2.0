#!/usr/bin/perl -w
#
# $Id: cutbyname.pl,v 1.3 2003/03/18 00:02:31 cavet Exp $
#
# Pulls column(s) by name or index from character-delimited files.
# Outputs columns in the order specified.
# If the first line of the file has trailing null fields then the 
# correct number of columns must be provided with the -n option..
#
# OPTIONS:
# -c <column name> - specify names of columns to output
# -d <delimiter>   - specify input field delimiter
# -r <delimiter>   - specify output field delimiter
# -i               - use case-insensitive name matching
# -f 2,1,5-10      - specify numbers of columns to output 
# -h               - keep header line when using -c option (default)
# -H               - no header line when using -c option
# -n <NCOLS>       - explicitly set number of columns to <NCOLS>

use strict;

my($indelimiter,
   $outdelimiter,
   @columnnames,
   @found,
   @indices,
   @fields,
   $option,
   $i,
   $j,
   $caseinsensitiveflag,
   $keepheaderflag,
   $field,
   $ncols,
   $firstflag,
   );

$indelimiter = "\t";
$outdelimiter = "\t";
$caseinsensitiveflag = 0;
$keepheaderflag = 1;
$firstflag = 1;

if ($#ARGV == -1) { helpfulexit(); }

while ((scalar(@ARGV)) && ($ARGV[0] =~ /^-(\S+)/)) {

  $option = $1;
  shift(@ARGV);

  if ($option =~ /c/) { push(@columnnames, shift(@ARGV)); }
  if ($option =~ /d/) { $indelimiter = shift(@ARGV);      }
  if ($option =~ /r/) { $outdelimiter = shift(@ARGV);     }
  if ($option =~ /i/) { $caseinsensitiveflag = 1;         }
  if ($option =~ /h/) { $keepheaderflag = 1;              }
  if ($option =~ /H/) { $keepheaderflag = 0;              }
  if ($option =~ /n/) { $ncols = shift(@ARGV);            }
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

# IF CUTTING BY COLUMN HEADINGS
if (scalar(@columnnames)) {
  
  $firstflag = 0;

  if (!defined($_ = <>)) {
    die "$0: error, can't read header line\n";
  }
  chomp;
  @fields = split(/$indelimiter/, $_);
  $ncols  = $#fields + 1; # TO MAKE SURE WE ALWAYS SPLIT TO THE RIGHT NUMBER OF COLS

  for ($j = 0; $j <= $#columnnames; $j++) {
    
    for ($i = 0; $i <= $#fields; $i++) {
      
      if (($fields[$i] eq $columnnames[$j]) ||
	  ($caseinsensitiveflag && (lc($fields[$i]) eq lc($columnnames[$j])))) {
	
	$found[$j] = 1;
	push(@indices, $i);
      }
    }
  }
  # CHECK THAT ALL COLUMN NAMES WERE FOUND
  for ($j = 0; $j <= $#columnnames; $j++) {
    
    if (!defined($found[$j])) {
      die "$0: error, column name '$columnnames[$j]' not found\n";
    }
  }
  # DON'T PRINT HEADER LINE UNLESS ASKED
  if ($keepheaderflag) { 
    print join($outdelimiter, @fields[@indices]), "\n";
  }
}

# PROCESS INPUT DATA
while (<>) {

  chomp;

  # IF FIRST LINE, NEED TO COUNT COLUMNS
  if (!defined($ncols)) {

    @fields = split(/$indelimiter/, $_);
    $ncols = $#fields + 1;
    $firstflag = 0;

  } else {

    @fields = split(/$indelimiter/, $_, $ncols);
  }

  # OUTPUT APPROPRIATE COLUMNS
  print join($outdelimiter, @fields[@indices]), "\n";
}

sub helpfulexit {

  print <<EOS;

$0 pulls column(s) by name or index from 
character-delimited files.  It outputs columns in the order specified.
If the first line of the file has trailing null fields then the 
correct number of columns must be provided with the -n option..

OPTIONS:
-c <column name> - specify names of columns to output
-d <delimiter>   - specify input field delimiter
-r <delimiter>   - specify output field delimiter
-i               - use case-insensitive name matching
-f 2,1,5-10      - specify numbers of columns to output 
-h               - keep header line when using -c option (default)
-H               - no header line when using -c option
-n <NCOLS>       - explicitly set number of columns to <NCOLS>

EOS

  exit;
}
