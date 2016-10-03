#!/usr/bin/perl -w
#
# Does vector arithmetic and logical filtering on text streams

$indelimiter = "\t";
@data = ();
$cmdlineflag = 0; # whether we're evaluating a single expression from the command line
$filterflag = 0;  # whether we're filtering the data rather than modifying it

# PROCESS OPTIONS
while (($#ARGV >= 0) && ($ARGV[0] =~ /^-(\w+)/)) {

  $option = $1;
  shift(@ARGV);

  if ($option =~ /d/) { $indelimiter = shift(@ARGV); }
  if ($option =~ /e/) { $cmdlineflag = 1; }
  if ($option =~ /f/) { $filterflag  = 1; }
  if ($option =~ /n/) { $ncols       = shift(@ARGV); }
}

# PROCESS EXPRESSION
if ($#ARGV < 0) { helpfulexit(); }
$expr = shift(@ARGV);
$expr =~ s/\b(c(\d+))/'$data['.($2 - 1).']'/ge;

# PROCESS INPUT FROM COMMAND LINE AND EXIT IF DESIRED
if ($cmdlineflag) {
  
  @data = @ARGV;
  print eval($expr), "\n";
  exit;
}

# PROCESS INPUT FROM STDIN
while (<>) {

  chomp;
  if (defined($ncols)) # USER HAS EXPLICITLY SET NUMBER OF INPUT COLUMNS
  {
    @data = split(/$indelimiter/, $_, $ncols);
  }
  else
  {
    @data = split(/$indelimiter/, $_);
  }

  if ($filterflag) {

    if (eval($expr)) { print $_, "\n"; }

  } else {    
    print eval($expr), "\n";
  }
}

exit;

sub helpfulexit {

  print <<EOS;
usage: $0 [ OPTIONS ] EXPRESSION [ FILE1 [ FILE2... ] ]

$0 manipulates (by default) or filters tabular text data one line 
at a time.  It can also provide one-off evaluation of expressions.
It breaks up each input line into fields delimited by a 
user-defined string (default is a tab) and allows the user to 
access these fields.  Successive fields may be referred to by the
strings c1, c2...cN in EXPRESSION.
When manipulating data, the value of EXPRESSION is output for 
each line in the input.
When filtering, each line in the input is output unchanged if
and only if EXPRESSION evaluates to true.
When performing one-off evaluation, EXPRESSION is evaluated and
the result output.

Options:
-f              : perform filtering (default is to manipulate data, not filter)
-d <DELIMITER>  : use input field delimiter <DELIMITER>
-e              : output evaluation of EXPRESSION directly, ignoring other input
-n <NCOLS>      : assume input has NCOLS columns (saves trailing null fields)

For example:

> more test
1       2
3       4
> vector.pl 'c2/c1' test
2
1.33333333333333
> vector.pl -f 'c1 + c2 == 3' test
1       2
> vector.pl -e '1 / 3'
0.333333333333333
> vector.pl '((c2 > (c1**2))? join("\t", c1, c2) : join("\t", c2, c1))' test
1       2
4       3

EOS

  exit;
}
