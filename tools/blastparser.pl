#!/usr/bin/perl

# Version 1.1
# BLASTALL output parser
# Extracts the query  names, matching sequence names, scores,
# E-values, lengths from BLASTALL output files.
# Works for blast 2.2 14 as well.

$S_cutoff = 50;  # default bitscore cutoff
$E_cutoff = 1.0; # default evalue cutoff
$L_cutoff = 0.5; # default length cutoff
$N_cutoff = 5;  # default number of hits to be kept, -1 represents that keep all the hits
$chaining = 0;   # indicate whether consider multiple HSPs or not 

open(IN, '<'.$ARGV[0]) ||
    die "$0: error, can't open input file ", $ARGV[0], "\n";

$S_cutoff = $ARGV[1];
$E_cutoff = $ARGV[2];
$L_cutoff = $ARGV[3];
$N_cutoff = $ARGV[4];
$chaining = $ARGV[5];

while (<IN>) {
  
  $ifprint = 0;  
  QUERY:
  if (/^Query=/) { 
     if( $ifprint) {# print previour query
       #print "Print", $_,"\n";
       $ifprint = 0; # reset the indicator
       for $a(1..$j){
         for $b(1..($HSP_count[$a]-1)){
           while($qstart[$a][$b] > $qstart[$a][$b+1]){
             $qs = $qstart[$a][$b];
             $qe = $qend[$a][$b];
             $ms = $mstart[$a][$b];
             $me = $mend[$a][$b];
             $qstart[$a][$b] = $qstart[$a][$b+1];
             $qend[$a][$b]   = $qend[$a][$b+1];
             $mstart[$a][$b] = $mstart[$a][$b+1];
             $mend[$a][$b]   = $mend[$a][$b+1];
             $qstart[$a][$b+1] = $qs;
             $qend[$a][$b+1]   = $qe;
             $mstart[$a][$b+1] = $ms;
             $mend[$a][$b+1]   = $me;
             --$b if ($b > 1);
           }
         }
       }
       for $k(1..($j-1)){
         while($score[$k] < $score[$k+1]){
           $tempM   = $match[$k];
           $tempML  = $HSP_length[$k];
           $tempTL  = $match_length[$k];
           $tempS   = $score[$k];
           $tempE   = $E_value[$k];
           $tempID  = $id[$k];
           $tempSIM = $pos[$k];
           $tempGAP = $ngaps[$k];
           $tempHC  = $HSP_count[$k];
           @tempQS  = @qstart[$k];
           @tempQE  = @qend[$k];
           @tempMS  = @mstart[$k];
           @tempME  = @mend[$k];
           $match[$k] = $match[$k+1];
           $HSP_length[$k] = $HSP_length[$k+1];
           $match_length[$k] = $match_length[$k+1];
           $score[$k] = $score[$k+1];
           $E_value[$k] = $E_value[$k+1]; 
           $id[$k] = $id[$k+1];
           $pos[$k] = $pos[$k+1];
           $ngaps[$k] = $ngaps[$k+1];
           $HSP_count[$k] = $HSP_count[$k+1];
           @qstart[$k] = @qstart[$k+1];
           @qend[$k] = @qend[$k+1];
           @mstart[$k] = @mstart[$k+1];
           @mend[$k] = @mend[$k+1];

           $match[$k+1] = $tempM;
           $HSP_length[$k+1] = $tempML;
           $match_length[$k+1] = $tempTL;
           $score[$k+1] = $tempS;
           $E_value[$k+1] = $tempE;
           $id[$k+1] = $tempID;
           $pos[$k+1] = $tempSIM;
           $ngaps[$k+1] = $tempGAP;
           $HSP_count[$k+1] = $tempHC;
           @qstart[$k+1] = @tempQS;
           @qend[$k+1] = @tempQE;
           @mstart[$k+1] = @tempMS;
           @mend[$k+1] = @tempME;
           --$k if ($k > 1);
         }
       }
       $i = $N_cutoff if ($N_cutoff != -1 and $i > $N_cutoff);
       for $k(1..$i){
         next if ($score[$k] < $S_cutoff);
         next if ($E_value[$k] > $E_cutoff);
         next if ($HSP_length[$k] < $L_cutoff*$query_length and $HSP_length[$k] < $L_cutoff*$match_length[$k]);

         print  $query,"\t";
         print  $match[$k],"\t";
         print  $score[$k],"\t";
         $score[$k] = 25 if ($score[$k] < 25); # To avoid underflow of the E-value
         #$E = $M*$N*(2**-$score[$k]);
         print $E_value[$k],"\t";
         print $query_length,"\t";

         # Alignments were found, print match_length, HSP_length,
         # identity %, match % and query-target ratio
         if ($HSP_length[$k]){
           print  $match_length[$k];
           $HSP_length[$k] -= $ngaps[$k] if ($ngaps[$k]);
           printf ("\t%.0f", $HSP_length[$k]);  # The length of actually aligned amino  acids
           for $hsp (1..$HSP_count[$k]){ # Query HSPs
             print "\t$qstart[$k][$hsp]-$qend[$k][$hsp]" if ($qstart[$k][$hsp]);
           }
#          print "\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
#          for $hsp (1..$HSP_count[$k]){ # Match HSPs
#            print "\t$mstart[$k][$hsp]-$mend[$k][$hsp]" if ($mstart[$k][$hsp]);
#          }
         }
         print  "\n";
       }

     }
      ############New query##################
      # Get the name of query  
      @tmp = split(/\s+/);
      $query = $tmp[1];
      # Get length of query
      until(/letters/){
         $_ = <IN>; # Take next line
         m/\s+\((\d+)/;
         $query_length = $1;
      }
   }

  if (/^Sequences producing/) { # Has some hits
    $ifprint = 1;
    $_ = <IN>; # Now on a blank line
    $i = 0; # reset match counter
    while (<IN>) { # Loops until the list of matches
      if (/^\s*$/){# What if blank line? Let's look one more line:
        $_ = <IN>; # Read next line after the end of match list
        last if /^\s*$/;
        if (/^>/){ # Now on alignment line
          $j=0; # alignment count
          goto ALIGNMENTS;
        }
      }
      ++$i;
      chomp;
      @tmp=split(/\s+/,$_);    # Take matching name, score and E
      $match[$i] = $tmp[0];   # match name
      $score[$i] = $tmp[1];
      $E_value[$i] = $tmp[2];
      @tmp = split(//,$E_value[$i]);# Check for correct writing of E-values
      $E_value[$i] = "1" . $E_value[$i] if $tmp[0] eq "e";
    }# End of match list, continue analyzing alignments

    # Now continue to find length of each matching protein/NA
    $j=0;
    while (<IN>){ # Loops for each query sequence
      goto QUERY if (/^Query=/);       
      ALIGNMENTS:
      if(/^>/){ # Found lines with alignments
        ++$j;$HSP_count[$j] = 0;
        @qstart[$j] = @qend[$j] = ();
        @mstart[$j] = @mend[$j] = ();
        chomp;
        s/\>//;
        @tmp = split(/\s+/);
        print "Error1",$match[$j],"\n tmp[0]=",$tmp[0]," \n" if($match[$j] ne $tmp[0]);
        until(/Length =/){ # repeat until target length is found
          $_ = <IN>; # Next line
          @tmp = split(/\s+/);
          $match_length[$j] = $tmp[3];
          print "Error2 \n" if($tmp[1] ne "Length");
          print "Error3" if($match_length[$j]==NULL or $match_length[$j] =~ /\D/);
        }
        ##############
        $_ = <IN>;$_ = <IN>; # Score line
        chomp;
        s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
        @tmp = split(/\/|\s+/);
        $score[$j] = $tmp[2]; # in bits. Overrides the value that was read from initial list
        #$E_value[$j] = $tmp[6]; # in Evalue. Overrides the value that was read from initial list
        
        $_ = <IN>; # Now on identity line
        chomp;
        s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
        @tmp = split(/\/|\s+/);
        $HSP_length[$j] = $tmp[3];
        $id[$j] = $tmp[2];     # Number of identities in match region
        $pos[$j] = $tmp[6];   # Number similarities in match region
        $ngaps[$j] = $tmp[10]; # Number of gaps

        $n = $HSP_count[$j] = 1;  #segment number
        $_ = <IN>;$_ = <IN>; # Get query start
        chomp;
        @tmp = split(/\s+/);
        $qstart[$j][$n] = $tmp[1]; # Query start
        $qend[$j][$n] = $tmp[3];   # Query end
        $_ = <IN>;$_ = <IN>;
        chomp;
        @tmp = split(/\s+/);
        $mstart[$j][$n] = $tmp[1]; # Match start
        $mend[$j][$n] = $tmp[3];   # Match end
        $_ = <IN>;$_ = <IN>;
        for(;;){
          last if /^\s*$/; # empty line after alignment
          chomp;
          @tmp = split(/\s+/);
          $qend[$j][$n] = $tmp[3];
          $_ = <IN>;$_ = <IN>;
          chomp;@tmp = split(/\s+/);
          $mend[$j][$n] = $tmp[3];
          $_ = <IN>;$_ = <IN>; 
        }
        #print $j,"\t",$n,"\t",$HSP_count[$j],"\t",$match_length[$j],"\t",$score[$j],"\t",$HSP_length[$j],"\t",$id[$j],"\t",$pos[$j] ,"\t",$ngaps[$j],"\n";
	#print $qstart[$j][$n],"\t",$qend[$j][$n],"\t",$mstart[$j][$n],"\t",$mend[$j][$n],"\n";
      }

      if( $chaining and /^ Score =/ and ($n < 5)){ # This happens if we have more than 1 segment reported
        ++$n; ++$HSP_count[$j];# HSP numbering
        chomp;
        s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
        @score_line = split(/\/|\s+/);
        $_ = <IN>; # Now on identity line
        chomp;
        s/\(|\)|,|=//g; # Remove brackets, commas and equalsigns
        @id_line = split(/\/|\s+/);
        $_ = <IN>;$_ = <IN>; # Get query start
        chomp;@tmp = split(/\s+/);
        $qstart[$j][$n] = $tmp[1];
        $qend[$j][$n] = $tmp[3];
        $_ = <IN>;$_ = <IN>;
        chomp;@tmp = split(/\s+/);
        $mstart[$j][$n] = $tmp[1]; # Match start
        $mend[$j][$n] = $tmp[3];   # Match end
        $_ = <IN>;$_ = <IN>;
        for(;;){
          if (/^\s*$/){ # 2 empty lines after alignment
            #Decide if we want to accept this HSP
            # Check for overlap
            $accepted = 1;
            for $c(1..($n-1)){
              # Which HSP is N-terminal?
              if ($qstart[$j][$n] < $qstart[$j][$c]){
                $qstart1 = $qstart[$j][$n];
                $qend1   = $qend[$j][$n];
                $mstart1 = $mstart[$j][$n];
                $mend1   = $mend[$j][$n];
                $qstart2 = $qstart[$j][$c];
                $qend2   = $qend[$j][$c];
                $mstart2 = $mstart[$j][$c];
                $mend2   = $mend[$j][$c];
              }
              else{
                #if ($qstart[$j][$n]+$mstart[$j][$n] > $qstart[$j][$c]+$mstart[$j][$c]){
                $qstart2 = $qstart[$j][$n];
                $qend2   = $qend[$j][$n];
                $mstart2 = $mstart[$j][$n];
                $mend2   = $mend[$j][$n];
                $qstart1 = $qstart[$j][$c];
                $qend1   = $qend[$j][$c];
                $mstart1 = $mstart[$j][$c];
                $mend1   = $mend[$j][$c];
              }
              #print "Working with $match[$j], HSP $n, checking against HSP $c\n";
              #print "$qstart1-$qend1 $mstart1-$mend1   $qstart2-$qend2 $mstart2-$mend2\n";
              $qlength1 = $qend1 - $qstart1;
              $qlength2 = $qend2 - $qstart2;
              $qlength  = ($qlength1 < $qlength2)?$qlength1:$qlength2;

              $mlength1 = $mend1 - $mstart1;
              $mlength2 = $mend2 - $mstart2;
              $mlength  = ($mlength1 < $mlength2)?$mlength1:$mlength2;

              $qoverlap = $qstart2 - $qend1;
              $moverlap = $mstart2 - $mend1;

             #print "$qoverlap\t$qlength\t$moverlap\t$mlength\n";
              if (!($qlength * $qlength)){ # Zero length - why?
                --$n; # Discard this HSP
                --$HSP_count[$j];
                $accepted = 0;
                last;
              }
              # 5% overlap is allowed:
              if (($qoverlap/$qlength < -0.05) or ($moverlap/$mlength < -0.05)){
                --$n; # Discard this HSP
                --$HSP_count[$j];
               #print "########### Rejected #########\n";
                $accepted = 0;
                last;
              }
            }
            $score[$j]      += $score_line[2] if($accepted); # add score of this HSP
	    $E_value[$j]    *= $score_line[6] if($accepted);
	    #print  $j,"\t", $score_line[6],"\t",$E_value[$j],"\n";
            $HSP_length[$j] += $id_line[3]    if($accepted);
            $id[$j]         += $id_line[2]    if($accepted);
            $pos[$j]        += $id_line[6]    if($accepted);
            $ngaps[$j]      += $id_line[10]   if($accepted);
            last;
          }
          chomp;@tmp = split(/\s+/);
          $qend[$j][$n] = $tmp[3];
          $_ = <IN>;$_ = <IN>;
          chomp;@tmp = split(/\s+/);
          $mend[$j][$n] = $tmp[3];
          $_ = <IN>;$_ = <IN>; 
        }
      }
    }
  }
}
#########print out the last query#############
if( $ifprint) {# print previour query
  $ifprint = 0; # reset the indicator
  for $a(1..$j){
    for $b(1..($HSP_count[$a]-1)){
      while($qstart[$a][$b] > $qstart[$a][$b+1]){
        $qs = $qstart[$a][$b];
        $qe = $qend[$a][$b];
        $ms = $mstart[$a][$b];
        $me = $mend[$a][$b];
        $qstart[$a][$b] = $qstart[$a][$b+1];
        $qend[$a][$b]   = $qend[$a][$b+1];
        $mstart[$a][$b] = $mstart[$a][$b+1];
        $mend[$a][$b]   = $mend[$a][$b+1];
        $qstart[$a][$b+1] = $qs;
        $qend[$a][$b+1]   = $qe;
        $mstart[$a][$b+1] = $ms;
        $mend[$a][$b+1]   = $me;
        --$b if ($b > 1);
      }
    }
  }
  for $k(1..($j-1)){
    while($score[$k] < $score[$k+1]){
      $tempM   = $match[$k];
      $tempML  = $HSP_length[$k];
      $tempTL  = $match_length[$k];
      $tempS   = $score[$k];
      $tempE   = $E_value[$k];
      $tempID  = $id[$k];
      $tempSIM = $pos[$k];
      $tempGAP = $ngaps[$k];
      $tempHC  = $HSP_count[$k];
      @tempQS  = @qstart[$k];
      @tempQE  = @qend[$k];
      @tempMS  = @mstart[$k];
      @tempME  = @mend[$k];

      $match[$k] = $match[$k+1];
      $HSP_length[$k] = $HSP_length[$k+1];
      $match_length[$k] = $match_length[$k+1];
      $score[$k] = $score[$k+1];
      $E_value[$k] = $E_value[$k+1];
      $id[$k] = $id[$k+1];
      $pos[$k] = $pos[$k+1];
      $ngaps[$k] = $ngaps[$k+1];
      $HSP_count[$k] = $HSP_count[$k+1];
      @qstart[$k] = @qstart[$k+1];
      @qend[$k] = @qend[$k+1];
      @mstart[$k] = @mstart[$k+1];
      @mend[$k] = @mend[$k+1];

      $match[$k+1] = $tempM;
      $HSP_length[$k+1] = $tempML;
      $match_length[$k+1] = $tempTL;
      $score[$k+1] = $tempS;
      $E_value[$k+1] = $tempE;
      $id[$k+1] = $tempID;
      $pos[$k+1] = $tempSIM;
      $ngaps[$k+1] = $tempGAP;
      $HSP_count[$k+1] = $tempHC;
      @qstart[$k+1] = @tempQS;
      @qend[$k+1] = @tempQE;
      @mstart[$k+1] = @tempMS;
      @mend[$k+1] = @tempME;
      --$k if ($k > 1);
    }
  }
  $i = $N_cutoff if ($N_cutoff != -1 and $i > $N_cutoff);
  for $k(1..$i){
    next if ($score[$k] < $S_cutoff);
    next if ($E_value[$k] > $E_cutoff);
    next if ($HSP_length[$k] < $L_cutoff*$query_length);
    next if ($HSP_length[$k] < $L_cutoff*$match_length[$k]);

    printf  ("%-25.25s", $query);
    printf  ("\t%-25.25s" , $match[$k]);
    printf  ("\t%.1f", $score[$k]);
    $score[$k] = 25 if ($score[$k] < 25); # To avoid underflow of the E-value
    #$E = $M*$N*(2**-$score[$k]);
    printf ("\t%3.1g", $E_value[$k]);
    print  "\t$query_length";

    # Alignments were found, print match_length, HSP_length,
    # identity %, match % and query-target ratio
    if ($HSP_length[$k]){
      print  "\t$match_length[$k]";
      $HSP_length[$k] -= $ngaps[$k] if ($ngaps[$k]);
      printf ("\t%.0f", $HSP_length[$k]);  # The length of actually aligned amino  acids
      for $hsp (1..$HSP_count[$k]){ # Query HSPs
        print "\t$qstart[$k][$hsp]-$qend[$k][$hsp]" if ($qstart[$k][$hsp]);
      }
#     print "\n\t\t\t\t\t\t\t\t\t\t\t\t\t\t";
#     for $hsp (1..$HSP_count[$k]){ # Match HSPs
#       print "\t$mstart[$k][$hsp]-$mend[$k][$hsp]" if ($mstart[$k][$hsp]);
#     }
    }
    print  "\n";
  }
}


