#!/usr/bin/perl
# 
# Copyright 2013, Gary W. Stuart <gstuart@indstate.edu>
# Rename this file as phyppa.pl
# Compare Phyper Subsamples: Phyppa3
#
# This is a free script: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
#
# This is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# A copy of the GNU General Public License can be found at <http://www.gnu.org/licenses/>.
#

#
# Compares fragments of a genome with another genome fragments to find identical fragments and similar fragments. 
# 
#
# By Gary W. Stuart <gstuart@indstate.edu>
#

$version = "Phyppa3";
print "\nPhylogenomic Profile Pairwise Analysis: $version\n\n";

#input first file name for analysis...
print "Enter filename #1 (i.e. GenSpeSubEnzXm.txt): ";
$infilename1 = <>;
chomp $infilename1;

open (INFILE1, "<$infilename1");
while (<INFILE1>) {
  chomp $_;
  push (@infile1, $_);
}
print "\n",$#infile1+1," fragments";

#input second filename for analysis...
print "\n\nEnter filename #2 (i.e. GenSpeSubEnzXm.txt): ";
$infilename2 = <>;
chomp $infilename2;

push (@namesort, substr($infilename1, 3, 6), substr($infilename2, 3, 6));
@namesort = sort(@namesort);

$outfile1 = substr($infilename1, 9, 5) . "nomatcha" . $namesort[0] . $namesort[1] . ".txt";
#$outfile2 = substr($infilename1, 9, 5) . "cnv" . $namesort[0] . $namesort[1] . ".txt";
$outfile3 = substr($infilename1, 9, 5) . "pairs" . $namesort[0] . $namesort[1] . ".txt";
$outfile4 = substr($infilename1, 9, 5) . "nomatchb" . $namesort[0] . $namesort[1] . ".txt";
$outfile5 = substr($infilename1, 9, 5) . "sum" . $namesort[0] . $namesort[1] . ".txt";

$Enz = substr($infilename1, 9, 3);
$sym = 0;
if ($Enz == "Bpl"){$sym = 1;}
if ($Enz == "Fal"){$sym = 1;}
if ($Enz == "Hae"){$sym = 1;}
if ($Enz == "Hin"){$sym = 1;}
if ($Enz == "Alf"){$sym = 1;}

open (INFILE2, "<$infilename2");
open (NOMATCH1, ">$outfile1");
#open (CNV, ">$outfile2");
open (PAIRS, ">$outfile3");
open (NOMATCH2, ">$outfile4");
open (SUM, ">$outfile5");

push (@varlist, -1);
$family = 1;

while (<INFILE2>) {
  chomp $_;
  push (@infile2, $_);
}
print "\n", $#infile2+1," fragments\n\n";

$size = length($infile2[1])-15;

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$StarTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 

print "\n\nFinding identical fragments...\n\n";


$loops1 = $#infile1+1;
print "Listing Identical Pairs in $outfile3...CNV's in $outfile2...\n\n";

for ($i = 0; $i < $#infile1+1; $i++){
  if ($i % 200 == 0){
    print " $i\n";
  }
  $test1 = $infile1[$i];
  for ($j = 0; $j < $#infile2+1; $j++){
    $test2 = $infile2[$j];
	if ($sym == 1){#if symmetric target site, then create rc frag...
		 $test2c = substr($test2, 0 , $size);
	     $test2c =~ tr/ACTG/TGAC/;
         $test2rc = reverse $test2c;
	}
    if((substr($test1, 0, $size) eq substr($test2, 0, $size))&&(length(substr($test2, 0, $size))) > 10){
      print $test1 . "=" . substr($test2, $size+1),"\n";
	  #if ($j % 10 == 0){print "$j-"}
      if(substr($test1, $size+10, 6) ne substr($test2, $size+10, 6)){
        push(@cnv, $test1 . " = " . substr($test2, $size+1));
        print PAIRS "CNV:\n";
		#print PAIRS $test1 . " = " . substr($test2, $size+1) . "\n";
        #print CNV $test1 . " = " . substr($test2, $size+1) . "\n";
      }
      push(@dups, $test1 . " = " . substr($test2, $size+1));
      print PAIRS $test1 . " = " . substr($test2, $size+1) . "\n";
      splice (@infile1,$i,1);
      splice (@infile2,$j,1);
      $i = $i-1;
      $j = $j-1;
    }elsif((substr($test1, 0, $size) eq substr($test2rc, 0, $size))&&(length(substr($test2, 0, $size))) > 10){
      print $test1 . "=" . substr($test2, $size+1),"\n";
	  #if ($j % 10 == 0){print "$j."}
      if(substr($test1, $size+10, 6) ne substr($test2, $size+10, 6)){
        push(@cnv, $test1 . " = " . substr($test2, $size+1));
        print PAIRS "CNV:\n";
		#print PAIRS $test1 . " = " . substr($test2, $size+1) . "\n";
        #print CNV $test1 . " = " . substr($test2, $size+1) . "\n";
      }
      push(@dups, $test1 . " = " . substr($test2, $size+1));
      print PAIRS $test1 . " = " . substr($test2, $size+1) . "\n";
      splice (@infile1,$i,1);
      splice (@infile2,$j,1);
      $i = $i-1;
      $j = $j-1;
    }
  }
}

print "\n\nListing unique fragments in $outfile1 and $outfile4 ...\n\n";

for ($i = 0; $i < $#infile1+1; $i++){
  push (@nomatch1, $infile1[$i]);
}

for ($j = 0; $j < $#infile2+1; $j++){
  push (@nomatch2, $infile2[$j]);
}

@nomatch1 = sort (@nomatch1);
@nomatch2 = sort (@nomatch2);

for ($k = 0; $k < $#nomatch1+1; $k++){
  print NOMATCH1 "$nomatch1[$k]\n";
}

for ($k = 0; $k < $#nomatch2+1; $k++){
  print NOMATCH2 "$nomatch2[$k]\n";
}

print "\n ",$#dups+1," identical fragments including ",$#cnv+1," with copy number variation \n";
print " ",$i," unique fragments in ", substr($infilename1, 0, 13),"\n";
print " ",$j," unique fragments in ", substr($infilename2, 0, 13),"\n";
print "\n (", $#dups+1," x 2) + $i + $j = ", $#infile2+$#infile1+2," total fragments\n";

#end of identical pairings and CNV discovery
#*******************************************
#begin variant discovery

$polys = 5;
$loops1 = 0;
$family = 1;
close NOMATCH1;
open (NOMATCH1, "<$outfile1");

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$VarTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 

@freq = (0, 0, 0, 0, 0, 0);
if ($polys > 0){
   $allvars = 0;
   print PAIRS "\n\n";
   print "\n\n","finding variants with $polys or less mismatches...\n\n";
   while (<NOMATCH1>) {
      chomp $_;
      if ($loops1 % 200 == 0){
        print "$loops1\n";
      }
      $test = $_;
      $variants = 1;
      $weakhits = 0;
      $loops2 = 0;
      $loops3 = $#CspCIt+1;
      $d = -1;

      for ($i = 0; $i < $#nomatch2+1; $i++){
         $test2 = $nomatch2[$i];
         $mismatch = -1;
		 if ($sym == 1){#if symmetric target site, then create rc frag...
		 	$test2c = substr($test2, 0 , $size);
	        $test2c =~ tr/ACTG/TGAC/;
            $test2rc = reverse $test2c;
		 }
		 
         #if (substr($test, 0, $size) eq substr($test2, 0 , $size)){
         #   $mismatch = 0;
         #   splice (@nomatch2, $i, 1);
         #   $i = $i-1;
         #}

		 if (substr($test, 0, $size) eq substr($test2, 0 , $size)){
            $mismatch = 0;
            splice (@nomatch2, $i, 1);
            $i = $i-1;
		 }elsif (substr($test, 0, $size) eq substr($test2rc, 0 , $size)){
			$mismatch = 0;
			splice (@nomatch2, $i, 1);
			$i = $i-1;
	     }
		 
         #if ($mismatch < 0){
         #   $k = 0;
         #   $mismatch = 0;
         #   while (($k < $size)&&($mismatch <= $polys)){
         #      if(substr($test, $k, 1) ne substr($test2, $k, 1)){
         #         $mismatch++;
         #      }
         #      $k++;
         #   }   
         #}

		 if (($mismatch < 0)&&($sym == 1)){
            $k = 0;
            $mismatch = 0;
            while (($k < $size)&&($mismatch <= $polys)){
                if(substr($test, $k, 1) ne substr($test2rc, $k, 1)){
                  $mismatch++;
                }
               $k++;
			   #if ($k % 20 == 0){print "$k ", substr($test2, 0, $size), " $test2rc\n";}
            }
         }
		 unless (($mismatch >= 0)&&($mismatch <= $polys)){
            $k = 0;
            $mismatch = 0;
            while (($k < $size)&&($mismatch <= $polys)){
                if(substr($test, $k, 1) ne substr($test2, $k, 1)){
                  $mismatch++;
                }
               $k++;
			   #if ($k % 12 == 0){print $k,"-";}
            }
         }
		 
         if(($mismatch <= $polys)&&($mismatch > 0)){
		    $freq[$mismatch-1] = $freq[$mismatch-1]+1;
            print $test2," - $family.$variants.$mismatch\n";
			print PAIRS $test2," - $family.$variants.$mismatch\n";
            #print DUOMATCHV $test2," - $family.$variants.$mismatch\n";
            push (@duomatchv, $test2 . " - $family.$variants.$mismatch");
            $variants++;
            splice (@nomatch2, $i, 1);
            $i = $i-1;
			$allvars++;
			
         }

      }
      #end of "for" loop

      if ($variants > 1){
         print $test," - $family.$variants\n\n";
		 print PAIRS $test," - $family.$variants\n\n";
		 #print DUOMATCHV $test," - $family.$variants\n\n";
         push (@duomatchv, $test . " - $family.$variants");
         $family++;

      }
      $loops1++;
   }
   #end of first "while" loop
}else{
   print "\n no search for variants...";
}

#end of variant discovery
#***************************
#begin creating final output

print SUM "\nPhylogenomic Profile Pairwise Analysis: $version\n";
print SUM "\n ",$#dups+1," identical fragments including ",$#cnv+1," with copy number variation \n";
print SUM " ",$loops1," unique fragments in ", substr($infilename1, 0, 14),"\n";
print SUM " ",$j," unique fragments in ", substr($infilename2, 0, 14),"\n";
print SUM "\n (", $#dups+1," x 2) + $loops1 + $j = ", 2*($#dups+1)+$#infile2+$#infile1+2," total fragments\n\n";

print SUM $#duomatchv+1, " variants with $polys or less mismatches in ", $family-1, " families\n";
print SUM "   variant distributions (1-6):","@freq", "\n\n";

print "\nPhylogenomic Profile Pairwise Analysis: $version\n";
print "\n ",$#dups+1," identical fragments including ",$#cnv+1," with copy number variation \n";
print " ",$loops1," unique fragments in ", substr($infilename1, 0, 14),"\n";
print " ",$j," unique fragments in ", substr($infilename2, 0, 14),"\n";
print "\n (", $#dups+1," x 2) + $loops1 + $j = ", 2*($#dups+1)+$#infile2+$#infile1+2," total fragments\n\n";

print $#duomatchv+1, " variants with $polys or less mismatches in ", $family-1, " families\n";
print "   variant distributions (1-6): ","@freq", "\n\n";
print "fragment size = $size\n\n";

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$DoneTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 

print SUM "Start $StarTime\n";
print SUM "Var   $VarTime\n";
print SUM "Done  $DoneTime\n";
print "Start $StarTime\n";
print "Var   $VarTime\n";
print "Done  $DoneTime\n"; 


close SUM;
close INFILE1;
close INFILE2;
close NOMATCH1;
close NOMATCH2;
#close CNV;
close PAIRS;


#_______________________________________
# leaves window open for double-click starts
print "\nEND - press return to exit.";
$input = <>; 

