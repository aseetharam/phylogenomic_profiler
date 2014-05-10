#!/usr/bin/perl
# 
# Copyright 2013, Gary W. Stuart <gstuart@indstate.edu>
# Rename this file as phyper.pl
# Subsample Phylogenomic Profiler: Phyper3
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
# Simulate restriction digestion for all 16 Type IIB endonuclease enzymes
# 
#
# By Gary W. Stuart <gstuart@indstate.edu>
#

$version = "Phyper3";
print "\nPhylogenomic Profiler: $version\n\n";

#default number of mismatches allowed during variant analysis...
$dpolys = 6;
#default enzyme used for profile...
$dEnz = "Cjp";
$sym = 0;

#input file name for analysis...
print "Enter fasta filename (i.e. GenSpeSub.fasta): ";
$infile = <>;
chomp $infile;
$name = "<$infile";
$GenSpeSub = substr($infile,0,9);

#input typeIIB enzyme name for analysis...
print "\nEnter typeIIB enzyme (default = $dEnz): ";

$Enz = <>;
chomp $Enz;

if (length($Enz) < 1){#Check whether Enz was left undefined...
	$Enz = $dEnz;
}
if ($Enz eq "Csp"){#CspCI (11/13)CAANNNNNGTGG(12/10) 8192 33
	$L = 13;$LS = "CAA";$M = 5;$RS = "GTGG";$R = 12; $size = 37;
}elsif ($Enz eq "Alo"){#AloI (7/12)GAACNNNNNNTCC(12/7) 8192 27
	$L = 12;$LS = "GAAC";$M = 6;$RS = "TCC";$R = 12; $size = 37;
}elsif ($Enz eq "Ppi"){#PpiI (7/12)GAACNNNNNCTC(13/8) 8192 27 
	$L = 12;$LS = "GAAC";$M = 5;$RS = "CTC";$R = 13; $size = 37;
}elsif ($Enz eq "Psr"){#PsrI (7/12)GAACNNNNNNTAC(12/7) 8192 27 
	$L = 12;$LS = "GAAC";$M = 6;$RS = "TAC";$R = 12; $size = 37;
}elsif ($Enz eq "Bpl"){#BplI (8/13)GAGNNNNNCTC(13/8) 4096 27 
	$L = 13;$LS = "GAG";$M = 5;$RS = "CTC";$R = 13; $size = 37; $sym = 1;
}elsif ($Enz eq "Fal"){#FalI (8/13)AAGNNNNNCTT(13/8) 4096 27 
	$L = 13;$LS = "AAG";$M = 5;$RS = "CTT";$R = 13; $size = 37; $sym = 1;
}elsif ($Enz eq "Bsp"){#Bsp24I (8/13)GACNNNNNNTGG(12/7) 2048 27 
	$L = 13;$LS = "GAC";$M = 6;$RS = "TGG";$R = 12; $size = 37;
}elsif ($Enz eq "Bsa"){#BsaXI (9/12)ACNNNNNCTCC(10/7) 2048 27 
	$L = 12;$LS = "AC";$M = 5;$RS = "CTCC";$R = 10; $size = 33;
}elsif ($Enz eq "Hae"){#HaeIV (7/13)GAYNNNNNRTC(14/9) 1024 27 
	$L = 13;$LS = "GA[CT]";$M = 5;$RS = "[AG]TC";$R = 14; $size = 38; $sym = 1;
}elsif ($Enz eq "Cje"){#CjeI (8/14)CCANNNNNNGT(15/9) 512 28 
	$L = 14;$LS = "CCA";$M = 6;$RS = "GT";$R = 15; $size = 40;
}elsif ($Enz eq "Cjp"){#CjePI (7/13)CCANNNNNNNTC(14/8) 512 27 
	$L = 13;$LS = "CCA";$M = 7;$RS = "TC";$R = 14; $size = 39;
}elsif ($Enz eq "Hin"){#Hin4I (8/13)GAYNNNNNVTC(13/8) 512 27 
	$L = 13;$LS = "GA[CT]";$M = 5;$RS = "[CAG]TC";$R = 13; $size = 37; $sym = 1;
}elsif ($Enz eq "Bae"){#BaeI (10/15)ACNNNNGTAYC(12/7) 4096 28 
	$L = 15;$LS = "AC";$M = 4;$RS = "GTA[CT]C";$R = 12; $size = 38;
}elsif ($Enz eq "Alf"){#AlfI (10/12)GCANNNNNNTGC(12/10) 4096 32 
	$L = 12;$LS = "GCA";$M = 6;$RS = "TGC";$R = 12; $size = 36; $sym = 1;
}elsif ($Enz eq "Bcg"){#BcgI (10/12)CGANNNNNNTGC(12/10) 2048 32 
	$L = 12;$LS = "CGA";$M = 6;$RS = "TGC";$R = 12; $size = 36;
}elsif ($Enz eq "Bsl"){#BslFI (6/10)GGGAC(10/14) 512 21 
	$L = 10;$LS = "GGG";$M = 0;$RS = "AC";$R = 14; $size = 29;
}elsif ($Enz ne "Csp"){#CspCI (11/13)CAANNNNNGTGG(12/10) 8192 33
	$L = 13;$LS = "CAA";$M = 5;$RS = "GTGG";$R = 12; $size = 37;
	print "Can't recognize enzyme, using Csp instead...";
	$Enz = "Csp";
}
$cRS = $RS; $cLS = $LS;
$cRS =~ tr/ACTG/TGAC/;
$cLS =~ tr/ACTG/TGAC/;
$rcRS = reverse $cRS;
$rcLS = reverse $cLS;
$rcRSHin = "GAG";

#Enz ($L,$LS,$M,$RS,$R) ($R,$rcRS,$M,$rcLS,$L)
#Hin4I enzyme requires special non-symmetric site search using $rcRSHin due to "V" reduncancy

#input maximum mismatches allowed for fragment matching...
print "\nEnter maximum mismatches allowed for fragment matching (1-9, default = $dpolys): ";
$polys = <>;
chomp $polys;
unless ($polys >= 1 && $polys <= 9){#if out of range, substitute default...
	$polys = $dpolys;
}
print "\nEnzyme = $Enz & mismatches = $polys\n";

#end of enzyme and input file identification
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#begin data preprocessing

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$StartTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 
print "\n$StartTime\n"; 

#remove line feeds from fasta file inputs..
print "\n processing file...\n\n";
open (INPUT, $name);
open (WORK, ">>rdelw.txt");
$chunks = 0;
$lines = 0;
while (<INPUT>) {
 if ((substr($_, 0, 1) eq ">")&&($lines == 0)){
    $_ = ">";
    $chunks++;
    }elsif (substr($_, 0, 1) eq ">"){
    $_ = "\n>";
    $chunks++;
 }
 chomp $_;
 $lines++;
    #if ($lines % 1000 == 0){
    #   print " ",$lines/1000,"k ";
    #}
print WORK $_;
}
close INPUT;
close WORK;

print " Loaded ",$chunks," clones/chunks/contigs";
print "\n   from ",$lines," lines of sequence...\n";

#end of data preprocessing
#:::::::::::::::::::::::::::::::::::::::::::
#begin "while" loop for fragment extraction

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$FragTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 
print "\n$FragTime\n";

$outfile = $GenSpeSub . "$Enz" . "Ux.txt";
$outfile1 = $GenSpeSub . "$Enz" . "Tx.txt";
$outfile2 = $GenSpeSub . "$Enz" . "V$polys.txt";
$outfile3 = $GenSpeSub . "$Enz" . "F$polys.txt";
$outfile4 = $GenSpeSub . "$Enz" . "Ax.txt";
$outfile5 = $GenSpeSub . "$Enz" . "R$polys.txt";
$outfile6 = $GenSpeSub . "$Enz" . "S$polys.txt";
$outfile7 = $GenSpeSub . "$Enz" . "Sum$polys.txt";

open (WORK, "<rdelw.txt");
open (MATCHU, ">$outfile");
open (MATCHT, ">$outfile1");
open (MATCHV, ">$outfile2");
open (MATCHF, ">$outfile3");
open (MATCHA, ">$outfile4");
open (MATCHR, ">$outfile5");
open (MATCHS, ">$outfile6");
open (SUM, ">$outfile7");

print "\n Loading Match-Work file...\n";

$chunknum = 0;
$wcount = 0; $crccount = 0; $wcount2 = 0; $crccount2 = 0;

#find matching patterns in data file and copy to "all" array...
while (<WORK>) {
  $chunknum++;
  print "\n Finding $Enz fragments in chunk $chunknum: ";
  @CspCIw = $_ =~ m/[ACGT]{$L}$LS[ACGT]{$M}$RS[ACGT]{$R}/g;
  @CspCIw2 = $_ =~ m/[ACGT]$LS[ACGT]{$M}$RS[ACGT]{0,$R - 1}$LS[ACGT]{$M}$RS[ACGT]{$R}/g;
  print $#CspCIw+1,"+",$#CspCIw2+1," $Enz fragments, ";
    if ($Enz eq "Hin"){#test for Hin4I site...
	@CspCIc = $_ =~ m/[ACGT]{$R}$rcRSHin[ACGT]{$M}$rcLS[ACGT]{$L}/g;
	@CspCIc2 = $_ =~ m/$rcRS[ACGT]{$M}$rcLS[ACGT]{0,$L - 1}$rcRS[ACGT]{$M}$rcLS[ACGT]{$L}/g;
	print $#CspCIc+1,"+",$#CspCIc2+1," rc-Match fragments...\n";
    }elsif ($LS ne $rcRS){#test for symmetric target site...
	@CspCIc = $_ =~ m/[ACGT]{$R}$rcRS[ACGT]{$M}$rcLS[ACGT]{$L}/g;
	@CspCIc2 = $_ =~ m/$rcRS[ACGT]{$M}$rcLS[ACGT]{0,$L - 1}$rcRS[ACGT]{$M}$rcLS[ACGT]{$L}/g;
	print $#CspCIc+1,"+",$#CspCIc2+1," rc-Match fragments...\n";
	}
  for ($i = 0; $i < $#CspCIc+1; $i++){
    if (length ($CspCIc[$i]) == $size){
      $CspCIc[$i] =~ tr/ACTG/TGAC/;
      $CspCIcrc[$i] = reverse $CspCIc[$i];
    }
  }

  for ($i = 0; $i < $#CspCIc2+1; $i++){
    if (length ($CspCIc2[$i]) == $size){
       $CspCIc2[$i] =~ tr/ACTG/TGAC/;
       $CspCIcrc2[$i] = reverse $CspCIc2[$i];
       #delete $CspCIc2[$i];
    }
  }

  for ($i = 0; $i < $#CspCIw+1; $i++){
    if (length ($CspCIw[$i]) == $size){
      push (@CspCIa, $CspCIw[$i]);
      print MATCHA "$CspCIw[$i]\n";
      $wcount++;
    }
  }

  print "  one strand done...";

  for ($i = 0; $i < $#CspCIcrc+1; $i++){
    if (length ($CspCIcrc[$i]) == $size){
      push (@CspCIa, $CspCIcrc[$i]);
      print MATCHA "$CspCIcrc[$i]\n";
      $crccount++;
    }
  }

#must clear "crc" arrays using splice to reset counts!!
  splice (@CspCIcrc, 0);
  splice (@CspCIcrc2, 0);

  print "other strand done...";

  for ($i = 0; $i < $#CspCIw2+1; $i++){
    if (length ($CspCIw2[$i]) == $size){
      $CspCIw2[$i] = substr($CspCIw2[$i],$size);
      push (@CspCIa, $CspCIw2[$i]);
      print MATCHA "$CspCIw2[$i]\n";
      $wcount2++;
    }
  }

  for ($i = 0; $i < $#CspCIcrc2+1; $i++){
    if (length ($CspCIcrc2[$i]) == $size){
      $CspCIcrc2[$i] = substr($CspCIcrc2[$i],$size);
      push (@CspCIa, $CspCIcrc2[$i]);
      print MATCHA "$CspCIcrc2[$i]\n";
      $crccount2++;
    }
  }
  print "overlaps resolved...\n";
}
print "\n";
print $wcount,"+",$wcount2," (+)strand matches\n",'N' x $L,$LS,'N' x $M,$RS,'N' x $R,"\n\n";
print $crccount,"+",$crccount2," (-)strand matches\n","N" x $L,$cLS,"N" x $M,$cRS,'N' x $R,"\n\n";
  
#end of "while" loop for fragment extraction
#:::::::::::::::::::::::::::::::::::::::::::
#begin copy number determination

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$CopyTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 
print "\n$CopyTime\n"; 
  
print " determining fragment copy numbers...\n\n";
$loops = $#CspCIa+1;
for ($i = 0; $i < $loops; $i++){
  if ($i % 200 == 0){print "$i\n";}
  $copies = 1;
  $test = pop(@CspCIa);
  for ($j = 0; $j < $#CspCIa+1; $j++){
    #print "?",$CspCIa[$j],"\n";
    if("$test" eq "$CspCIa[$j]"){
      $copies++;
      $CspCIa[$j] = "$CspCIa[$j]" . " $copies";
      #print $CspCIa[$j],"\n";
	}elsif ($sym > 0){#test for symmetric target site...
	  $CspCIc[$j] = $CspCIa[$j];
	  $CspCIc[$j] =~ tr/ACTG/TGAC/;
      $CspCIrc[$j] = reverse $CspCIc[$j];
	  if("$test" eq "$CspCIrc[$j]"){
        $copies++;
        $CspCIa[$j] = "$CspCIa[$j]" . " $copies";
	  }
    }
  }

  $test1 = "$test" . " $copies";
  if($copies != 1){
    #print $copies," ";

    print substr($test1, 0, $size) . " $GenSpeSub$Enz.$copies\n";
    #$dups{$test} = "$copies";
    push(@CspCId, substr($test1, 0, $size+1) . "$GenSpeSub$Enz.$copies");
    push(@CspCIt, substr($test1, 0, $size+1) . "$GenSpeSub$Enz.$copies");
  }elsif (length ($test) == $size){
    push(@CspCIu, substr($test1, 0, $size+1) . "$GenSpeSub$Enz.$copies");
    push(@CspCIt, substr($test1, 0, $size+1) . "$GenSpeSub$Enz.$copies");
  }
}

#create summary output with copy number info

@CspCIt = sort @CspCIt;
@CspCItt = @CspCIt;
@CspCIu = sort @CspCIu;

for ($i = 0; $i < $#CspCIt+1; $i++){
  #print $CspCIt[$i]," ",$i+1,"\n";
  print MATCHT $CspCIt[$i],"\n";
  #print MATCHX $CspCIt[$i],"\n";
}

close MATCHA;
close MATCHT;
#close MATCHX;

#end of copy number determination
#*********************************
#begin variant discovery

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$VarTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 
print "\n$VarTime\n"; 

open (MATCHT, "<$outfile1");
#open (MATCHX, "<$outfile5");

$loops1 = 0;
$family = 1;
if ($polys > 0){
   print "\n\n","finding variants with $polys or less mismatches...\n\n";
   while (<MATCHT>) {
      #print "$loops1:$loops2  ";
      chomp $_;
      if ($loops1 % 200 == 0){
        print "$loops1\n";
      }
      $test = $_;
      #print "$test,$test2\n";
      $variants = 1;
      $weakhits = 0;
      $loops2 = 0;
      $loops3 = $#CspCIt+1;
      $d = -1;

      for ($i = 0; $i < $#CspCIt+1; $i++){
         $test2 = $CspCIt[$i];
         #$test2 = pop(@CspCIt);
         #print "$test2\n";
         $mismatch = -1;
         $maxhits = $polys + 1; $tries = 0; $hits = 0;
		 if ($sym == 1){#if symmetric target site, then create rc frag...
		 	$test2c = substr($test2, 0 , $size);
	        $test2c =~ tr/ACTG/TGAC/;
            $test2rc = reverse $test2c;
		 }
		 if (substr($test, 0, $size) eq substr($test2, 0 , $size)){
            $mismatch = 0;
            splice (@CspCIt, $i, 1);
            $i = $i-1;
		 }elsif (substr($test, 0, $size) eq substr($test2rc, 0 , $size)){
			$mismatch = 0;
			splice (@CspCIt, $i, 1);
			$i = $i-1;
	     }
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
            $copies = substr($test2, $size+14, length($test2)+1);
            #$tcopies = $copies =~ m/\d{1,5}/;
            $test2 = substr($test2, 0, $size);
            print $test2," $GenSpeSub$Enz.$copies.$family.$variants.$mismatch\n";
            print MATCHV $test2," $GenSpeSub$Enz.$copies.$family.$variants.$mismatch\n";
            push (@CspCIv, $test2 . " $GenSpeSub$Enz.$copies.$family.$variants.$mismatch");
            $variants++;
            splice (@CspCIt, $i, 1);
            $i = $i-1;
         }

      }
      #end of "for" loop

      if ($variants > 1){
         #$variants++;
         $copies = substr($test, $size+14, length($test)+1);
         #$first = substr($first, 0, $size);
         $test = substr($test, 0, $size);
         print "$test $GenSpeSub$Enz.$copies.$family.$variants\n";
         print MATCHV "$test $GenSpeSub$Enz.$copies.$family.$variants\n";
         push (@CspCIv, "$test $GenSpeSub$Enz.$copies.$family.$variants");
         $family++;
         #splice (@CspCIt, $d, 1);
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

for ($i = 0; $i < $#CspCIu+1; $i++){
  $CspCIu[$i] = substr($CspCIu[$i], 0, $size);
  #print $CspCIu[$i]," ",$i+1,"\n";
  print MATCHU $CspCIu[$i]," $GenSpeSub$Enz.1","\n";
  $stop1 = 0;
  for ($j = 0; $j < $#CspCIv+1; $j++){
    if(substr($CspCIu[$i], 0, $size) eq substr($CspCIv[$j], 0, $size)){
       $stop1++;
    }
  }
  if ($stop1 == 0){
     push (@CspCIs,$CspCIu[$i]);
     print MATCHS $CspCIu[$i]," $GenSpeSub$Enz.1","\n";
  }
}

print "\nPhylogenomic Profiler: $version\n\n";


$rcount = 0;

for ($i = 0; $i < $#CspCItt+1; $i++){
  $skip = 0;
  for ($j = 0; $j < $#CspCIv+1; $j++){
      if(substr($CspCIv[$j], 0, $size) eq substr($CspCItt[$i], 0, $size)){
        #print $CspCIv[$j],"\n";
        #print substr($CspCIv[$j], 0, $size),"\n";
        push (@CspCIf,$CspCIv[$j]);
        print MATCHF $CspCIv[$j],"\n";
        $skip = 1;
      }
  }
  if ($skip == 0) {
    $copies = substr($CspCItt[$i], $size+14, length($CspCItt[$i])+1);
    #$tcopies = $copies =~ m/\d{1,5}/;
    $CspCItt[$i] = substr($CspCItt[$i], 0, $size);
    push (@CspCIf,$CspCItt[$i]);
    print MATCHF $CspCItt[$i]," $GenSpeSub$Enz.$copies","\n";
    print MATCHR $CspCItt[$i]," $GenSpeSub$Enz.$copies","\n";
    $rcount++
    }
}
print SUM "\nPhylogenomic Profiler: $version\n\n";

print SUM "\n",$#CspCIu+1,"/",$loops1,"/",$loops," (unique/distinct/total) fragments listed in...\n";
print SUM $outfile3,"\n";

print SUM "\n",$#CspCIs+1," singular fragments listed in...\n";
print SUM $outfile6,"\n";

print SUM "\n",$rcount," representative fragments listed in...\n";
print SUM $outfile5,"\n";

print SUM "\n",$#CspCIv+1," total variants in ",$family-1," family clusters listed in...\n";
print SUM $outfile2,"\n\n";


print "\n",$#CspCIu+1,"/",$loops1,"/",$loops," (unique/distinct/total) fragments listed in...\n";
print $outfile3,"\n";

print "\n",$#CspCIs+1," singular fragments listed in...\n";
print $outfile6,"\n";

print "\n",$rcount," representative fragments listed in...\n";
print $outfile5,"\n";

print "\n",$#CspCIv+1," total variants in ",$family-1," family clusters listed in...\n";
print $outfile2,"\n\n";

@months = qw(01 02 03 04 05 06 07 08 09 10 11 12);
#@weekDays = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($sec, $min, $hr, $date, $mon, $y, $day, $doy, $dsav) = localtime();
$year = 1900 + $y;
$months[$mon];
$DoneTime = "$year:$months[$mon]:$date:$hr:$min:$sec"; 
#$weekDays[$day];
print SUM "Start $StartTime\n";
print SUM "Frag  $FragTime\n";
print SUM "Copy  $CopyTime\n";
print SUM "Var   $VarTime\n";
print SUM "Done  $DoneTime\n"; 

print "Start $StartTime\n";
print "Frag  $FragTime\n";
print "Copy  $CopyTime\n";
print "Var   $VarTime\n";
print "Done  $DoneTime\n"; 

close WORK;
open (WORK, ">rdelw.txt");
print WORK "";
close WORK;
close MATCHU;
close MATCHT;
close MATCHV;
close MATCHF;
#close MATCHX;
close MATCHS;
close MATCHR;
close SUM;


#_______________________________________
# leaves window open for double-click starts
print "\nEND - press return to exit.";
$input = <>; 

