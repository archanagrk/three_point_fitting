#!/usr/bin/perl

use List::Util qw[min max];
use POSIX qw(ceil floor);

die "plot_three_point_fit.pl <filename_header>\n" unless $#ARGV == 0;


$filename = $ARGV[0];

$plotfilename = "${filename}";

# make a gnuplot input file
$random = int(rand(1000));
$tmpfile = "/tmp/plot_three_point_fit.${random}.gnu";

$ymin = 2.0;
$ymax = -2.0;

open(FILE , " < ${plotfilename}");
$corrname = <FILE>;
$firstline = <FILE>;
$F = <FILE>;
$chisq = <FILE>;
<FILE>; <FILE>; <FILE>; <FILE>;

#now into active data
while($line = <FILE>){

  if($line =~ m/# ensem fit/){last;}
  if($line =~ m/# inactive data/){last;}
  if($line =~ m/# active data/){next;}
  if($line =~ m/# F values/){last;}
  if($line eq /^#/){next;}
  if($line eq "\n"){next;}

  chomp $line;
  ($dt, $x, $y, $err) = split(' ', $line);
  
  if( ($y + $err) > $ymax){ $ymax = $y + $err; }
  if( ($y - $err) < $ymin){ $ymin = $y - $err; }

}
close(FILE);

$ymax += 0.005*abs($ymax - $ymin);
$ymin -= 0.005*abs($ymax - $ymin);

chomp $firstline; chomp $F; chomp $chisq; chomp $corrname;
$corrname = substr($corrname, 8);
($x, $x, $tmin, $x, $DT, $x, $nfits) = split(' ', $firstline);
$F = substr($F, 3);
$chisq = substr($chisq, 3);

my @Dt_val = split(',', $DT);

$Dt = @Dt_val;

open(OUT, "> $tmpfile");
# print OUT "set terminal pdf\n";  #add options for x11, or pdf
# print OUT "set output '.pdf'";

$end = $Dt - 1;


print OUT "set lmargin 5\nset rmargin 0\n";
print OUT "set multiplot\n";
print OUT "unset key\n\n\n";


$dim = ceil( sqrt($Dt) );
$s = 0.95/$dim;
$size = "$s, $s";

## time dep plot 
foreach my $i (0..$end) {

  $tmax = $Dt_val[$i];

  $x = $s*($i % $dim);
  $y = 1.0 - $s - $s * int($i / $dim);  

  print OUT "set origin $x,$y\n";
  print OUT "set size $size\n";

  print OUT "set xrange[0:${tmax}]\n set yrange[${ymin}:${ymax}]\n ";


  #$xpos = 0.25*$tmax ;
  $ypos1 = $ymax - 0.10 * abs($ymax - $ymin);
  $ypos2 = $ymax - 0.15 * abs($ymax - $ymin);
  $ypos3 = $ymax - 0.20 * abs($ymax - $ymin);

  print OUT "set label 1 \"${corrname}\" at ${tmin},${ypos1} noenhanced textcolor rgb \"#000000\" front \n";
  print OUT "set label 2 \"$F\" at ${tmin},${ypos2} noenhanced textcolor rgb \"#000000\" front \n";
  print OUT "set label 3 \"$chisq\" at ${tmin},${ypos3} noenhanced textcolor rgb \"#000000\" front \n";

  $idx_one   =  $i;
  $idx_two   =  $Dt + $i;
  $idx_three =  $Dt + $Dt + $i;
  $idx_four =  $Dt + $Dt + $Dt + $i;



  #print OUT "plot   \'${plotfilename}\' index $idx_three using 2:3:5 w filledcu fs solid 0.15 fc rgb \"#C0272D\",\\\n";
  print OUT "plot   \'${plotfilename}\' index $idx_three using 2:5 w lines lw 2 lc rgb \"#C0272D\",\\\n"; #don't fill in aquaterm to see better
  print OUT "       \'${plotfilename}\' index $idx_three using 2:3 w lines lw 2 lc rgb \"#C0272D\",\\\n"; #don't fill in aquaterm to see better
  print OUT "       \'${plotfilename}\' index $idx_three using 2:4 w lines lw 2 lc rgb \"#C0272D\",\\\n";
  print OUT "       \'${plotfilename}\' index $idx_four using 2:5 w lines lw 2 lc rgb \"#4B0082\",\\\n"; #don't fill in aquaterm to see better
  print OUT "       \'${plotfilename}\' index $idx_four using 2:3 w lines lw 2 lc rgb \"#4B0082\",\\\n"; #don't fill in aquaterm to see better
  print OUT "       \'${plotfilename}\' index $idx_four using 2:4 w lines lw 2 lc rgb \"#4B0082\",\\\n";
  print OUT "       \'${plotfilename}\' index $idx_one using 2:3:4 with yerr pt 70 lc rgb \"#000000\",\\\n";
  print OUT "       \'${plotfilename}\' index $idx_two using 2:3:4 with yerr pt 69 lc rgb \"#2F7A79\"\n";

}


print OUT "set nomultiplot\n";
print OUT "pause -1\n";

close(OUT);

system("gnuplot -geometry 1500x800 -persist $tmpfile");

#clean up
system("rm $tmpfile");


exit(0);
