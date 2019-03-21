#!/usr/bin/perl

use List::Util qw[min max];

die "plot_three_point_fit.pl <filename_header> <#Dt>\n" unless $#ARGV == 1;
# should probably also have options to use pdf and save the file ...

$filename = $ARGV[0];
$Dt       = $ARGV[1];
$plotfilename = "${filename}.plot";

# make a gnuplot input file
$random = int(rand(1000));
$tmpfile = "/tmp/plot_three_point_fit.${random}.gnu";

$ymin = 100.0;
$ymax = -100.0;

open(FILE , " < ${plotfilename}");
$firstline = <FILE>;
$F = <FILE>;
$chisq = <FILE>;
<FILE>; <FILE>; <FILE>;
#now into active data

while($line = <FILE>){
  if($line eq "\n"){last;}
  chomp $line;
  ($dt, $x, $y, $err) = split(' ', $line);
  
  if( ($y + $err) > $ymax){ $ymax = $y + $err; }
  if( ($y - $err) < $ymin){ $ymin = $y - $err; }
}
close(FILE);

$ymax += 0.05;
$ymin -= 0.05;

chomp $firstline; chomp $F; chomp $chisq;
($x, $x, $tmin, $x, $tmax, $x, $nfits) = split(' ', $firstline);
$F = substr($F, 3);
$chisq = substr($chisq, 3);

open(OUT, "> $tmpfile");
#print OUT "set term aqua noenhanced\n";  #add options for x11, or pdf
$end = $Dt - 1;

$size = $Dt*500;
print OUT "set term x11 noenhanced size 800,$size\n set multiplot\n";

## time dep plot 
foreach my $i (0..$end) {

  $top = 0.95+($i*0.85)/$Dt;
  $bottom = 0.1+($i*0.85)/$Dt;
  

  print OUT "set lmargin at screen 0.05\n set rmargin at screen 0.65\n set bmargin at screen $bottom\n set tmargin at screen $top\n";

  print OUT "unset label 1\n unset label 2\n unset label 3\n";
  print OUT "set border 3\n";

  print OUT "set xrange[0:${tmax}]\n set yrange[${ymin}:${ymax}]\n set size ratio 1.0\n unset key\n";
  print OUT "set xtics nomirror\n";
  print OUT "set ytics nomirror\n";

  #$xpos = max($tmin, $tmax - 10);
  $ypos1 = $ymax - 0.1*($ymax - $ymin) + ($i*0.85);
  $ypos2 = $ymin + 0.15*($ymax - $ymin) + ($i*0.85);
  $ypos3 = $ymin + 0.1*($ymax - $ymin) + ($i*0.85);

  print OUT "set label 1 \"${filename}\" at ${tmin},${ypos1} noenhanced textcolor rgb \"#000000\"\n";
  print OUT "set label 2 \"$F\" at ${tmin},${ypos2} noenhanced textcolor rgb \"#000000\"\n";
  print OUT "set label 3 \"$chisq\" at ${tmin},${ypos3} noenhanced textcolor rgb \"#000000\"\n";

  $idx_one   =  $i;
  $idx_two   =  $Dt + $i;
  $idx_three =  $Dt + $Dt + $i;



  print OUT "plot   \'${plotfilename}\' index $idx_three using 2:3:4 w filledcu fs solid 0.15 fc rgb \"#C0272D\",\\\n";
  print OUT "       \'${plotfilename}\' index $idx_three using 2:4 w lines lw 2 lc rgb \"#C0272D\",\\\n";
  print OUT "       \'${plotfilename}\' index $idx_one using 2:3:4 with yerr pt 70 lc rgb \"#000000\",\\\n";
  print OUT "       \'${plotfilename}\' index $idx_two using 2:3:4 with yerr pt 69 lc rgb \"#2F7A79\"\n";

}


################################################################################################################
## now the syst plot
$systfilename = "${filename}.syst";

open(FILE , " < ${systfilename}");
$firstline = <FILE>;
$lastline = $_ while <FILE>;
close(FILE);

chomp $firstline;
($a, $ensem, $chisq) = split(/\|/, $firstline);
($x, $mean, $err) = split(' ', $a);

chomp $lastline;
($ymax, $rest) = split(' ', $lastline); $ymax++;

$mean_plus = $mean + $err;
$mean_minus = $mean - $err;

# convert to comma separated
$datafile = "/tmp/plot_param_syst.${random}.dat";
open(DAT, "> $datafile");

$mmin = 1000.;
$mmax = -1000;

open(FILE , " < ${systfilename}");
$line = <FILE>; $line = <FILE>; $line = <FILE>; 
while($line = <FILE>){
  ($a, $fit, $chisq) = split(/\|/, $line);
  ($x, $z, $zz) = split(' ', $a);
  print DAT "${x},${z},${zz},${fit},${chisq}\n";
  if($z + $zz > $mmax){ $mmax = $z + $zz; }
  if($z - $zz < $mmin){ $mmin = $z - $zz; }
    
};
close(FILE);

$top = $mmax * 1.05;
$bot = 0.975 * $mmin;

print OUT "set lmargin at screen 0.65\n set rmargin at screen 0.95\n set bmargin at screen 0.15\n set tmargin at screen 0.95\n";

print OUT "unset label 1\n unset label 2\n unset label 3\n";
print OUT "set border 5\n";
print OUT "set xrange[${bot}:${top}]\n set yrange[0:${ymax}]\n set size ratio 2.2\n set grid xtics\n unset ytics\n unset key\n";
print OUT "set xtics rotate\n";
print OUT "set datafile separator \",\"\n";

print OUT "set obj 20 rect from ${mean_minus},0 to ${mean_plus},(${ymax} - 0.75) fs solid 0.15 noborder fc rgb \"#C0272D\" behind\n";
print OUT "set arrow from ${mean},0 to ${mean},(${ymax} - 0.75) nohead lc rgb \"#C0272D\" lw 2\n";
print OUT "set label \"${mean} +/- ${err}\" at ${mean},${ymax} center offset 0, char -1 noenhanced textcolor rgb \"#C0272D\"\n";

print OUT "plot   \'${datafile}\' index 0 using 2:1:3 with xerr lc rgb \"#2F7A79\",\\\n";
print OUT "       \'${datafile}\' index 0 using (${mmax} + 0*\$2):1:4 with labels left offset char 2, 0 ,\\\n";
print OUT "       \'${datafile}\' index 0 using (${mmax} + 0*\$2):1:5 with labels left offset 0, char -0.75 textcolor rgb \"#2F7A79\"\n";

print OUT "replot\n pause -1\n";

close(OUT);

system("gnuplot $tmpfile");
#system("gnuplot $tmpfile");

#clean up
system("rm $tmpfile $datafile");


exit(0);
