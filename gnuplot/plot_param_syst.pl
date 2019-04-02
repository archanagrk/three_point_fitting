#!/usr/bin/perl

use List::Util qw[min max];

die "plot_three_point_fit.pl <filename_header>\n" unless $#ARGV == 0;
# should probably also have options to use pdf and save the file ...

$filename = $ARGV[0];
$plotfilename = "${filename}.plot";

# make a gnuplot input file
$random = int(rand(1000));
$tmpfile = "/tmp/plot_three_point_fit.${random}.gnu";

open(OUT, "> $tmpfile");

#print OUT "set term aqua noenhanced\n";  #add options for x11, or pdf


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

$mmin = 10;
$mmax = -10;

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
