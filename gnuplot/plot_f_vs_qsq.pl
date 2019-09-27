#!/usr/bin/perl

use List::Util qw[min max];
use POSIX qw(ceil floor);

die "plot_f_vs_qsq.pl <filename>\n" unless $#ARGV == 0;


$plotfilename = $ARGV[0];


# make a gnuplot input file
$random = int(rand(1000));
$tmpfile = "/tmp/plot_f_vs_qsq.${random}.gnu";

$ymin = 0.01;
$ymax = -0.01;

$xmin = 0.01;
$xmax = -0.01;

my @irreps;

open(FILE , " < ${plotfilename}");

#now into active data
while($line = <FILE>){

  if($line =~ m/## irrep=/){($tmp, $irrep) = split('=', $line); push (@irreps, $irrep);}
  if($line eq /^#/){next;}
  if($line eq "\n"){next;}

  chomp $line;
  ($x, $xerr, $y, $yerr) = split(' ', $line);
  
  if( ($y + $yerr) > $ymax){ $ymax = $y + $yerr; }
  if( ($y - $yerr) < $ymin){ $ymin = $y - $yerr; }

  if( ($x + $xerr) > $xmax){ $xmax = $x + $xerr; }
  if( ($x - $xerr) < $xmin){ $xmin = $x - $xerr; }

}
close(FILE);

$ymax += 0.05*abs($ymax - $ymin);
$ymin -= 0.05*abs($ymax - $ymin);

$xmax += 0.05*abs($xmax - $xmin);
$xmin -= 0.05*abs($xmax - $xmin);


## loop over indeces
$end = scalar(@irreps);
$i = 0;


#foreach my $i(0 .. $#irreps) {

  open(OUT, "> $tmpfile");
  print OUT "set term aqua noenhanced\n";  #add options for x11, or pdf
  print OUT "set border 3\n";

  # print OUT "set xrange[${xmin}:${xmax}]\n set yrange[-10:10]\n ";
  print OUT "set key default\n";
  
  ## do the plot 
  $irrep = $irreps[$i];
  print OUT "plot \'${plotfilename}\' index $i u 1:3:4  w yerr  lc rgb 'black' pt 70 title \"${irrep}\" \\\n";

  close(OUT);

  #system("gnuplot -geometry 1500x800 -persist /tmp/${tmp_file}");
  system("gnuplot $tmpfile");

  #clean up
  system("rm $tmpfile");


  #$i++;

#}

exit(0);