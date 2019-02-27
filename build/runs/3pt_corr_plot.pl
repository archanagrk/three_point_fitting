#!/usr/bin/perl

use POSIX qw(ceil floor);

die "3pt_corr_fns_plot.pl <ini.file> \n" unless $#ARGV >= 0;

print "three point corr fn plots \n";


$random = int(rand(1000));

open(OUT, "> ./3pt_corr_fns_plot_${random}.gnu");

print OUT "set lmargin 5\nset rmargin 0\n";
print OUT "set 3pt_corr_fns_plot\n";
print OUT "unset key\n\n\n";


$filename = shift(@ARGV);

    
($text, $xlow, $xhigh, $ylow, $yhigh, $chisq) = &make_gnu_from_plot($filename);

print OUT "set yrange[$ylow:$yhigh]\n";
print OUT "set xrange[$xlow:$xhigh]\n";
    
$xpos = 0.25*$xhigh ;
$ypos = $yhigh - 0.15 * abs($yhigh - $ylow);
print OUT "set label 1 \"chisq=$chisq\" at $xpos,$ypos \n";
$ypos = $yhigh - 0.30 * abs($yhigh - $ylow);

print OUT "$text\n\n";



print OUT "set no3pt_corr_fns_plot\n";

print OUT "pause -1\n";

# active this to keep a postscript file - doesn't really work
#    print OUT "set term postscript enhanced solid color eps\n";
#    print OUT "set out \"princorr.ps\"\n";
#    print OUT "replot\n";


close(OUT);

#make the window a sensible size
#$read = `xdpyinfo  | grep \'dimensions:\'`; chomp $read;
#($a, $dims, $a, $a,$a) = split(' ', $read);
#($width, $height) = split('x', $dims);
#$width = int( 0.9* $width);
#$height = int( 0.9* $height);

#system("gnuplot -geometry ${width}x${height} -persist /tmp/multiplot_${random}.gnu");
system("gnuplot -geometry 1500x800 -persist ./3pt_corr_fns_plot_${random}.gnu");

system("rm ./3pt_corr_fns_plot_${random}.gnu");

exit(0);


sub make_gnu_from_plot{
    local($filename) = @_ ;

    my $f = "\'${filename}\'";

    #below fit region
    my $text = "plot $f index 0 using 1:2:3 with lines ls 3, \\\n";
    $text = $text . "$f index 1 using 1:2:3 with lines ls 3, \\\n";
    $text = $text . "$f index 2 using 1:2:3 with lines ls 3, \\\n";

    #in fit region
    $text = $text . "$f index 3 using 1:2:3 with lines ls 1, \\\n";
    $text = $text . "$f index 4 using 1:2:3 with lines ls 1, \\\n";
    $text = $text . "$f index 5 using 1:2:3 with lines ls 1, \\\n";

    #above fit region
    $text = $text . "$f index 6 using 1:2:3 with lines ls 3, \\\n";
    $text = $text . "$f index 7 using 1:2:3 with lines ls 3, \\\n";
    $text = $text . "$f index 8 using 1:2:3 with lines ls 3, \\\n";

    #data
    #in fit region
    $text = $text . "$f index 9 using 1:2:3:4 with yerr ls 6, \\\n";
    #out of fit region
    $text = $text . "$f index 10 using 1:2:3:4 with yerr ls 3\n";

    #label text
    my $line = `cat $filename | grep 'chisq'`;

    my ($junk, $junk,$tmp) = split('=', $chisq);
    $chisq = $tmp;

    #get range
    $line =  `cat $filename | grep 'x '`;
    my ($a, $xlow, $xhigh) = split(' ', $line);

    $line = `cat $filename | grep 'y'`;
    my ($a, $ylow, $yhigh) = split(' ', $line);

    
    @out = ($text, $xlow, $xhigh, $ylow, $yhigh, $chisq);

    return @out;
}
