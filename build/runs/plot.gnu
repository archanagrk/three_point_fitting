gnuplot
set table 'lower.dat'
plot 't2_test32_three_pt_exp_fit.plot' using 2:($3-$4) smooth cspline
set table 'upper.dat'
plot 't2_test32_three_pt_exp_fit.plot' using 2:($3+$4) smooth cspline
unset table
set terminal pngcairo
set output 'data.png'

set style fill transparent solid 0.2 noborder
plot '< paste lower.dat upper.dat' using 1:2:5 with filledcurves title '95% confidence', \
     't2_test32_three_pt_exp_fit.plot' using 2:3 with lines lt 1 smooth cspline title 'mean value',\
     '' using 2:3 with points lt 1 pt 7 ps 1.5 lw 3 title 'data points'
exit
