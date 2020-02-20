#!/bin/bash
while IFS='' read -r line || [[ -n "$line" ]]; do
    #echo "calcbc ' real( $line )' | axis -e -c \\sq | xplot" | bash -
     echo "$line"
     end='_three_pt_fit' 
     fit="$line$end"
     echo $fit
     echo "cd $line; ../../../../gnuplot/plot_timeslice_three_pt_fit.pl $fit; cd .." | bash -
done < "$1"
