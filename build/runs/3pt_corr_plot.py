from numpy import *
import numpy as np


import math
import sys
import re
import os
from os import system
import Gnuplot, Gnuplot.funcutils
import subprocess
import os, glob
from scipy import *

if len(sys.argv)!=2:
    print "3pt_corr_fns_plot: <ini_file>"
    sys.exit(1)

EXT = sys.argv[1]

EXT='*.x_y'
FILELIST=glob.glob(EXT)
SIZEX = len(io.array_import.read_array(FILELIST[0]))
DATAMATRIX = zeros((SIZEX,len(FILELIST)), Float)
TWOTHETA=io.array_import.read_array(FILELIST[0])[:,0]
TIMESTEP=150

for y in range(len(FILELIST)):
		DATAMATRIX[:,y]=sqrt(io.array_import.read_array(FILELIST[y])[:,1])

file = open("3ddata.dat", "w")

for y in range(len(FILELIST)):
	for x in range(1126,1968):
		file.write(repr(TIMESTEP*y)+" "\
		+repr(TWOTHETA[x])+" "+repr(DATAMATRIX[x,y]))
		file.write("\n")
	file.write("\n")
file.close()

f=os.popen('gnuplot' ,'w')
print >>f, "set ticslevel 0.0 ;set xlabel 'Time [s]'; set ylabel 'Diffraction angle'"
print >>f, "set pm3d; unset surface; set view 60,75; splot '3ddata.dat' notitle"
print >>f, "set terminal png large transparent size 600,400; set out '3ddata_1.png'"
print >>f, "replot"
f.flush()
