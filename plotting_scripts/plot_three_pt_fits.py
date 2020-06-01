#!/usr/local/bin/python3

import sys
import numpy as np
import math 

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import os
import glob

import matplotlib.cm as cm

pi = 3.141592653589

#------------------------------------------------------------------------------------#

def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
    return i + 1

#------------------------------------------------------------------------------------#
def plot(fit_plot_file, ax):

    n_lines = file_len(fit_plot_file)
    print("reading ", fit_plot_file, ", found n =", n_lines," lines")

    # ax.set(ylabel=r'$F$')
    # ax.set(xlabel='t')
    # ax.labelsize : 4

    f = open(fit_plot_file, 'r')

    data_in = f.read()

    ff_qsq_data = data_in.split("\n\n")

    for idx,dat in enumerate(ff_qsq_data, 1):
        if idx == 1:
            name    = dat.split("\n")[0].rstrip('\n').split()[2]
            frmf = float(dat.split("\n")[2].rstrip('\n').split()[2])
            frmf_err = float(dat.split("\n")[2].rstrip('\n').split()[4])
            chisq   = float(dat.split("\n")[3].rstrip('\n').split()[8])

        if idx == 2:
            active = dat

        elif idx == 3:
            inactive = dat

        elif idx == 4:
            ensem = dat
        else: continue


        plotLabel = r"$F$ = " + ("{:7.5f}".format(frmf)) + "(" + ("{:7.5f}".format(frmf_err)) + ")" + '\n' + r"$chisq/dof$ = " + ("{:6.4f}".format(chisq))


    active_dat = active.split("\n")


    for elem in active_dat:

        this_line = (elem.rstrip('\n')).split()

        if ( len(this_line) == 0):
            continue
        if (this_line[0]=="#"):
            #commented out lines begin with "# "
            continue

        t = float(this_line[1])
        ff  = float(this_line[2])
        err = float(this_line[3])

        ax.errorbar( t, ff, yerr=err, markersize=1.5,fmt='.',color='black', mfc='black',mec='black',elinewidth=0.5, capsize=1, mew=0.7,zorder =10, label='_nolegend_')

    inactive_dat = inactive.split("\n")

    for elem in inactive_dat:

        this_line = (elem.rstrip('\n')).split()

        if ( len(this_line) == 0 ):
            continue
        if (this_line[0]=="#"):
            #commented out lines begin with "# "
            continue

        t = float(this_line[1])
        ff  = float(this_line[2])
        err = float(this_line[3])

        ax.errorbar( t, ff, yerr=err, markersize=1.5,fmt='.',color='darkgreen', mfc='darkgreen',mec='darkgreen',elinewidth=0.5, capsize=1, mew=0.7,zorder =5, label='_nolegend_')


    ensem_dat = ensem.split("\n")

    t = []
    ff  = []
    formfac = []
    ffu = []
    ffb = []

    for elem in ensem_dat:

    
        this_line = (elem.rstrip('\n')).split()

        if ( len(this_line) == 0 ):
            continue
        if (this_line[0]=="#"):
            #commented out lines begin with "# "
            continue

        t.append(float(this_line[1]))
        ff.append(float(this_line[3]))
        ffb.append(float(this_line[2]))
        ffu.append(float(this_line[4]))
        formfac.append(frmf)



    ax.plot(t, ff, color='steelblue', lw = 1,zorder =5, label='_nolegend_')
    ax.fill_between(t, ffb, ffu, color='steelblue', alpha = 0.5,zorder =0 )

    ax.plot(t, formfac, color='grey', lw = 1, label= plotLabel)
    ax.fill_between(t, frmf - frmf_err, frmf + frmf_err, color='grey', alpha = 0.3,zorder =0 )

    ax.set_xlim(min(t), max(t))
    # ax.set_aspect(aspect=1/ax.get_data_ratio())
    ax.legend(loc='best',fontsize='xx-small') 
    ax.set_title('{}'.format(name), fontsize='xx-small')
    ax.tick_params(axis='both', which='major', labelsize='xx-small')



#------------------------------------------------------------------------------------#
# ----- begin main program ----- #

input = len(sys.argv)

if (input<2):
    print("usage: ", str(sys.argv[0]), "<Q2_xxx> ..")
    exit(1)

qsq_files = []

if(str(sys.argv[1]) == "all"): qsq_files = glob.glob(os.getcwd()+'/Q2_*')
    
else:
    for qsq in range(input-1):
        qsq_files.append(os.getcwd()+'/'+str(sys.argv[qsq+1]))
    

for qsq_file in qsq_files:

    print("plotting ", qsq_file)

    if os.path.exists(qsq_file):

        os.chdir(qsq_file)

        num = len(glob.glob(os.getcwd()+'/*/Real/*.plot'))
        if num == 0: 
            print("could not find a plot file")
            continue
        row = int(math.sqrt(num))
        col = int(math.ceil(float(num)/float(row)))

        count = 0
        fig = plt.figure(figsize=(20, 10))
        #fig = plt.figure()
        outer = gridspec.GridSpec(row, col, wspace=0.1, hspace=0.2)

        for plot_file in glob.glob(os.getcwd()+'/*/Real/*.plot'):

            inner = gridspec.GridSpecFromSubplotSpec(3, 3,
                    subplot_spec=outer[count], wspace=0.0, hspace=0.0)

            ax = plt.Subplot(fig, inner[:,:])
            plot(plot_file, ax)
            fig.add_subplot(ax)

            count += 1

        fig.patch.set_facecolor('white')

	plt.suptitle(os.path.basename(os.getcwd()))
        plt.show()
        os.chdir(os.getcwd()+'/../')        
    
    else: print(qsq_file, " does not exist")



exit(0)



# End of plotter
