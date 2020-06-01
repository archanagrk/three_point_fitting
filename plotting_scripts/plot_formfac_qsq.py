#!/usr/local/bin/python3

import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import matplotlib.patches as mpatches

import matplotlib.cm as cm


pi = 3.141592653589

#------------------------------------------------------------------------------------#

def file_len(fname):
  with open(fname) as f:
    for i, l in enumerate(f):
      pass
    return i + 1

#------------------------------------------------------------------------------------#
def plot_ff_qsq(fit_qsq_plot_file):

    n_lines = file_len(fit_qsq_plot_file)
    print("reading ", fit_qsq_plot_file, ", found n =", n_lines," lines")

    plt.figure()
    plt.ylabel(r'$F(Q^{2})$')
    plt.xlabel(r'$(a_tQ)^{2}$')

    f = open(fit_qsq_plot_file, 'r')

    data_in = f.read()

    ff_qsq_data = data_in.split("\n\n")

    for idx,dat in enumerate(ff_qsq_data, 1):
        if idx == 1:
            name    = dat.split("\n")[0].rstrip('\n').split()[2]
            chisq   = dat.split("\n")[4].rstrip('\n').split()[8]
            plotLabel = "fit to " + name  + '\n' + r"$chisq/dof$ = " + ("{:6.4f}".format(float(chisq))) 
        if idx == 2:
            active = dat
        elif idx == 3:
            inactive = dat
        elif idx == 4:
            ensem = dat
        else: continue


    active_dat = active.split("\n")


    for elem in active_dat:

        this_line = (elem.rstrip('\n')).split()

        if ( len(this_line) == 0):
            continue
        if (this_line[0]=="#"):
            #commented out lines begin with "# "
            continue

        qsq = float(this_line[0])
        ff  = float(this_line[1])
        err = float(this_line[2])

        plt.errorbar( qsq, ff, yerr=err, markersize=0.2,fmt='.',color='black', mfc='black',mec='black',elinewidth=0.5, capsize=1, mew=0.7,zorder =10, label='_nolegend_')

    inactive_dat = inactive.split("\n")

    for elem in inactive_dat:

        this_line = (elem.rstrip('\n')).split()

        if ( len(this_line) == 0 ):
            continue
        if (this_line[0]=="#"):
            #commented out lines begin with "# "
            continue

        qsq = float(this_line[0])
        ff  = float(this_line[1])
        err = float(this_line[2])

        plt.errorbar( qsq, ff, yerr=err, markersize=0.2,fmt='.',color='red', mfc='red',mec='red',elinewidth=0.5, capsize=1, mew=0.7,zorder =5, label='_nolegend_')


    ensem_dat = ensem.split("\n")

    qsq = []
    ff  = []
    ffu = []
    ffb = []

    for elem in ensem_dat:

    
        this_line = (elem.rstrip('\n')).split()

        if ( len(this_line) == 0 ):
            continue
        if (this_line[0]=="#"):
            #commented out lines begin with "# "
            continue

        qsq.append(float(this_line[0]))
        ff.append(float(this_line[2]))
        ffb.append(float(this_line[1]))
        ffu.append(float(this_line[3]))



    plt.plot(qsq, ff, color='red', lw = 1,zorder =5, label=plotLabel)
    plt.fill_between(qsq, ffb, ffu, color='orange', alpha = 0.7,zorder =0 )


    #hactive = mpatches.Patch(color='black', label='active')
    #hinactive  = mpatches.Patch(color='red', label='inactive')
    #hensem  = mpatches.Patch(color='orange', label='fit to active',alpha = 0.7)

    #plt.legend(handles=[hactive,hinactive,hensem])
    plt.legend(loc='best',fontsize='x-small')
    plt.title('K' r'$\pi$ ' f'{name} matrix elements for ' r'$\rho$' ' current')



    fontP = FontProperties()
    fontP.set_size('small')
    plt.tight_layout()


#------------------------------------------------------------------------------------#
# ----- begin main program ----- #

if (len(sys.argv)!=2):
    print("usage: ", str(sys.argv[0]), "fit_qsq_plot_file")
    exit(1)

fit_qsq_plot_file = str(sys.argv[1])

plot_ff_qsq(fit_qsq_plot_file)


plt.savefig(fit_qsq_plot_file+".ff_qsq.pdf", bbox_inches='tight')
plt.show()
exit(0)

# End of plotter
