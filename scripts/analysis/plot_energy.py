import numpy as np
import matplotlib.pyplot as plt
import getopt
import sys
from glob import glob
from HiDEM import Vtu

def usage():
    print("Usage: " + sys.argv[0] + " -i in_glob")

try:
	opts, args = getopt.getopt(sys.argv[1:],"i:")

except getopt.GetoptError:
	usage()
	sys.exit(2)

for opt, arg in opts:
	if(opt =="-i"):
	    in_glob = arg
	    got_infile = True
	else:
	    usage()
	    sys.exit(2)
            

#all data
labels = ["Time","Collis","Stiff","Kin","Pot","BeamDamp","Drag","Press","BedInt","BedDamp"]
#which to plot?
idx = [1,2,3,4,6,8]
colours = ['blue','green','red','orange','yellow','purple']
alpha = [1.0, 1.0, 1.0, 1.0, 1.0, 0.5]
infiles = glob(in_glob+"*dtop00")
infiles.sort(key=Vtu.vtu_sort_fun)
infiles = [f.split("dtop00")[0] for f in infiles]

def read_dtop(run_name):
    # dtopr: t, sum?, sum?
    # dtopr not actually used at present
    dtopr = np.genfromtxt(run_name+"dtopr")

    # dtop00: t, WENS (collis), ENMS (beam stiff), KINS (translational kinetic), MGHS (potential)
    dtop00 = np.genfromtxt(run_name+"dtop00")

    # dtop01: t, DPES (beam damping), DMPENS (drag), PSUMS (pressure), GSUMS (bed int), BDES (bed damp)
    dtop01 = np.genfromtxt(run_name+"dtop01")

    dtopall = np.concatenate((dtop00, dtop01[:,1:]),1)

    return dtopall


def plot_run(run_name, data, tens, nineties):

    pltcount = len(idx)
    fig, ax1 = plt.subplots(figsize=(10,10))
    axes = [ax1]
    [axes.append(ax1.twinx()) for i in range(pltcount-1)]

    lns = []
    for i in range(pltcount):
        ln = axes[i].plot(data[:,0], data[:,idx[i]],label=labels[idx[i]],color=colours[i], alpha=alpha[i
])
        lns+=ln
        axes[i].set_autoscale_on(False)

        if i % 2:
            axes[i].yaxis.tick_left()
        else:
            axes[i].yaxis.tick_right()

    for i in range(pltcount):
        axes[i].set_ylim(ymin=tens[idx[i]], ymax=nineties[idx[i]])
        #axes[i].set_yticklabels(axes[i].get_yticks())
        axes[i].tick_params(axis='y', labelcolor=colours[i])

    # lns = l1+l2+l3+l4+l5+l6
    labs = [l.get_label() for l in lns]

    leg = ax1.legend(lns,labs,loc=0)


    axes[0].set_xlabel('time (s)')
    plt.savefig(run_name+"_dtopr.png")
    plt.close()

#get all the data
data = [read_dtop(f) for f in infiles]

#combine all to 1 array to compute stats for lims
full_series = np.concatenate(data, 0)

#percentile calc for plot range
nineties = np.percentile(full_series, 90, 0)
tens = np.percentile(full_series, 10, 0)
rangies = nineties - tens
nineties+= rangies*0.5
tens-= rangies*0.5

mins = np.min(full_series, 0)
maxs = np.max(full_series, 0)

for i,f in enumerate(data):
    plot_run(infiles[i], f, tens, nineties)
