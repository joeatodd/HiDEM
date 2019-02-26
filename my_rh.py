##Based on rh.py from main HiDEM repo, but improved (hopefully)

import numpy as np
import getopt
import sys
import matplotlib.pyplot as plt
import math
from glob import glob
import re
from timeit import default_timer as timer
from HiDEM import Strain

# get input files from user

def usage():
    print ("Usage: " + sys.argv[0] + " -i input_file_glob -b buffer_dist (m, optional) -n interval (process every nth) -v max_str -d dx")

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:b:d:n:lv:")
except getopt.GetoptError:
    usage()
    sys.exit(2)

#default buffer distance (added to max_x and max_y
buff = 500.0
dx = 60.0
interval = 1
legacy_input = False
vmax = 0.01

for opt, arg in opts:
    if(opt =="-i"):
        inglob = arg
    elif(opt =='-b'):
        buff = float(arg)
    elif(opt == '-d'):
		dx = float(arg)
    elif(opt =='-n'):
        interval == int(arg)
    elif(opt =='-l'):
        legacy_input = True
    elif(opt =='-v'):
        vmax = float(arg)
    else:
        print "help"
        usage()
        sys.exit(2)

#Find the input files
if legacy_input:
    infiles = glob("*"+inglob+"*STR*.csv")
else:
    infiles = glob("*"+inglob+"*STR*.bin")
    
infiles.sort()

#Read the initial state
x0,y0,z0,str0 = Strain.str_from_file(infiles[0],legacy_input)

xx,yy,zz = Strain.grid_gen(x0,y0,z0,buff=buff,dx=dx)

for j,f in enumerate(infiles[:]):

    if j % interval != 0: continue
    print f

    x,y,z,str1 = Strain.str_from_file(f,legacy_input)

    strate = abs(str1[:])

    sx_plan, sx_side = Strain.grid_strain(x,y,z,strate,xx,yy,zz,buff=buff, dx=dx)

    plt.matshow(np.flipud(sx_side.T),vmin=0.0, vmax=vmax,cmap='jet')
    plt.savefig(f+"_side.png")
    plt.close()

    plt.matshow(sx_plan,vmin=0.0, vmax=vmax,cmap='jet')
    plt.savefig(f+"_plan.png")
    plt.close()
