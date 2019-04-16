##Based on rh.py from main HiDEM repo, but improved (hopefully)

import numpy as np
import getopt
import sys
import matplotlib.pyplot as plt
import math
from glob import glob
import re
from timeit import default_timer as timer
from HiDEM import Strain,Vtu

# get input files from user

def usage():
    print ("Usage: " + sys.argv[0] + " -i input_file_glob -b buffer_dist (m, optional) -n interval (process every nth) -v max_str -d dx")

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:b:d:n:v:")
except getopt.GetoptError:
    usage()
    sys.exit(2)

#default buffer distance (added to max_x and max_y
buff = 500.0
dx = 60.0
interval = 1
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
    elif(opt =='-v'):
        vmax = float(arg)
    else:
        print "help"
        usage()
        sys.exit(2)

#Find the input files
str_infiles = glob("*"+inglob+"*STR*.bin")
vtu_infiles = glob("*"+inglob+"*JYR*.vtu")
    
str_infiles.sort()
vtu_infiles.sort()

#x0,y0,z0,str0 = Strain.str_from_file_old(infiles[0],legacy_input)
#points0,str0 = Strain.pointstr_from_files(vtu_infiles[0],str_infiles[0])

#Get NN info (this doesn't change)
nn, efs0 = Strain.bonds_from_file(str_infiles[0])
points = Vtu.points_from_vtu(vtu_infiles[0])

fname_sth = Strain.fname_sth_from_str(str_infiles[0])
SCL = Strain.get_scl(fname_sth)
LNN = SCL * 1.1225

xx,yy,zz = Strain.grid_gen(points[:,0],points[:,1],points[:,2],buff=buff,dx=dx)

for j,f in enumerate(str_infiles[:]):

    if j % interval != 0: continue
    print f

    points = Vtu.points_from_vtu(vtu_infiles[j])

    lens = Strain.get_bond_lengths(points, nn)
    strain = (lens-LNN)/LNN
    cpoints = Strain.get_bond_centrepoints(points,nn)


    strate = abs(strain)

    print cpoints.shape
    print strate.shape
    sx_plan, sx_side = Strain.grid_strain(cpoints[:,0],cpoints[:,1],cpoints[:,2],strate,xx,yy,zz,buff=buff, dx=dx)

    plt.matshow(np.flipud(sx_side.T),vmin=0.0, vmax=vmax,cmap='jet')
    plt.savefig(f+"_side.png")
    plt.close()

    plt.matshow(sx_plan,vmin=0.0, vmax=vmax,cmap='jet')
    plt.savefig(f+"_plan.png")
    plt.close()
