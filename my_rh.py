##Based on rh.py from main HiDEM repo, but improved (hopefully)

import numpy as np
import getopt
import sys
import matplotlib.pyplot as plt
import math
from glob import glob
import re
from timeit import default_timer as timer


# get input files from user

def usage():
    print ("Usage: " + sys.argv[0] + " -i input_file_glob -b buffer_dist (m, optional) -n interval (process every nth) -v max_str -d dx")

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:b:d:n:lv:")
except getopt.GetoptError:
    usage()
    sys.exit(2)

def str_from_file(fname,legacy_input):
    if(legacy_input):
        the_data = np.fromfile(fname,sep=" ")
    else:
        count_re = re.compile("Count: ([0-9]*)")
        type_re = re.compile('Type: ([a-zA-Z0-9]*)')

        indata = open(fname, 'rb')
        header = indata.readline()
        num_count = int(count_re.search(header).group(1))
        float_type = type_re.search(header).group(1).lower()
        the_data = np.fromfile(indata, dtype=np.dtype(float_type), count=num_count*4)

    return the_data.reshape((-1,4))

#default buffer distance (added to max_x and max_y
buff = 500.0
dx = 60.0
buff_boxes = math.floor(buff/dx)
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
str0 = str_from_file(infiles[0],legacy_input)

extent_min = str0[:,:-1].min(0)
extent_max = str0[:,:-1].max(0)
extent = [extent_min[0]-buff,extent_max[0]+buff,extent_min[1],extent_max[1]+buff,extent_min[2],extent_max[2]]

xmin = int(math.floor(extent[0] / dx)) * dx
xmax = int(math.ceil(extent[1] / dx)) * dx
ymin = int(math.floor(extent[2] / dx)) * dx
ymax = int(math.ceil(extent[3] / dx)) * dx
zmin = int(math.floor(extent[4] / dx)) * dx
zmax = int(math.ceil(extent[5] / dx)) * dx

xx = np.arange(xmin,xmax+dx,dx)
yy = np.arange(ymin,ymax+dx,dx)
zz = np.arange(zmin,zmax+dx,dx)

for j,f in enumerate(infiles[1:]):

    if j % interval != 0: continue
    print f

    str1 = str_from_file(f,legacy_input)

    #get the strain 'rate'
    strate = abs(str1[:,3] - str0[:,3])


    sx_plan = np.zeros([len(xx), len(yy)])
    nx_plan = np.zeros([len(xx), len(yy)])
    sx_side = np.zeros([len(yy), len(zz)])
    nx_side = np.zeros([len(yy), len(zz)])

    #find the voxel in which each particle belongs
    boxes = np.rint((str1[:,:-1] / np.array((dx,dx,dx)))).astype(int)
    boxes[:,0] += int(math.floor(buff_boxes))

    #check box bounds
    boxes[:,0] = boxes[:,0].clip(0,sx_plan.shape[0]-1)
    boxes[:,1] = boxes[:,1].clip(0,sx_plan.shape[1]-1)
    boxes[:,2] = boxes[:,2].clip(0,sx_side.shape[1]-1)

    np.add.at(sx_plan, (boxes[:,0],boxes[:,1]),strate)
    np.add.at(sx_side, (boxes[:,1],boxes[:,2]),strate)
    np.add.at(nx_plan, (boxes[:,0],boxes[:,1]),np.ones(strate.size,dtype=int))
    np.add.at(nx_side, (boxes[:,1],boxes[:,2]),np.ones(strate.size,dtype=int))

    sx_plan /= nx_plan
    sx_side /= nx_side

    sx_side[np.isnan(sx_side)] = 0.0
    sx_plan[np.isnan(sx_plan)] = 0.0

    #stdev = np.std(strate)

    plt.matshow(np.flipud(sx_side.T),vmin=0.0, vmax=vmax,cmap='jet')
    plt.savefig(f+"_side.png")
    plt.close()

    plt.matshow(sx_plan,vmin=0.0, vmax=vmax,cmap='jet')
    plt.savefig(f+"_plan.png")
    plt.close()
