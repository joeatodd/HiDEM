import numpy as np
import getopt
import sys
import matplotlib.pyplot as plt
import math

# get input files from user

def usage():
    print ("Usage: " + sys.argv[0] + " -1 first_strain.csv -2 second_strain.csv")

try:
    opts, args = getopt.getopt(sys.argv[1:],"1:2:")
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if(opt =="-1"):
        str1_file = arg
    elif(opt =="-2"):
        str2_file = arg
    else:
        print "help"
        usage()
        sys.exit(2)

#read files
str1 = np.fromfile(str1_file,sep=" ").reshape((-1,4))
str2 = np.fromfile(str2_file,sep=" ").reshape((-1,4))

#get the strain 'rate'
strate = str2[:,3] - str1[:,3]

extent_min = str1[:,:-1].min(0)
extent_max = str1[:,:-1].max(0)
extent = [extent_min[0],extent_max[0],extent_min[1],extent_max[1],extent_min[2],extent_max[2]]

dx = 60.0
dz = 60.0

xmin = int(math.floor(extent[0] / dx)) * dx
xmax = int(math.ceil(extent[1] / dx)) * dx
ymin = int(math.floor(extent[2] / dx)) * dx
ymax = int(math.ceil(extent[3] / dx)) * dx
zmin = int(math.floor(extent[4] / dz)) * dz
zmax = int(math.ceil(extent[5] / dz)) * dz

xx = np.arange(0,xmax+dx,dx)
yy = np.arange(0,ymax+dx,dx)
zz = np.arange(0,zmax+dz,dz)

sx_plan = np.zeros([len(xx), len(yy)])
nx_plan = np.zeros([len(xx), len(yy)])
sx_side = np.zeros([len(yy), len(zz)])
nx_side = np.zeros([len(yy), len(zz)])

for i in range(len(strate)):

    xbox = int(round(str1[i,0] / dx))
    ybox = int(round(str1[i,1] / dx))
    zbox = int(round(str1[i,2] / dz))

    sx_plan[xbox,ybox] += strate[i]
    nx_plan[xbox,ybox] += 1
    sx_side[ybox,zbox] += strate[i]
    nx_side[ybox,zbox] += 1

sx_plan /= nx_plan
sx_side /= nx_side

sx_side[np.isnan(sx_side)] = 0.0
sx_plan[np.isnan(sx_plan)] = 0.0

stdev = np.std(strate)

plt.matshow(np.flipud(sx_side.T),vmax=0.01)
plt.savefig(str2_file+".png")

plt.matshow(sx_plan,vmax=0.01)
plt.savefig(str2_file+"_plan.png")
