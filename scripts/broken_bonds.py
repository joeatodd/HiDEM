import sys
import getopt
import numpy as np
from glob import glob
import re
from HiDEM import Strain

def usage():
    print ("Usage: " + sys.argv[0] + " -i input_file_glob -v strain_threshold -l legacy_input (switch)")

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:lv:")
except getopt.GetoptError:
    usage()
    sys.exit(2)

legacy_input = False
thresh = 0.5
for opt, arg in opts:
    if(opt =="-i"):
        inglob = arg
    elif(opt =='-l'):
        legacy_input = True
    elif(opt =='-v'):
        thresh = float(arg)
    else:
        print "help"
        usage()
        sys.exit(2)

if legacy_input:
    infiles = glob("*"+inglob+"*STR*.csv")
else:
    infiles = glob("*"+inglob+"*STR*.bin")
infiles.sort()

for f in infiles:
    x,y,z,strain = Strain.str_from_file(f,legacy_input)

    broken = strain > thresh

    xout = x[broken]
    yout = y[broken]
    zout = z[broken]
    strout = strain[broken]

    outfile = f.rsplit(".",1)[0] + "_broken.bin"
    print outfile
    Strain.str_to_file(xout,yout,zout,strout,outfile)

# output =  np.stack((xout,yout,zout,strout),axis=-1)
# print output.shape
# np.savetxt("test_out.csv",output,header="X,Y,Z,Strain", comments="", delimiter=",")
