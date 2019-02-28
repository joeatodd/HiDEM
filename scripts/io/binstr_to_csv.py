import numpy as np
import getopt
import sys
from glob import glob
import re

def usage():
    print ("Usage: " + sys.argv[0] + " -i input_file_glob -l [legacy input?]")

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:l")
except getopt.GetoptError:
    usage()
    sys.exit(2)

legacy_input = False

for opt, arg in opts:
    if(opt =="-i"):
        inglob = arg
    elif(opt =='-l'):
        legacy_input = True

#Find the input files
if legacy_input:
    infiles = glob("*"+inglob+"*STR*.csv")
else:
    infiles = glob("*"+inglob+"*STR*.bin")

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
        return the_data.reshape((4,-1)).T

for f in infiles:

    print("Processing: "+f)
    strain = str_from_file(f,legacy_input)
    np.savetxt(f+'.csv',strain,header="X,Y,Z,Strain",comments="",delimiter=",", fmt="%.12e")
