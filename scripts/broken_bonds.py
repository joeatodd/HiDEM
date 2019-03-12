"""
Extract list of broken bonds from a STR.bin file and
write to new .bin
Idea is that this information could be extracted server-side and
then transferred.
NOTE - this is dumb - better to take the .vtu files and single STH file and forget the STRs.
"""
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
low_thresh = 0.5
for opt, arg in opts:
    if(opt =="-i"):
        inglob = arg
    elif(opt =='-l'):
        legacy_input = True
    elif(opt =='-v'):
        low_thresh = float(arg)
    else:
        print "help"
        usage()
        sys.exit(2)

high_thresh = 10.0

if legacy_input:
    str_infiles = glob("*"+inglob+"*STR*.csv")
else:
    str_infiles = glob("*"+inglob+"*STR*.bin")
    sth_infile = Strain.fname_sth_from_str(str_infiles[0])
    with open(sth_infile,'rb') as sth:
        sth_header = sth.readline()

vtu_infiles = glob("*"+inglob+"*JYR*.vtu")

str_infiles.sort()
vtu_infiles.sort()

for i in range(len(str_infiles)):
    f = str_infiles[i]
    v = vtu_infiles[i]

    nn, efs, strain = Strain.strain_from_file(f,v)

    #x,y,z,strain = Strain.str_from_file_old(f,legacy_input)

    broken = (strain > low_thresh) & (strain < high_thresh)

    nn_out = nn[broken,:]
    efs_out = efs[broken]
    # str_out = strain[broken]

    tstep_str = f.rsplit("STR",1)[-1].split(".")[0]
    sth_name = Strain.fname_sth_from_str(f)

    sth_outfile = sth_name.rsplit(".",1)[0] + tstep_str + "_broken.bin"
    str_outfile = f.rsplit(".",1)[0] + "_broken.bin"

    print sth_outfile

    Strain.sth_to_file(nn_out,sth_header,sth_outfile)
    Strain.str_to_file(efs_out, str_outfile)

# output =  np.stack((xout,yout,zout,strout),axis=-1)
# print output.shape
# np.savetxt("test_out.csv",output,header="X,Y,Z,Strain", comments="", delimiter=",")
