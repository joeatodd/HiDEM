"""
Script to fix the cray-nonsense output in HiDEM: "2*0.0", "3*T" to "0.0 0.0" and "T T T" respectively
NB: This is potentially destructive, though ".BAK"ups will be made
"""

from glob import glob
import argparse
import os
import re

#Set up args
parser = argparse.ArgumentParser()
parser.add_argument("infile_glob")
args = parser.parse_args()

#regex for cray nonsense
cray_re = re.compile("([0-9]+)\*([TtFfEe\-0-9\.]+)")

#Get glob pattern (prefix, without asterisks!) for files
infile_glob = args.infile_glob
infiles = glob(infile_glob+"*")

#cycle input files
for f in infiles:

    #move to backup
    bak_name = f+".BAK"
    os.rename(f,bak_name)

    #cycle lines in now-backup
    with open(bak_name,'r') as source:
        lines = source.readlines()
        with open(f,'w') as dest:
            for line in lines:
                #find all cray nastiness
                reps = cray_re.finditer(line)
                for rep in reps:
                    old_str = rep.group(0)
                    cnt = int(rep.group(1))
                    num = rep.group(2)
                    #create string replacement (and remove last ', ')
                    new_str = ((num+", ") * cnt)[:-2]
                    print old_str
                    print new_str
                    line = line.replace(old_str, new_str, 1)
                    print line

                dest.write(line)

