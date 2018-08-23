import struct
import re
import getopt
import sys
import numpy as np

def usage():
    print ("Usage: " + sys.argv[0] + " -i input_file")

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:")
except getopt.GetoptError:
    usage()
    sys.exit(2)

for opt, arg in opts:
    if(opt =="-i"):
        infile = arg
    else:
        print "help"
        usage()
        sys.exit(2)

#output file name
outfile = infile.split(".")[0]+'.csv'

#regexps for parsing vtu file
type_re = re.compile(r"type=\"([A-Za-z0-9]*)\"")
offset_re = re.compile(r"offset=\"([0-9]*)\"")
endian_re = re.compile(r"byte_order=\"([A-Za-z]*)\"")


indata = open(infile,'rb')

#Read info from vtu xml header
while True:
    l = indata.readline()

    if 'Position' in l:
        float_type = type_re.search(l).group(1)
        point_offset = offset_re.search(l).group(1)
        print float_type
        print point_offset
    if '<VTKFile' in l:
        endianness = endian_re.search(l).group(1)
    if '<Appended' in l: break

print 'Reading data in '+float_type+' format from '+infile

indata.read(1) # the '_' symbol which starts VTK appended data

bytecount = struct.unpack('i',indata.read(4))[0]

#define the expected structure

if 'little' in endianness.lower():
    float_format = '<'
elif 'big' in endianness.lower():
    float_format = '>'
else:
    raise Exception("Didn't understand endianness of input: "+endiannness)


if '32' in float_type:
    float_format += 'f'
    num_count = bytecount / 4
elif '64' in float_type:
    float_format += 'd'
    num_count = bytecount / 8
else:
    raise Exception("Didn't understand float format for points: "+float_type)

np_data = np.fromfile(indata, dtype=np.dtype(float_format), count=num_count).reshape((-1,3))

indata.close()

np.savetxt(outfile, np_data, fmt='%16.6f')
