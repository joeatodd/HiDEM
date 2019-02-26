import numpy as np
from glob import glob
import re
from HiDEM import Strain

x,y,z,strain = Strain.str_from_file("ElmerHiDEM_store_summer_bg1pl1_uc032_r01_STR0050.bin")
thresh = np.arange(0,2.1,0.1)

print strain.shape
for t in thresh:
    print t, (np.sum(strain > t) / float(strain.size))


thresh = 0.5

broken = strain > thresh

xout = x[broken]
yout = y[broken]
zout = z[broken]
strout = strain[broken]

output =  np.stack((xout,yout,zout,strout),axis=-1)
print output.shape
np.savetxt("test_out.csv",output,header="X,Y,Z,Strain", comments="", delimiter=",")
