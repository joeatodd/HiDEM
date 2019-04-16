import sys
import numpy as np
import matplotlib.pyplot as plt

indata_file = sys.argv[1]

indata = np.genfromtxt(indata_file, skip_header=True)

x_siz = np.unique(indata[:,0]).size
y_siz = np.unique(indata[:,1]).size

plt.matshow(np.reshape(indata[:,-1],(y_siz, x_siz)))

plt.savefig("geom_mask.png")
