import numpy as np
import argparse
from paraview import simple as PV
from paraview.numpy_support import vtk_to_numpy
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

#Set up args
parser = argparse.ArgumentParser()
parser.add_argument("vtu_in")
args = parser.parse_args()

vtu_in = args.vtu_in
# vtu_in = "ElmerHiDEM_store_melange_genmel_7_8_JYR1312.vtu"

box_pixel = 50.0
SCL = 20.0
EFS = 1.0e9

LNN = SCL * 1.1225
Thresh = LNN * 0.9 #filter noise

#polygon for points to check
polygon = np.zeros((4,2))
polygon[0,:] = (-400,4500.0)
polygon[1,:] = (-400,7500.0)
polygon[2,:] = (5200.0,7500.0)
polygon[3,:] = (5200.0,4500.0)

bboxx = (np.min(polygon[:,0]),np.max(polygon[:,0]))
bboxy = (np.min(polygon[:,1]),np.max(polygon[:,1]))

#read the vtu input
reader = PV.XMLUnstructuredGridReader( FileName=vtu_in )
reader.UpdatePipeline()

data = PV.servermanager.Fetch(reader)
GlacierSource = PV.GetActiveSource()

points = vtk_to_numpy(data.GetPoints().GetData())

#clear points beyond poly
points = points[points[:,0] > bboxx[0],:]
points = points[points[:,0] < bboxx[1],:]
points = points[points[:,1] > bboxy[0],:]
points = points[points[:,1] < bboxy[1],:]

npoints = points.shape[0]

#TODO - Point In Polygon

#Identify close points & compute dists

#construct kd tree & find pairs closer than LNN
tree = cKDTree(points)
close = np.asarray(list(tree.query_pairs(LNN)))

#compute dists of above
rc = np.linalg.norm((points[close[:,0],:] - points[close[:,1],:]),axis=1)

#filter out those too close
close = close[rc > Thresh]
rc = rc[rc > Thresh]

#xyz components of particle distance
rcx = (points[close[:,0],0] - points[close[:,1],0])
rcy = (points[close[:,0],1] - points[close[:,1],1])
rcz = (points[close[:,0],2] - points[close[:,1],2])

#centrepoints of interaction
cpoints = (points[close[:,0],:] + points[close[:,1],:])/2.0

#force (missing youngs modulus)
#pulled this from circ.f90
force = (rc * EFS * (LNN - rc)**1.5) / 1.0e6 #meganewtons
frx = (abs(rcx * EFS * (LNN - rc)**1.5)) / 1.0e6
fry = (abs(rcy * EFS * (LNN - rc)**1.5)) / 1.0e6
frz = (abs(rcz * EFS * (LNN - rc)**1.5)) / 1.0e6


#Plot the collisions w/ stress - this is not very compelling, try something else
# plt.scatter(collis_centre[:,0],collis_centre[:, 1], c=force)
# ax = plt.gca()
# ax.set_aspect(1)

xx = np.arange(bboxx[0],bboxx[1],box_pixel)
yy = np.arange(bboxy[0],bboxy[1],box_pixel)

nx = xx.size
ny = yy.size

# boxed_force = np.zeros((nx,ny))
# for i in range(nx-1):
#     for j in range(ny-1):
#         #gather stress in box
#         boxforce = force[(cpoints[:,0] > xx[i]) & (cpoints[:,0] < xx[i+1]) & (cpoints[:,1] > yy[j]) & (cpoints[:,1] < yy[j+1])]
#         if(boxforce.size > 0):
#             boxed_force[i,j] = np.sum(boxforce)#/boxforce.size
#         else:
#             boxed_force[i,j] = 0.0

# plt.matshow(boxed_force,vmax=5.0)


#compute pixel values for force in Y direction
boxed_force = np.zeros((nx,ny))
for i in range(nx-1):
    for j in range(ny-1):
        #gather stress in box
        boxforce = fry[(cpoints[:,0] > xx[i]) & (cpoints[:,0] < xx[i+1]) & (cpoints[:,1] > yy[j]) & (cpoints[:,1] < yy[j+1])]
        if(boxforce.size > 0):
            boxed_force[i,j] = np.sum(boxforce)#/boxforce.size
        else:
            boxed_force[i,j] = 0.0

plt.matshow(boxed_force, vmax=2000.0)
plt.savefig(vtu_in+".png")
