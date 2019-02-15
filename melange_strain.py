import numpy as np
import argparse
from paraview import simple as PV
from paraview.numpy_support import vtk_to_numpy
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

#Set up args
parser = argparse.ArgumentParser()
parser.add_argument("vtu_in")
args = parser.parse_args()

vtu_in = args.vtu_in

vtu_in = "ElmerHiDEM_store_melange_genmel_5_2_JYR0352.vtu"

SCL = 20.0
LNN = SCL * 1.1225

#polygon for points to check
polygon = np.zeros((4,2))
polygon[0,:] = (-400,5050.0)
polygon[1,:] = (-400,5300.0)
polygon[2,:] = (5200.0,5300.0)
polygon[3,:] = (5200.0,5050.0)

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


dists = cdist(points, points)
close = np.where((dists < LNN) & (dists > LNN*0.99))

collis = (points[close[0],:] + points[close[1],:]) / 2
collis_dist = dists[close]

collis_stress = LNN - collis_dist 

plt.scatter(collis[collis_stress > 0.01,0],collis[collis_stress > 0.01, 1], c=collis_stress[collis_stress > 0.01])
