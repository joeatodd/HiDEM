import numpy as np
import argparse
from paraview.simple import servermanager, XMLUnstructuredGridReader, GetActiveSource
from paraview.numpy_support import vtk_to_numpy
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt
from glob import glob

#Set up args
parser = argparse.ArgumentParser()
parser.add_argument("vtu_inglob", type=str, help="e.g. HiDEM_Run*.vtu")
args = parser.parse_args()

#glob pattern for files to operate on
vtu_inglob = args.vtu_inglob
vtu_inglob = vtu_inglob.split("vtu")[0] #strip off 'vtu' if specified (re-added below)

#user specifiable
box_pixel = 50.0

#These two should match simulation
SCL = 20.0 #particle size
EFS = 1.0e9 #youngs modulus (1e9 is HiDEM default)

LNN = SCL * 1.1225
Thresh = LNN * 0.9 #filter noise

#polygon for points to check (assume rectangle atm)
polygon = np.zeros((4,2))
polygon[0,:] = (-400,4500.0)
polygon[1,:] = (-400,7500.0)
polygon[2,:] = (5200.0,7500.0)
polygon[3,:] = (5200.0,4500.0)

bboxx = (np.min(polygon[:,0]),np.max(polygon[:,0]))
bboxy = (np.min(polygon[:,1]),np.max(polygon[:,1]))

#get the files
vtu_infiles = glob(vtu_inglob+"vtu")
vtu_infiles.sort()


def compression_plot(vtu_in):

    #read the vtu input
    reader = XMLUnstructuredGridReader( FileName=vtu_in )
    reader.UpdatePipeline()

    data = servermanager.Fetch(reader)
    GlacierSource = GetActiveSource()

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


#loop over files, making figs
for f in vtu_infiles:
    print "Processing: "+f
    compression_plot(f)
