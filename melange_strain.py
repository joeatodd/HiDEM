import argparse

#Set up args - do this before expensive paraview import in case of errors - fail quickly!
parser = argparse.ArgumentParser()
parser.add_argument("vtu_inglob", type=str, help="e.g. HiDEM_Run*.vtu")
parser.add_argument("--no-mpi", dest='mpi', default=True, action="store_false", help="If you get nasty mpi errors trying to run this in serial, use this flag.")
parser.add_argument("--redo", default=False, action="store_true", help="If image file already exists, by default processing is not redone")
args = parser.parse_args()

import numpy as np
from paraview.simple import servermanager, XMLUnstructuredGridReader, GetActiveSource
from paraview.numpy_support import vtk_to_numpy
from scipy.spatial import cKDTree
import re

#weird stuff required for matplotlib (temp file for matplotlib workspace I guess?)
import os
import tempfile
os.environ['MPLCONFIGDIR'] = tempfile.mkdtemp()
import matplotlib.pyplot as plt

from glob import glob


#Check for serial or MPI (see notes on --no-mpi if errors occur)
if args.mpi:
    try:
        from mpi4py import MPI
        mpi_job = True
    except:
        print("Failed to load mpi - falling back to serial")
        mpi_job = False
else:
    mpi_job = False

if mpi_job:
    mpi_size = MPI.COMM_WORLD.Get_size()
    mpi_rank = MPI.COMM_WORLD.Get_rank()

redo = args.redo

#glob pattern for files to operate on
vtu_inglob = args.vtu_inglob
vtu_inglob = vtu_inglob.rsplit("vtu",1)[0] #strip off 'vtu' if specified (re-added below)

#user specifiable
box_pixel = 50.0

#These two should match simulation
SCL = 20.0 #particle size
EFS = 1.0e9 #youngs modulus (1e9 is HiDEM default)

LNN = SCL * 1.1225
Thresh = LNN * 0.9 #filter noise

#regex for extracting run name
title_re = re.compile("_([0-9]*_[0-9]*)_JYR([0-9]*)")

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

assert len(vtu_infiles) > 0, 'Bad glob!'

#take a subset of files or don't (embarrassingly parallel)
if mpi_job:
    vtu_infiles = vtu_infiles[mpi_rank::mpi_size]

def compression_plot(vtu_in):

    outfile = vtu_in+".png"
    if os.path.isfile(outfile) and not redo:
        print("Skipping existing output: "+outfile)
        return

    #generate figure title
    figtitle = "Run: "+title_re.search(vtu_in).group(1) + " Step: " + title_re.search(vtu_in).group(2)

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
    rcx = (points[close[:,0],0] - points[close[:,1],0])/rc
    rcy = (points[close[:,0],1] - points[close[:,1],1])/rc
    rcz = (points[close[:,0],2] - points[close[:,1],2])/rc

    #centrepoints of interaction
    cpoints = (points[close[:,0],:] + points[close[:,1],:])/2.0

    #pulled this from circ.f90
    force = (SCL * EFS * (LNN - rc)**1.5) / 1.0e6 #meganewtons
    frx = abs(rcx * SCL * EFS * (LNN - rc)**1.5) / 1.0e6
    fry = abs(rcy * SCL * EFS * (LNN - rc)**1.5) / 1.0e6
    frz = abs(rcz * SCL * EFS * (LNN - rc)**1.5) / 1.0e6

    #Cycle through points in box, adding up forces from 'force'
    particle_force = np.zeros(points.shape)

    np.add.at(particle_force[:,0], close[:,:], np.stack((frx,frx)).T)
    np.add.at(particle_force[:,1], close[:,:], np.stack((fry,fry)).T)
    np.add.at(particle_force[:,2], close[:,:], np.stack((frz,frz)).T)
    #force magnitude
    particle_fmag = np.linalg.norm(particle_force, axis=1)


    xx = np.arange(bboxx[0],bboxx[1],box_pixel)
    yy = np.arange(bboxy[0],bboxy[1],box_pixel)

    nx = xx.size
    ny = yy.size

    #find xy bin for each particle
    xbins = np.digitize(points[:,0], xx) - 1
    ybins = np.digitize(points[:,1], yy) - 1

    #compute pixel values for total force (vectorised!)
    boxed_force = np.zeros((nx,ny,3))
    np.add.at(boxed_force[:,:,0], (xbins, ybins), particle_force[:,0])
    np.add.at(boxed_force[:,:,1], (xbins, ybins), particle_force[:,1])
    np.add.at(boxed_force[:,:,2], (xbins, ybins), particle_force[:,2])

    boxed_forcemag = np.linalg.norm(boxed_force,axis=2) / (box_pixel**2.0)

    boxed_force /= (box_pixel**2.0)
    mat = plt.matshow(boxed_force[:,:,1], vmax=2.0, cmap='jet')
    fig = plt.gcf()
    ax = plt.gca()
    fig.colorbar(mat,fraction=0.05)
    ax.set_title(figtitle)
    plt.savefig(outfile)


#loop over files, making figs
for f in vtu_infiles:
    print "Processing: "+f
    compression_plot(f)
