"""
Script to analyse the mélange particle interactions from a HiDEM simulation.
Set up for Store Glacier but would work with any, given modification of bounding box.

Particle stress tensor is defined using the Love-Weber formula, described most succintly here:
https://journals.aps.org/pre/pdf/10.1103/PhysRevE.72.041307

but laid out in more detail here:
https://www.sciencedirect.com/science/article/pii/S0020768313001492

Script is able to run in parallel using MPI, allowing processing of a whole series of files at once.

Modifiable variables which are not yet exposed as arguments:

SCL - the particle size (diameter)
EFS - Youngs modulus of simulation
box_pixel - the resolution of output image
polygon - rectangle defining the region of interest

Output:
_melthick.png - thickness of melange
_melstress.png - depth integrated stress in pixel (but this is potentially incorrectly formulated!)
_melfchains.png - force chain analysis - the most compressive pixels marked by principal compressive direction

"""
import argparse

#Set up args - do this before expensive paraview import in case of errors - fail quickly!
#=====================================
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

    outfile = vtu_in+"_melstress.png"
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

    #filter out those too close - TODO investigate further
    # close = close[rc > Thresh]
    # rc = rc[rc > Thresh]

    #xyz components of particle distance
    rcx = (points[close[:,0],0] - points[close[:,1],0])/rc
    rcy = (points[close[:,0],1] - points[close[:,1],1])/rc
    rcz = (points[close[:,0],2] - points[close[:,1],2])/rc

    #centrepoints of interaction
    cpoints = (points[close[:,0],:] + points[close[:,1],:])/2.0

    #Compute force of each collision
    #==============================

    #pulled this from circ.f90
    force = (SCL * EFS * (LNN - rc)**1.5) / 1.0e6 #meganewtons
    frx = rcx * SCL * EFS * (LNN - rc)**1.5 / 1.0e6
    fry = rcy * SCL * EFS * (LNN - rc)**1.5 / 1.0e6
    frz = rcz * SCL * EFS * (LNN - rc)**1.5 / 1.0e6

    force_vec = np.stack((frx,fry,frz)).T

    #Limit all stresses to a 180 degree range (invert those which have -ve Y)
    #Prevents equivalent compressions from cancelling each other out
    fry_dir = np.sign(fry)
    # frx *= fry_dir
    # fry *= fry_dir
    # frz *= fry_dir

    #the normalised force vector
    force_direction = (np.stack((frx,fry,frz))/force).T

    #Collect force contributions on each particle
    #=========================================

    particle_force = np.zeros(points.shape)

    np.add.at(particle_force[:,0], close[:,:], np.abs(np.stack((frx,frx))).T)
    np.add.at(particle_force[:,1], close[:,:], np.abs(np.stack((fry,fry))).T)
    np.add.at(particle_force[:,2], close[:,:], np.abs(np.stack((frz,frz))).T)
    #force magnitude
    particle_fmag = np.linalg.norm(particle_force, axis=1)

    #the force direction...
    #this causes issues - how does one 'average' a direction...
    particle_forcedir = np.zeros(points.shape)

    # pf_cnt = np.zeros(points.shape[0])
    # pf_ones = np.ones(force_direction[:,0].shape)
    # np.add.at(pf_cnt, close[:,:], np.stack((pf_ones, pf_ones)).T)

    np.add.at(particle_forcedir[:,0], close[:,:], np.stack((force_direction[:,0], 
                                                            force_direction[:,0])).T)
    np.add.at(particle_forcedir[:,1], close[:,:], np.stack((force_direction[:,1], 
                                                            force_direction[:,1])).T)
    np.add.at(particle_forcedir[:,2], close[:,:], np.stack((force_direction[:,2], 
                                                            force_direction[:,2])).T)

    particle_forcedir /= np.linalg.norm(particle_forcedir, axis=1)[:,None]


    #Particle based full force vector for force chains
    #=======================================
    # Eq.1 from https://doi.org/10.1103/PhysRevE.72.041307 
    # TODO -  What is the representative volume of these particles?
    # Their mass in wave2.f90 is given by a cube but they interact as
    # spheres....
    # 
    # https://www.sciencedirect.com/science/article/pii/S0020768313001492
    # ^- good discussion of particle stress vectors (this approach is Love-Weber)


    #find the radius vectors for each particle collision
    #actually these are symmetric because of same particle size
    #so just compute one and reverse it
    rvecs = cpoints - points[close[:,0]]

    #the stress contribution (3x3) for both sides of each collision
    #dyadic product
    p0_stress = rvecs[:,:,None] * force_vec[:,None,:]
    p1_stress = -rvecs[:,:,None] * -force_vec[:,None,:]

    #sum up contributions from collisions
    particle_stress = np.zeros((points.shape[0],3,3))
    np.add.at(particle_stress, close[:,0], p0_stress)
    np.add.at(particle_stress, close[:,1], p1_stress)
    
    # * (1/V) to give a stress
    particle_stress /= ((4/3.0) * np.pi * ((SCL/2.0) ** 3.0))

    eigv, eigd = np.linalg.eigh(particle_stress)

    #compressive principal stress direction
    eigd_compr = eigd[:,:,0]
    #limit eigd to 180 degrees
    eigd_compr *= np.sign(eigd_compr[:,1])[:,None]

    #Threshold first on z component (not interested in vertical stress for force chains)
    bigforce_thresh = -0.4
    bigforce = (np.abs(eigd_compr[:,-1]) < (1/(2.0**0.5)))
    #Then thresh on value of stress - NB arbitrary factor 0.5 here
#    bigforce = (eigv[:,0] < 0.5*np.mean(eigv[bigforce,0])) & bigforce #(np.abs(eigd_compr[:,-1]) < (1/(2.0**0.5)))
    bigforce = (eigv[:,0] < bigforce_thresh) & bigforce #(np.abs(eigd_compr[:,-1]) < (1/(2.0**0.5)))

    M = -eigd_compr[bigforce,1]


    #Rasterise/bin the particle data for display
    #====================================
    #Cycle through points in box, adding up forces from 'force'

    xx = np.arange(bboxx[0],bboxx[1],box_pixel)
    yy = np.arange(bboxy[0],bboxy[1],box_pixel)

    nx = xx.size
    ny = yy.size

    #find xy bin for each particle
    xbins = np.digitize(points[:,0], xx) - 1
    ybins = np.digitize(points[:,1], yy) - 1

    #compute pixel values for total force (vectorised!)
    #These have already been 'abs'd
    boxed_force = np.zeros((nx,ny,3))

    np.add.at(boxed_force[:,:,0], (xbins, ybins), particle_force[:,0])
    np.add.at(boxed_force[:,:,1], (xbins, ybins), particle_force[:,1])
    np.add.at(boxed_force[:,:,2], (xbins, ybins), particle_force[:,2])

    boxed_forcemag = np.linalg.norm(boxed_force,axis=2) / (box_pixel**2.0)
    boxed_force /= (box_pixel**2.0)


    #Mean force direction in boxes
    boxed_forcedir = np.zeros((nx,ny,3))
    np.add.at(boxed_forcedir, (xbins, ybins), particle_forcedir)
    boxed_forcedir /= np.linalg.norm(boxed_forcedir, axis=2)[:,:,None]

    #Melange thickness per box
    #===========================

    #Count particles in each box
    boxed_pcount = np.zeros((nx,ny))
    particle_counters = np.ones(npoints)
    np.add.at(boxed_pcount, (xbins, ybins), particle_counters)

    pvol = SCL**3.0 #although they interact as spheres, particles have cubic mass
    #Thickness is (number of particles * pvol) / box_area
    box_melthick = (boxed_pcount * pvol)/(box_pixel ** 2.0)

    #Figures
    #===========================

    #Melange thickness
    mat = plt.matshow(np.flipud(box_melthick.T), vmax=300.0, cmap='terrain')
    fig = plt.gcf()
    ax = plt.gca()
    cbar=fig.colorbar(mat,fraction=0.05)
    cbar.set_label("Melange Thickness (m)")
    ax.set_title(figtitle)
    plt.savefig(vtu_in+"_melthick.png",dpi=600)
    plt.close()

    #Pressure magnitude figure
    mat = plt.matshow(np.flipud(boxed_force[:,:,1].T), vmax=2.0, cmap='jet')
    fig = plt.gcf()
    ax = plt.gca()
    cbar=fig.colorbar(mat,fraction=0.05)
    cbar.set_label("Melange stress (MPa)")
    ax.set_title(figtitle)
    plt.savefig(outfile,dpi=600)
    plt.close()

    # #Just plot the x part of the boxed force direction
    # mat = plt.matshow(boxed_forcedir[:,:,0], vmin=-0.2, vmax=0.2, cmap='seismic')

    #Particle forces for above average force
    plt.quiver(points[bigforce,0],points[bigforce,1],eigd_compr[bigforce,0],
               eigd_compr[bigforce,1],M,pivot='mid',scale_units="width", 
               scale=100.0)
    ax = plt.gca()
    ax.set_aspect(1)
    plt.savefig(vtu_in+"_melfchains.png",dpi=600)
    plt.close()


#loop over files, making figs
for f in vtu_infiles:
    print "Processing: "+f
    compression_plot(f)

