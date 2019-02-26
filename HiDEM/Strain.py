"""
Module for strain/bond breakage related stuff. i.e. interaction with *STR* files
"""

import re
import numpy as np
import math

def str_from_file(fname,legacy_input=False):
    """
    Returns the bond strain rate from given file
    Strain is change_in_length / original_length
    """
    if(legacy_input):
        the_data = np.fromfile(fname,sep=" ")
    else:
        count_re = re.compile("Count: ([0-9]*)")
        type_re = re.compile('Type: ([a-zA-Z0-9]*)')

        indata = open(fname, 'rb')
        header = indata.readline()
        num_count = int(count_re.search(header).group(1))
        float_type = type_re.search(header).group(1).lower()
        the_data = np.fromfile(indata, dtype=np.dtype(float_type), count=num_count*4)
    
    the_data = the_data.reshape((-1,4))
    x,y,z,strain = np.hsplit(the_data,4)
    x = x[:,0]
    y = y[:,0]
    z = z[:,0]
    strain = strain[:,0]
    return x,y,z,strain

def str_to_file(x,y,z,strain,fname):
    """
    Spits out a .bin file like HiDEM would produce. Use case: extract a subset of
    data on server for download
    """
    data = np.stack((x,y,z,strain))
    count = strain.size
    with open(fname,'wb') as fout:
        fout.write("Count: "+str(count)+" Type: Float32\n")
        data.tofile(fout)

def grid_gen(x,y,z,buff=500.0,dx=60.0):
    """
    Generate an xx,yy,zz grid for gridding strain.
    This is seperate from the gridding process so the same
    grid can be easily reused
    """
    extent = [x.min()-buff,x.max()+buff,y.min(),y.max()+buff,z.min(),z.max()]

    xmin = int(math.floor(extent[0] / dx)) * dx
    xmax = int(math.ceil(extent[1] / dx)) * dx
    ymin = int(math.floor(extent[2] / dx)) * dx
    ymax = int(math.ceil(extent[3] / dx)) * dx
    zmin = int(math.floor(extent[4] / dx)) * dx
    zmax = int(math.ceil(extent[5] / dx)) * dx

    xx = np.arange(xmin,xmax+dx,dx)
    yy = np.arange(ymin,ymax+dx,dx)
    zz = np.arange(zmin,zmax+dx,dx)

    return xx,yy,zz

def grid_strain(x,y,z,strain,xx,yy,zz,buff=500.0,dx=60.0):
    """
    Given a strain at points, grid it and return gridded values
    """

    buff_boxes = math.floor(buff/dx)

    sx_plan = np.zeros([len(xx), len(yy)])
    nx_plan = np.zeros([len(xx), len(yy)])
    sx_side = np.zeros([len(yy), len(zz)])
    nx_side = np.zeros([len(yy), len(zz)])

    #find the voxel in which each particle belongs
    boxes = np.rint(np.stack((x,y,z)).T / np.array((dx,dx,dx))).astype(int)
    boxes[:,0] += int(math.floor(buff_boxes))

    #check box bounds
    boxes[:,0] = boxes[:,0].clip(0,sx_plan.shape[0]-1)
    boxes[:,1] = boxes[:,1].clip(0,sx_plan.shape[1]-1)
    boxes[:,2] = boxes[:,2].clip(0,sx_side.shape[1]-1)

    np.add.at(sx_plan, (boxes[:,0],boxes[:,1]),strain)
    np.add.at(sx_side, (boxes[:,1],boxes[:,2]),strain)
    np.add.at(nx_plan, (boxes[:,0],boxes[:,1]),np.ones(strain.size,dtype=int))
    np.add.at(nx_side, (boxes[:,1],boxes[:,2]),np.ones(strain.size,dtype=int))

    sx_plan /= nx_plan
    sx_side /= nx_side

    sx_side[np.isnan(sx_side)] = 0.0
    sx_plan[np.isnan(sx_plan)] = 0.0

    return sx_plan, sx_side
