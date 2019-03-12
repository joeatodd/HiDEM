"""
Module for strain/bond breakage related stuff. i.e. interaction with *STR* files
"""

import re
import numpy as np
import math
import os
from glob import glob
from HiDEM import Vtu

float_re = re.compile('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *[-\+]?\ *[0-9]+)?')

def bonds_from_file(fname_bin, fname_sth=None):
    """
    Returns the bond str info from given file
    """

    #Get vtu and sth filenames if not provided
    if not fname_sth:
        fname_sth = fname_sth_from_str(fname_bin)

    count_re = re.compile("Count: ([0-9]*)")
    type_re = re.compile('Type: ([a-zA-Z0-9]*)')

    indata_str = open(fname_bin, 'rb')
    header = indata_str.readline()
    num_count = int(count_re.search(header).group(1))
    float_type = type_re.search(header).group(1).lower()

    indata_sth = open(fname_sth, 'rb')
    indata_sth.readline()

    nn_data = np.fromfile(indata_sth, dtype=np.dtype(np.int32), count=num_count*2)
    efs_data = np.fromfile(indata_str, dtype=np.dtype(float_type), count=num_count)

    nn_data = nn_data.reshape((-1,2))

    return nn_data, efs_data

def strain_from_file(fname_str, fname_vtu=None, fname_sth=None):
    """
    Returns strain rates info from given file
    """

    #Get vtu and sth filenames if not provided
    if not fname_sth:
        fname_sth = fname_sth_from_str(fname_str)
    if not fname_vtu:
        fname_vtu = fname_vtu_from_str(fname_str)

    nn_data, efs = bonds_from_file(fname_str, fname_sth)
    SCL = get_scl(fname_sth)
    LNN = SCL * 1.1225

    #Get points from vtu and compute bond centrepoints
    points = Vtu.points_from_vtu(fname_vtu)

    #Do we need these?
    #cpoints = get_bond_centrepoints(points,nn_data)
    lens = get_bond_lengths(points,nn_data)
    strain = (lens-LNN)/LNN

    return nn_data, efs, strain


def sth_to_file(nn,sth_header_template,fname_sth):
    """
    Spits out a STH file like HiDEM would produce. See also str_to_file

    nn - the NAN data to write
    sth_header_template - first line of previous STH file (count will be replaced)
    fname_sth - the name of the file produced
    """
    count = nn.shape[0]

    beamcnt_re = re.compile("(?<=BeamCount: )([0-9]*)")
    sth_header = re.sub(beamcnt_re,str(count),sth_header_template)

    with open(fname_sth,'wb') as fout:
        fout.write(sth_header)
        nn.tofile(fout)

def str_to_file(efs, fname_str):
    """
    Spits out STR .bin files like HiDEM would produce. Use case: extract a subset of
    data on server for download
    """
    count = efs.size
    with open(fname_str,'wb') as fout:
        fout.write("Count: "+str(count)+" Type: Float32\n")
        efs.tofile(fout)


def get_scl(sth_file):
    """
    Return simulation SCL from sth file
    """
    scl_re = re.compile("SCL:  ("+float_re.pattern+")")

    with open(sth_file, 'rb') as f:
        header = f.readline()
        try:
            SCL = float(scl_re.search(header).group(1))
        except:
            raise Exception("Failed to get SCL info from STH header")
    return SCL

def fname_sth_from_str(fname_in, check=True):
    """
    Given a STR000*.bin filename, returns equivalent STH filename
    Fuzzy matching for restarted runs
    """
    fname_out = fname_in.rsplit("STR",1)[0]+"STH.bin"

    #If that file doesn't exist, try fuzzy matching on
    #text characters from STR.bin filename
    if (not os.path.isfile(fname_out)) and check:
        text_re = re.compile("([a-zA-Z_-]+)")
        strings = text_re.findall(fname_in)
        search_glob = "*".join(strings)
        search_glob = search_glob.replace("STR*","STH.")
        infiles = glob(search_glob)
        if len(infiles) > 1:
            raise Exception("Fuzzy matching found more than 1 possible STH file")
        else:
            return infiles[0]
    else:
        return fname_out

def fname_vtu_from_str(fname_in):
    return fname_in.replace("STR","JYR").replace(".bin",".vtu")

def pointstr_from_files(vtuf,strf,sthf=None):

    from HiDEM import Vtu

    points = Vtu.points_from_vtu(vtuf)
    nn, efs, strain = strain_from_file(strf,vtuf,sthf)
    cpoints = get_bond_centrepoints(points,nn)

    return cpoints,strain

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

def get_bond_centrepoints(points,nn):
    """
    Returns the x,y,z centrepoint of each bond specified in NN
    """

    #-1 because NN refers to particle global IDs, which start at 1
    return np.mean(points[nn[:,:]-1,:],1)

def get_bond_lengths(points, nn):
    """
    Returns the length of bonds given points, nn
    """
    bond_points = points[nn[:,:]-1,:]
    diffs = bond_points[:,0,:] - bond_points[:,1,:]
    return np.linalg.norm(diffs, axis=1)


##########################
##### LEGACY CODE ########
##########################

def str_from_file_old(fname,legacy_input=False):
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

def str_to_file_old(x,y,z,strain,fname):
    """
    Spits out a .bin file like HiDEM would produce. Use case: extract a subset of
    data on server for download
    """
    data = np.stack((x,y,z,strain))
    count = strain.size
    with open(fname,'wb') as fout:
        fout.write("Count: "+str(count)+" Type: Float32\n")
        data.tofile(fout)

