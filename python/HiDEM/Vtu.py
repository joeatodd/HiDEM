"""
Module for interaction with .vtu files
"""

import re
import numpy as np
import math
import vtk
from vtk import VTK_VERSION
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy, numpy_to_vtkIdTypeArray

def points_from_vtu(fname):
    """
    Returns particle points from a HiDEM .vtu file
    """

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()
    data = reader.GetOutput()

    points = vtk_to_numpy(data.GetPoints().GetData())
    return points


def data_from_vtu(fname):
    """
    Returns all particle data from a HiDEM .vtu file
    """

    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(fname)
    reader.Update()
    data = reader.GetOutput()

    #The point coordinates
    points = vtk_to_numpy(data.GetPoints().GetData())

    #The point data array
    pd = data.GetPointData()
    nvars = pd.GetNumberOfArrays()
    point_data = {} #dict to hold variables

    for i in range(nvars):
        name = pd.GetArray(i).GetName()
        point_data[name] = vtk_to_numpy(pd.GetArray(i))

    return points, point_data

def data_to_vtu(points, point_data, fname, dp=False):
    """
    Writes out a (presumably updated) .vtu file from HiDEM results
    points - the xyz coords
    point_data - dictionary of var names & data arrays
    """

    if dp:
        datatype = "float64"
    else:
        datatype = "float32"

    #Create the 'mesh' object
    mesh = vtk.vtkUnstructuredGrid()

    #Create & add points to mesh
    pointsVtk = vtk.vtkPoints()
    pointsVtk.SetData(numpy_to_vtk(points))
    mesh.SetPoints(pointsVtk)

    #Add all variables
    for varname, values in point_data.iteritems():
        dataVtk = numpy_to_vtk(values.real.astype(datatype), deep=True)
        dataVtk.SetName(varname)
        mesh.GetPointData().AddArray(dataVtk)

    #Set up the writer object
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(fname)
    writer.SetDataModeToAppended()
    if float(VTK_VERSION.split('.')[0]) >= 6:
        writer.SetInputData(mesh)
    else:
        writer.SetInput(mesh)

    #write it!
    writer.Write()

    #TODO - presumably might need to clean up some data after this? C memory leaks?

def vtu_sort_fun(fname):
    """
    Sorting key for non-zero padded filenames
    """
    #substitute the _integer__ with a zero padded equivalent string
    #return re.sub(fname_re, "_"+str(int(fname_re.search(fname).group(1))).zfill(4),fname)
    fname_re = re.compile("[0-9]+")
    return re.sub(fname_re, lambda x: x.group(0).zfill(4), fname)
