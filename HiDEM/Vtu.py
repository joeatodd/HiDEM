"""
Module for interaction with .vtu files
"""

import re
import numpy as np
import math
import paraview.simple as PV
from paraview.numpy_support import vtk_to_numpy

def points_from_vtu(fname):
    """
    Returns all particle positions from a HiDEM .vtu file
    """
    
    reader = PV.XMLUnstructuredGridReader( FileName=fname )
    reader.UpdatePipeline()

    data = PV.servermanager.Fetch(reader)
    GlacierSource = PV.GetActiveSource()

    points = vtk_to_numpy(data.GetPoints().GetData())
    #TODO - clear out paraview datastructures here
    #TODO - optionally return point data too
    return points
