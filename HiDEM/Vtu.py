"""
Module for interaction with .vtu files
"""

import re
import numpy as np
import math
import vtk
import paraview.simple as PV
from paraview.numpy_support import vtk_to_numpy

def data_from_vtu(fname):
    """
    Returns all particle data from a HiDEM .vtu file
    """
    
    reader = PV.XMLUnstructuredGridReader( FileName=vtu_in )
    reader.UpdatePipeline()

    data = PV.servermanager.Fetch(reader)
    GlacierSource = PV.GetActiveSource()

    points = vtk_to_numpy(data.GetPoints().GetData())
