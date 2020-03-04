import numpy as np

#modified from Icebergs.py for Jan's coordinate system
def compute_rotation_matrix(plane_vector):
    """
    Given a 3D vector describing the normal direction of a glacier front, return
    the rotation matrix which will convert this glacier domain into HiDEM format 
    (flow purely in y direction)
    
    Assuming that the z component of input vector is zero, should always return 
    matrix for which new z component is (0,0,1) i.e. untransformed in this 
    direction
    """
    plane_vector = np.asarray(plane_vector)
    assert plane_vector.size == 3, "Input not a 3D vector!"

    ex = np.empty(3)
    ey = np.empty(3)
    ez = np.empty(3)

    #normalise
    plane_vector = plane_vector / np.linalg.norm(plane_vector)

    ey = np.asarray(plane_vector)

    MaxIndex = np.argmax(abs(ey))

    #Odd approach here to catch an edge case (0,1,0), (1,0,0)
    #Basically, we do argmin but take last occurrence (z) if two are equal
    #MinIndex = np.argmin(abs(ey))
    MinIndex = 2 - np.argmin(abs(ey)[::-1])
    
    for i in range(3):
        if (i == MaxIndex or i == MinIndex):
            continue
        MidIndex = i

    ez[MinIndex] = 1.0
    ez[MidIndex] = 0.0

    ez[MaxIndex] = -ey[MinIndex]/ey[MaxIndex]
    ez = ez / (sum(ez ** 2) ** 0.5)

    #The new y-axis is orthogonal to new x and z axes
    ex = -np.cross(ez, ey)
    ex = ex / (sum(ex ** 2)**0.5)

    rot_mat = np.empty((3,3))
    rot_mat[0,:] = ex
    rot_mat[1,:] = ey
    rot_mat[2,:] = ez
    rot_mat = rot_mat.T

    return np.matrix(rot_mat)

class HiDEMTransformation:
    def __init__(self, filename=None):
        if not filename:
            raise Exception("Attempted to intialise a HiDEMTransformation object without file!")

        self.rot_mat = np.empty((3,3))
        self.grid_xrange = np.empty(2)
        self.grid_yrange = np.empty(2)

        with open(filename, 'r') as transformin:
            transformin.readline()
            transformin.readline()
            transformin.readline()
            for i in range(3):
                rotline = transformin.readline().split(" ")[:-1]
                self.rot_mat[i,:] = [float(j) for j in rotline]

            transformin.readline()
            transformin.readline()
            translate_line = transformin.readline().split(" ")[:-1]
            self.translation = np.asarray(translate_line,dtype=np.float64)

            ## Don't really need from here...

            transformin.readline()
            transformin.readline()
            self.bed_nanval = float(transformin.readline())
            transformin.readline()
            transformin.readline()
            grid_line = transformin.readline().split(" ")

            self.grid_xrange[0], self.grid_xrange[1], \
                self.grid_res, self.grid_yrange[0], \
                self.grid_yrange[1], self.grid_res = [float(j) for j in grid_line]


        # xsteps = int((grid_xrange[-1] - grid_xrange[0]) / grid_res)
        # ysteps = int((grid_yrange[-1] - grid_yrange[0]) / grid_res)
        # mat_shape = (ysteps+1, xsteps+1)

        self.rot_mat = np.matrix(self.rot_mat)
        self.unrot_mat = self.rot_mat.T
