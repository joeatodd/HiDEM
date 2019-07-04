import numpy as np

#modified from Icebergs.py for Jan's coordinate system
def compute_rotation_matrix(plane_vector):

	ex = np.empty(3)
	ey = np.empty(3)
	ez = np.empty(3)

	#normalise
	plane_vector = plane_vector / np.linalg.norm(plane_vector)

	ey = np.asarray(plane_vector)

	MaxIndex = np.argmax(abs(ey))
	MinIndex = np.argmin(abs(ey))

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
