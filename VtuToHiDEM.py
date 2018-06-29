import numpy as np
import math
import sys
import getopt
from paraview.simple import *
#from vtk.util.numpy_support import vtk_to_numpy
from paraview.numpy_support import vtk_to_numpy
import paraview.vtk as vtk
from scipy import interpolate as interp
from scipy.interpolate import Rbf
import shapefile
from shapely.geometry import shape, Point
import shapely.affinity as shaf
import timeit

######################
### Read Arguments ###
######################

def usage():
    print("Usage: " + sys.argv[0] + " -i input.pvtu -b bed_dem_file  -p param_file -t transform_file (optional)")

try:
	opts, args = getopt.getopt(sys.argv[1:],"i:f:b:p:t:h")
except getopt.GetoptError:
	usage()
	sys.exit(2)

got_infile = False
got_beddem = False
got_transform = False
got_param = False
got_fjordfile = False

for opt, arg in opts:
	if(opt =="-i"):
		result_file = arg
		got_infile = True
	elif(opt =="-f"):
		fjord_file = arg
		got_fjordfile = True
	elif(opt =="-b"):
		bed_dem_file = arg
		got_beddem = True
	elif(opt =="-t"):
		transform_file = arg
		got_transform = True
	elif(opt =="-p"):
		param_file = arg
		got_param = True
	else:
		usage()
		sys.exit(2)

if(not(got_infile) or not(got_beddem) or not(got_param)):
	usage()
	sys.exit(2)

determine_transform = not got_transform

###################
### Parameters ####
###################


# maps label to attribute name and types
label_attr_map = {
       "front_orientation:": ["front_orientation", float, float, float],
        "upstream_extent:": [ "upstream_extent", float],
        "downstream_extent:": [ "downstream_extent", float],
        "grid_res:": [ "grid_res", float],
    "front_id:": [ "front_id", int],
    "bed_id:": [ "bed_id", int],
    "surf_id:": [ "surf_id", int]
}

class Params(object):
    def __init__(self, input_file_name):
        with open(input_file_name, 'r') as input_file:
            for line in input_file:
                if "#" in line: continue
                if line == "\n": continue
                row = line.split()
                label = row[0]
                data = row[1:]  # rest of row is data list

                attr = label_attr_map[label][0]
                datatypes = label_attr_map[label][1:]

                values = [(datatypes[i](data[i])) for i in range(len(data))]
                self.__dict__[attr] = values if len(values) > 1 else values[0]

params = Params(param_file)
print params.downstream_extent
print params.upstream_extent
print params.front_orientation
print params.front_id
print params.bed_id
print params.surf_id
print params.grid_res

#result_file = sys.argv[1]
# bed_dem_file = "BedOut.xyz"
# transform_file = "StoreCalving3D_Sik1_PlaM1_BasM1_PluM1_10__precalve0021.transform"
# result_file = "StoreCalving3D_Sik1_PlaM1_BasM1_PluM1_10__precalve0021.pvtu"


#vector normal to calving front
#front_orientation = [-0.9008556490535132, -0.43411876205523897, 0.0]
front_orientation = params.front_orientation

#the new origin
upstream_extent = params.upstream_extent
downstream_extent = params.downstream_extent
grid_res = params.grid_res

#Define the grid for HiDEM (Store front)

# #how far upstream behind the front do we want to go?
# upstream_extent = 5000.0
# #how far into the fjord do we want bathy?
# downstream_extent = 1000.0
# #required resolution
# grid_res = 50.0

#the geometry IDs of the BCs
front_id = params.front_id
bed_id = params.bed_id
surf_id = params.surf_id

#the value to replace unfound values - is this used??
bed_nanval = 0.0


###################
### Functions ####
###################

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

#stuff to compute the rotation angle
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in degrees between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    radians =  np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    return radians * (180.0/math.pi)

##############################
#### Initial Computations ####
##############################

#compute the rotation matrix
rot_mat = compute_rotation_matrix(front_orientation)

#check the vtu file extension
extension = result_file.split('.')[-1]
assert extension[1:] == 'vtu'

#set the output filename
outfile = result_file.replace(extension,"dat")
transform_outfile = result_file.replace(extension,"transform")

##############################
#### Do Paraview Things! #####
##############################

#read the vtu file(s)
if extension == 'pvtu':
    reader = XMLPartitionedUnstructuredGridReader(FileName=result_file)
else:
    reader = XMLUnstructuredGridReader(FileName=result_file)

reader.UpdatePipeline()
data = servermanager.Fetch(reader)

#transform the elmer mesh
#Paraview does rotation then translation
rotate_angle = angle_between(front_orientation, np.array((0.0,1.0,0.0)))
if(front_orientation[0] < 0): rotate_angle *= -1

if(determine_transform):
	print "determining transformation"

	transform1 = Transform(Input=reader)
	transform1.Transform = 'Transform'
	transform1.Transform.Rotate = [0.0,0.0,rotate_angle]
	trans_data = servermanager.Fetch(transform1)
	trans_points = vtk_to_numpy(trans_data.GetPoints().GetData())

	#this won't be constant between domains - is that a problem?
	y_translate = -(np.max(trans_points[:,1]) - upstream_extent)

	### Need this somewhere above...
	clip1 = Clip(Input=transform1)
	clip1.ClipType = 'Plane'
	clip1.ClipType.Origin = [0.0, -y_translate, 0.0]
	clip1.ClipType.Normal = [0.0, 1.0, 0.0]
	clipdata = servermanager.Fetch(clip1)

	clip_points = vtk_to_numpy(clipdata.GetPoints().GetData())
	x_translate = -clip_points[:,0].min()
	x_max = clip_points[:,0].max() + x_translate

	transform1.Transform.Translate = [x_translate, y_translate, 0.0]

	#determine the shape of the HiDEM grid
	grid_xrange = [0.0, x_max]
	grid_yrange = [0.0, upstream_extent + downstream_extent]

	xsteps = int(np.ceil((grid_xrange[1] - grid_xrange[0])/grid_res))
	ysteps = int(np.ceil((grid_yrange[1] - grid_yrange[0])/grid_res))
	grid_xrange[1] = grid_xrange[0] + xsteps*grid_res
	grid_yrange[1] = grid_yrange[0] + ysteps*grid_res

	mat_shape = (ysteps+1, xsteps+1)


else:
	print "using existing transformation"

	#already determined transformation?
	#not implemented yet but I suspect I might need it...
	#Potential to write out the translation/rotation and 
	#then read it back in and use it here
	rot_mat = np.empty((3,3))
	grid_xrange = np.empty(2)
	grid_yrange = np.empty(2)

	with open(transform_file, 'r') as transformin:
		transformin.readline()
		transformin.readline()
		for i in range(3):
			rotline = transformin.readline().split(" ")[:-1]
			rot_mat[i,:] = [float(j) for j in rotline]

		transformin.readline()
		transformin.readline()
		translate_line = transformin.readline().split(" ")[:-1]
		x_translate, y_translate, z_translate = [float(j) for j in translate_line]
		transformin.readline()
		transformin.readline()
		bed_nanval = float(transformin.readline())
		transformin.readline()
		transformin.readline()
		grid_line = transformin.readline().split(" ")
		grid_xrange[0], grid_xrange[1], grid_res, grid_yrange[0], grid_yrange[1], grid_res = [float(j) for j in grid_line]


	xsteps = int((grid_xrange[-1] - grid_xrange[0]) / grid_res)
	ysteps = int((grid_yrange[-1] - grid_yrange[0]) / grid_res)
	mat_shape = (ysteps+1, xsteps+1)

	rot_mat = np.matrix(rot_mat)

	transform1 = Transform(Input=reader)
	transform1.Transform = 'Transform'
	transform1.Transform.Rotate = [0.0,0.0,rotate_angle]
	transform1.Transform.Translate = [x_translate, y_translate, 0.0]


#create the plane mesh for sampling
SamplePlane = Plane()
SamplePlane.Origin = [grid_xrange[0],grid_yrange[0],0.0]
SamplePlane.Point1 = [grid_xrange[1],grid_yrange[0],0.0]
SamplePlane.Point2 = [grid_xrange[0],grid_yrange[1],0.0]
SamplePlane.XResolution = xsteps
SamplePlane.YResolution = ysteps

#read in the bed dem and filter -9999 (nodata)
buf = 500.0
bed_dem = np.genfromtxt(bed_dem_file)
#rotate bed dem
bed_dem = np.asarray(bed_dem * rot_mat)
#translate bed dem
bed_dem[:,0] += x_translate
bed_dem[:,1] += y_translate

bed_dem = bed_dem[bed_dem[:,2] != -9999]
bed_dem = bed_dem[bed_dem[:,0] > (grid_xrange[0] - buf)]
bed_dem = bed_dem[bed_dem[:,0] < (grid_xrange[1] + buf)]
bed_dem = bed_dem[bed_dem[:,1] > (grid_yrange[0] - buf)]
bed_dem = bed_dem[bed_dem[:,1] < (grid_yrange[1] + buf)]

#extract the boundary dataf
frontthresh = Threshold(Input=transform1)
frontthresh.Scalars = ['CELLS', 'GeometryIds']
frontthresh.ThresholdRange = [100+front_id, 100+front_id]
frontdata = servermanager.Fetch(frontthresh)

calculatorf = Calculator(Input=frontthresh)
calculatorf.Function = '-1'
calculatorf.ResultArrayName = 'gl'

bedthresh = Threshold(Input=transform1)
bedthresh.Scalars = ['CELLS', 'GeometryIds']
bedthresh.ThresholdRange = [100+bed_id, 100+bed_id]
beddata = servermanager.Fetch(bedthresh)

calculatorb = Calculator(Input=bedthresh)
calculatorb.Function = '((groundedmask>-0.5)*2) - 1'
calculatorb.ResultArrayName = 'gl'

surfthresh = Threshold(Input=transform1)
surfthresh.Scalars = ['CELLS', 'GeometryIds']
surfthresh.ThresholdRange = [100+surf_id, 100+surf_id]
surfdata = servermanager.Fetch(surfthresh)

#group the bed and front for ice surface extraction
bedfront = GroupDatasets(Input=[calculatorf, calculatorb])
bedfrontdata = servermanager.Fetch(bedfront)

#flatten base/front group
calculator1 = Calculator(Input=bedfront)
calculator1.Function = 'iHat + jHat + kHat*coordsZ'
calculator1.ResultArrayName = 'zvec'

flatbase = WarpByVector(Input=calculator1)
flatbase.Vectors = ['POINTS', 'zvec']
flatbase.ScaleFactor = -1.0

#and resample onto plane
print "Resampling base onto sample plane"
baseresamp = ResampleWithDataset(Input=flatbase,
    Source=SamplePlane)

#flatten ice surface
calculator2 = Calculator(Input=surfthresh)
calculator2.Function = 'kHat*coordsZ'
calculator2.ResultArrayName = 'zvec'

flatsurf = WarpByVector(Input=calculator2)
flatsurf.Vectors = ['POINTS', 'zvec']
flatsurf.ScaleFactor = -1.0

#and resample onto plane
print "Resampling surf onto sample plane"
surfresamp = ResampleWithDataset(Input=flatsurf,
    Source=SamplePlane)

#fetch base data
baseresampdata = servermanager.Fetch(baseresamp)
baseplane_points = vtk_to_numpy(baseresampdata.GetPoints().GetData())
baseplane_z = vtk_to_numpy(baseresampdata.GetPointData().GetArray("zvec"))[:,-1]
baseplane_slip = vtk_to_numpy(baseresampdata.GetPointData().GetArray("slip coefficient 2"))
baseplane_mask = vtk_to_numpy(baseresampdata.GetPointData().GetArray("vtkValidPointMask"))
baseplane_grounded = vtk_to_numpy(baseresampdata.GetPointData().GetArray("gl"))

#fetch surf data
surfresampdata = servermanager.Fetch(surfresamp)
surfplane_points = vtk_to_numpy(surfresampdata.GetPoints().GetData())
surfplane_z = vtk_to_numpy(surfresampdata.GetPointData().GetArray("zvec"))[:,-1]
surfplane_mask = vtk_to_numpy(surfresampdata.GetPointData().GetArray("vtkValidPointMask"))

print "Finished resampling onto plane."

#the mask for valid glacier points
glac_mask = (baseplane_mask == 1) & (surfplane_mask == 1)


# # Bed DEM strategy:

# 1) Just set bed directly from DEM interp
# 2) Then, at grounded ice points, set bed = base
# 3) Finally ensuring that no ice sticks thru bed

#MODIFIED FOR FJORD

##############################
#Read Polygon defining fjord #
##############################

if got_fjordfile:
    fjord_shapefile = shapefile.Reader(fjord_file)
    fjord_poly = fjord_shapefile.shapeRecords()[0].shape.__geo_interface__

    fjord_shapely = shape(fjord_poly) #convert to shapely format

    #rotate
    fjord_shapely = shaf.rotate(fjord_shapely,rotate_angle,origin=(0,0)) 
    #translate
    fjord_shapely = shaf.translate(fjord_shapely,x_translate, y_translate)

    print 'fjord centroid: ',fjord_shapely.centroid.xy

    #Use a np vectorized function to find points in poly
    def shapely_pip(x, y, poly):
        return Point((x,y)).within(poly)
    shapely_pipv = np.vectorize(shapely_pip)
    fjord_mask = shapely_pipv(baseplane_points[:,0],baseplane_points[:,1],fjord_shapely)



    #ensure no fjord mask where glac mask:
    fjord_mask[glac_mask] = False

    outside_mask = (~glac_mask) & (~fjord_mask)

else:

    outside_mask = ~glac_mask

## Interpolate the bed from DEM
gridbedinterp = interp.griddata(bed_dem[:,:-1], bed_dem[:,-1], (baseplane_points[:,0], baseplane_points[:,1]))
final_bed = gridbedinterp.copy()
#2) grounded points get ice height for bed
final_bed[baseplane_grounded > 0] = baseplane_z[baseplane_grounded > 0]
#3) ensure no ice sticks through bed
final_bed[(final_bed > baseplane_z) & (baseplane_mask == 1)] = baseplane_z[(final_bed > baseplane_z) & (baseplane_mask == 1)]


#Set bed outside glacier/fjord to high value
#BAD STRATEGY - to be replaced
if False:
    high_bed = 1.0E6
    final_bed[outside_mask] = high_bed

#the minimum bed height
z_translate = -np.nanmin(final_bed)

#Where there's no ice, set surf=base=bed
final_base = baseplane_z.copy()
final_base[~glac_mask] = final_bed[~glac_mask]

final_surf = surfplane_z.copy()
final_surf[~glac_mask] = final_bed[~glac_mask]

#in the fjord, set melange if desired
do_melange = False
if do_melange and got_fjordfile:
    final_base[fjord_mask] = -50.0
    final_surf[fjord_mask] = 5.0

final_slip = baseplane_slip.copy()
#scale the slip coeff from MPa-metres-years to SI
final_slip = final_slip / (1.0E-6/(365.25*24*60*60))

final_slip[final_slip < 0.0] = -final_slip[final_slip < 0.0]

if False:
    #set friction outside glacier/fjord area to high value
    #Part of bad strategy from before
    final_slip[outside_mask] = 1.0E10

#translate Z
final_bed += z_translate
final_base += z_translate
final_surf += z_translate

#fill nan with no data
final_bed[np.isnan(final_bed)] = bed_nanval
final_base[np.isnan(final_base)] = bed_nanval
final_surf[np.isnan(final_surf)] = bed_nanval
final_surf[np.isnan(final_slip)] = 0.0

#check the data
assert((final_surf >= final_base).all())
assert((final_base >= final_bed).all())
assert((final_slip >= 0.0).all())


###Output: x, y, ice_surf, ice_base, bed, slip

with open(outfile, 'wb') as output:
	output.write(str(len(baseplane_points))+"\n")
	for i in range(len(baseplane_points)):
		outstr = "%.9f" % baseplane_points[i,0] + "\t"
		outstr += "%.9f" % baseplane_points[i,1] + "\t"
		outstr += "%.9f" % final_surf[i] + "\t"
		outstr += "%.9f" % final_base[i] + "\t"
		outstr += "%.9f" % final_bed[i] + "\t"
		outstr += "%.9f" % final_slip[i]
		outstr += "\n"
		output.write(outstr)

#write out the transformation info too
if(determine_transform):
	with open(transform_outfile, 'wb') as output:
		output.write("Transformation Info\n")
		output.write("Rotation Matrix:\n")
		for i in range(3):
			for j in range(3):
				output.write("%.9f " % rot_mat[i,j])
			output.write("\n")
		output.write("\n")
		output.write("X, Y, Z translation:\n")
		output.write("%.9f" % x_translate + " " + "%.9f" % y_translate + " " + "%.9f" % z_translate + " \n\n")
		output.write("No Data Value: \n")
		output.write("%.9f" % bed_nanval)
		output.write("\n\n")
		output.write("HiDEM Grid Spec: xmin,xmax,xstep,ymin,ymax,ystep\n")
		output.write("%.6f" % grid_xrange[0] + " " + "%.6f" % grid_xrange[1] + " " + "%.6f" % grid_res + " " + "%.6f" % grid_yrange[0] + " " + "%.6f" % grid_yrange[1] + " " + "%.6f" % grid_res + "\n\n")
