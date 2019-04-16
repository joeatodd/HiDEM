"""
Script for cluster analysis of HiDEM results. Appends cluster data to point data in .vtu (optionally in-place)

Needs the igraph module, installable with pip:
   pip install python-igraph
The package depends on libxml2-dev

"""
import sys
import getopt
import numpy as np
from HiDEM import Vtu, Strain
import igraph
from glob import glob

def usage():
    print ("Usage: " + sys.argv[0] + " -i input_file_glob -v strain_thresh -o [overwrite existing files]")

try:
    opts, args = getopt.getopt(sys.argv[1:],"i:v:o")
except getopt.GetoptError:
    usage()
    sys.exit(2)

strain_thresh = 1.0e-4
inplace = False

for opt, arg in opts:
    if(opt =="-i"):
        inglob = arg
    elif(opt =="-v"):
        strain_thresh = float(arg)
    elif(opt=="-o"):
        inplace = True
    else:
        print "help"
        usage()
        sys.exit(2)

#Find the vtu files
vtu_inglob = inglob+"*.vtu"
vtu_infiles = glob(vtu_inglob)
vtu_infiles.sort()

#Get the STR & STH filenames
str_infiles = [Strain.fname_str_from_vtu(f) for f in vtu_infiles]
sth_infile = Strain.fname_sth_from_str(str_infiles[0])

SCL = Strain.get_scl(sth_infile)
LNN = SCL * 1.1225

nn, efs0 = Strain.bonds_from_file(str_infiles[0])

for i in range(len(vtu_infiles)):

    vtu_in = vtu_infiles[i]
    str_in = str_infiles[i]

    print("Processing clusters for "+vtu_in)

    points, pointdata = Vtu.data_from_vtu(vtu_in)
    npoints = points.shape[0]
    efs = Strain.efs_from_file(str_in)

    lens = Strain.get_bond_lengths(points,nn)
    strain = (lens-LNN)/LNN

    #Create the graph
    #Start by filtering nn where EFS < 1.0

    broken = (efs < 0.01) & (strain > strain_thresh)
    print("broken: "+str(np.sum(broken)))
    nn_joined = nn[~broken]

    graph = igraph.Graph()
    graph.add_vertices(npoints)
    graph.add_edges(nn_joined-1) #0 indexed

    #Have confirmed by manual inspection of efs & nn that this is working
    cluster_output = graph.clusters()
    clust = np.asarray(cluster_output.membership)
    pointdata["Cluster"] = clust

    #TODO - what other/alternative output would be useful?
    #- Should we mark the cluster size (i.e. iceberg size) for
    #  visualisation?
    #
    #- Should the clusters be size ordered? This would lend some
    #  temporal consistency to the results
    #
    #- Is colouring possible? I doubt it because the graph/cluster
    #  has no knowledge of the cartesian space
    #
    #- What stats does Jan want?

    if(inplace):
        outfname = vtu_in
    else:
        outfname = "_clust_".join(vtu_in.split("_JYR"))

    Vtu.data_to_vtu(points, pointdata, outfname) #TODO
