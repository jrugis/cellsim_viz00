# -*- coding: utf-8 -*-
#
# bin2vtk-all.py
# Convert calcium concentration simulation output to vtk time series for ParaView.
#
# J Rugis
# 16.05.16
#

import numpy as np
import subprocess
from evtk.hl import pointsToVTK
import cs

##################################################################
# main program
##################################################################

dist_names = subprocess.check_output("ls *.bin", shell=True).split()
print dist_names

#d = {}
#for i in range(len(dist_names)):
#  fname = dist_names[i]
#  dist = cs.get_data(fname)[:,0]
#  d[fname.split('.')[0]] = dist

# set the simulation time scale
#i_start = 828
#i_end = 1811

# get the mesh coordinates
#x,y,z = cs.get_mesh_coords(meshfname + ".msh")

# get calcium data
#c = cs.get_data(cfname + ".bin")

# write vtk time series files
#for i in xrange(i_start, i_end, 1):
#  fname = fdir + "cell" + str(i)
#  d = {}
#  d["c"] = c[:, i]
#  pointsToVTK(fname, x, y, z, data = d)

