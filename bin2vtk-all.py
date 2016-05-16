# -*- coding: utf-8 -*-
#
# bin2vtk-all.py
# Convert calcium concentration simulation output to vtk time series for ParaView.
#
# J Rugis
# 16.05.16
#

import numpy as np
import os
import subprocess
from evtk.hl import pointsToVTK
import cs

##################################################################
# main program
##################################################################

dist_names = subprocess.check_output("ls *.bin", shell=True).split()

for i in range(len(dist_names)):
  dist_name = dist_names[i]
  cell_name = dist_name.split('_')[0]
  mesh_name =  cell_name + 'm_HARMONIC_100p.msh'
  print '*** ' + cell_name + ' ***'
  print 'reading mesh file: ' + mesh_name 
  x,y,z = cs.get_mesh_coords(mesh_name)
  print 'reading data file: ' + dist_name 
  dist = cs.get_data(dist_name)

  # get start and finish indices for a single cycle
  i_start = 100
  i_finish = 120

  # write vtk time series files
  if os.path.isdir(cell_name):
    os.system("rm -rf " + cell_name)
  os.mkdir(cell_name)
  print 'creating vtk files', 
  for j in xrange(i_start, i_finish, 1):
    fname = cell_name + '/' + cell_name + '_' + str(j)   
    print '.',
    d = {}
    d["c"] = dist[:, j]
    pointsToVTK(fname, x, y, z, data = d)
  print 

