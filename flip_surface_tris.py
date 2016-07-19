# -*- coding: utf-8 -*-
#
# curvature.py
# Calculate surface curvature at each surface node.
#
# J Rugis
# 04.07.16
#

import numpy as np
import matplotlib.pyplot as plt
from pyevtk.hl import pointsToVTK

##################################################################
# flip surface triangles
def flip_tris(fname):
  f1 = open(fname + "R.msh", 'r')
  f2 = open(fname + ".msh", 'w')
  flag = 0
  for line in f1: 
    if flag == 2:            # surface element?
      v = line.split()
      if int(v[1]) == 2:     # triangle element?
        line = v[0]+" "+v[1]+" "+v[2]+" "+v[3]+" "+v[4]+" "+v[5]+" "
        line += v[7]+" "+v[6]+"\n" # reverse the triangle node order
      else: flag = 3         # no more triangles
    elif flag == 1: flag = 2 # skip the element count
    elif line.startswith("$Elements"): flag = 1 # elements?
    f2.write(line)
  f2.close
  f1.close
  return

##################################################################
# main program
##################################################################

cell_name = "cell04"
mesh_name =  cell_name + "m_HARMONIC_100p"

print mesh_name
flip_tris(mesh_name)

