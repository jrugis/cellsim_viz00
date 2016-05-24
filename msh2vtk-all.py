import numpy as np
import subprocess
from evtk.hl import pointsToVTK

##################################################################
# functions
##################################################################

def convert_mesh(fname):
  f1 = open(fname + '.msh', 'r') # open the mesh file

  # get the mesh coordinates
  for line in f1: 
    if line.startswith("$Nodes"): break
  pcount = int(f1.next())
  xyz = np.empty([3, pcount])
  for t in range(pcount):
    v = f1.next().split()
    xyz[0,t] = float(v[1])
    xyz[1,t] = float(v[2])
    xyz[2,t] = float(v[3])

  # get "distance to nearest lumen" data
  for line in f1: 
    if line.startswith('"distance to nearest lumen"'): break
  for t in range(6): # skip 6 lines
    f1.next()
  dnl = np.empty(pcount)
  for t in range(pcount):
    v = f1.next().split()
    dnl[t] = float(v[1])

  f1.close # close the mesh file 

  # write out to vtk file
  d = {}
  d["dnl"] = dnl
  pointsToVTK(fname, xyz[0,:], xyz[1,:], xyz[2,:], data = d) # write out vtk file

##################################################################
# main program
##################################################################

mesh_names = subprocess.check_output("ls *.msh", shell=True).split()
for mesh in mesh_names:
  fname = mesh.split('.')[0]
  print fname
  convert_mesh(fname)

