import numpy as np
import struct
import subprocess

##################################################################
# functions
##################################################################

##################################################################
# get the mesh node coordinates
def get_mesh_coords(fname):
  f1 = open(fname, 'r')
  for line in f1: 
    if line.startswith("$Nodes"): break
  ncount = int(f1.next())
  xyz = np.empty((ncount, 3))
  t = 0
  for line in f1:
    v = line.split()
    xyz[t] = [float(v[1]), float(v[2]), float(v[3])]
    t += 1
    if t == ncount: break
  f1.close
  return xyz
##################################################################
# get the mesh elements (i.e. surface triangles and tetrahedrons)
def get_mesh_elements(fname):
  f1 = open(fname, 'r')
  for line in f1: 
    if line.startswith("$Elements"): break
  ecount = int(f1.next())
  tt = np.zeros((ecount, 4), dtype = np.uint32)
  t = 0
  tricount = 0
  tetcount = 0
  for line in f1:       # assumes all tris listed before tets
    v = line.split()
    if int(v[1]) == 2:   
      tt[t] = [int(v[5])-1, int(v[6])-1, int(v[7])-1, 0] # change to zero indexing
      tricount += 1
    elif int(v[1]) == 4:
      tt[t] = [int(v[5])-1, int(v[6])-1, int(v[7])-1, int(v[8])-1] # change to zero indexing
      tetcount += 1
    t += 1
    if t == ecount: break
  f1.close
  return tt[0:tricount, 0:3], tt[tricount:, :] 
##################################################################
# get the distance to nearest lumen (dnl) data
def get_mesh_dnl(fname):
  f1 = open(fname, 'r')
  for line in f1: 
    if line.startswith('"distance to nearest lumen"'): break
  for t in range(5): # skip 5 lines
    f1.next()
  ncount = int(f1.next())
  dnl = np.empty(ncount)
  for t in range(ncount):
    v = f1.next().split()
    dnl[t] = float(v[1])
  f1.close
  return dnl
##################################################################
# read in a bin data file
def get_data(fname):
  f1 = open(fname, "rb")
  rows = struct.unpack('l', f1.read(8))[0]
  cols = struct.unpack('l', f1.read(8))[0]
  data = np.zeros((rows, cols), dtype=np.float32)
  for j in range(0, cols):     # the data is in column order
    for i in range(0, rows):
      t = f1.read(4)
      data[i,j] = struct.unpack('f', t)[0]
  f1.close()
  return data
##################################################################

