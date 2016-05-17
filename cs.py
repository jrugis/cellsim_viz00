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
  for line in f1:
    pcount = int(line)
    break
  x = np.empty((pcount))
  y = np.empty((pcount))
  z = np.empty((pcount))
  t = 0
  for line in f1:
    v = line.split()
    x[t] = float(v[1])
    y[t] = float(v[2])
    z[t] = float(v[3])
    t += 1
    if t == pcount: break
  f1.close
  return x, y, z
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


