import numpy as np
import subprocess
import os

import cs

##################################################################
# functions
##################################################################

def get_membrane_data(fname):
  xyz = cs.get_mesh_coords(fname + '.msh')          # get the node data
  tris, tets = cs.get_mesh_elements(fname + '.msh') # get the element data
  dnl = cs.get_mesh_dnl(fname + '.msh')             # get "distance to nearest lumen" data
  ncount = xyz.shape[0]

  # calculate surface area
  area = 0.0
  for tri in tris:
    v = np.vstack((xyz[tri[0]], xyz[tri[1]], xyz[tri[2]]))
    area += np.abs(np.linalg.det(np.array((v[0]-v[2],v[1]-v[2],(1.0,1.0,1.0)))))/2.0
 
  # calculate volume
  volume = 0.0
  for tet in tets:
    v = np.vstack((xyz[tet[0]], xyz[tet[1]], xyz[tet[2]], xyz[tet[3]]))
    volume += np.abs(np.linalg.det((np.array((v[0]-v[3],v[1]-v[3],v[2]-v[3])))))/6.0

  # identify the surface nodes
  sn = np.zeros(ncount, dtype = np.uint8)
  for v in tris: # mark the surface element nodes
    for i in range(3):
      sn[v[i]] = 1
  sncount = 0 # surface node count
  for n in sn:
    if n == 1: sncount += 1

  # create surface node dnl array
  sndnl = np.zeros(sncount)
  i = 0;
  for n in range(ncount):
    if sn[n] == 1:
      sndnl[i] = dnl[n]
      i += 1

  extents = np.vstack((np.amin(xyz, axis=0), np.amax(xyz, axis=0)))
  return ncount, area, volume, sndnl, extents

##################################################################
# main program
##################################################################

import matplotlib.pyplot as plt
np.set_printoptions(precision=2)

mesh_names = subprocess.check_output("ls *.msh", shell=True).split()
plt.rcParams['axes.titlesize'] = 'small'
plt.rcParams['figure.figsize'] = 5, 15
fig, plots = plt.subplots(len(mesh_names), sharex='col', sharey='col')
fig.subplots_adjust(hspace = 0.4)
fig.suptitle('Surface membrane dnl histograms')

fmin = np.finfo(np.float64).max
vmin = (fmin,fmin,fmin)
fmax = np.finfo(np.float64).min
vmax = (fmax,fmax,fmax)
bbox = np.array((vmin,vmax))

for i in range(len(mesh_names)):
  fname = mesh_names[i].split('.')[0]
  cell_name = fname.split('_')[0]
  print cell_name
  plots[i].set_title(cell_name)
  ncount, area, volume, sndnl, extents = get_membrane_data(fname)
  sncount = len(sndnl)
  bbox[0] = np.minimum(bbox[0], extents[0])
  bbox[1] = np.maximum(bbox[1], extents[1])

  print " ", ncount, "nodes"
  print "  %.1f volume (cubic um3)" % volume
  print "  %.1f nodes per volume (nodes / cubic um3)" % (ncount / volume)
  print " ", sncount, "surface nodes"
  print "  %.1f surface area (square um)" % area
  print "  %.1f nodes per surface area (nodes / square um)" % (sncount / area)
  print "  %.2f surface node dnl minimum (um)" % np.amin(sndnl)
  print "  %.2f surface node dnl maximum (um)" % np.amax(sndnl)
  print "  %.2f surface node dnl average (um)" % np.mean(sndnl)
  print "  %.1f surface node percentage of total nodes" % (100.0 * sncount / ncount)
  app = 100.0 * np.histogram(sndnl, bins=8, range=(0.0, 8.0))[0][0] / sncount
  print "  %.1f apical percentage of surface nodes (1um dnl cutoff)" % app

  #plots[i].hist(sndnl, bins=16, range=(0.0, 8.0), normed=True, rwidth=0.9)
  plots[i].hist(sndnl, bins=16, range=(0.0, 8.0), rwidth=0.9)

print "coordinates bounding box (um):" 
print bbox
bboxc = (bbox[1] + bbox[0]) / 2.0
print "coordinates bounding box center (um):" 
print bboxc

plots[len(mesh_names) - 1].set_xlabel('dnl')
open('temp.pdf', 'w').close()
plt.savefig('temp.pdf')
os.rename('temp.pdf', 'seven-cells_membrane-dist.pdf')
plt.show()


