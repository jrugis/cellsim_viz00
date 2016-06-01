import numpy as np
import subprocess
import os

##################################################################
# functions
##################################################################

def get_membrane_data(fname):
  f1 = open(fname + '.msh', 'r') # open the mesh file

  # get the node count
  for line in f1: 
    if line.startswith("$Nodes"): break
  ncount = int(f1.next())
  print " ", ncount, "nodes"

  # identify the surface nodes
  sn = np.zeros(ncount, dtype = np.uint8)
  for line in f1: 
    if line.startswith("$Elements"): break
  f1.next()
  for line in f1:  
    v = line.split()
    if int(v[1]) != 2: break # finished with surface elements?
    for i in range(5, 8):    # mark the surface element nodes
      sn[int(v[i])] = 1
  sncount = 0 # surface node count
  for n in sn:
    if n == 1: sncount += 1
  print " ", sncount, "surface nodes"

  # get "distance to nearest lumen" data
  for line in f1: 
    if line.startswith('"distance to nearest lumen"'): break
  for t in range(6): # skip 6 lines
    f1.next()
  dnl = np.empty(ncount)
  for t in range(ncount):
    v = f1.next().split()
    dnl[t] = float(v[1])

  # create surface node dnl array
  sndnl = np.zeros(sncount)
  i = 0;
  for n in range(ncount):
    if sn[n] == 1:
      sndnl[i] = dnl[n]
      i += 1
  print " ", np.amin(sndnl), "minimum"
  print " ", np.amax(sndnl), "maximum"

  f1.close # close the mesh file 
  return sndnl

##################################################################
# main program
##################################################################

import matplotlib.pyplot as plt

mesh_names = subprocess.check_output("ls *.msh", shell=True).split()
plt.rcParams['axes.titlesize'] = 'small'
plt.rcParams['figure.figsize'] = 5, 15
fig, plots = plt.subplots(len(mesh_names), sharex='col', sharey='col')
fig.subplots_adjust(hspace = 0.4)
fig.suptitle('Surface membrane dnl histograms')

for i in range(len(mesh_names)):
  fname = mesh_names[i].split('.')[0]
  cell_name = fname.split('_')[0]
  print cell_name
  plots[i].set_title(cell_name)
  sndnl = get_membrane_data(fname)
  #plots[i].hist(sndnl, bins=16, range=(0.0, 8.0), normed=True, rwidth=0.9)
  plots[i].hist(sndnl, bins=16, range=(0.0, 8.0), rwidth=0.9)

plots[len(mesh_names) - 1].set_xlabel('dnl')
open('temp.pdf', 'w').close()
plt.savefig('temp.pdf')
os.rename('temp.pdf', 'seven-cells_membrane-dist.pdf')
plt.show()


