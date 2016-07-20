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
import os

##################################################################
# get surface nodes and triangles
def get_mesh_tris(fname):

  # get all of the nodes
  f1 = open(fname, 'r')
  for line in f1: 
    if line.startswith("$Nodes"): break
  mncount = int(f1.next())
  mnodes = np.empty((mncount, 3), dtype=np.float32)
  i = 0
  for line in f1:
    v = line.split()
    mnodes[i] = [float(v[1]), float(v[2]), float(v[3])]
    i += 1
    if i == mncount: break

  # get the triangle elements
  for line in f1: 
    if line.startswith("$Elements"): break
  ecount = int(f1.next())
  tris = np.zeros((ecount, 3), dtype = np.int32) # more than needed
  tcount = 0
  for line in f1:       # assumes all tris listed before tets
    v = line.split()
    if int(v[1]) == 2:
      tris[tcount] = [int(v[5])-1, int(v[6])-1, int(v[7])-1] # change to zero indexing
      tcount += 1
      if tcount == ecount: break
    else: break
  f1.close
  tris = tris[0:tcount, :] # trim to actual

  # copy out the triangle nodes...
  flags = np.full((mncount, 1), -2, dtype=np.int32) # "-2" = unused
  for t in tris:
    for ni in t:
      flags[ni] = -1                                # "-1" = used
  ncount = np.count_nonzero(flags==-1)
  nodes = np.empty((ncount, 3), dtype=np.float32)
  i = 0
  for j in range(mncount):
    if flags[j] == -1:
      nodes[i] = mnodes[j] # copy used triangle vertices
      flags[j] = i                                  # the new node index
      i += 1
    j += 1

  #  and update the node indices
  for i in range(tcount):
    for j in range(3):
      tris[i, j] = flags[tris[i, j]]

  return nodes, tris

##################################################################
# create normals array
def create_normals(nodes, tris):
  tcount = tris.shape[0]
  normals = np.empty((tcount, 3), dtype=np.float32)
  for i in range(tcount):
    v1 = nodes[tris[i][1]] - nodes[tris[i][0]]
    v2 = nodes[tris[i][2]] - nodes[tris[i][0]]
    temp = np.cross(v1, v2) # assumes counter-clockwise ordering
    normals[i] = temp / np.sqrt(temp.dot(temp)) # normalize
  return normals

##################################################################
# get edge and triangle fans at a node
def get_fan(node, tris):

  tri_ind = np.where(np.any(tris == node, 1))[0] # get the triangle indicies
  tri_nodes = tris[tri_ind]

  fan_size = tri_ind.shape[0]
  node_fan = np.zeros((fan_size), dtype=np.int32)
  tri_fan = np.zeros((fan_size), dtype=np.int32)

  # counter-clockwise fan identification
  tri_node = 0
  tflags = np.zeros((fan_size, 1), dtype=np.int32) # to flag triangle hits
  for i in range(fan_size):
    # store and flag the triangle 
    tri_fan[i] = tri_ind[tri_node]
    tflags[tri_node] = True

    cnodes = tri_nodes[tri_node] # get the current triangle nodes
    # find the triangle node that is clockwise from the center node...
    cwnode = cnodes[np.mod(np.where(cnodes == node)[0][0] - 1, 3)]
    node_fan[i] = cwnode #  and add it to the node fan

    # find the triangle that is counter-clockwise from the current triangle
    tpair = np.where(np.any(tri_nodes == cwnode, 1))[0] 
    if tflags[tpair[0]] == True:  # two triangles share each fan node...
      tri_node = tpair[1]         #   so pick the unused triangle
    else: 
      tri_node = tpair[0]

  return node_fan, tri_fan

##################################################################
# calculate area of triangle
def calc_area(v):
  s1 = v[1] - v[0] # edge vectors
  s2 = v[2] - v[0]
  temp = np.cross(s1, s2)
  area = np.sqrt(temp.dot(temp)) / 2.0
  return area

##################################################################
# calculate signed, weighted angle between normals
def calc_wangle(n0, n1, v0, v1, e0, e1):

  # calculate the angle between normals
  temp = np.dot(n0, n1)
  if temp > 1.0: temp = 1.0 # check precision error bound
  angle = np.arccos(temp)  

  # negative angle for concave
  temp = v1 - v0
  temp = temp / np.sqrt(temp.dot(temp))
  temp = np.dot(n1, temp)
  if temp > 1.0: temp = 1.0 # check precision error bound
  temp = np.arccos(temp)
  if temp < (np.pi / 2.0): # concave?
    angle *= -1.0

  # calculate edge weight
  temp = e1 - e0
  weight = np.sqrt(temp.dot(temp)) # length of edge vector

  return angle * weight

##################################################################
# calculate curvatures
def calc_curvature(nodes, tris, normals):
  ncount = nodes.shape[0]
  curvatures = np.empty((ncount), dtype=np.float32)

  for node in range(ncount):
    node_fan, tri_fan = get_fan(node, tris)
    fan_size = node_fan.shape[0]
    wangle_sum = 0.0
    area_sum = 0.0
    for i in range(fan_size):
      n0 = normals[tri_fan[np.mod(i - 1, fan_size)]]
      n1 = normals[tri_fan[i]] # normals
      v0 = nodes[node_fan[np.mod(i - 1, fan_size)]]  # flap verticies   
      v1 = nodes[node_fan[np.mod(i + 1, fan_size)]]      
      e0 = nodes[node]                               # edge verticies
      e1 = nodes[node_fan[i]]
      wangle_sum += calc_wangle(n0, n1, v0, v1, e0, e1)
      area_sum += calc_area(nodes[tris[tri_fan[i]]])
    curvatures[node] = 0.75 * wangle_sum / area_sum
  return curvatures

##################################################################
# write out a mesh file with node attributes 
def save_mesh(cell_name, nodes, tris, curvatures):
  fname = cell_name + "_curvature.msh"
  f1 = open(fname, 'w')
  f1.write("$MeshFormat\n")
  f1.write("2.2 0 8\n")
  f1.write("$EndMeshFormat\n")
  f1.write("$Nodes\n")
  ncount = len(nodes)
  f1.write(str(ncount) + "\n")
  for i in range(ncount):
    f1.write(str(i+1)+" "+str(nodes[i, 0])+" "+str(nodes[i, 1])+" "+str(nodes[i, 2])+"\n")
  f1.write("$EndNodes\n")
  f1.write("$Elements\n")
  ecount = len(tris)
  f1.write(str(ecount) + "\n")
  for i in range(ecount):
    f1.write(str(i+1)+" 2 2 101 100 "+str(tris[i, 0]+1)+" "+str(tris[i, 1]+1)+" "+str(tris[i, 2]+1)+"\n")
  f1.write("$EndElements\n")
  f1.write('$NodeData\n1\n"curvature"\n1\n0.0\n3\n0\n1\n')
  ccount = len(curvatures)
  f1.write(str(ccount) + "\n")
  for i in range(ccount):
    f1.write(str(i+1)+" "+str(curvatures[i])+"\n")
  f1.write("$EndNodeData\n")
  f1.close()
  return

##################################################################
# main program
##################################################################

plt.rcParams['figure.figsize'] = 6, 14
fig, plots = plt.subplots(7, sharex='col')
fig.suptitle("Surface Curvature Histograms")
plots[6].set_xlabel("surface curvature")

for i in range(1, 8):
  cell_name = "cell0" + str(i)
  mesh_name =  cell_name + "m_HARMONIC_100p.msh"

  print mesh_name
  nodes, tris = get_mesh_tris(mesh_name)
  #print "      surface node count: ", nodes.shape[0]
  #print "  surface triangle count: ", tris.shape[0]

  normals = create_normals(nodes, tris)
  curvatures = calc_curvature(nodes, tris, normals)
  print "           curvature min: ", curvatures.min()
  print "          curvature mean: ", curvatures.mean()
  print "           curvature max: ", curvatures.max()
  print "        curvature median: ", np.median(curvatures)

  plots[7-i].set_ylim([0, 450])
  plots[7-i].set_ylabel(cell_name)
  plots[7-i].hist(curvatures, bins=61, range=(-1.0, 1.0))

  save_mesh(cell_name, nodes, tris, curvatures)

  # write out to vtk file
  fname = cell_name
  d = {}
  d["curvature"] = curvatures
  nodes = np.transpose(nodes)
  pointsToVTK(fname, nodes[0,:], nodes[1,:], nodes[2,:], data = d) # write out vtk file

open('temp.pdf', 'w').close()
plt.savefig('temp.pdf')
os.rename('temp.pdf', 'curvature_histograms.pdf')
plt.show()

