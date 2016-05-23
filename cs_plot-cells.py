import matplotlib.pyplot as plt
import numpy as np
import os
import subprocess

import cs

##################################################################
# main program
##################################################################

# plot the calcium data
dist_names = subprocess.check_output("ls *.bin", shell=True).split()
plt.rcParams['axes.color_cycle'] = ['r', 'g', 'b']
fig, plots = plt.subplots(len(dist_names), 1, sharex='col')

for i in range(len(dist_names)):
  cell_name = dist_names[i].split('_')[0]
  print cell_name
  data = cs.get_data(dist_names[i])
  max_per_row = np.amax(data, axis=1)
  rows = [np.argmax(max_per_row), np.argmin(max_per_row)]
  plots[i].set_ylabel(cell_name)
  plots[i].set_ylim([0.0, 0.8])
  plots[i].plot(np.transpose(data[rows, :]), lw=0.5)

fig.set_size_inches(5, 13)
open('temp.pdf', 'w').close()
plt.savefig('temp.pdf')
os.rename('temp.pdf', 'seven_cells.pdf')
plt.show()

