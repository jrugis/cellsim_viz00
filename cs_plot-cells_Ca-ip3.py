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
plt.rcParams['axes.titlesize'] = 'small'
plt.rcParams['figure.figsize'] = 12, 15
fig, plots = plt.subplots(len(dist_names)/2, 2, sharex='col')
fig.subplots_adjust(hspace = 0.4)

for i in range(len(dist_names)/2):
  cfname = dist_names[2*i]
  ip3fname = dist_names[2*i + 1]
  vals = cfname.split('_')[2] + '_' + cfname.split('_')[3].split('.')[0]
  cell_name = cfname.split('_')[0]
  print cell_name

  plots[i, 0].set_title(cell_name + '_' + vals)

  cdata = cs.get_data(cfname)
  plots[i, 0].set_ylim([0.0, 1.0])
  plots[i, 0].set_ylabel('Ca')
  plots[i, 0].plot(np.transpose(cdata), lw=0.5)

  ip3data = cs.get_data(ip3fname)
  plots[i, 1].set_ylim([0.15, 0.35])
  plots[i, 1].set_ylabel('ip3')
  plots[i, 1].plot(np.transpose(ip3data), lw=0.5)

open('temp.pdf', 'w').close()
plt.savefig('temp.pdf')
os.rename('temp.pdf', 'seven-cells_Ca-ip3.pdf')
plt.show()

