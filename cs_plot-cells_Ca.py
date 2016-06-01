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
plt.rcParams['figure.figsize'] = 6, 15
fig, plots = plt.subplots(len(dist_names)/2, sharex='col')
fig.subplots_adjust(hspace = 0.4)

for i in range(len(dist_names)/2):
  cfname = dist_names[2*i]
  vals = cfname.split('_')[2] + '_' + cfname.split('_')[3].split('.')[0]
  cell_name = cfname.split('_')[0]
  print cell_name

  plots[i].set_title(cell_name + '_' + vals)

  cdata = cs.get_data(cfname)
  plots[i].set_ylim([0.0, 0.8])
  plots[i].set_ylabel('Ca')
  plots[i].plot(np.transpose(cdata), lw=0.5)

open('temp.pdf', 'w').close()
plt.savefig('temp.pdf')
os.rename('temp.pdf', 'seven-cells_Ca.pdf')
plt.show()

