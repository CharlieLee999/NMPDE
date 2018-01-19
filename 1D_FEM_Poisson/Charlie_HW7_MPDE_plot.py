# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 15:37:22 2018

@author: charlie
"""
### Plot 

import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
 
N = [100, 1000]

line_color = ['r', 'b']
title_file = 'data_N_'

fig, ax = plt.subplots(2, sharex = True)
plt.xlabel('x')
plt.ion()

def Plot(N_in):
    flag = int(np.log10(N_in)-2)
    df  = pd.read_csv(title_file + str(N_in) + '.csv', sep = ',', header = None, error_bad_lines=False)
    data = df.values[1:, 1:]
    ax[0].plot(data[0, :], data[1, :], line_color[flag] + '--', label = 'Numerical N = ' + str(N_in))
    ax[0].plot(data[0, :], data[2, :], line_color[flag] + '-', label = 'Theoretical N = ' + str(N_in))
    ax[1].plot(data[0, :], data[3, :], line_color[flag], label = 'Approxi Errors N = ' + str(N_in))
    
Plot(N[0])
Plot(N[1])

ax[0].legend(loc = 'best')
ax[1].legend(loc = 'best')
ax[0].set_ylabel('$u(x)$')
ax[1].set_ylabel('$e(x)$')
fig.text(0.5, 1, '$u(x)$ and $e(x)$ under Different grids', \
         verticalalignment = 'top', \
         horizontalalignment = 'center', fontsize=16)

fig.savefig('u_x_e_x_2_grids.png', format = 'png')

n_in = 1000
df  = pd.read_csv(title_file + str(n_in) + '.csv', sep = ',', header = None, error_bad_lines=False)
data = df.values[1:, 1:]
plt.figure()
plt.plot(data[0, :], data[3, :] )
plt.xlabel('$x$')
plt.ylabel('$e(x)$')
plt.title('Approximation errors N = ' + str(n_in))
plt.savefig('Approximation_errors_N_' + str(n_in) + '.png', format = 'png')
