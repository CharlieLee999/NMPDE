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
max_inter = np.zeros( 2 )
max_error = np.zeros( 2 )

line_color = ['r', 'b']
title_file = 'data_N_'

fig, ax = plt.subplots(2, sharex = True)
plt.xlabel('x')
plt.ion()


def Plot_Get_error(N_in):
    flag = int(np.log10(N_in)-2)
    df  = pd.read_csv(title_file + str(N_in) + '.csv', sep = ',', header = None, error_bad_lines=False)
    data = df.values[1:, 1:]
    
    error = data[3, 1:]
    inter = data[0, 1:] - data[0, 0:-1]
    error_norm = error**2
    error_norm_l0 = error_norm.dot(inter.T)
    ax[0].plot(np.log10(data[0, :]), np.log10(data[1, :]), line_color[flag] + '--', label = 'Numerical N = ' + str(N_in))
    ax[0].plot(np.log10(data[0, :]), np.log10(data[2, :]), line_color[flag] + '-', label = 'Theoretical N = ' + str(N_in))
    ax[1].plot(np.log10(data[0, :]), np.log10(data[3, :]), line_color[flag], label = 'Approxi Errors N = ' + str(N_in))
    max_interval = max( inter )
#    max_error 	 = max( error )
    return max_interval, error_norm_l0
    
max_inter[0], max_error[0] = Plot_Get_error(N[0])
max_inter[1], max_error[1] = Plot_Get_error(N[1])

ax[0].legend(loc = 'best')
ax[1].legend(loc = 'best')
ax[0].set_ylabel('$u(x)$')
ax[1].set_ylabel('$e(x)$')
fig.text(0.5, 1, '$u(x)$ and $e(x)$ under Different grids', \
         verticalalignment = 'top', \
         horizontalalignment = 'center', fontsize=16)
#
#fig.savefig('u_x_e_x_2_grids.png', format = 'png')

fig_log = plt.figure()
ax = fig_log.add_subplot(111)
plt.plot(np.log10(max_inter), np.log10(max_error))
plt.scatter(np.log10(max_inter), np.log10(max_error), c = 'k')
plt.ylabel('$log_{10} (||e(x)||_{0, \Omega})$')
plt.xlabel('$log_{10} (h) $ ')
plt.title( '$log_{10} (h) - log_{10} (||e(x)||_{0, \Omega})$' )

ax.text(-2, -9.5, u'$N = 1000$')
ax.text(-1, -5, u'$N = 100$')
plt.savefig('log-log.eps', format = 'eps')

#n_in = 1000
#df  = pd.read_csv(title_file + str(n_in) + '.csv', sep = ',', header = None, error_bad_lines=False)
#data = df.values[1:, 1:]
#plt.figure()
#plt.plot(data[0, :], data[3, :] )
#plt.xlabel('$x$')
#plt.ylabel('$e(x)$')
#plt.title('Approxi Errors N = ' + str(n_in))
#plt.savefig('Approximation_errors_N_' + str(n_in) + '.png', format = 'png')
