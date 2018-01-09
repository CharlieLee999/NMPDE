# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 20:09:35 2018

@author: charlie
"""

### Plotting code 

import pandas as pd  
import numpy as np
import matplotlib.pyplot as plt

flag_flux_scalar = [1, 2]
N_th = 1000
N_hun  = 100

label_th = ', N=1000'
label_hun = ', N=100'

line_types =  ['-.', ':', '-', '--']

fig, ax = plt.subplots(3, sharex = True)
plt.xlabel('x')

### read data from csv files to Dataframes
def Plot(flag_flux_scalar):
    flag = flag_flux_scalar
    title_th = []
    for s in ['density', 'velocity', 'mach_number', 'pressure', 'time']:
        title_temp = 'x-t-' + s + '_flux_' + str(flag_flux_scalar) + '_N_'+ str(N_th)    
        title_th.append(title_temp)
    
    title_hun = []
    for s in ['density', 'velocity', 'mach_number', 'pressure', 'time']:
        title_temp = 'x-t-' + s + '_flux_' + str(flag_flux_scalar) + '_N_'+ str(N_hun)    
        title_hun.append(title_temp)
        
    df_den_th   = pd.read_csv(title_th[0] + '.csv', sep = ',' ,header = None)
    df_vel_th   = pd.read_csv(title_th[1] + '.csv', sep = ',' ,header = None)
#    df_mach_th  = pd.read_csv(title_th[2] + '.csv', sep = ',' ,header = None)
    df_press_th = pd.read_csv(title_th[3] + '.csv', sep = ',' ,header = None)
    
    df_den_hun   = pd.read_csv(title_hun[0] + '.csv', sep = ',' ,header = None)
    df_vel_hun   = pd.read_csv(title_hun[1] + '.csv', sep = ',' ,header = None)
#    df_mach_hun  = pd.read_csv(title_hun[2] + '.csv', sep = ',' ,header = None)
    df_press_hun = pd.read_csv(title_hun[3] + '.csv', sep = ',' ,header = None)
    
    
    ### Extrate time, x, density, velocity, mach number and pressure from Dataframes 
#    time_th = df_den_th.values[1:, 0]
    x_th    = df_den_th.values[0, 1:]
    
    data_den_th     = df_den_th.values[1:, 1:]
    data_vel_th     = df_vel_th.values[1:, 1:]
#    data_mach_th    = df_mach_th.values[1:, 1:]
    data_press_th   = df_press_th.values[1:, 1:]
    
    
#    time_hun = df_den_hun.values[1:, 0]
    x_hun    = df_den_hun.values[0, 1:]
    
    data_den_hun     = df_den_hun.values[1:, 1:]
    data_vel_hun     = df_vel_hun.values[1:, 1:]
#    data_mach_hun    = df_mach_hun.values[1:, 1:]
    data_press_hun   = df_press_hun.values[1:, 1:]
    
#    return time_th, x_th, data_den_th, data_vel_th, data_mach_th, data_press_th, \
#        time_hun, x_hun, data_den_hun, data_vel_hun, data_mach_hun, data_press_hun

    ### plot
    
    ax[0].plot(x_th, data_den_th[-1, :], 'r'+line_types[2*flag-2], label = 'flux=' + str(flag) + label_th)
    ax[1].plot(x_th, data_vel_th[-1, :], 'b'+line_types[2*flag-2], label = 'flux=' + str(flag) + label_th)
    ax[2].plot(x_th, data_press_th[-1, :], 'k'+line_types[2*flag-2], label = 'flux=' + str(flag) + label_th)
    
    ax[0].plot(x_hun, data_den_hun[-1, :], 'r'+line_types[2*flag-1], label = 'flux=' + str(flag) + label_hun)
    ax[1].plot(x_hun, data_vel_hun[-1, :], 'b'+line_types[2*flag-1], label = 'flux=' + str(flag) + label_hun)
    ax[2].plot(x_hun, data_press_hun[-1, :], 'k'+line_types[2*flag-1], label = 'flux=' + str(flag) + label_hun)
    
    ax[0].set_ylabel('density')
    ax[1].set_ylabel('velocity')
    ax[2].set_ylabel('pressure')
    
    ax[0].set_title('x - final density')
    ax[1].set_title('x - final velocity')
    ax[2].set_title('x - final pressure')
    
    ax[0].legend( loc = 'best')
    ax[1].legend( loc = 'best')
    ax[2].legend( loc = 'best')


Plot( flag_flux_scalar[0] )
Plot( flag_flux_scalar[1] )

title_png = 'x-t-den-vel-pre' #+ str(flag_flux_scalar)
fig.text(0.5, 1, title_png, verticalalignment = 'top', \
    horizontalalignment = 'center', fontsize=16)
fig.savefig(title_png + '.png', format = 'png')


