# -*- coding: utf-8 -*-
"""
Created on Wed Jan 10 14:14:31 2018

@author: charlie
"""

### FEM of 1D Poisson-Equation

import csv
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

### Input data ###
N = 1000  # 100, 1000

x_start = 0.0
x_end   = 1.0

alpha = 10.0
x = (np.exp( alpha * np.arange( N+1 ) / N) - 1.0) / \
        (np.exp(alpha) - 1.0)

f = np.pi**2 * np.sin( np.pi * x)

### Initialization

phi_mat = np.eye( N-1, N+1, k = 1 )
Z_vec   = np.zeros( N-1 )
A_mat   = np.zeros( (N-1, N-1) )
f_inter = np.zeros( N-1 )
### Calculate the A matrix 

sp_mesh = 1.0 / ( x[1:] - x[0:-1])

for i in range(N-1):
    A_mat[i, i] = sp_mesh[i] + sp_mesh[i+1]
    
    if i == N-2:
        pass
    else:
        A_mat[i, i+1] = -sp_mesh[i+1]
    
    if i == 0:
        pass
    else: 
        A_mat[i, i-1] = -sp_mesh[i]

f_inter = 0.5 * f[1:-1] * (x[2:] - x[0:-2])

A_mat_inv = np.linalg.inv(A_mat)
Z_vec = A_mat_inv.dot(f_inter)

title = ['Theoretical and numerical u_x n ' + str(N), 'Approximation errors n ' + str(N)]

u_vec_num = np.append(np.asarray([0.0]), np.append(Z_vec, 0.0))
u_vec_the = np.sin(np.pi * x)


#plt.figure()
#plt.plot(x, u_vec_num, label = 'numerical')
#plt.plot(x, u_vec_the, label = 'theoretical')
#plt.legend()
#plt.xlabel('x')
#plt.ylabel('$u(x)$')
#plt.title('Theoretical and numerical $u(x)$ (N = ' + str(N) + ')')
#plt.savefig(title[0] + '.png', format = 'png')
#
errors = u_vec_num - u_vec_the
#plt.figure()
#plt.plot(x, errors)
#plt.xlabel('x')
#plt.ylabel('$e(x)$')
#plt.title('Approximation errors (N =' + str(N) + ')')
#plt.savefig(title[1] + '.png', format = 'png')


data_out = np.zeros( (4, len(x)) )
data_out[0, :] = x.copy()
data_out[1, :] = u_vec_num.copy()
data_out[2, :] = u_vec_the.copy()
data_out[3, :] = errors.copy()

df = pd.DataFrame(data_out)
df.to_csv('data_N_' + str(N) +'.csv')