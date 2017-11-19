# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 13:52:26 2017

@author: charlie
"""

import numpy as np
import matplotlib.pyplot as plt

### Input grid sizes
N = input('The size of spatial mesh:')
#N = int(N)

length = 1.0
h = length/N
mu = 0.9
a = 1.0
tho = mu * h / a
M = int(N / mu)
T = M * tho

size_A_mat = N + 1

x = np.linspace(0, length, N+1)
U_mat_upwind = np.zeros((N+1, M+1))
U_mat_lw = np.zeros((N+1, M+1))
U_mat_lw[:, 0] = U_mat_upwind[:, 0] = np.sin(2 * np.pi * x)

s_upwind = abs(mu)
s_lw = mu**2


#u_temp  = np.zeros(size(x))
#u_1 = 3

### Function of generating the general A matrix with input parameter s
def Get_A_matrix(s):
    #A_mat = np.zeros((size_A_mat,size_A_mat))
    middle_diag  = (1.0 - s) * np.ones(size_A_mat)
    upper_diag   = (s - mu) / 2.0 * np.ones(size_A_mat - 1)
    lower_diag   = (s + mu) / 2.0 * np.ones(size_A_mat - 1)
    A_mat        = np.diag(middle_diag) + np.diag(upper_diag, k = 1) + np.diag(lower_diag, k = -1)
    A_mat[0, -2] = (s + mu) / 2.0
    A_mat[-1, 1] = (s - mu) / 2.0
    return A_mat


### Generate both A matrix of upwind and LW scheme
A_mat_upwind = Get_A_matrix(s_upwind)
A_mat_lw = Get_A_matrix(s_lw)

### U iterate with time 
for i in range(1, M+1):
    U_mat_upwind[:, i] = A_mat_upwind.dot(U_mat_upwind[:, i-1])
    U_mat_lw[:, i] = A_mat_lw.dot(U_mat_lw[:, i-1])

U_exact_final = np.sin( 2.0 * np.pi * (x - T))

error_final_upwind = U_mat_upwind[:, -1] - U_exact_final
error_final_lw     = U_mat_lw[:, -1] - U_exact_final

max_error_final_upwind = max(error_final_upwind)
max_error_final_lw     = max(error_final_lw)

### print the max error of Upwind and Lax-Wendroff scheme respectively
print ('Max final appro. error of Upwind scheme is')
print max_error_final_upwind
print ('Max final appro. error of Lax-Wendroff scheme is')
print max_error_final_lw

### Plot the final U at each position of both scheme
plt.figure()
plt.plot(x, U_mat_upwind[:, -1], label = 'Upwind scheme')
plt.plot(x, U_mat_lw[:, -1], label = 'Lax-Wendroff scheme')
plt.plot(x, U_exact_final, label = 'Exact solution')
plt.xlabel('Position')
plt.ylabel('$U(x, T)$')
plt.title('$Final- U(x, Y)$')
plt.legend(loc = 'best')
plt.show()

### plot final approximation errors at each position of both scheme
plt.figure()
plt.plot(x, error_final_upwind, label = 'Upwind scheme')
plt.plot(x, error_final_lw, label = 'Lax-Wendroff scheme')
plt.xlabel('Position x')
plt.ylabel('Appr. Error')
plt.title('Approximation errors of Upwind and L-W schemes')
plt.legend(loc = 'best')
plt.show()