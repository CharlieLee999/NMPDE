# -*- coding: utf-8 -*-
"""
Created on Sun Nov 19 13:52:26 2017

@author: charlie
"""

import numpy as np
import matplotlib.pyplot as plt

### Input grid sizes
N = input('The size of spatial mesh:')
N = int(N)


### Parameters to choose initial condition
### flag_scheme == 1 --> modified equation of upwind scheme
### flag_scheme == 0 --> normal equation of upwind scheme
### flag_init_condition == 1 --> discontinuous initial conditionS
### flag_init_condition == 0 --> sin wave initial conditions
flag_scheme = 1
flag_init_condition = 0

length = 1.0
h = length/N
nu = 0.9
a = 1.0
tho = nu * h / a
M = int(N / nu)
T = M * tho

size_A_mat = N + 1

x = np.linspace(0, length, N+1)
U_mat_upwind = np.zeros((N+1, M+1))
U_mat_lw = np.zeros((N+1, M+1))
U_mat_upwind_modified = np.zeros((N+1, M+1))

### Initial conditions
if flag_scheme == 1:
    U_mat_upwind_modified[:, 0] = U_mat_lw[:, 0] = U_mat_upwind[:, 0] = np.sin(2 * np.pi * x)
elif flag_init_condition == 0:
    U_mat_lw[:, 0] = U_mat_upwind[:, 0] = np.sin(2 * np.pi * x)
else:
    half_M = int(M / 2.0)    
    U_mat_lw[0 : half_M, 0] = U_mat_upwind[0 : half_M, 0] = 1.0


### Calculate s of different schemes
s_upwind = abs(nu)
s_lw = nu**2
mu_modified = 1.0 * abs(a) * h * (1 - abs(nu)) / 2
s_upwind_modified = mu_modified * tho * 2 / h**2 + s_upwind


### Function of generating the general A matrix with input parameter s
def Get_A_matrix(s):
    #A_mat = np.zeros((size_A_mat,size_A_mat))
    middle_diag  = (1.0 - s) * np.ones(size_A_mat)
    upper_diag   = (s - nu) / 2.0 * np.ones(size_A_mat - 1)
    lower_diag   = (s + nu) / 2.0 * np.ones(size_A_mat - 1)
    A_mat        = np.diag(middle_diag) + np.diag(upper_diag, k = 1) + np.diag(lower_diag, k = -1)
    A_mat[0, -2] = (s + nu) / 2.0
    A_mat[-1, 1] = (s - nu) / 2.0
    return A_mat


### Generate both A matrix of upwind and LW scheme
A_mat_upwind = Get_A_matrix(s_upwind)
A_mat_lw     = Get_A_matrix(s_lw)
A_mat_upwind_modified = Get_A_matrix(s_upwind_modified)


### U iterates with time 
for i in range(1, M+1):
    U_mat_upwind[:, i] = A_mat_upwind.dot(U_mat_upwind[:, i-1])
    U_mat_lw[:, i] = A_mat_lw.dot(U_mat_lw[:, i-1])
    U_mat_upwind_modified[:, i] = A_mat_upwind_modified.dot(U_mat_upwind_modified[:,i-1])


### Plot final approximated U and errors
if flag_scheme == 1:
    U_exact_final_modified = np.sin(2.0 * np.pi * (x - T)) * np.exp(-mu_modified * 4 * np.pi**2 * T)
    
    error_final_upwind_modified = U_mat_upwind_modified[:,-1] - U_exact_final_modified
    max_error_final_upwind_modified = max(error_final_upwind_modified)
    print ('The max error of modified equation is:')
    print max_error_final_upwind_modified
    
    plt.figure()
    plt.plot(x, U_mat_upwind_modified[:, -1], label = 'Upwind modified scheme')
    plt.plot(x, U_exact_final_modified, label= 'Exact solution')
    plt.xlabel('Position')
    plt.ylabel('$U(x, T)$')
    plt.title('Final approximated $U(x, T)$ of modified Upwind scheme')
    plt.legend(loc = 'best')
    plt.savefig('Appr_U_Modified_upwind_schemes')
    plt.show()

else:
    ## when flag_scheme is 0, normal advection equation, or modified equation of Upwind scheme
    
    if flag_init_condition == 0:
        U_exact_final = np.sin( 2.0 * np.pi * (x - T))
    else:
        U_exact_final = U_mat_upwind[:, 0]
    #    U_exact_final = x - a * T
    #    U_exact_final[(U_exact_final<0.5) & (U_exact_final>=0)] = 1.0
    #    U_exact_final[(U_exact_final>=0.5) | (U_exact_final<0)] = 0.0
    
    error_final_upwind = U_mat_upwind[:, -1] - U_exact_final
    error_final_lw     = U_mat_lw[:, -1] - U_exact_final
    
    max_error_final_upwind = max(error_final_upwind)
    max_error_final_lw     = max(error_final_lw)
    
    ### print the max error of Upwind and Lax-Wendroff scheme respectively
    print ('Max final appro. error of Upwind scheme is')
    print max_error_final_upwind
    print ('Max final appro. error of Lax-Wendroff scheme is')
    print max_error_final_lw
    
    ### Plot the final U at each position of the Upwind scheme
    plt.figure()
    plt.plot(x, U_mat_upwind[:, -1], label = 'Upwind scheme')
    plt.plot(x, U_mat_lw[:, -1], label = 'Lax-Wendroff scheme')
    plt.plot(x, U_exact_final, label = 'Exact solution')
    plt.xlabel('Position')
    plt.ylabel('$U(x, T)$')
    plt.title('Final approximated $U(x, T)$ of Upwind and L-W schemes')
    plt.legend(loc = 'best')
    plt.savefig('Appr_U_upwind_lw_schemes')
    plt.show()
    
    ### plot final approximation errors at each position of the LW scheme
    plt.figure()
    plt.plot(x, error_final_upwind, label = 'Upwind scheme')
    plt.plot(x, error_final_lw, label = 'Lax-Wendroff scheme')
    plt.xlabel('Position x')
    plt.ylabel('Appr. Error')
    plt.title('Approximation errors of Upwind and L-W schemes')
    plt.legend(loc = 'best')
    plt.savefig('Appr_err_upwind_lw_schemes')
    plt.show()
