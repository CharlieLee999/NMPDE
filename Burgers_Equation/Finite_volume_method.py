# -*- coding: utf-8 -*-
"""
Created on Sun Dec 10 16:08:21 2017

@author: charlie
"""

import numpy as np
import matplotlib.pyplot as plt

N           = 100        # 100 or 1000
alpha_coeff = 2.0       # 1.0 or 2.0
flag_init   = 0         # if flag_init = 1 --> sin
                        # if flag_init = 0 --> step

T         = 1
start_p   = 0.0
end_p     = 1.0
h         = (end_p - start_p) / N

u_equal_epsl = 1e-5
x            = np.linspace(start_p, end_p, N, endpoint = False)
len_x        = len(x)

num_iter   = 10
count_iter = 0

a_i_plus   = np.zeros( N )
a_i_minus  = np.zeros( N )
c_i_plus   = np.zeros( N )
c_i_minus  = np.zeros( N ) 
### Initialization of U_n (U_0 initial condition) 

if flag_init == 1:
    u_n     = 1 + 0.5 * np.sin(2 * np.pi * x)
else:
    u_n     = np.zeros(len_x)
    u_n[int(0.1*len_x): int(0.3*len_x+1)] = 0.3 

u_n_new = np.zeros(len_x)
tau_sum   = 0

fig = plt.figure()
plt.ion()
plt.xlabel('$x$')
plt.ylabel('$t$')
plt.title('$u(x,t)$')

while tau_sum <= T:
    
    ### Calculate f_u 
    f_u = 0.5 * u_n**2
    
    ### Calculate a(i + 1/2)
    for i in range(N-1):  ### calculate i = 0 ,...., N-2
        if(abs(u_n[i+1] - u_n[i]) > u_equal_epsl):
            a_i_plus[i] = (f_u[i+1] - f_u[i]) / (u_n[i+1] - u_n[i])
        else:
            a_i_plus[i] = u_n[i]  
            
    ### Calculate a(i+ 1/2)  when i = N - 1,*** u_n[N-1+1] = u_n[N] = u_n[0]
    if(abs(u_n[0] - u_n[N-1]) > u_equal_epsl):
        a_i_plus[-1] = (f_u[0] - f_u[-1]) / (u_n[0] - u_n[-1])
    else:
        a_i_plus[-1] = u_n[-1]
    
    ### Calculate a(i - 1/2) with a(i + i/2)
    a_i_minus[1:] = a_i_plus[0:-1]
    a_i_minus[0]  = a_i_plus[-1]
    
    
    ### Calculate alpha(i + 1/2) and alpha(i - 1/2)
    alpha_i_plus  = abs(a_i_plus)
    alpha_i_minus = abs(a_i_minus)
    
    c_i_plus    = 0.5 * (alpha_i_plus - a_i_plus)
    c_i_minus   = 0.5 * (alpha_i_minus + a_i_minus)
    c_temp = c_i_plus + c_i_minus    
                                                ### Calculate c_temp = c^n (i, i-1) + c^n (i, i+1)
#    c_temp = 0.5 * (alpha_i_plus - a_i_plus) + 0.5 * (a_i_minus + alpha_i_plus)
    
    ### Calculate tau = tau_n = min(h/c_temp)
    
    tau = min(h/c_temp)
    tau_sum += tau

    ### Calculate u_n_new when i = 1, ..., N-2
#    for i in range(1, N - 1):
#        f_n_temp   = 0.5 *(f_u[i+1] - f_u[i-1] - alpha_i_plus[i] *(u_n[i+1] -u_n[i]) \
#                   + alpha_i_minus[i] * (u_n[i] - u_n[i-1]))
#        u_n_new[i] = u_n[i] - 1.0 * tau * N * f_n_temp  # h = 1 / size_mesh   
#    
#    ### Calculate u_n_new when i = 0
#    f_n_temp    = 0.5 * (f_u[1] - f_u[-1] - alpha_i_plus[0] * (u_n[1] - u_n[0])  \
#                + alpha_i_minus[0] * (u_n[0] - u_n[-1]))
#    u_n_new[0]  = u_n[0] - 1.0 * tau * N * f_n_temp
#    
#    ### Calculate u_n_new when i = N-1
#    f_n_temp   = 0.5 *(f_u[0] - f_u[N-2] - alpha_i_plus[N-1] *(u_n[0] -u_n[N-1]) \
#                   + alpha_i_minus[N-1] * (u_n[N-1] - u_n[N-2]))
#    u_n_new[-1] = u_n[N-1] - 1.0 * tau * N * f_n_temp



    for i in range(1, N-1):
        u_n_new[i] = (1.0 - tau/h * (c_i_plus[i] + c_i_minus[i])) * u_n[i] + tau/h * (c_i_plus[i] * u_n[i+1] + c_i_minus[i] * u_n[i-1])
    
    u_n_new[0] = (1.0 - tau/h * (c_i_plus[0] + c_i_minus[0])) * u_n[0] + tau/h * (c_i_plus[0] * u_n[1] + c_i_minus[0] * u_n[N-1])
    
    u_n_new[N-1] = (1.0 - tau/h * (c_i_plus[N-1] + c_i_minus[N-1])) * u_n[N-1] + tau/h * (c_i_plus[N-1] * u_n[0] + c_i_minus[N-1] * u_n[N-2])
    u_n = u_n_new.copy()
    
    count_iter += 1
    
#    fig.clear()
    plt.scatter(x, tau_sum*np.ones(len_x), c = u_n)
#    plt.plot(x, u_n)
    plt.legend(str(count_iter))
    plt.pause(0.5)

str1 = 'T = ' + str(T)
title_ = 'final u(x,' + str(T)+'),' + '  N = ' + str(N) + ', coeff_a  = '+str(alpha_coeff)
plt.figure()
plt.plot(x, u_n)
plt.xlabel('$x$')
plt.ylabel('final u(x,' + str(T)+')')
plt.legend(str(count_iter))
plt.title(title_)
#plt.legend('')
plt.savefig(title_ + '.eps')







