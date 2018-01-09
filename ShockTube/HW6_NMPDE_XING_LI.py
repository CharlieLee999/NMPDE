
### HW6: The Shock Tube problem
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

### Input data ###

flag_flux_scalar = 1 # if 2 -> A matrix or 1 -> alpha coefficiency 
flag_continuous = 0 #  if 0 -> continuous initial, or 1 -> discrete initial 

N = 1000 # 1000, 100 

x_start = 0.0
x_end   = 1.0
T = 0.2

gama = 1.4 
C = 1.0 
C_roe = 1.0
A_lambda = 1.0
CFL = 1.0

### output data ###
Time = []
Density = []
Velocity = []
Mach_num = []
Pressure = []

### Initial conditions ###

size_mesh = N + 1

delta_x = (x_end - x_start) / N
x = np.linspace( x_start, x_end, size_mesh ) # spatial mesh 

den = np.ones( size_mesh )    # u[0]
den[ x >=0.5 ] = 0.125
E   = 2.5 * np.ones( size_mesh )    # u[2]
E[ x >=0.5 ] = 0.25
velocity = np.zeros( size_mesh )

### Initialization ###

tau_sum = 0.0
counter_iter = 0
den_times_vel = den * velocity  # u[1], f[0]

u_mat = np.zeros( (3, size_mesh) )
u_mat_temp = np.zeros( (3, N-1 ) )

den_vel_sqr = np.zeros( size_mesh )         # rho * u ^ 2 
pressure    = np.zeros( size_mesh )         # pressure p

mach_num = np.zeros( size_mesh )

### Second formulation of discrete flux function
R = np.ones( (9, N) )
R_inv = np.zeros( (9, N) )
delta_mat_abs = np.zeros( (9, N) )
A_mat = np.zeros( (9, N) )
A_mat_temp = np.zeros( (3,3) )
lambda_mat = np.zeros( (3, N) )

#f_2 = np.zeros( size_mesh )     # f[1]
#f_3 = np.zeros( size_mesh )     # f[2]

f_mat = np.zeros( (3, size_mesh) )

alpha = np.zeros( N )           # alpha_{i+ 1/2}
f_hat = np.zeros( (3, N) )           # much easier to calculate f_{i+1/2} and f_{i-1/2}
f_i_plus_half  = np.zeros( (3, N-1) ) #f_{i+1/2}, i = 1 ... N-1 
f_i_minus_half = np.zeros( (3, N-1) ) #f_{i-1/2}, i = 1 ... N-1

enthalpy  = np.zeros( size_mesh )
spd_sound = np.zeros( size_mesh )
# lambda = np.zeros( size_mesh  )
lambda_max = 0.0 

### Definition of functions ###
def Get_den_vel_sqr(density, velocity):
    return density * velocity**2 
    
def Get_pressure(gama, Energy, den_vel_sqr):
    return (gama - 1) * ( Energy - 0.5 * den_vel_sqr )


### Calculate u_vector

u_mat[0, :] = den.copy()
u_mat[1, :] = den_times_vel.copy()
u_mat[2, :] = E.copy()

### initial output data ###
#Time.append( tau_sum )
#Density.append( den )

title_state = 'final state' + ' flux ' + str(flag_flux_scalar) + ' N '+ str(N)

f, ax = plt.subplots(3, sharex = True)
plt.ion()
plt.xlabel('x')
ax[0].plot(x, u_mat[0, :], c = 'r')
ax[1].plot(x, u_mat[1, :], c = 'b')
ax[2].plot(x, u_mat[2, :], c = 'k')
ax[0].set_title('density')
ax[1].set_title('den_times_vel')
ax[2].set_title('energy')

#plt.xlim( (0.4, 0.6) )
plt.pause(0.5)


while tau_sum < T :   # tau_sum < T  # counter_iter < 30
    ### Calculate u_vector and f_vecotr ###
    den_vel_sqr = Get_den_vel_sqr(u_mat[0, :], velocity)
    pressure  = Get_pressure(gama, u_mat[2, :], den_vel_sqr)
        
    
    ### Calculate f_mat 
    f_mat[0, :] = u_mat[1, :].copy() 
    f_mat[1, :] = den_vel_sqr + pressure
    f_mat[2, :] = velocity * (u_mat[2, :] + pressure)
    
#    print min(u_mat[0,:])
    
#    print pressure[N/2]
    ### Calculate the speed of sound ###
    spd_sound = np.sqrt( gama * pressure / u_mat[0, :]) 
    mach_num = velocity / spd_sound 
#    print max(mach_num)
    vel_abs   = abs(velocity)
    
    ### 
    ### Calculate the alpha_{i+1/2} and f_hat
    vel_abs_plus_c = vel_abs + spd_sound
#    for i in range(len(velocity)):
#        if vel_abs_plus_c[i] < spd_sound[i] : 
#            vel_abs_plus_c[i] = 0.5 * (spd_sound[i] + vel_abs_plus_c[i]**2 / spd_sound[i])
#        vel_abs_plus_c[vel_abs_plus_c < spd_sound]  = 0.5 * (spd_sound + vel_abs_plus_c**2 / spd_sound)   
    
    if flag_flux_scalar == 2:
        enthalpy = ( u_mat[2,:] + pressure ) / u_mat[0, :]
        
        den_sqrt = np.sqrt( u_mat[0, :])    # N+1
        den_sqrt_vel = den_sqrt * velocity  # N+1
        u_hat = ( den_sqrt_vel[0:-1] + den_sqrt_vel[1:] ) / ( den_sqrt[0:-1] + den_sqrt[1:] ) # N
        
        den_sqrt_h = den_sqrt * enthalpy
        h_hat = ( den_sqrt_h[0: -1] + den_sqrt_h[1:] ) / ( den_sqrt[0:-1] + den_sqrt[1:] )
        
        c_hat = np.sqrt( ( gama-1 ) * (h_hat  - 0.5 * u_hat**2 ) )    
        
        
        Gamma = ( gama - 1.0 ) / c_hat**2
        for i in range(3,6):
            R[i, :] = u_hat + (i - 4.0) * c_hat
        R[6, :] = h_hat - u_hat * c_hat 
        R[7, :] = 0.5 * u_hat**2
        R[8, :] = h_hat + u_hat * c_hat
        
        Gamma_u_sqr_h = Gamma * ( u_hat**2 - h_hat)
        
        ### R_Inv = Get_R_inv(Gamma, c_hat, h_hat, u_hat)
        R_inv[0, :] = 0.5 * ( 1 + Gamma_u_sqr_h + u_hat/c_hat)
        R_inv[1, :] = -0.5 * ( Gamma * u_hat + 1.0/c_hat ) 
        R_inv[2, :] = 0.5 * Gamma * np.ones( R_inv.shape[1] )
        R_inv[3, :] = -Gamma_u_sqr_h
        R_inv[4, :] = Gamma * u_hat
        R_inv[5, :] = -Gamma * np.ones( R_inv.shape[1] )
        R_inv[6, :] = R_inv[0, :] - u_hat/c_hat
        R_inv[7, :] = R_inv[1, :] + 1.0/c_hat
        R_inv[8, :] = R_inv[2, :].copy()
        
        for i in range( 3 ):
            lambda_mat[i, :] = u_hat + (i-1) * c_hat
            ### lambda_mat[i, :] = Get_disc_eigenvalue(A_lambda_ss, lambda_mat[i,:])
        
        A_lambda_ss = A_lambda * c_hat
        
        ### change lambda when initial conditions are discrete
        if flag_continuous == 0:
            for i in range( 3 ):
#                for j in range ( N ):
#                    if lambda_mat[i,j] < A_lambda_ss[j] :
#                        lambda_mat[i, j] = 0.5 * (A_lambda_ss[j] + lambda_mat[i,j]**2 / A_lambda_ss[j])
                lambda_greater = lambda_mat[i, :] < A_lambda_ss        
                lambda_mat[i, lambda_greater] = \
                 0.5 * ( A_lambda_ss[lambda_greater] + \
                 lambda_mat[i, lambda_greater]**2 / A_lambda_ss[lambda_greater] )
        
        ### delta_mat_abs 
        delta_mat_abs[0, :] = np.abs(lambda_mat[0, :])
        for i in range(2, 4):
            delta_mat_abs[2**i, :] = np.abs(lambda_mat[i-1, :])
            
        for j in range( N ): 
            A_mat_temp = R[:, j].reshape((3,3)).dot( delta_mat_abs[:, j].reshape((3,3)).dot( \
                           R_inv[:, j].reshape((3,3)) ) )
            f_hat[:, j] = 0.5 * (f_mat[:, j+1] + f_mat[:, j]) \
                        - 0.5 * C_roe * A_mat_temp.dot(u_mat[:, j+1] - u_mat[:, j]) 

    else: 
        for i in range(len(alpha)):
            alpha[i] = C * max(vel_abs_plus_c[i], vel_abs_plus_c[i+1])
        
    #    for i in range (0, 3):
        f_hat = 0.5 * (f_mat[:, 1:] + f_mat[:, 0:-1]) \
                         - 0.5 * alpha * (u_mat[:, 1:] - u_mat[:, 0:-1])
    
    
    f_i_plus_half  = f_hat[:, 1:]
    f_i_minus_half= f_hat[:, 0:-1]
    
    ### Calculate the tau^n
    lambda_max = np.max(  vel_abs_plus_c )
    tau = CFL * delta_x / lambda_max
    
    ### Update u_vector and velocity
    temp = f_i_plus_half - f_i_minus_half
    u_mat[: ,1:-1] = u_mat[:, 1:-1]  - tau / delta_x * (f_i_plus_half - f_i_minus_half)
    
    velocity = u_mat[1, :] / u_mat[0, :]    
    
    ### update tau_sum and iteration time 
    tau_sum += tau
    counter_iter += 1
    
    ### Output data ###
#    print u_mat[0][50]
    Time.append( tau_sum )
    Density.append( u_mat[0,:].copy() )
    Velocity.append( velocity )
    Mach_num.append( mach_num )
    Pressure.append( pressure )   
#    print Density[counter_iter-1][N/2]
    
    ### plot
    ax[0].clear()
    ax[1].clear()
    ax[2].clear()   
    
    ax[0].plot(x, Density[counter_iter-1][:], c = 'r')  #u_mat[1,:]
    ax[1].plot(x, u_mat[1, :], c = 'b')
    ax[2].plot(x, u_mat[2, :], c = 'k')
    ax[0].set_title('density')
    ax[1].set_title('den_times_vel')
    ax[2].set_title('energy')
#    plt.xlim( (0.4, 0.6) )
#    plt.pause(0.05)

f.savefig(title_state, format = 'png') 


### plot states with respect to x and time 
title = []
for s in ['density', 'velocity', 'mach_number', 'pressure', 'time']:
    title_temp = 'x-t-' + s + '_flux_' + str(flag_flux_scalar) + '_N_'+ str(N)    
    title.append(title_temp)


fig_den     = plt.figure()
fig_vel     = plt.figure()
fig_mach    = plt.figure()
fig_press   = plt.figure()


axes_den    = fig_den.add_subplot(111, projection = '3d')
axes_vel    = fig_vel.add_subplot(111, projection = '3d')
axes_mach   = fig_mach.add_subplot(111, projection = '3d')
axes_press  = fig_press.add_subplot(111, projection = '3d')

for i in range(len(Time)):
    axes_den.plot(x, Time[i] * np.ones(len(x)), Density[i][:])
    axes_vel.plot(x, Time[i] * np.ones(len(x)), Velocity[i][:])
    axes_mach.plot(x, Time[i] * np.ones(len(x)), Mach_num[i][:])
    axes_press.plot(x, Time[i] * np.ones(len(x)), Pressure[i][:])


axes_den.set_xlabel('x')
axes_den.set_ylabel('time')
axes_den.set_zlabel('Density')
#axes_den.set_title(title[0], verticalalignment = 'top')
fig_den.text(0.5, 0.9, title[0], verticalalignment = 'top', \
    horizontalalignment = 'center', fontsize=16)
fig_den.savefig(title[0]+'.png', format = 'png')

axes_vel.set_xlabel('x')
axes_vel.set_ylabel('time')
axes_vel.set_zlabel('Velocity')
#axes_vel.set_title(title[1], verticalalignment = 'top')
fig_vel.text(0.5, 0.9, title[1], verticalalignment = 'top', \
    horizontalalignment = 'center', fontsize=16)
fig_vel.savefig(title[1]+'.png', format = 'png')

axes_mach.set_xlabel('x')
axes_mach.set_ylabel('time')
axes_mach.set_zlabel('Mach number')
#axes_mach.set_title(title[2], verticalalignment = 'top')
fig_mach.text(0.5, 0.9, title[2], verticalalignment = 'top', \
    horizontalalignment = 'center', fontsize=16)
fig_mach.savefig(title[2]+'.png', format = 'png')

axes_press.set_xlabel('x')
axes_press.set_ylabel('time')
axes_press.set_zlabel('Pressure')
#axes_press.set_title(title[3], verticalalignment = 'top')
fig_press.text(0.5, 0.9, title[3], verticalalignment = 'top', \
    horizontalalignment = 'center', fontsize=16)
fig_press.savefig(title[3]+'.png', format = 'png')





#
#plt.figure()
#plt.plot(x, Density[-1])
##plt.xlim( (0.4, 0.6) )
#plt.xlabel('x')
#plt.ylabel('final density')
#plt.title('x-final density')
#plt.show()
#
#plt.figure()
#plt.plot(x, Mach_num[-1])
##plt.xlim( (0.4, 0.6) )
#plt.xlabel('x')
#plt.ylabel('final Mach_num')
#plt.title('x-final Mach_num')
#plt.show()
#
#plt.figure()
#plt.plot(x, Pressure[-1])
##plt.xlim( (0.4, 0.6) )
#plt.xlabel('x')
#plt.ylabel('final pressure')
#plt.title('x-final pressure')
#plt.show()
#

#plt.figure()
#for i in range( len(Density) ):
#    print Density[i][N/2]
#    plt.plot(x, Density[i][:])

x_str = []
for xx in x:
    x_str.append( str(xx) )


with open(title[0]+'.csv', 'wb') as den:
    writer = csv.writer(den, delimiter = ',')
    writer.writerow( ['time,value,x'] + x_str )
    for i in range( len(Density) ):
        writer.writerow( [Time[i]] + Density[i].tolist() ) 

with open(title[1]+'.csv', 'wb') as mach:
    writer = csv.writer(mach, delimiter = ',')
    writer.writerow( ['time,alue,x'] + x_str )
    for i in range( len(Mach_num) ):
        writer.writerow( [Time[i]] + Velocity[i].tolist() )

with open(title[2]+'.csv', 'wb') as mach:
    writer = csv.writer(mach, delimiter = ',')
    writer.writerow( ['time,alue,x'] + x_str )
    for i in range( len(Mach_num) ):
        writer.writerow( [Time[i]] + Mach_num[i].tolist() ) 

with open(title[3]+'.csv', 'wb') as pressure:
    writer = csv.writer(pressure, delimiter = ',')
    writer.writerow( ['time,value,x'] + x_str )
    for i in range( len(Pressure) ):
        writer.writerow( [Time[i]] + Pressure[i].tolist() ) 

with open(title[4]+'.csv', 'wb') as time:
    writer = csv.writer(time, delimiter = ',')
    writer.writerow(['x', 'time'])
    writer.writerow( x )
    writer.writerow( Time )
