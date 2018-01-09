
### HW6: The Shock Tube problem
import csv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

### Input data ###

flag_flux = 0 # if 0 -> A matrix or 1 -> alpha coefficiency 
N = 400 # 1000, 100
T = 0.2 

x_start = 0.0
x_end   = 1.0

gama = 1.4 
C = 1.0 
A_lamda = 1.0
CFL = 1.0

### output data ###
Time = []
Density = []
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

#f_2 = np.zeros( size_mesh )     # f[1]
#f_3 = np.zeros( size_mesh )     # f[2]

f_mat = np.zeros( (3, size_mesh) )

alpha = np.zeros( N )           # alpha_{i+ 1/2}
f_hat = np.zeros( (3, N) )           # much easier to calculate f_{i+1/2} and f_{i-1/2}
f_i_plus_half  = np.zeros( (3, N-1) ) #f_{i+1/2}, i = 1 ... N-1 
f_i_minus_half = np.zeros( (3, N-1) ) #f_{i-1/2}, i = 1 ... N-1

enthalpy  = np.zeros( size_mesh )
spd_sound = np.zeros( size_mesh )
# lamda = np.zeros( size_mesh  )
lamda_max = 0.0 

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

fig = plt.figure()
axes = fig.add_subplot(111, projection = '3d')


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
    
    for i in range(len(alpha)):
        alpha[i] = C * max(vel_abs_plus_c[i], vel_abs_plus_c[i+1])
    
#    for i in range (0, 3):
    f_hat[0: ,:] = 0.5 * (f_mat[0: , 1:] + f_mat[0: , 0:-1]) \
                     - 0.5 * alpha * (u_mat[0: , 1:] - u_mat[0: , 0:-1])
    f_i_plus_half[0:, :]  = f_hat[0:, 1:]
    f_i_minus_half[0:, :] = f_hat[0:, 0:-1]
    
    ### Calculate the tau^n
    lamda_max = np.max(  vel_abs_plus_c )
    tau = CFL * delta_x / lamda_max
    
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
    Mach_num.append( mach_num )
    Pressure.append( pressure )   
    print Density[0][N/2]
    
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
    plt.pause(0.05)
    
    axes.plot(x, tau_sum * np.ones(len(x)), Density[counter_iter-1][:], c='k')
    axes.set_xlabel('x')
    axes.set_ylabel('time')
    axes.set_zlabel('density')
    axes.set_title('x-t-density')
    

f.savefig('final_state.png',format = 'png') 
fig.savefig('x_t_density.png', format = 'png') 
#
fig2 = plt.figure()
axes2 = fig2.add_subplot(111, projection = '3d')
for i in range(len(Time)):
    axes2.plot(x, Time[i] * np.ones(len(x)), Density[i][:])

axes2.set_xlabel('x')
axes2.set_ylabel('time')
axes2.set_zlabel('density')
plt.title('x-t-density')

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
for i in range( len(Density) ):
    print Density[i][N/2]
#    plt.plot(x, Density[i][:])

x_str = []
for xx in x:
    x_str.append( str(xx) )

 
with open('density.csv', 'wb') as den:
    writer = csv.writer(den, delimiter = ',')
    writer.writerow( ['time,value,x'] + x_str )
    for i in range( len(Density) ):
        writer.writerow( [Time[i]] + Density[i].tolist() ) 

with open('mach_num.csv', 'wb') as mach:
    writer = csv.writer(mach, delimiter = ',')
    writer.writerow( ['time,alue,x'] + x_str )
    for i in range( len(Mach_num) ):
        writer.writerow( [Time[i]] + Mach_num[i].tolist() ) 

with open('pressure.csv', 'wb') as pressure:
    writer = csv.writer(pressure, delimiter = ',')
    writer.writerow( ['time,value,x'] + x_str )
    for i in range( len(Pressure) ):
        writer.writerow( [Time[i]] + Pressure[i].tolist() ) 

with open('time.csv', 'wb') as time:
    writer = csv.writer(time, delimiter = ',')
    writer.writerow(['x', 'time'])
    writer.writerow( x )
    writer.writerow( Time )