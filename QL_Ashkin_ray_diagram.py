# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 20:48:44 2020

@author: liuqi
"""



import scipy as sp
import numpy as np

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pylab as plt
import scipy.integrate as spi
from scipy.integrate import quad
import seaborn
from scipy.integrate import odeint
from scipy.integrate import dblquad
import Module_Trialplot_Q as TQ
import Module_table_parameter as MTP
import Module_integration_manually as MIM


import time

##########################################################
# integration method control: stiffness plot of Qx, Qz
##########################################################



integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 100

plt.close('all')

###########################

t_start = time.time()


#Our sphere

g = 9.8 #gravitational acceleration
c = 3 * 10**8
#m = 4/3 * np.pi * Rs**3 * (sig_s - sig_0)

#TEM01* reflective target Table 6 matching

w_0 = 10* 10 ** (-6)


Lambda = 1.064 * 10**(-6)
#P = 1
z_R = np.pi* w_0 ** 2 / Lambda
rho = 30 * 10 ** (-6)
 
weight = 3.4 * 10**(-10)

n_0 = 1

n_s_n = 0.04

k = 7.6097

n_s = n_s_n - k*1j

sig_s = 10.49 * 10 ** 3 * ( 3 ** 3 - 2.9 ** 3 / 3 ** 3 ) #density of sphere in kg/m^3

sig_0 = 0 #density of medium in kg/m^3



m = 4/3 * np.pi * rho ** 3 * ( sig_s - sig_0 )

Permittivity = 8.85 * 10**(-12)

P = 0.5 * c * n_0 * Permittivity    #total power of the LG01 beam



theta = np.linspace(0, np.pi/2, grid_size)

#phi = np.linspace(0, np.pi * 2, grid_size)

#theta,phi = np.meshgrid(theta, phi)


theta_c = np.arcsin( ( n_s.real ) / n_0 )

rlist = []
Qz = []
Qy = []
    
for i in theta:
    print (i)
    if i > theta_c:                         #total internal reflection
        
        theta_2 = np.pi/2
        
    else:                                       #for non TIR
        
        theta_2 = np.arcsin(n_0*np.sin(i)/n_s.real) 
    


    Reflection_coe = TQ.rfcoe_sm(i,theta_2, n_0, n_s)
    
    
    R = np.around(Reflection_coe, decimals= 15) 

    Fz = n_0/ c * (1 + R * np.cos(2 * i))
    print (Fz)

    Q_z = Fz * c / ( n_0 * P )
    print (Q_z)

    Fy = n_0/c * R * np.sin(2 * i)

    Q_y = Fy * c / ( n_0 * P )
    print (Q_y)

    r = rho * np.sin(i)

    rlist.append(r)
    Qz.append(Q_z)
    Qy.append(Q_y)

plt.plot(np.asarray(rlist) * 10 ** 6, np.asarray(Qz) / max(Qz), lw=2, c="c", label = 'Fz')
plt.plot(np.asarray(rlist) * 10 ** 6, np.asarray(Qy) / max(Qy), lw=2, c="r", label = 'Fr')




new_ticks1 = np.linspace(0, 30, 4) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(0, 1, 3),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('X(um)',fontsize=20)
plt.ylabel('nomalized force
           ',fontsize=20)

#plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


    
