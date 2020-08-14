# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 11:35:58 2020

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
import Will_Module_Trialplot_Q as TQ
import Module_table_parameter as MTP
import Module_integration_manually as MIM

import time

##########################################################
# integration method control: Qz vs w plots
##########################################################



integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 200

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


###################################
#6 plot of Q_z vs displacement Rs/w
###################################

#6 plot of Q_z vs displacement Rs/w


rho_0 = [0,0]   #no offset

w = np.sqrt(2) * 30 * 10 ** (-6)

rho = np.linspace(10** (-6), rho, 100)


Axial_flist_vs_Rs =  np.asarray(TQ.Fz_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z = Axial_flist_vs_Rs * c / ( n_0 * P )

plt.figure(6)


plt.plot(rho * 10 ** 6, Q_z , lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 30, 4) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.3, 0, 4),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
#plt.title('Qz vs radius/beam waist w0 at w = w0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30', 'related to w', '0','0', '1 to 30')