# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 17:53:28 2020

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
# integration method control: Merging plots
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

'''
#################################################################
#3 plot of Q_z, Qy, Qx vs displacement w at offset = 2.5um, rho = 30um
#################################################################

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*rho, 100)




rho_0 = [0, 0]

rho_0[0] = 0.25 * w_0#[0.25 * w_0, 0.5 * w_0, w_0 / np.sqrt(2)]



radial_flistx0 = np.asarray(TQ.Fx_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0  * c / ( n_0 * P )



Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )




plt.figure(1)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="Qz")
plt.plot( w * 10 ** (6), Q_x0 , lw=2, c="r", label="Qx")
#plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho_0x = 5um")
#plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho_0x = 10/sqrt(2)um")


new_ticks1 = np.linspace(0, 60, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0.5, 6),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=16)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Q',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','2.5', '30')

t_finish = time.time()

print(f"script executed in {t_finish-t_start:.2f} seconds")
'''


'''
#################################################################
#3 plot of Q_z, Qy, Qx vs displacement w at offset = 7.07um, rho = 30um
#################################################################

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*rho, 100)




rho_0 = [0, 0]

rho_0[0] =  w_0 / np.sqrt(2)#[0.25 * w_0, 0.5 * w_0, w_0 / np.sqrt(2)]



radial_flistx0 = np.asarray(TQ.Fx_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0  * c / ( n_0 * P )



Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )




plt.figure(1)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="Qz")
plt.plot( w * 10 ** (6), Q_x0 , lw=2, c="r", label="Qx")
#plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho_0x = 5um")
#plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho_0x = 10/sqrt(2)um")


new_ticks1 = np.linspace(0, 60, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0.5, 6),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=16)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Q',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','7.07', '30')

t_finish = time.time()

print(f"script executed in {t_finish-t_start:.2f} seconds")

'''

'''
########################################################################
#3 plot of Q_z, Qy, Qx vs displacement w at offset = 7.07um, rho = 15um
########################################################################

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*rho, 100)


rho_0 = [0, 0]

rho_0[0] = w_0 / np.sqrt(2)

rho = 0.5 * rho #[0.5 * a, 0.75 * a, a]



radial_flistx0 = np.asarray(TQ.Fx_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0  * c / ( n_0 * P )

Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )


plt.figure(3)


plt.plot( w * 10 ** (6), Q_x0 , lw=2, c="c", label="Qx")
plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="r", label="Qz")
#plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho_0x = 5um")
#plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho_0x = 10/sqrt(2)um")


new_ticks1 = np.linspace(0, 60, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0.5, 6),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=16)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Q',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','7.07', '15')

t_finish = time.time()

print(f"script executed in {t_finish-t_start:.2f} seconds")
'''

########################################################################
#3 plot of Q_z, Qy, Qx vs displacement w at offset = 7.07um, rho = 22.5um
########################################################################

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*rho, 100)


rho_0 = [0, 0]

rho_0[0] = w_0 / np.sqrt(2)

rho = 0.75 * rho #[0.5 * a, 0.75 * a, a]



radial_flistx0 = np.asarray(TQ.Fx_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0  * c / ( n_0 * P )

Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )


plt.figure(3)


plt.plot( w * 10 ** (6), Q_x0 , lw=2, c="c", label="Qx")
plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="r", label="Qz")
#plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho_0x = 5um")
#plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho_0x = 10/sqrt(2)um")


new_ticks1 = np.linspace(0, 60, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0.5, 6),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=16)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Q',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','7.07', '22.5')

t_finish = time.time()

print(f"script executed in {t_finish-t_start:.2f} seconds")