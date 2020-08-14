# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 10:41:09 2020

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

############################
# integration method control
############################

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

w_0 = 10 * 10 ** (-6)


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
#############################################################
#3 plot of Q_z vs displacement w for various offset rho_0x
#############################################################

#3 plot of Q_z vs displacement w for various offset rho_0x/w0


rho_0 = [0,0]   #no offset


rho_0[0] = [0, 0.25*w_0, 0.5*w_0, w_0/np.sqrt(2)]
w = np.linspace(w_0, 60 * 10** (-6), 100)

#d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))



Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )


Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][1],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][2],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )
Axial_flist_vs_d3 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][3],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z3 = Axial_flist_vs_d3 * c / ( n_0 * P )


plt.figure(1)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="rho_0x = 0")
plt.plot( w * 10 ** (6), Q_z1 , lw=2, c="r", label="rho_0x = 2.5um")
plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho_0x = 5um")
plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho_0x = 10/sqrt(2)um")


new_ticks1 = np.linspace(0, 70, 8) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=13)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','0, 2.5, 5, 7.07', '30')

t_finish = time.time()

print(f"script executed in {t_finish-t_start:.2f} seconds")

'''




'''
##########################################################
#4 plot of Q_z vs displacement w for various radius rho
##########################################################


rho_0 = [0,0]   #no offset



w = np.linspace(w_0, 2* rho, 100)

rho = [0.5*rho, 0.75*rho, rho]



Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )


Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )
#Axial_flist_vs_d3 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[3], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
#Q_z3 = Axial_flist_vs_d3 * c / ( n_0 * P )


plt.figure(2)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="rho = 15um")
plt.plot( w * 10 ** (6), Q_z1 , lw=2, c="r", label="rho = 22.5um")
plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho = 30um")
#plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho = 11.2*w_0")


new_ticks1 = np.linspace(0, 60, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=13)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','0',  '15, 22.5, 30')

t_finish = time.time()

print(f"script executed in {t_finish-t_start:.2f} seconds")
'''


'''
#1 plot of Q_z vs displacement d


rho_0 = [0,0]   #no offset


w = np.linspace(w_0, 60 * 10** (-6), 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))


Axial_flist_vs_d =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z = Axial_flist_vs_d * c / ( n_0 * P )

plt.figure(3)


plt.plot( d * 10 ** (6), Q_z , lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 1800, 10) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0z(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('x-axis w0 to 60', 'related to w', '0','0',  '30')
'''