# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 10:58:43 2020

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

import Module_Fz_stiffness_at_parameter as MFSAP

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

sig_s = 10.49 * 10 ** 3 * (( 3 ** 3 - 2.25 ** 3 )/ 3 ** 3 ) #density of sphere in kg/m^3

sig_0 = 0 #density of medium in kg/m^3



m = 4/3 * np.pi * rho ** 3 * ( sig_s - sig_0 )

Permittivity = 8.85 * 10**(-12)

P = 0.5 * c * n_0 * Permittivity    #total power of the LG01 beam


'''
######################################
#5 plot of Q_z vs rho_0x for various w
######################################

rho = 30 * 10 ** (-6)

rho_0 = [0,0]   #no offset


w = [0.65*rho, 0.675*rho, 0.7* rho, 0.725* rho, 0.75*rho, 0.775*rho, 0.8* rho]
#w = [w_0, 0.25* np.sqrt(2)*rho, 0.5 * np.sqrt(2)*rho, 0.75* np.sqrt(2)*rho,np.sqrt(2)*rho, 1.25*np.sqrt(2)*rho, 2*np.sqrt(2)*rho]

rho_00 = np.linspace(-3*w[0], 3*w[0], 100)

rho_01 = np.linspace(-3*w[1], 3*w[1], 100)

rho_02 = np.linspace(-3*w[2], 3*w[2], 100)

rho_03 = np.linspace(-3*w[3], 3*w[3], 100)

rho_04 = np.linspace(-3*w[4], 3*w[4], 100)

rho_05 = np.linspace(-3*w[5], 3*w[5], 100)

rho_06 = np.linspace(-3*w[6], 3*w[6], 100)


Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_00,rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )

Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_01,rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_02,rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )

Axial_flist_vs_d3 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_03,rho_0[1], rho, n_0, n_s, w_0, w[3], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z3 = Axial_flist_vs_d3 * c / ( n_0 * P )

Axial_flist_vs_d4 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_04,rho_0[1], rho, n_0, n_s, w_0, w[4], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z4 = Axial_flist_vs_d4 * c / ( n_0 * P )

Axial_flist_vs_d5 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_05,rho_0[1], rho, n_0, n_s, w_0, w[5], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z5 = Axial_flist_vs_d5 * c / ( n_0 * P )

Axial_flist_vs_d6 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_06,rho_0[1], rho, n_0, n_s, w_0, w[6], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z6 = Axial_flist_vs_d6 * c / ( n_0 * P )

plt.figure(1)


plt.plot(rho_00 / w[0], Q_z0 , lw=2, c="c", label="w/(sqrt(2)rho) = 0.65")
plt.plot(rho_01 / w[1], Q_z1 , lw=2, c="r", label="w/(sqrt(2)rho) = 0.675")
plt.plot(rho_02 / w[2], Q_z2 , lw=2, c="g", label="w/(sqrt(2)rho) = 0.7")
plt.plot(rho_03 / w[3], Q_z3 , lw=2, c="y", label="w/(sqrt(2)rho) = 0.725")
plt.plot(rho_04 / w[4], Q_z4 , lw=2, c="b", label="w/(sqrt(2)rho) = 0.75")
plt.plot(rho_05 / w[5], Q_z5 , lw=2, c="k", label="w/(sqrt(2)rho) = 0.775")
plt.plot(rho_06 / w[6], Q_z6 , lw=2, c="m", label="w/(sqrt(2)rho) = 0.8")

print ((np.sqrt(2)*rho_03 / w[3])[np.argmin(Q_z3)])
print ((np.sqrt(2)*rho_04 / w[4])[np.argmin(Q_z4)])
print ((np.sqrt(2)*rho_05 / w[5])[np.argmin(Q_z5)])
print ((np.sqrt(2)*rho_06 / w[6])[np.argmin(Q_z6)])



new_ticks1 = np.linspace(-3, 3, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1, 0, 6),fontsize=20)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-3))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=15)
plt.title('rho = 30um, w0 = 10um',fontsize=20)

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('see legend', 'related to w', '0','x-aixs rho_0x/(w/sqrt(2)) = -3 to 3', '30')
'''



################################################
#Investigation of axial force change vs offsets
################################################

rho = 30 * 10 ** (-6)

rho_0 = [0,0]   #no offset

w = np.linspace(0.4* np.sqrt(2)*rho, 2*np.sqrt(2)*rho, 50)

rho_00 = np.linspace(-6*rho, 6*rho, 300)

ratio = []

axialmin = []
axialmax = []
diffe = []
w_minloc = []
for we in w:
    axial_force_array = np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_00,rho_0[1], rho, n_0, n_s, w_0, we, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
    
    Q_zmin = axial_force_array  * c / ( n_0 * P )

    axial_min = np.min(Q_zmin)
    
    axialmin.append(axial_min)
    
    min_location = (np.sqrt(2)*rho_00 / we)[np.argmin(Q_zmin)]
    
    w_minloc.append(min_location)
    
    axial_centre_max = TQ.F_total_manual_integration(0, rho_0[1], rho, n_0, n_s, w_0, we, z_R, P, target = "reflective", coordinate = 'z', grid_size = grid_size)['force_total']
    
    Q_zmax = axial_centre_max  * c / ( n_0 * P )
    
    axialmax.append(Q_zmax)
    
    diff = Q_zmax - axial_min
    
    diffe.append(diff)
    fluct_ratio = (Q_zmax - axial_min)/abs(Q_zmax)
    
    ratio.append(fluct_ratio)

print(axialmin)
print(axialmax)
print(diffe)
print(w_minloc)
plt.plot(w/(np.sqrt(2)*rho), ratio , lw=2, c="c", label="w = w_0")

plt.title('rho = 30um, w0 = 10um',fontsize=20)

plt.xlabel('w/(sqrt(2)*rho)',fontsize=20)
plt.ylabel('ratio',fontsize=20)
plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
plt.grid()
plt.show()
 

#####################################
#gradient of Qz vs rho_0x
#####################################
rho = 30 * 10 ** (-6)


rho_0 = [0 , 0]

resolution = 0.01 * 10 ** (-6)

w = np.linspace(0.6*rho, 1.5*rho, 100)

grad_Fzlist = []

for we in w:

    rho_0x_in = 0.5 * we

    grad_Fz = np.mean(MFSAP.Fz_rho_0x_stiffness_at_parameter(rho_0x_in, rho_0[1], rho, n_0, n_s, w_0, we, z_R, P, resolution, target = 'reflective', integration_method = integration_method, grid_size = grid_size))

    grad_Fzlist.append(grad_Fz*10**8)

plt.figure(13)
plt.plot(w/( rho), grad_Fzlist, lw=2, c="c", label="rho_0x/w = 0.5")


new_ticks1 = np.linspace(0.6, 1.5, 10) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-5, 2, 8),fontsize=20)

plt.rc('xtick',labelsize=20)
plt.rc('ytick',labelsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data', 0.6))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=15)

plt.xlabel('w/rho',fontsize=20)
plt.ylabel('grad_Fz(stiffness)10^(-8)',fontsize=20)

#plt.title('grad_Fx vs x offset waist w0 at w = w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

'''
########################################
#6 plot of Q_z vs rho_0x for various rho
########################################

a= 30 * 10 ** (-6)



w = np.sqrt(2)*rho

rho_0 = [0,0]   #no offset

rho_0[0] = np.linspace(-a, a, 100)

rho = [0.5*a, 0.75*a, a]




Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )

Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )

plt.figure(2)


plt.plot(rho_0[0]/a, Q_z0 , lw=2, c="c", label="rho = 15um")
plt.plot(rho_0[0]/a, Q_z1 , lw=2, c="r", label="rho = 22.5um")
plt.plot(rho_0[0]/a, Q_z2 , lw=2, c="g", label="rho = 30um")


new_ticks1 = np.linspace(-1, 1, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.4, 0, 5),fontsize=20)
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

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30', 'related to w', '0','x-aixs -30 to 30', '15, 22.5, 30')
'''


