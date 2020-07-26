# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 21:39:12 2020

@author: liuqi
"""


import scipy as sp
import numpy as np
import matplotlib.pylab as plt
import scipy.integrate as spi
from scipy.integrate import quad
import seaborn
from scipy.integrate import odeint
from scipy.integrate import dblquad
import Module_Roosen_function as RF


#TEM01* reflective target Table 6 matching

w_0 = 2.5* 10 ** (-6)
w = 20 * 10 ** (-6)

Lambda = 0.5145 * 10**(-6)
P = 1
z_R = np.pi* w_0 ** 2 / Lambda
Rs = 15 * 10 ** (-6)
 
weight = 3.4 * 10**(-10)

n_0 = 1

n_s_n = 0.05

k =  2.87 

n_s = n_s_n - k*1j


'''
#laser beam parameters for Gauthier tranparent sphere

w_0 = 2.5* 10 ** (-6)
w = 20 * 10 ** (-6)

Lambda = 0.5145 * 10**(-6)
P = 1
z_R = np.pi* w_0 ** 2 / Lambda
Rs = 15 * 10 ** (-6)
 
weight = 3.4 * 10**(-10)

#trapping item parameters

#Rs = 8.5/2 * w_0

#sig_s = 1100 #density of sphere in kg/m^3

#sig_0 = 1000 #density of medium in kg/m^3

n_0 = 1

n_s = 1.5#sphere density in Gauthier


g = 9.8 #gravitational acceleration
c = 3 * 10**8
#m = 4/3 * np.pi * Rs**3 * (sig_s - sig_0)
'''




a  = 8* 10 ** (-6)
#d = 55.95 * 10**(-6) 
axial_forcenet = RF.Axial_force_total(a, Rs, n_0, n_s, w_0, w, z_R, P)

#axial_force_Ashkin = RF.axial_force_Ashkin(theta,phi, d, a, Rs, n_0, w_0, z_R, P)

print (axial_forcenet[0])
#F = [abs(number) for number in axial_forcenet]
#print(d[#

radial_forcenet = RF.Radial_force_total(a, Rs, n_0, n_s, w_0, w, z_R, P)

print (radial_forcenet[0])

'''
a = 0
Rs = np.linspace(0,2*w,100)
axial_flist = RF.Axial_force_total_plot(a,Rs, n_0, n_s, w_0, w, z_R, P)



plt.figure(1)
plt.plot(Rs/w,np.log10(axial_flist), lw=2, c="c", label="plot figure")

new_ticks1 = np.linspace(0, 2, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1)
plt.yticks(np.linspace(-3, 1, 5))
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

ax.legend(['n_s = 1.5 offset rho_0 = 0'], fontsize=13)

plt.xlabel('Rs/w',fontsize=20)
plt.ylabel('log10(F_axial(nN))',fontsize=20)

plt.grid()
plt.show()





a = np.linspace(0.1*w,2.6*w,100)
w = 10 * 10 ** (-6)
Rs = 1.5* w
radial_flist = RF.Radial_force_total_plot(a, Rs, n_0, n_s, w_0, w, z_R, P)
#plot out radial force vs offset a

plt.figure(2)
plt.plot(a/w,np.log10(radial_flist), lw=2, c="c", label="plot figure")
new_ticks1 = [0.1, 1, 2] # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-3, 1, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

ax.legend(['n_s = 1.5 Rs/w = 1.5'], fontsize=13)

plt.xlabel('offset rho_0/w',fontsize=20)
plt.ylabel('log10(F_radial(nN))',fontsize=20)
plt.grid()
plt.show()
'''

