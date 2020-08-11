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

'''
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



#laser beam parameters for Roosen table 1.B tranparent sphere

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
#individual force term reproduction of roosen 1979
a  = 4* 10 ** (-6)
#d = 55.95 * 10**(-6) 
axial_forcenet = RF.Axial_force_total(a, Rs, n_0, n_s, w_0, w, z_R, P)

#axial_force_Ashkin = RF.axial_force_Ashkin(theta,phi, d, a, Rs, n_0, w_0, z_R, P)

print (axial_forcenet[0])
#F = [abs(number) for number in axial_forcenet]
#print(d[#

radial_forcenet = RF.Radial_force_total(a, Rs, n_0, n_s, w_0, w, z_R, P)

print (radial_forcenet[0])
'''


#plot reproduction of Fz vs offset rho0 Roosen 1979
a = np.linspace(0, 20 * 10** (-6),100)
axial_flist = RF.Axial_force_total_vs_offset_plot(a,Rs, n_0, n_s, w_0, w, z_R, P)

#Roosen reflective sphere 1979 table 6 cylindrical beam

a_roosen = np.linspace(0,16,9)
Fz_roosen = [-8.31, -8.42, -8.74, -9.20, -9.69, -10.10, -10.33, -10.29, -9.94]

plt.figure(1)
p1 = plt.scatter(a_roosen, Fz_roosen, label = "Roosen 1979" )
p2 = plt.plot(a * 10** (6),axial_flist, lw=2, c="c", label="My result")

new_ticks1 = np.linspace(0, 20, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-8, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('offset rho_0(um)',fontsize=20)
plt.ylabel('Fz(10^(-10)N)',fontsize=20)

plt.grid()
plt.show()



#plot reproduction of Fy vs offset rho0 Roosen 1979
radial_flist = RF.Radial_force_total_vs_offset_plot(a, Rs, n_0, n_s, w_0, w, z_R, P)
#plot out radial force vs offset a

#Roosen reflective sphere 1979 table 6 cylindrical beam

a_roosen = np.linspace(0,16,9)
Fy_roosen = [0, 0.66,1.18, 1.45, 1.41, 1.06, 0.46, -0.28, -1.04]


plt.figure(2)
plt.plot( a*10** (6), radial_flist, lw=2, c="c", label="My result")
plt.scatter(a_roosen, Fy_roosen,label="Roosen 1979")

new_ticks1 = np.linspace(0, 20, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(0, 4, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('offset rho_0(um)',fontsize=20)
plt.ylabel('Fy(10^(-10)N)',fontsize=20)
plt.grid()
plt.show()


#plot reproduction of Fx vs offset rho0 Roosen 1979
radial_flistx = RF.Radial_force_total_vs_offsetx_plot(a, Rs, n_0, n_s, w_0, w, z_R, P)
#plot out radial force vs offset a

#Roosen reflective sphere 1979 table 6 cylindrical beam

a_roosen = np.linspace(0,16,9)
Fy_roosen = [0, 0.66,1.18, 1.45, 1.41, 1.06, 0.46, -0.28, -1.04]


plt.figure(3)
plt.plot( a*10** (6), radial_flistx, lw=2, c="c", label="My result")
plt.scatter(a_roosen, Fy_roosen,label="Roosen 1979")

new_ticks1 = np.linspace(0, 20, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(0, 4, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('offset rho_0x(um)',fontsize=20)
plt.ylabel('Fx(10^(-10)N)',fontsize=20)
plt.grid()
plt.show()

