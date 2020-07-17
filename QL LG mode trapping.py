# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:46:02 2020

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
import Module_LGmode_function as GLF

#laser beam parameters for Gauthier tranparent sphere

w_0 = 2*10 ** (-6)

Lambda = 0.514 * 10**(-6)
P = 20 * 10**(-3)
z_R = np.pi* w_0 ** 2 / Lambda


#trapping item parameters

Rs = 15/2*10**(-6)

sig_s = 1100 #density of sphere in kg/m^3

sig_0 = 1000

n_0 = 1.333

#n_s = 1 #sphere density in Gauthier


g = 9.8 #gravitational acceleration
c = 3 * 10**8
m = 4/3 * np.pi * Rs**3 * (sig_s - sig_0)


d = np.linspace(0, 1600 * 10**(-6), 100)
axial_forcenet = GLF.Axial_force_total(d, 0, m, Rs, n_0, w_0, z_R)

print (axial_forcenet)
F = [abs(number) for number in axial_forcenet]
print(d[F.index(min(F))])


plt.figure(1)
plt.plot(d * 10**6,axial_forcenet, lw=2, c="c", label="plot figure")

new_ticks1 = np.linspace(0, 1600, 9) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1)
plt.yticks(np.linspace(-5, 1, 7))
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.xlabel('d(um)',fontsize=20)
plt.ylabel('F_axial(pN)',fontsize=20)
plt.grid()
plt.show()


a = np.linspace(0,20*10**(-6),100)
radial_forcenet = GLF.Radial_force_total(0, a, m, Rs, n_0, w_0, z_R)

print (radial_forcenet)
Fr = [abs(number) for number in radial_forcenet]
print(a[Fr.index(min(Fr))])

#plot out radial force vs offset a

plt.figure(2)
plt.plot(a*10**(6),radial_forcenet, lw=2, c="c", label="plot figure")
new_ticks1 = np.linspace(0, 20, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-5, 1, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.xlabel('a(um)',fontsize=20)
plt.ylabel('F_radial(pN)',fontsize=20)
plt.grid()
plt.show()



a = 0
y0 = [0.0, 0.0] #initial displacement d0 and velocity vz0

t = np.linspace(0, 20, 301)

d = 0
r0 = [0.0, 0.0] #initial displacement a0 and velocity vr0




axial_disp, axial_vel, radial_disp, radial_vel = GLF.disp_velo(d,a,m, Rs, n_0, w_0, z_R, t, y0, r0)

#plot trajectory
plt.figure(5)  
plt.plot(axial_disp, radial_disp)
plt.xlabel('a(m)',fontsize = 20)
plt.ylabel('d(m)',fontsize = 20)
plt.grid()