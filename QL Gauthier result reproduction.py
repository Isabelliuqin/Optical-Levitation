# -*- coding: utf-8 -*-
"""
Created on Fri Jul 10 15:41:44 2020

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
import Module_Gauthier_result_functions as GRF



#laser beam parameters for Gauthier tranparent sphere

w_0 = 2*10 ** (-6)

Lambda = 0.514 * 10**(-6)
P = 20 * 10**(-3)
z_R = np.pi* w_0 ** 2 / Lambda


#trapping item parameters

Rs = 15/2 * w_0

sig_s = 1100 #density of sphere in kg/m^3

sig_0 = 1000 #density of medium in kg/m^3

n_0 = 1.333

n_s = 1.5468 #sphere density in Gauthier


g = 9.8 #gravitational acceleration
c = 3 * 10**8
m = 4/3 * np.pi * Rs**3 * (sig_s - sig_0)



d = np.linspace(0, 1600 * 10**(-6), 100)
axial_forcenet = GRF.Axial_force_total(d, 0, m, Rs, n_0, n_s,w_0, z_R, P)

print (axial_forcenet)
F = [abs(number) for number in axial_forcenet]
print(d[F.index(min(F))])


plt.figure(1)
plt.plot(d * 10**6,axial_forcenet, lw=2, c="c", label="plot figure")

new_ticks1 = np.linspace(0, 1600, 9) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1)
plt.yticks(np.linspace(-20, 180, 11))
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
radial_forcenet = GRF.Radial_force_total(0, a, m, Rs, n_0, n_s,w_0, z_R, P)

print (radial_forcenet)
Fr = [abs(number) for number in radial_forcenet]
print(a[Fr.index(min(Fr))])

#plot out radial force vs offset a

plt.figure(2)
plt.plot(a*10**(6),radial_forcenet, lw=2, c="c", label="plot figure")
new_ticks1 = np.linspace(0, 20, 5) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-20, 140, 9),fontsize=20)
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









#Axial displacement
'''
def target_displacment(y, t, b):
    
    d, vz = y
    dydt = [vz, b/m]
    return dydt
'''
'''
a = 0
y0 = [0.0, 0.0] #initial displacement d0 and velocity vz0

t = np.linspace(0, 20, 301)
d_list = []
vz_list = []
ya = [1,1]
for i in t: 
    if i == 0:
        F_1rz = dblquad(integrand, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (y0[0],a)) #returns a tuple with first element the integral result and second element = upper bound error
        F_znet = F_1rz[0] - m * g
        b = F_znet
        disp = 0.5 * b/m * i**2 + y0[1] * i
        vz = y0[1] + b/m * i
    else:
        F_1rz = dblquad(integrand, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (ya[0],a)) #returns a tuple with first element the integral result and second element = upper bound error
        F_znet = F_1rz[0] - m * g
        b = F_znet
        disp = 0.5 * b/m * i**2 + ya[1] * i
        vz = ya[1] + b/m * i
    
    ya[0] = disp
    ya[1] = vz
    d_list.append(disp)
    vz_list.append(vz)
    
plt.figure(3)    
plt.plot(t, d_list, 'b', label='displacement')
plt.plot(t, vz_list, 'g', label='velocity')
plt.legend(loc='best')
plt.xlabel('t(sec)', fontsize = 20)
plt.ylabel('axial disp(m) and velocity(m/s)', fontsize = 20)
plt.grid()
plt.show()
'''













'''
#Radial displacement

d = 0
r0 = [0.0, 0.0] #initial displacement a0 and velocity vr0

t = np.linspace(0, 20, 301)
a_list = []
vr_list = []
ra = [1,1]
for i in t: 
    if i == 0:
        F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,r0[0])) #returns a tuple with first element the integral result and second element = upper bound error #returns a tuple with first element the integral result and second element = upper bound error
        b = F_1rr[0]
        disp = 0.5 * b/m * i**2 + r0[1] * i
        vr = r0[1] + b/m * i
    else:
        F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,ra[0])) #returns a tuple with first element the integral result and second element = upper bound error #returns a tuple with first element the integral result and second element = upper bound error
        b = F_1rr[0]
        disp = 0.5 * b/m * i**2 + ra[1] * i
        vr = ra[1] + b/m * i
    
    ra[0] = disp
    ra[1] = vr
    a_list.append(disp)
    vr_list.append(vr)
    
plt.figure(4)    
plt.plot(t, a_list, 'b', label='displacement')
plt.plot(t, vr_list, 'g', label='velocity')
plt.legend(loc='best')
plt.xlabel('t(sec)',fontsize = 20)
plt.ylabel('radial disp(m) and veolocity(m/s)',fontsize = 20)
plt.grid()
plt.show()
'''

'''
#plot trajectory
plt.figure(5)  
plt.plot(a_list, d_list)
plt.xlabel('a(m)',fontsize = 20)
plt.ylabel('d(m)',fontsize = 20)
plt.grid()
'''
