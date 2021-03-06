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


seaborn.set(style='ticks')

#laser beam parameters
w_0 = 2 * 10 ** (-6)
Lambda = 0.514 * 10**(-6)
P = 20 * 10**(-3)
z_R = np.pi* w_0 ** 2 / Lambda


#trapping item parameters
Rs = 0.75 * w_0
sig_s = 1100 #density of sphere in kg/m^3
sig_0 = 1000 #density of medium in kg/m^3
n_s = 1.5468
n_0 = 1.333
g = 9.81 #gravitational acceleration
c = 3 * 10**8
m = 4/3 * np.pi * Rs**3 * (sig_s - sig_0)

def r_s(theta, theta2):
    '''Fresnel coefficient of s-polarisation'''
    return (n_0*np.cos(theta) - n_s*np.cos(theta2))/(n_0*np.cos(theta) + n_s*np.cos(theta2))

def r_p(theta, theta2):
    '''Fresnel coefficient of p-polarisation'''
    return (n_0 * np.cos(theta2) - n_s * np.cos(theta))/ (n_s * np.cos(theta) + n_0 * np.cos(theta2))

def rfcoe_sm(theta,theta2):
    '''Average reflectance'''
    return ((r_s(theta,theta2) + r_p(theta,theta2))/2)**2

def rfcoe(theta, theta2):
    '''reflectance given in paper'''
    return (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta2))**2)**2/ \
    (n_0*n_s*((np.cos(theta))**2 + (np.cos(theta2))**2) + (n_0**2 + n_s**2) * np.cos(theta)*np.cos(theta2))**2


def I_GB(rou, z):
    '''GB intensity profile'''
    #Gaussian beam propagation parameters
    w = w_0 * np.sqrt(1 + (z/z_R)**2)
    #R = z + z_R**2/z
    return (2 * P/(np.pi * w**2)) * np.exp( - 2 * rou ** 2 / w ** 2)

def I_LG01(rou,z):
    '''LG01 intensity profile'''
    w = w_0 * np.sqrt(1 + (z/z_R)**2)
    return 2/(np.pi * w**2) * 2 * rou**2/w**2* np.exp(-2*rou**2/w**2)

#Axial Force calculation
def integrand(theta,d):
    
    '''integrand of lower surface reflection'''
    
    rou = Rs * np.sin(theta) #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta) * np.cos(theta_2))**2
    
    return (np.pi/c) * n_0 * (1 + np.cos(2*theta)) * \
        I_GB(rou, z) * Reflection_coe * Rs**2 * np.sin(2*theta)

def integrand_1tz(theta,d):
    
    '''integrand of lower surface transmission'''
    
    rou = Rs * np.sin(theta) #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta)*np.cos(theta_2))**2
    
    Transmission_coe = 1 - Reflection_coe
    
    return (np.pi/c) * (n_0 - n_s*np.cos(theta - theta_2)) * \
        I_GB(rou,z) * Transmission_coe *Rs**2 * np.sin(2*theta)

def integrand_2rz(theta,d):
    
    '''integrand of upper surface reflection'''
    
    rou = Rs * np.sin(theta) #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta)*np.cos(theta_2))**2
    
    Transmission_coe = 1 - Reflection_coe
    
    return (np.pi/c) * n_s * (np.cos(theta - theta_2) + np.cos(3*theta_2 - theta)) * \
        I_GB(rou,z) * Reflection_coe * Transmission_coe * Rs**2 * np.sin(2*theta)

def integrand_2tz(theta,d):
    
    '''integrand of upper surface transmission'''
    
    rou = Rs * np.sin(theta) #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta)*np.cos(theta_2))**2
    
    Transmission_coe = 1 - Reflection_coe
    
    return (np.pi/c) * (n_s * np.cos(theta - theta_2) - n_0 * np.cos(2*(theta - theta_2)))* \
    I_GB(rou,z) * Transmission_coe * Transmission_coe * Rs**2 * np.sin(2*theta)

'''
#Calculation of axial force integral
d = np.linspace(0, 835 * 10**(-6), 400) 
forcenet = []
for d_e in d:
    F_z = quad(integrand, 0, np.pi/2, args = (d_e)) #returns a tuple with first element the integral result and second element = upper bound error
    F_1tz = quad(integrand_1tz, 0, np.pi/2, args = (d_e))
    F_2rz = quad(integrand_2rz, 0, np.pi/2, args = (d_e))
    F_2tz = quad(integrand_2tz, 0, np.pi/2, args = (d_e))
    
    
    F_znet = F_z[0] + F_1tz[0] + F_2rz[0] + F_2tz[0] - m * g
    forcenet.append(F_znet * 10 ** 12)

print (forcenet)
'''
'''
plt.figure(1)
plt.plot(d * 10**6,forcenet, lw=2, c="c", label="plot figure")

new_ticks1 = np.linspace(0, 900, 10) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1)
plt.yticks(np.linspace(-2, 5, 8))
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.xlabel('d(um)',fontsize=15)
plt.ylabel('F_axial(pN)',fontsize=15)
plt.show()
'''


'''
d = 825 * 10 **(-6)
F_z = quad(integrand, 0, np.pi/2, args = (d)) #returns a tuple with first element the integral result and second element = upper bound error
F_1tz = quad(integrand_1tz, 0, np.pi/2, args = (d))
F_2rz = quad(integrand_2rz, 0, np.pi/2, args = (d))
F_2tz = quad(integrand_2tz, 0, np.pi/2, args = (d))
    
    
F_znet = F_z[0] + F_1tz[0] + F_2rz[0] + F_2tz[0] - 4/3 * np.pi * Rs**3 * (sig_s - sig_0) * g

#F_znet = F_z[0] - 4/3 * np.pi * Rs**3 * (sig_s - sig_0) * g
print (F_znet)
'''  


#Axial displacement
'''
def target_displacment(y, t, b):
    
    d, vz = y
    dydt = [vz, b/m]
    return dydt
'''

'''
y0 = [0.0, 0.0] #initial displacement d0 and velocity vz0

t = np.linspace(0, 10, 101)
d_list = []
vz_list = []
for i in t: 
    F_z = quad(integrand, 0, np.pi/2, args = (y0[0])) #returns a tuple with first element the integral result and second element = upper bound error
    F_znet = F_z[0] - m * g

    b = F_znet
    disp = 0.5 * b/m * i**2 + y0[1] * i
    vz = y0[1] + b/m * i
    y0[0] = disp
    y0[1] = vz
    d_list.append(disp)
    vz_list.append(vz)
    
    
plt.plot(t, d_list, 'b', label='displacement')
#plt.plot(t, vz_list, 'g', label='velocity')
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
'''



#Radial force

def integrand_1rr(theta,phi,d, a):
    
    '''integrand of lower surface reflection'''
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta) * np.cos(theta_2))**2
    
    return (- n_0/(2*c) * np.sin(2*theta) * I_GB(rou,z) * Reflection_coe * Rs**2 * np.cos(phi) * np.sin(2*theta))

def integrand_1tr(theta,phi, d, a):
    
    '''integrand of lower surface transmission'''
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta)*np.cos(theta_2))**2
    
    Transmission_coe = 1 - Reflection_coe
    
    return n_s/(2*c) * np.sin(theta - theta_2) * I_GB(rou,z) * Transmission_coe * Rs**2 * np.cos(phi) * np.sin(2*theta)

def integrand_2rr(theta, phi, d, a):
    
    '''integrand of upper surface reflection'''
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta)*np.cos(theta_2))**2
    
    Transmission_coe = 1 - Reflection_coe
    
    return n_s/(2*c) * (np.sin(3*theta_2 - theta) - np.sin(theta - theta_2)) * I_GB(rou,z) *\
        Reflection_coe * Transmission_coe * Rs**2 * np.cos(phi) * np.sin(2*theta)

def integrand_2tr(theta, phi, d, a):
    
    '''integrand of upper surface transmission'''
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    #Reflection_coe = (n_0 * n_s)**2 * ((np.cos(theta))**2 - (np.cos(theta_2))**2)**2/ \
    #(n_0*n_s*((np.cos(theta))**2 + (np.cos(theta_2))**2) + (n_0**2 + n_s**2) * np.cos(theta)*np.cos(theta_2))**2
    
    Transmission_coe = 1 - Reflection_coe
    
    return 1/(2*c) * (n_0 * np.sin(2 * (theta - theta_2)) - n_s * np.sin(theta - theta_2)) * I_GB(rou,z) * \
        Transmission_coe * Transmission_coe * Rs**2 * np.cos(phi) * np.sin(2*theta)

#Calculation of radial force integral
d = 0


a = np.linspace(0, 16*10**(-6), 100) 
forcenet_r = []
for a_e in a:
    F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e)) #returns a tuple with first element the integral result and second element = upper bound error
    F_1tr = dblquad(integrand_1tr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e))
    F_2rr = dblquad(integrand_2rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e))
    F_2tr = dblquad(integrand_2tr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e))
    
    
    F_znet = F_1rr[0] + F_1tr[0] + F_2rr[0] + F_2tr[0]
    forcenet_r.append(F_znet * 10 ** 12)

#print (forcenet_r)
    
#plot out radial force vs offset a

plt.figure(2)
plt.plot(a*10**(6),forcenet_r, ls="-.", lw=2, c="c", label="plot figure")
new_ticks1 = np.linspace(0, 16, 9) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1)
plt.yticks(np.linspace(-10, 5, 16))
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.xlabel('a(um)',fontsize=15)
plt.ylabel('F_radial(pN)',fontsize=15)
plt.show()

'''
#Radial displacement
r0 = [0.0, 0.0] #initial displacement a0 and velocity vr0

t = np.linspace(0, 10, 101)
a_list = []
vr_list = []
for i in t: 
    F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,r0[0])) #returns a tuple with first element the integral result and second element = upper bound error #returns a tuple with first element the integral result and second element = upper bound error
    b = F_1rr[0]
    
    disp = 0.5 * b/m * i**2 + r0[1] * i
    vr = r0[1] + b/m * i
    r0[0] = disp
    r0[1] = vr
    a_list.append(disp)
    vr_list.append(vr)
    
    
#plt.plot(t, a_list, 'b', label='displacement')
#plt.plot(t, vr_list, 'g', label='velocity')
#plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()
'''

#plot trajectory
#plt.plot(a_list, d_list)


