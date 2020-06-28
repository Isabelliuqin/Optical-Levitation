# -*- coding: utf-8 -*-
"""
Created on Sat Jun 27 12:03:12 2020

@author: liuqi
"""

import numpy as np
import cmath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

#use the data from Vaity 2013
Lambda = 6.32*10**(-7)
k = 2*np.pi/Lambda
d = 0.05*10**(-3)
theta = 0.04
z = 0.2
w_0 = np.sqrt(Lambda*d/(np.pi*theta))
z_R = np.pi*w_0**2/Lambda


w = w_0 * np.sqrt(1 + (z/z_R)**2)

R = z + z_R**2/z




def factorial(n):
    if n == 0:
        return 1
    else:
        return n * factorial(n-1)

def Lagpolynomial(p,l,x):
    if p == 0:
        return 1
    if p == 1:
        return (l+1)-x
    if p == 2:
        return (l+1)*(l+2)/factorial(2) - (l+2)*x + x**2/factorial(2)
    else:
        return 0

def LG_pl(p,l,x,y):
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    Guoy_phase = (2*p + l + 1)*np.arctan2(z,z_R)
    Kpl = (-1)**p*np.sqrt(2/np.pi)*np.sqrt(factorial(p)/factorial(p+l))
    
    return Kpl/w * Lagpolynomial(p,l,2*r**2/w**2)*(np.sqrt(2)*r/w)**l*np.exp(-l*phi*1j)*np.exp(-r**2/w**2)*np.exp(k*r**2/(2*R)*1j)*np.exp(-Guoy_phase*1j)

def LG_phase(p,l,x,y):
    r = np.sqrt(x**2 + y**2)
    phi_m = np.arctan2(y, x)
    Guoy_phase = (2*p + l + 1)*np.arctan2(z,z_R)
    phase = l*phi_m + k*r**2/(2*R) - Guoy_phase
    #phase = k*r**2/(2*R)
    
    while (np.ptp(phase) > 2*np.pi):
    
        phase = np.where(phase <= np.pi, phase, phase - 2*np.pi)
        phase = np.where(phase >= -np.pi, phase, phase + 2*np.pi)
        
    
    else:
 
        return phase
def E_linear(x,theta_inc):
    return 20*np.exp((k*np.cos(theta_inc)*z + k*np.sin(theta_inc)*x)*1j)

def inteferogram(x,y,p,l,theta_inc):
    E_int = E_linear(x,theta_inc)+ LG_pl(p,l,x,y)
    Er = E_int.real
    Eim = E_int.imag
    
    return np.sqrt(Er**2+Eim**2)
           
    
p = 0
l = 1
theta_inc = 0

#inteferogram plot to see fork pattern
x = np.arange(-50*10**(-5), 50*10**(-5), 50*10**(-8))
y = np.arange(-50*10**(-5), 50*10**(-5), 50*10**(-8))
X, Y = np.meshgrid(x, y)
plt.imshow(inteferogram(X,Y,p,l,theta_inc))
plt.colorbar()


'''
#close look
x = np.arange(-50*10**(-6), 50*10**(-6), 50*10**(-8))
y = np.arange(-50*10**(-6), 50*10**(-6), 50*10**(-8))
X, Y = np.meshgrid(x, y)
plt.imshow(LG_phase(p, l, X, Y))
plt.colorbar()
'''

'''
#further away from centre

x = np.arange(-20*10**(-5), 20*10**(-5), 50*10**(-8))
y = np.arange(-20*10**(-5), 20*10**(-5), 50*10**(-8))
X, Y = np.meshgrid(x, y)
plt.imshow(LG_phase(p,l,X,Y))
plt.colorbar()
'''
'''
fig = plt.figure()
ax = Axes3D(fig)
x = np.arange(-50*10**(-4), 50*10**(-4), 50*10**(-7))
y = np.arange(-50*10**(-4), 50*10**(-4), 50*10**(-7))
X, Y = np.meshgrid(x, y)





#Kpl = (-1)**p*np.sqrt(2/np.pi)*np.sqrt(factorial(p)/factorial(p+l))

#LG_pl = Kpl/w * Lagpolynomial(p,l,2*r**2/w**2)*(np.sqrt(2)*r/w)**l*np.exp(-l*phi*1j)*np.exp(-r**2/w**2)*np.exp(k*r**2/(2*R)*1j)*np.exp(-Guoy_phase*1j)

#E = (w_0 * np.sqrt(2)/w**2)*(X + Y*1j)*np.exp(-r**2/w**2)*np.exp(k*r**2/(2*R)*1j)*np.exp(-2*np.arctan(z/z_R)*1j)

#phase = l*phi 
#+ k*r**2/(2*R) - Guoy_phase
Er = LG_pl(p,l,X,Y).real
Eim = LG_pl(p,l,X,Y).imag
I = np.sqrt(Er**2+Eim**2)


ax.plot_surface(X, Y, LG_phase(p,l,X,Y), cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
'''