# -*- coding: utf-8 -*-
"""
Created on Fri Jun 26 11:58:46 2020

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
d = 0.001*10**(-3)
theta = 0.04
z = 0.8
w_0 = np.sqrt(Lambda*d/(np.pi*theta))
z_R = np.pi*w_0**2/Lambda


w = w_0 * np.sqrt(1 + (z/z_R)**2)

R = z + z_R**2/z

Guoy_phase = np.arctan2(z,z_R)

def E_tilt(x,y,d,theta):
    return (w/w_0)*np.exp(-(1/w**2)*((x + z*theta)**2 + (y + d)**2)) * np.exp(-Guoy_phase*1j)\
        *np.exp(k*z*1j + k/(2*R)*((x + z*theta)**2 + (y + d)**2)*1j)*np.exp(-k*(x*theta + z*theta**2/2)*1j)

def I(x,y):
    E_vortex = E_tilt(X,Y,d,-theta) + E_tilt(X,Y,-d,theta)*np.exp(np.pi*1j)
    Er = E_vortex.real
    Eim = E_vortex.imag
    
    return np.sqrt(Er**2+Eim**2)

def I_analytic(x,y):

    h = k * z
    m = k * (x**2 + y**2) / (2 * R)
    E_g = (w/w_0) * np.exp(-(x**2 + y**2) / w**2) * np.exp(-(Guoy_phase)*1j) * np.exp(h*1j + m*1j)

    a = k * (x * theta * (1 - (z / R)) + d * y / R)
    b = (2 / w**2) * (X * z * theta - d * Y)
    c = k / (2*R)
    n = 1/ w**2 
    f = z**2 * theta**2 + d**2
    g = k * z * theta**2 / 2
    E = 2*1j * np.sin(a - b*1j) * E_g * np.exp((c*1j - n ) * f) * np.exp(-g*1j)
    Er = E.real
    Eim = E.imag
    
    return np.sqrt(Er**2+Eim**2)

def E_linear(x,theta_inc):
    return 20*np.exp((k*np.cos(theta_inc)*z + k*np.sin(theta_inc)*x)*1j)

def inteferogram(x,y,theta_inc):
    E_int = E_linear(x,theta_inc)+ E_tilt(x,y)
    Er = E_int.real
    Eim = E_int.imag
    
    return np.sqrt(Er**2+Eim**2)


theta_inc = 0.01

#Generation of Vortex beam
x = np.arange(-50*10**(-3), 50*10**(-3), 50*10**(-6))
y = np.arange(-50*10**(-3), 50*10**(-3), 50*10**(-6))
X, Y = np.meshgrid(x, y)
plt.imshow(I_analytic(X,Y))
plt.colorbar()


'''

fig = plt.figure()
ax = Axes3D(fig)
x = np.arange(-80*10**(-4), 80*10**(-4), 50*10**(-7))
y = np.arange(-80*10**(-4), 80*10**(-4), 50*10**(-7))
X, Y = np.meshgrid(x, y)
'''




#ax.plot_surface(X, Y, I(X,Y,z,d,theta), cmap=cm.coolwarm, linewidth=0, antialiased=False)
