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
d = 0.05*10**(-3)
theta = 0.04
z = 0.2
w_0 = np.sqrt(Lambda*d/(np.pi*theta))
z_R = np.pi*w_0**2/Lambda


w = w_0 * np.sqrt(1 + (z/z_R)**2)

R = z + z_R**2/z

Guoy_phase = np.arctan(z/z_R)


#Generation of Vortex beam

fig = plt.figure()
ax = Axes3D(fig)
x = np.arange(-50*10**(-3), 50*10**(-3), 50*10**(-6))
y = np.arange(-50*10**(-3), 50*10**(-3), 50*10**(-6))
X, Y = np.meshgrid(x, y)

h = k * z
m = k * (X**2 + Y**2) / (2 * R)
E_g = (w/w_0) * np.exp(-(X**2 + Y**2) / w**2) * np.exp(-(Guoy_phase)*1j) * np.exp(h*1j + m*1j)

a = k * (X * theta * (1 - (z / R)) + d * Y / R)
b = (2 / w**2) * (X * z * theta - d * Y)
c = k / (2*R)
n = 1/ w**2 
f = z**2 * theta**2 + d**2
g = k * z * theta**2 / 2
E = 2*1j * np.sin(a - b*1j) * E_g * np.exp((c*1j - n ) * f) * np.exp(-g*1j)


#E_g1 = np.exp(-(X**2 + Y**2) / w**2)
Er = E.real
Eim = E.imag
I = np.sqrt(Er**2+Eim**2)


ax.plot_surface(X, Y, I, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)