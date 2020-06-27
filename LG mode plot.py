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

Guoy_phase = np.arctan(z/z_R)

fig = plt.figure()
ax = Axes3D(fig)
x = np.arange(-50*10**(-4), 50*10**(-4), 50*10**(-7))
y = np.arange(-50*10**(-4), 50*10**(-4), 50*10**(-7))
X, Y = np.meshgrid(x, y)


r = np.sqrt(X**2 + Y**2)

E = (w_0 * np.sqrt(2)/w**2)*(X + Y*1j)*np.exp(-r**2/w**2)*np.exp(k*r**2/(2*R)*1j)*np.exp(-2*np.arctan(z/z_R)*1j)

Er = E.real
Eim = E.imag
I = np.sqrt(Er**2+Eim**2)


ax.plot_surface(X, Y, I, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
