# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 12:57:39 2020

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

import time

##########################################################
# Integration control: Fz 2D force surface plots
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

sig_s = 10.49 * 10 ** 3 * ( 3 ** 3 - 2.9 ** 3 / 3 ** 3 ) #density of sphere in kg/m^3

sig_0 = 0 #density of medium in kg/m^3



m = 4/3 * np.pi * rho ** 3 * ( sig_s - sig_0 )

Permittivity = 8.85 * 10**(-12)

P = 0.5 * c * n_0 * Permittivity    #total power of the LG01 beam


'''
####################################################
#Qz vs 2D w and offset full range
####################################################

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*np.sqrt(2)*rho, 40)


#rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0, 0]

rho_0[0] = np.linspace(-3*rho, 3*rho, 40)

w,rho_0x = np.meshgrid(w, rho_0[0])


Fzarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 
Qzarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 
        
for i in range(len(w)):
            
    for k in range(len(rho_0x)):        


        Fzarray[k][i] = TQ.F_total_manual_integration(rho_0x[k][i],rho_0[1], rho, n_0, n_s, w_0, w[k][i], z_R, P, target = 'reflective', coordinate = 'z', grid_size = grid_size)['force_total']
            
        

        Qzarray[k][i] = Fzarray[k][i] * c / ( n_0 * P )

fig = plt.figure() 
axes = fig.gca(projection ='3d') 
axes.plot_surface(w/(np.sqrt(2)*rho), rho_0x/(rho), np.array(Qzarray),cmap='viridis', edgecolor='none') 

plt.xlabel('w/(sqrt(2)rho)',fontsize=20)
plt.ylabel('rho_0x/rho',fontsize=20)
axes.set_zlabel('Qz', fontsize=20)

plt.title('rho = 30um, w0 = 10um',fontsize=20)
plt.grid()
plt.show() 




####################################################
#Qx vs 2D w and offset full range
####################################################

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*np.sqrt(2)*rho, 40)


#rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0, 0]

rho_0[0] = np.linspace(-3*rho, 3*rho, 40)

w,rho_0x = np.meshgrid(w, rho_0[0])


Fxarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 
Qxarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 


for i in range(len(w)):
            
    for k in range(len(rho_0x)):        


        Fxarray[k][i] = TQ.F_total_manual_integration(rho_0x[k][i],rho_0[1], rho, n_0, n_s, w_0, w[k][i], z_R, P, target = 'reflective', coordinate = 'x', grid_size = grid_size)['force_total']
            
        Qxarray[k][i] = Fxarray[k][i] * c / ( n_0 * P )



fig = plt.figure() 
axes = fig.gca(projection ='3d') 
axes.plot_surface(w/(np.sqrt(2)*rho), rho_0x/(rho), np.array(Qxarray),cmap='viridis', edgecolor='none') 

plt.xlabel('w/(sqrt(2)rho)',fontsize=20)
plt.ylabel('rho_0x/rho',fontsize=20)
axes.set_zlabel('Qz', fontsize=20)

plt.title('rho = 30um, w0 = 10um',fontsize=20)
plt.grid()
plt.show() 

'''



####################################################
#Qz vs 2D w and offset valid range
####################################################

rho = 30 * 10 ** (-6)
w = np.linspace(np.sqrt(2)*rho, 2*np.sqrt(2)*rho, 40)


#rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0, 0]

rho_0[0] = np.linspace(-np.sqrt(2)*rho, 3*rho, 40)

w,rho_0x = np.meshgrid(w, rho_0[0])


Fzarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 
Qzarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 
        
for i in range(len(w)):
            
    for k in range(len(rho_0x)):        


        Fzarray[k][i] = TQ.F_total_manual_integration(rho_0x[k][i],rho_0[1], rho, n_0, n_s, w_0, w[k][i], z_R, P, target = 'reflective', coordinate = 'z', grid_size = grid_size)['force_total']
            
        

        Qzarray[k][i] = Fzarray[k][i] * c / ( n_0 * P )

fig = plt.figure() 
axes = fig.gca(projection ='3d') 
axes.plot_surface(w/(np.sqrt(2)*rho), rho_0x/(rho), np.array(Qzarray),cmap='viridis', edgecolor='none') 

plt.xlabel('w/(sqrt(2)rho)',fontsize=20)
plt.ylabel('rho_0x/rho',fontsize=20)
axes.set_zlabel('Qz', fontsize=20)

plt.title('rho = 30um, w0 = 10um',fontsize=20)
plt.grid()
plt.show() 


####################################################
#Qx vs 2D w and offset valid range
####################################################

rho = 30 * 10 ** (-6)
w = np.linspace(np.sqrt(2)*rho, 2*np.sqrt(2)*rho, 40)


#rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0, 0]

rho_0[0] = np.linspace(-np.sqrt(2)*rho, np.sqrt(2)*rho, 40)

w,rho_0x = np.meshgrid(w, rho_0[0])


Fxarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 
Qxarray = [ [ None for p in range( len(w) ) ]
               for q in range( len(rho_0x) ) ] 


for i in range(len(w)):
            
    for k in range(len(rho_0x)):        


        Fxarray[k][i] = TQ.F_total_manual_integration(rho_0x[k][i],rho_0[1], rho, n_0, n_s, w_0, w[k][i], z_R, P, target = 'reflective', coordinate = 'x', grid_size = grid_size)['force_total']
            
        Qxarray[k][i] = Fxarray[k][i] * c / ( n_0 * P )



fig = plt.figure() 
axes = fig.gca(projection ='3d') 
axes.plot_surface(w/(np.sqrt(2)*rho), rho_0x/(rho), np.array(Qxarray),cmap='viridis', edgecolor='none') 

plt.xlabel('w/(sqrt(2)rho)',fontsize=20)
plt.ylabel('rho_0x/rho',fontsize=20)
axes.set_zlabel('Qz', fontsize=20)

plt.title('rho = 30um, w0 = 10um',fontsize=20)
plt.grid()
plt.show() 