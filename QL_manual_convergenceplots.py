# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 10:30:53 2020

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
# integration method control: Fx, Fy, Fz convergence plots
##########################################################



integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 1000

plt.close('all')

###########################

t_start = time.time()


#Our sphere

g = 9.8 #gravitational acceleration
c = 3 * 10**8
#m = 4/3 * np.pi * Rs**3 * (sig_s - sig_0)

#TEM01* reflective target Table 6 matching

w_0 = 2.67 * 10 ** (-6)


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
##########################################
#Fx convergence
##########################################

rho_0 = [w_0, 0]
w = w_0






Fx_list = []
timex_list = []

for i in range(10, grid_size, 20):
    
    start_time = time.time()
    
    integration_results = TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'x', grid_size = i)
            
    Fx = integration_results['force_total']
    
    #Fx =  TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'x', grid_size = i)
    
    Fx_list.append(Fx)
    
    
    
    timex_list.append(time.time() - start_time)
    print (i)
    
print(Fx_list)

print(timex_list)



number_of_grid = range(30, grid_size, 20)

plt.figure(1)

plt.scatter(number_of_grid, np.asarray(Fx_list[1:]) * 10 ** 12, label="Fx", marker = "^")

plt.axhline(y=TQ.Fx_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")[0] * 10 ** 12, color='r', linestyle='-', label = "quad built-in")
plt.legend(loc=4,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Fx(pN)',fontsize=20)
plt.grid()
plt.show()


plt.figure(2)
plt.scatter(number_of_grid, np.asarray(timex_list[1:]), label="Execution time",marker = "8")

plt.legend(loc=2,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('x Runningtime(sec)',fontsize=20)
plt.grid()
plt.show()






################################################
#Fy convergence
################################################

rho_0 = [0, w_0]
w = w_0

start_time = time.time()





Fy_list = []
timey_list = []

for i in range(10, grid_size, 20):
    
    start_time = time.time()
    
    integration_results = TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'y', grid_size = i)
            
    Fy = integration_results['force_total']
    
    Fy_list.append(Fy)
    
    
    
    timey_list.append(time.time() - start_time)
    print (i)
    
print(Fy_list)

print(timey_list)



number_of_grid = range(30, grid_size, 20)

plt.figure(3)

plt.scatter(number_of_grid, np.asarray(Fy_list[1:]) * 10 ** 12, label="Fy", marker = "^")

plt.axhline(y=TQ.Fy_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")[0] * 10 ** 12, color='r', linestyle='-', label = "quad built-in")
plt.legend(loc=4,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Fy(pN)',fontsize=20)
plt.grid()
plt.show()


plt.figure(4)

plt.scatter(number_of_grid, np.asarray(timey_list[1:]), label="Execution time",marker = "8")

plt.legend(loc=2,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('y Runningtime(sec)',fontsize=20)
plt.grid()
plt.show()





##################################################
#Fz convergence
##################################################


rho_0 = [0,0]
w = w_0

#gridsize = 1005



Fz_list = []
time_list = []

for i in range(10, grid_size, 20):
    
    start_time = time.time()
    
    integration_results =  TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'z', grid_size = i)
            
    Fz = integration_results['force_total']
    
    #Fz =  TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'z', grid_size = i)
    
    Fz_list.append(Fz)
    
    
    
    time_list.append(time.time() - start_time)
    print (i)
    
print(Fz_list)

print(time_list)



number_of_grid = range(30, grid_size, 20)

plt.figure(5)

plt.scatter(number_of_grid, np.asarray(Fz_list[1:]) * 10 ** 12, label="Fz", marker = "^")

plt.axhline(y=TQ.Fz_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")[0] * 10 ** 12, color='r', linestyle='-', label = "quad built-in")
plt.legend(loc=4,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Fz(pN)',fontsize=20)
plt.grid()
plt.show()


plt.figure(6)

plt.scatter(number_of_grid, np.asarray(time_list[1:]), label="Execution time",marker = "8")

plt.legend(loc=4,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('z Runningtime(sec)',fontsize=20)
plt.grid()
plt.show()


#print (MIM.Fz_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

#print("--- %s seconds ---" % (time.time() - start_time))



start_time = time.time()

print (TQ.Fz_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective"))

print("--- %s seconds ---" % (time.time() - start_time))
'''

'''
###########################
#check intensity looks nice
###########################

rho_0 = [0, 0]
w = 5 * w_0

gridsize = 100

theta = np.linspace(0, np.pi/2, gridsize)
phi = np.linspace(0, np.pi * 2, gridsize)


Theta,Phi = np.meshgrid(theta, phi)

x = rho * np.sin(Theta) * np.cos(Phi)#np.linspace(-5 * 10 ** (-6), 5 * 10 ** (-6), gridsize)

y = rho * np.sin(Theta) * np.sin(Phi)#np.linspace(-5 * 10 ** (-6), 5 * 10 ** (-6), gridsize)

#X, Y = np.meshgrid(x,y)

#x = rho * np.sin(Theta) * np.cos(Phi)

#y = rho * np.sin(theta) * np.sin(Phi)

r = np.sqrt( (rho * np.sin(Theta) * np.cos(Phi) - rho_0[0]) ** 2 + (rho * np.sin(Theta) * np.sin(Phi) - rho_0[1]) ** 2) #represent r(the parameter of laser beam) by theta

I = MIM.LG_01_Intensity(r,w)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(x * 10 ** 6, y * 10 ** 6, I, c='r', marker='o')

ax.set_xlabel('x(um)', fontsize = 15)
ax.set_ylabel('y(um)', fontsize = 15)
ax.set_zlabel('LG_01 intensity', fontsize = 15)

plt.show()

'''

###############################
#power of LG_01 beam
###############################

rho_0 = [w_0, 0]
w = 70*w_0

gridsize = 800
#rrange = 7 * 10 ** (-6)

P_list = []


for i in range(5, 50, 2):
    
    start_time = time.time()
    
    P =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, rrange = i* 10 ** (-6), target = "reflective")
    
    P_list.append(P)
    
    
    
    
    print (i)
    
print(P_list)



#print(MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, rrange, target = "reflective"))

P = 0.5 * c * n_0 * Permittivity


rmax = range(5, 50, 2)

plt.figure(7)

plt.scatter(rmax, np.asarray(P_list), label="P", marker = "^")

plt.axhline(y= P, color='r', linestyle='-', label = "analytic solution")

plt.title('Power vs the range of r taken, number of grids 800, rmax from 5um to 50um, w = 70w0, offset x = w0',fontsize=20)
plt.legend(loc=4,fontsize=20)

plt.xlabel('rmax(um)',fontsize=20)
plt.ylabel('P(pN)',fontsize=20)
plt.grid()
plt.show()
