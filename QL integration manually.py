# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 16:52:31 2020

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

import Module_integration_manually as MIM

import Module_Trialplot_Q as TQ

import time
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
rho_0 = [0,0]
w = w_0


#Fz convergence



gridsize = 5005



Fz_list = []
time_list = []

for i in range(5, gridsize, 100):
    
    start_time = time.time()
    
    Fz =  MIM.Fz_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize = i, target = "reflective")
    
    Fz_list.append(Fz)
    
    
    
    time_list.append(time.time() - start_time)
    print (i)
    
print(Fz_list)

print(time_list)



number_of_grid = range(105, gridsize, 100)

plt.figure(1)

plt.scatter(number_of_grid, np.asarray(Fz_list[1:]) * 10 ** 12, label="Fz", marker = "^")

plt.axhline(y=TQ.Fz_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")[0] * 10 ** 12, color='r', linestyle='-', label = "quad built-in")
plt.legend(loc=4,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Fz(pN)',fontsize=20)
plt.grid()
plt.show()


plt.figure(2)

plt.scatter(number_of_grid, np.asarray(time_list[1:]), label="Execution time",marker = "8")

plt.legend(loc=2,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Runningtime(sec)',fontsize=20)
plt.grid()
plt.show()


#print (MIM.Fz_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

#print("--- %s seconds ---" % (time.time() - start_time))



start_time = time.time()

print (TQ.Fz_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective"))

print("--- %s seconds ---" % (time.time() - start_time))
'''

'''
#Fy convergence


rho_0 = [0, w_0]
w = w_0

start_time = time.time()

gridsize = 800



Fy_list = []
timey_list = []

for i in range(10, gridsize, 20):
    
    Fy =  MIM.Fy_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize = i,target = "reflective")
    
    Fy_list.append(Fy)
    
    
    
    timey_list.append(time.time() - start_time)
    print (i)
    
print(Fy_list)

print(timey_list)



number_of_grid = range(10, gridsize, 20)

plt.figure(3)

plt.scatter(number_of_grid, np.asarray(Fy_list) * 10 ** 12, label="Fy", marker = "^")

plt.axhline(y=TQ.Fy_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")[0] * 10 ** 12, color='r', linestyle='-', label = "quad built-in")
plt.legend(loc=4,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Fy(pN)',fontsize=20)
plt.grid()
plt.show()


plt.figure(4)

plt.scatter(number_of_grid, np.asarray(timey_list), label="Execution time",marker = "8")

plt.legend(loc=2,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Runningtime(sec)',fontsize=20)
plt.grid()
plt.show()
'''

'''
#Fx convergence


rho_0 = [w_0, 0]
w = w_0



gridsize = 200



Fx_list = []
timex_list = []

for i in range(10, gridsize, 20):
    
    start_time = time.time()
    
    Fx =  MIM.Fx_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize = i,target = "reflective")
    
    Fx_list.append(Fx)
    
    
    
    timex_list.append(time.time() - start_time)
    print (i)
    
print(Fx_list)

print(timex_list)



number_of_grid = range(10, gridsize, 20)

plt.figure(5)

plt.scatter(number_of_grid, np.asarray(Fx_list) * 10 ** 12, label="Fx", marker = "^")

plt.axhline(y=TQ.Fx_total(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")[0] * 10 ** 12, color='r', linestyle='-', label = "quad built-in")
plt.legend(loc=4,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Fx(pN)',fontsize=20)
plt.grid()
plt.show()


plt.figure(6)
plt.scatter(number_of_grid, np.asarray(timex_list), label="Execution time",marker = "8")

plt.legend(loc=2,fontsize=20)

plt.xlabel('number of grids',fontsize=20)
plt.ylabel('Runningtime(sec)',fontsize=20)
plt.grid()
plt.show()
'''


'''
#check intensity looks nice

rho_0 = [0, 0]
w = w_0

gridsize = 300

theta = np.linspace(0, np.pi/2, gridsize)
phi = np.linspace(0, np.pi * 2, gridsize)


Theta,Phi = np.meshgrid(theta, phi)

r = np.sqrt( (rho * np.sin(theta) * 1 - rho_0[0]) ** 2 + (rho * np.sin(theta) * 0 - rho_0[1]) ** 2) #represent r(the parameter of laser beam) by theta

I = MIM.LG_01_Intensity(r,w)

x = rho * np.sin(theta) 



plt.scatter(x * 10 ** 6, I, label="LG_01 intensity",marker = "8")
plt.title('intensity plot with 300 numbers of grids',fontsize=20)

plt.legend(loc=1,fontsize=20)

plt.xlabel('X(um)',fontsize=20)
plt.ylabel('LG_01 intensity',fontsize=20)
plt.grid()
plt.show()
'''
'''
#power of LG_01 beam

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
'''

'''
#3 plot of Q_z vs displacement w for various offset rho_0x/w0

rho = 30 * 10 ** (-6)
rho_0 = [0,0]   #no offset


rho_0[0] = [0, 0.5*w_0, w_0]

gridsize = 400
w = np.linspace(w_0, 2*rho, 80)

#d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))



Axial_flist_vs_d0 =  np.asarray(MIM.Fz_total_vs_w_plot(rho_0[0][0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

Axial_flist_vs_d1 =  np.asarray(MIM.Fz_total_vs_w_plot(rho_0[0][1],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

Axial_flist_vs_d2 =  np.asarray(MIM.Fz_total_vs_w_plot(rho_0[0][2],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

P_d0 =  MIM.Power_Manual_integration(rho_0[0][0],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P_d0 )


P_d1 =  MIM.Power_Manual_integration(rho_0[0][1],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P_d1 )


P_d2 =  MIM.Power_Manual_integration(rho_0[0][2],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P_d2 )



plt.figure(3)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="rho_0x = 0")
plt.plot( w * 10 ** (6), Q_z1 , lw=2, c="r", label="rho_0x = 0.5w0")
plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho_0x = w0")



new_ticks1 = np.linspace(0, 70, 8) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=30)
plt.yticks(np.linspace(-2, 0, 5),fontsize=30)
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=15)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()
'''
'''

#4 plot of Q_z vs displacement w for various offset rho/w0


rho_0 = [0,0]   #no offset



w = np.linspace(w_0, 70 * 10** (-6), 100)

rho = [0.75*w_0, w_0, 5*w_0, 11.2*w_0]



Axial_flist_vs_d0 =  np.asarray(MIM.Fz_total_vs_w_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

Axial_flist_vs_d1 =  np.asarray(MIM.Fz_total_vs_w_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

Axial_flist_vs_d2 =  np.asarray(MIM.Fz_total_vs_w_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

Axial_flist_vs_d3 =  np.asarray(MIM.Fz_total_vs_w_plot(rho_0[0],rho_0[1], rho[3], n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

P_d0 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P_d0 )


P_d1 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P_d1 )


P_d2 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P_d2 )


P_d3 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho[3], n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z3 = Axial_flist_vs_d2 * c / ( n_0 * P_d3 )




plt.figure(4)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="rho = 0.75*w_0")
plt.plot( w * 10 ** (6), Q_z1 , lw=2, c="r", label="rho = w_0")
plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho = 5*w_0")
plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho = 11.2*w_0")


new_ticks1 = np.linspace(0, 70, 8) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()
'''


'''

#5 plot of Q_z vs rho_0x/w0
rho_0 = [0,0]   #no offset

rho_0[0] = np.linspace(0, w_0, 50)

rho = 30 * 10 ** (-6)

w = rho

gridsize = 500


Axial_flist_vs_d =  np.asarray(MIM.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))


P =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z = Axial_flist_vs_d * c / ( n_0 * P )

plt.figure(5)


plt.plot(rho_0[0]/w_0, Q_z , lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 5, 6) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x/w0',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()



#6 plot of Q_z vs displacement Rs/w


rho_0 = [0,0]   #no offset

w = 30 * 10 ** (-6)

gridsize = 500

rho = np.linspace(10** (-6), w, 100)


Axial_flist_vs_Rs =  np.asarray(MIM.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))


P =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_z = Axial_flist_vs_d * c / ( n_0 * P )

plt.figure(6)


plt.plot(rho/w_0, Q_z , lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 30, 4) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-2, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho/w_0',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.title('Qz vs radius/beam waist w0 at w = w0',fontsize=15)
plt.grid()
plt.show()



#9 plot Qx vs w
w = np.linspace(w_0, 70 * 10** (-6), 100)

rho_0 = [w_0 , 0]

rho = 30 * 10 ** (-6)

gridsize = 500

radial_flistx = np.asarray(TQ.Fx_total_vs_w_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

P =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x = radial_flistx * c / ( n_0 * P )





plt.figure(9)
plt.plot( w * 10 ** 6, Q_x, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 70, 8) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1, 0.5, 4),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('Qx vs w at rho_0x = 10*w_0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()




#10 plot Qx vs rho_0/w0 for various target radius rho
w = 30 * 10 ** (-6)

rho_0 = [ 0 , 0]

rho_0[0] = np.linspace(-w_0, w_0, 100)

rho = [0.5*w, 0.75 * w, 1*w]

gridsize = 500

radial_flistx0 = np.asarray(MIM.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

radial_flistx1 = np.asarray(MIM.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

radial_flistx2 = np.asarray(MIM.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))

P0 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x0 = radial_flistx0 * c / ( n_0 * P0 )


P1 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x1 = radial_flistx0 * c / ( n_0 * P1 )


P2 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x2 = radial_flistx0 * c / ( n_0 * P2 )




plt.figure(10)
plt.plot(rho_0[0]/w_0, Q_x0, lw=2, c="c", label="rho = 0.5w")
plt.plot(rho_0[0]/w_0, Q_x1, lw=2, c="r", label="rho = 0.75w")
plt.plot(rho_0[0]/w_0, Q_x2, lw=2, c="g", label="rho = w")



new_ticks1 = np.linspace(-25, 25, 11) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.5, 0.5, 3),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-25))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x/w_0',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


'''

#11 plot Qx vs rho_0x/w0 for various w

rho = 30 * 10 ** (-6)

w = [rho, 1.5* rho, 2 * rho]

rho_0 = [ 0 , 0]

rho_0[0] = np.linspace(-rho, rho, 100)

gridsize = 300

#rho = np.linspace(10** (-6), 150 * w_0, 100)

radial_flistx0 = np.asarray(MIM.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, gridsize, target = "reflective"))


radial_flistx1 = np.asarray(MIM.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, gridsize, target = "reflective"))



radial_flistx2 = np.asarray(MIM.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, gridsize, target = "reflective"))


P0 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x0 = radial_flistx0 * c / ( n_0 * P0 )


P1 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x1 = radial_flistx1 * c / ( n_0 * P1 )


P2 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x2 = radial_flistx2 * c / ( n_0 * P2 )


plt.figure(11)
plt.plot(rho_0[0]/w_0, Q_x0, lw=2, c="c", label="w = rho")
plt.plot(rho_0[0]/w_0, Q_x1, lw=2, c="r", label="w = 1.5rho")
plt.plot(rho_0[0]/w_0, Q_x2, lw=2, c="g", label="w = 2rho")

new_ticks1 = np.linspace(-12, 12, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.5, 0.5, 3),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-12))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x/w_0',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

'''

#16 compare
radial_flistx0 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective"))

Q_x0 = radial_flistx0 * c / ( n_0 * P )


radial_flistx1 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective"))

Q_x1 = radial_flistx1* c / ( n_0 * P )

radial_flistx2 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective"))

Q_x2 = radial_flistx2 * c / ( n_0 * P )


plt.figure(16)
plt.plot(rho_0[0]/w_0, Q_x0, lw=2, c="c", label="w = rho")
plt.plot(rho_0[0]/w_0, Q_x1, lw=2, c="r", label="w = 1.5rho")
plt.plot(rho_0[0]/w_0, Q_x2, lw=2, c="g", label="w = 2rho")

new_ticks1 = np.linspace(-12, 12, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.12, 0.12, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-12))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x/w_0',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()
'''
'''

#12 plot Qx vs Rs/w0
w =  30 * 10 ** (-6)

rho_0 = [w_0 , 0 ]

rho = np.linspace(10** (-6), w, 100)

gridsize = 500

radial_flistx = np.asarray(MIM.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, target = "reflective"))


P0 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_x = radial_flistx  * c / ( n_0 * P0 )





plt.figure(12)
plt.plot( rho/w_0, Q_x, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 150, 4) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1, 0.5, 4),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho/w_0',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('Qx vs radius/beam waist w0 at w = 10w0 and rho_0x =10 * w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()
'''



#14 plot Qy vs rho_0y/w0 for various w
'''
rho = 30 * 10 ** (-6)

w = [w_0, 2.5* w_0, 5 * w_0]

rho_0 = [ 0 , 0]

rho_0[1] = np.linspace(-25 * w_0, 25 *w_0,100)
'''
'''
rho = 30 * 10 ** (-6)

w = [rho, 1.5 * rho, 2 * rho]

rho_0 = [0, 0]

rho_0[1] = np.linspace(-rho, rho, 100)

gridsize = 300

#rho = np.linspace(10** (-6), 150 * w_0, 100)


radial_flisty0 = np.asarray(MIM.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w[0], z_R, P, gridsize, target = "reflective"))

radial_flisty1 = np.asarray(MIM.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w[1], z_R, P, gridsize, target = "reflective"))

radial_flisty2 = np.asarray(MIM.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, gridsize, target = "reflective"))


P0 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_y0 = radial_flisty0 * c / ( n_0 * P0 )


P1 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_y1 = radial_flisty1 * c / ( n_0 * P1 )


P2 =  MIM.Power_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w_0, z_R, P, gridsize, rrange = 50 * 10 ** (-6), target = "reflective")
Q_y2 = radial_flisty2 * c / ( n_0 * P2 )


plt.figure(14)
plt.plot(rho_0[1]/w_0, Q_y0, lw=2, c="c", label="w = rho")
plt.plot(rho_0[1]/w_0, Q_y1, lw=2, c="r", label="w = 1.5rho")
plt.plot(rho_0[1]/w_0, Q_y2, lw=2, c="g", label="w = 2rho")

new_ticks1 = np.linspace(-15, 15, 11) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.05, 0.05, 3),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-15))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0y/w_0',fontsize=20)
plt.ylabel('Qy',fontsize=20)

#plt.title('Qy vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()
'''
'''
#15 compare with builtin result 
radial_flistx0 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective"))

Q_x0 = radial_flistx0 * c / ( n_0 * P )


radial_flistx1 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective"))

Q_x1 = radial_flistx1* c / ( n_0 * P )

radial_flistx2 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective"))

Q_x2 = radial_flistx2 * c / ( n_0 * P )


plt.figure(15)
plt.plot(rho_0[1]/w_0, Q_x0, lw=2, c="c", label="w = rho")
plt.plot(rho_0[1]/w_0, Q_x1, lw=2, c="r", label="w = 1.5w0")
plt.plot(rho_0[1]/w_0, Q_x2, lw=2, c="g", label="w = 2w0")

new_ticks1 = np.linspace(-15, 15, 11) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.12, 0.12, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-15))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0y/w_0',fontsize=20)
plt.ylabel('Qy',fontsize=20)

#plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()'''
