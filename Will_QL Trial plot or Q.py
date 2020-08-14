# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 20:23:03 2020

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

############################
# integration method control
############################

integration_method = 'manual'   # 'manual' or 'integrated'
grid_size = 200

plt.close('all')

###########################

t_start = time.time()

'''
#laser beam parameters for Roosen table 1.B tranparent sphere

w_0 = 2.5* 10 ** (-6)


Lambda = 0.5145 * 10**(-6)
P = 1
z_R = np.pi* w_0 ** 2 / Lambda
rho = 15 * 10 ** (-6)
 
weight = 3.4 * 10**(-10)

c = 3 * 10**8

#trapping item parameters

#Rs = 8.5/2 * w_0

#sig_s = 1100 #density of sphere in kg/m^3

#sig_0 = 1000 #density of medium in kg/m^3

n_0 = 1

n_s = 1.5 #sphere density in Gauthier
'''


'''
#TEM01* reflective target Table 6 matching

c = 3 * 10**8

w_0 = 2.5* 10 ** (-6)
w = 20 * 10 ** (-6)

Lambda = 0.5145 * 10**(-6)
P = 1
z_R = np.pi* w_0 ** 2 / Lambda
rho = 15 * 10 ** (-6)
 
weight = 3.4 * 10**(-10)

n_0 = 1

n_s_n = 0.05

k =  2.87 

n_s = n_s_n - k*1j

'''


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


##################Qy vs rho_0y#############################


######################################################
#15 plot Qy vs rho_0y for various target radius rho
######################################################

w = w_0

rho_0 = [0 , 0]

rho_0[1] = np.linspace(-25 * w_0, 25*w_0, 100)

rho = [0.75 * w_0, w_0, 1.5*w_0, 11.2*w_0, 15*w_0]

#rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

radial_flisty0 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0], rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = 'integrated', grid_size = grid_size))

Q_y0 = radial_flisty0 * c / ( n_0 * P )


radial_flisty1 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0], rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = 'integrated', grid_size = grid_size))

Q_y1 = radial_flisty1* c / ( n_0 * P )

radial_flisty2 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0], rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = 'integrated', grid_size = grid_size))

Q_y2 = radial_flisty2 * c / ( n_0 * P )


radial_flisty3 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0], rho_0[1], rho[3], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = 'integrated', grid_size = grid_size))

Q_y3 = radial_flisty3 * c / ( n_0 * P )

radial_flisty4 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0], rho_0[1], rho[4], n_0, n_s, w_0, w, z_R, P, target = 'reflective', integration_method = 'integrated', grid_size = grid_size))

Q_y4 = radial_flisty4 * c / ( n_0 * P )


plt.figure(15)
plt.plot(rho_0[1] / w_0, Q_y0, lw=2, c="c", label="rho = 0.75w0")
plt.plot(rho_0[1] / w_0, Q_y1, lw=2, c="r", label="rho = w0")
plt.plot(rho_0[1] / w_0, Q_y2, lw=2, c="g", label="rho = 1.5w0")
plt.plot(rho_0[1] / w_0, Q_y3, lw=2, c="b", label="rho = 11.2w0")
plt.plot(rho_0[1] / w_0, Q_y4, lw=2, c="k", label="rho = 15w0")




new_ticks1 = np.linspace(-25, 25, 11) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1, 1, 5),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-25))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0y/w_0',fontsize=20)
plt.ylabel('Qy',fontsize=20)

#plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('2.67', 'related to w', 'x-axis -25w0 to 25w0', rho_0[0]* 10 ** 6, '0.75w0, w0, 1.5w0, 11.2w0, 15w0')


'''




############# test beam overlap and sphere plot ####################

w = 2*rho

rho_0x = w/10
rho_0y = w/5


integration_result = TQ.F_total_manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'x', grid_size = grid_size)

theta = integration_result['theta']
phi = integration_result['phi']

TQ.plot_intensity_sphere_overlap(theta, phi, rho_0x, rho_0y, rho, w)
####################




##################################################

'''
'''
#Fy convergence


rho_0 = [0, w_0]
w = w_0

start_time = time.time()

gridsize = 200



Fy_list = []
timey_list = []

for i in range(10, gridsize, 20):
    
    start_time = time.time()
    
    integration_results = TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'y', grid_size = grid_size)
            
    Fy = integration_results['force_total']
    
    Fy_list.append(Fy)
    
    
    
    timey_list.append(time.time() - start_time)
    print (i)
    
print(Fy_list)

print(timey_list)



number_of_grid = range(10, gridsize, 20)

plt.figure(50)

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

for i in range(1, gridsize, 2):
    
    start_time = time.time()
    
    integration_results = TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'x', grid_size = grid_size)
            
    Fx = integration_results['force_total']
    
    #Fx =  TQ.F_total_manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", coordinate = 'x', grid_size = i)
    
    Fx_list.append(Fx)
    
    
    
    timex_list.append(time.time() - start_time)
    print (i)
    
print(Fx_list)

print(timex_list)



number_of_grid = range(1, gridsize, 2)

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
#1 plot of Q_z vs displacement d


rho_0 = [0,0]   #no offset


w = np.linspace(w_0, 70 * 10** (-6), 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))


Axial_flist_vs_d =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective"))



Q_z = Axial_flist_vs_d * c / ( n_0 * P )

plt.figure(1)


plt.plot( d * 10 ** (6), Q_z , lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 550, 12) # plot axis
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

plt.xlabel('d(um)',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()




#2 plot of Q_z vs displacement d/w0
rho_0 = [0,0]   #no offset


w = np.linspace(w_0, 70 * 10** (-6), 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))


Axial_flist_vs_d =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z = Axial_flist_vs_d * c / ( n_0 * P )

plt.figure(2)


plt.plot( d/w_0, Q_z , lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 250, 6) # plot axis
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

plt.xlabel('d/w0',fontsize=20)
plt.ylabel('Qz',fontsize=20)
plt.grid()
plt.show()
'''



'''


#3 plot of Q_z vs displacement w for various offset rho_0x/w0


rho_0 = [0,0]   #no offset


rho_0[0] = [0, 0.5*w_0, w_0, 2*w_0]
w = np.linspace(w_0, 70 * 10** (-6), 100)

#d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))



Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )


Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][1],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][2],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )
Axial_flist_vs_d3 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0][3],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z3 = Axial_flist_vs_d3 * c / ( n_0 * P )


plt.figure(3)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="rho_0x = 0")
plt.plot( w * 10 ** (6), Q_z1 , lw=2, c="r", label="rho_0x = 0.5w0")
plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho_0x = w0")
plt.plot( w * 10 ** (6), Q_z3 , lw=2, c="y", label="rho_0x = 2w0")


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
#4 plot of Q_z vs displacement w for various offset rho/w0


rho_0 = [0,0]   #no offset



w = np.linspace(w_0, 70 * 10** (-6), 100)

rho = [0.75*w_0, w_0, 1.5*w_0, 11.2*w_0]



Axial_flist_vs_d0 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



Q_z0 = Axial_flist_vs_d0 * c / ( n_0 * P )


Axial_flist_vs_d1 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z1 = Axial_flist_vs_d1 * c / ( n_0 * P )

Axial_flist_vs_d2 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z2 = Axial_flist_vs_d2 * c / ( n_0 * P )
Axial_flist_vs_d3 =  np.asarray(TQ.Fz_total_vs_d_plot(rho_0[0],rho_0[1], rho[3], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))
Q_z3 = Axial_flist_vs_d3 * c / ( n_0 * P )


plt.figure(4)


plt.plot( w * 10 ** (6), Q_z0 , lw=2, c="c", label="rho = 0.75*w_0")
plt.plot( w * 10 ** (6), Q_z1 , lw=2, c="r", label="rho = w_0")
plt.plot( w * 10 ** (6), Q_z2 , lw=2, c="g", label="rho = 1.5*w_0")
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



#5 plot of Q_z vs rho_0x/w0
rho_0 = [0,0]   #no offset

rho_0[0] = np.linspace(0, 5*w_0, 100)

w = w_0


Axial_flist_vs_d =  np.asarray(TQ.Fz_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))



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


'''
#6 plot of Q_z vs displacement Rs/w


rho_0 = [0,0]   #no offset

w = w_0

rho = np.linspace(10** (-6), 30 *w_0, 100)


Axial_flist_vs_Rs =  np.asarray(TQ.Fz_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z = Axial_flist_vs_Rs * c / ( n_0 * P )

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
'''



'''

#7 plot Qx vs d

w = np.linspace(w_0, 70 * 10** (-6), 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))


rho_0 = [10*w_0 , 0 ]



radial_flistx = np.asarray(TQ.Fx_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x = radial_flistx  * c / ( n_0 * P )





plt.figure(7)
plt.plot( d * 10 ** 6, Q_x, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 600, 7) # plot axis
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

plt.xlabel('d(um)',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('Qx vs d at rho_0x = 10*w_0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()




#8 plot Qx vs d/w0
w = np.linspace(w_0, 70 * 10** (-6), 100)

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))


rho_0 = [10*w_0 , 0 ]



radial_flistx = np.asarray(TQ.Fx_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x = radial_flistx  * c / ( n_0 * P )





plt.figure(8)
plt.plot( d/w_0, Q_x, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 250, 6) # plot axis
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

plt.xlabel('d/w_0',fontsize=20)
plt.ylabel('Qx',fontsize=20)

plt.title('Qx vs d/w0 at rho_0x = 10*w_0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()
'''

'''
#9 plot Qx vs w

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*rho, 100)


rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [w_0 , 0 ]



radial_flistx = np.asarray(TQ.Fx_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x = radial_flistx  * c / ( n_0 * P )





plt.figure(9)
plt.plot( w * 10 ** 6, Q_x, lw=2, c="c", label="My result")



new_ticks1 = np.linspace(0, 60, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.2, 0.1, 4),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
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

#plt.title('Qx vs w at rho_0x = 10*w_0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('x-axis', 'x-axis', rho_0[1]* 10 ** 6, rho_0[0]* 10 ** 6, 30)

'''


'''
#10 plot Qx vs rho_0x/w0 for various target radius rho

a = 30 * 10 ** (-6)
w = a * np.sqrt(2)

rho_0 = [0 , 0]

rho_0[0] = np.linspace(-a, a, 100)

rho = [0.5*a, 0.75*a, a]

rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

radial_flistx0 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0 * c / ( n_0 * P )


radial_flistx1 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x1 = radial_flistx1* c / ( n_0 * P )

radial_flistx2 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x2 = radial_flistx2 * c / ( n_0 * P )





plt.figure(10)
plt.plot(rho_0[0] * np.sqrt(2)/w, Q_x0, lw=2, c="c", label="rho = 15um")
plt.plot(rho_0[0] * np.sqrt(2)/w, Q_x1, lw=2, c="r", label="rho = 22.5um")
plt.plot(rho_0[0] * np.sqrt(2)/w, Q_x2, lw=2, c="g", label="rho = 30um")



new_ticks1 = np.linspace(-1, 1, 3) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.08, 0.08, 5),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('Qx',fontsize=20)

#plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30um', '334', rho_0[1]* 10 ** 6, 'x-axis -30um to 30um', '15um, 22.5um, 30um')

'''


'''

#11 plot Qx vs rho_0x/w0 for various w
rho = 30 * 10 ** (-6)

w = [np.sqrt(2) * rho, 2*rho, 2.5*rho]

rho_0 = [ 0 , 0]

rho_0[0] = np.linspace(-rho, rho, 100)


radial_flistx0 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x0 = radial_flistx0 * c / ( n_0 * P )


radial_flistx1 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x1 = radial_flistx1* c / ( n_0 * P )

radial_flistx2 = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x2 = radial_flistx2 * c / ( n_0 * P )


plt.figure(11)
plt.plot(rho_0[0] * 10 ** 6, Q_x0, lw=2, c="c", label="w = sqrt(2)*30um")
plt.plot(rho_0[0] * 10 ** 6, Q_x1, lw=2, c="r", label="w = 60um")
plt.plot(rho_0[0] * 10 ** 6, Q_x2, lw=2, c="g", label="w = 75um")

new_ticks1 = np.linspace(-30, 30, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.08, 0.08, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-30))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x(um)',fontsize=20)
plt.ylabel('Qx',fontsize=20)

#plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('sqrt(2)*30um, 60um, 75um', 'related to w', rho_0[1]* 10 ** 6, 'x-axis -30um to 30um', '30um')


'''



'''

#12 plot Qx vs Rs/w0

a = 30 * 10 ** (-6)

w =  a * np.sqrt(2)

rho_0 = [a, 0 ]

rho = np.linspace(15 * 10 ** (-6), a, 100)

radial_flistx = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x = radial_flistx  * c / ( n_0 * P )





plt.figure(12)
plt.plot( rho * np.sqrt(2)/w, Q_x, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0.5, 1, 6) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(0, -0.008, 9),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0.5))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('sqrt(2)rho/w',fontsize=20)
plt.ylabel('Qx',fontsize=20)

#plt.title('Qx vs radius/beam waist w0 at w = 10w0 and rho_0x =10 * w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30um', 'related to w', rho_0[1]* 10 ** 6, '30um', 'x-axis 15um to 30um')

'''



'''

#13 plot gradient Fx vs rho_0x

a = 30 * 10 ** (-6)
w = a*np.sqrt(2)

rho_0 = [0 , 0]

rho_0[0] = np.linspace(-a, a, 100)

rho = [0.5*a, 0.75*a, a]


#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_x0 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective")[0])


grad_x1 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective")[0])

grad_x2 = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective")[0])


rho_0x = np.asarray(TQ.Fx_total_gradient(rho_0[0], rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective")[1])

#rho_0xratio = list(map(lambda x: x/w_0, rho_0x))


plt.figure(13)
plt.plot( rho_0x, grad_x0*10**8, lw=2, c="c", label="rho = 15um")

plt.plot( rho_0x, grad_x1*10**8, lw=2, c="r", label="rho = 22.5um")

plt.plot( rho_0x, grad_x2*10**8, lw=2, c="g", label="rho = 30um")


new_ticks1 = np.linspace(-1, 1, 3) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-3, 3, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=15)

plt.xlabel('rho_0x/(w/sqrt(2))',fontsize=20)
plt.ylabel('grad_Fx(stiffness)10^(-8)',fontsize=20)

#plt.title('grad_Fx vs x offset waist w0 at w = w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('sqrt(2)*30um', 'related to w', rho_0[1]* 10 ** 6, 'x-axis -30um to 30um', '15um, 22.5um, 30um')
'''



'''
#14 plot Qy vs w

rho = 30 * 10 ** (-6)
w = np.linspace(w_0, 2*rho, 100)


#rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

rho_0 = [0, w_0 ]



radial_flisty = np.asarray(TQ.Fy_total_vs_d_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y = radial_flisty  * c / ( n_0 * P )





plt.figure(14)
plt.plot( w * 10 ** 6, Q_y, lw=2, c="c", label="My result")



new_ticks1 = np.linspace(0, 60, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.2, 0.1, 4),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('w(um)',fontsize=20)
plt.ylabel('Qy',fontsize=20)

#plt.title('Qx vs w at rho_0x = 10*w_0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('x-axis', 'x-axis', rho_0[1]* 10 ** 6, rho_0[0]* 10 ** 6, 30)

'''



'''
#15 plot Qy vs rho_0y for various target radius rho

a = 30 * 10 ** (-6)
w = a * np.sqrt(2)

rho_0 = [0 , 0]

rho_0[1] = np.linspace(-a, a, 100)

rho = [0.5*a, 0.75*a, a]

#rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

radial_flisty0 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y0 = radial_flisty0 * c / ( n_0 * P )


radial_flisty1 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y1 = radial_flisty1* c / ( n_0 * P )

radial_flisty2 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y2 = radial_flisty2 * c / ( n_0 * P )





plt.figure(15)
plt.plot(rho_0[1] * np.sqrt(2)/w, Q_y0, lw=2, c="c", label="rho = 15um")
plt.plot(rho_0[1] * np.sqrt(2)/w, Q_y1, lw=2, c="r", label="rho = 22.5um")
plt.plot(rho_0[1] * np.sqrt(2)/w, Q_y2, lw=2, c="g", label="rho = 30um")



new_ticks1 = np.linspace(-1, 1, 3) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.08, 0.08, 5),fontsize=20)

plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0y/(w/sqrt(2))',fontsize=20)
plt.ylabel('Qy',fontsize=20)

#plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30um', '334', 'x-axis -30um to 30um', rho_0[0]* 10 ** 6, '15um, 22.5um, 30um')

'''


'''

#16 plot Qy vs rho_0y for various w
rho = 30 * 10 ** (-6)

w = [np.sqrt(2) * rho, 2*rho, 2.5*rho]

rho_0 = [ 0 , 0]

rho_0[1] = np.linspace(-rho, rho, 100)



radial_flisty0 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y0 = radial_flisty0 * c / ( n_0 * P )


radial_flisty1 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y1 = radial_flisty1* c / ( n_0 * P )

radial_flisty2 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y2 = radial_flisty2 * c / ( n_0 * P )


plt.figure(16)
plt.plot(rho_0[1] * 10 ** 6, Q_y0, lw=2, c="c", label="w = sqrt(2)*30um")
plt.plot(rho_0[1] * 10 ** 6, Q_y1, lw=2, c="r", label="w = 60um")
plt.plot(rho_0[1] * 10 ** 6, Q_y2, lw=2, c="g", label="w = 75um")

new_ticks1 = np.linspace(-30, 30, 7) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.08, 0.08, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-30))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0y(um)',fontsize=20)
plt.ylabel('Qy',fontsize=20)

#plt.title('Qx vs rho_0x/beam waist w0 rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('sqrt(2)*30um, 60um, 75um', 'related to w', 'x-axis -30um to 30um', rho_0[0]* 10 ** 6, '30um')
'''





'''

#17 plot Qx vs Rs/w0

a = 30 * 10 ** (-6)

w =  a * np.sqrt(2)

rho_0 = [0, a]

rho = np.linspace(15 * 10 ** (-6), a, 100)

radial_flisty = np.asarray(TQ.Fy_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y = radial_flisty  * c / ( n_0 * P )





plt.figure(17)
plt.plot( rho * np.sqrt(2)/w, Q_y, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0.5, 1, 6) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(0, -0.008, 9),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0.5))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('sqrt(2)rho/w',fontsize=20)
plt.ylabel('Qy',fontsize=20)

#plt.title('Qx vs radius/beam waist w0 at w = 10w0 and rho_0x =10 * w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()

MTP.table_parameter('sqrt(2)*30um', 'related to w', '30um', rho_0[0]* 10 ** 6, 'x-axis 15um to 30um')

'''



'''

#18 plot gradient Fy vs rho_0y

a = 30 * 10 ** (-6)
w = a*np.sqrt(2)

rho_0 = [0, 0]

rho_0[1] = np.linspace(-a, a, 100)

rho = [0.5*a, 0.75*a, a]


#F = TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective")

grad_y0 = np.asarray(TQ.Fy_total_gradient(rho_0[0], rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective")[0])

grad_y1 = np.asarray(TQ.Fy_total_gradient(rho_0[0], rho_0[1], rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective")[0])

grad_y2 = np.asarray(TQ.Fy_total_gradient(rho_0[0], rho_0[1], rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective")[0])




rho_0y = np.asarray(TQ.Fy_total_gradient(rho_0[0], rho_0[1], rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective")[1])

#rho_0xratio = list(map(lambda x: x/w_0, rho_0x))


plt.figure(18)
plt.plot( rho_0y, grad_y0*10**8, lw=2, c="c", label="rho = 15um")

plt.plot( rho_0y, grad_y1*10**8, lw=2, c="r", label="rho = 22.5um")

plt.plot( rho_0y, grad_y2*10**8, lw=2, c="g", label="rho = 30um")


new_ticks1 = np.linspace(-1, 1, 3) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-3, 3, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-1))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0y/(w/sqrt(2))',fontsize=20)
plt.ylabel('grad_Fy(stiffness)10^(-8)',fontsize=20)

#plt.title('grad_Fx vs x offset waist w0 at w = w0, rho_0y = 0',fontsize=15)
plt.grid()
plt.show()


MTP.table_parameter('sqrt(2)*30um', 'related to w', 'x-axis -30um to 30um',rho_0[0]* 10 ** 6,  '30um')

'''

'''
#Q vs x
rho = 30 * 10 ** (-6)

rho_0 = [0, 0]
w = w_0

gridsize = 30

theta = np.linspace(0, np.pi/2, gridsize)


dFz = MIM.dFz_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, theta, target = "reflective")[0]

x =  MIM.dFz_Manual_integration(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, gridsize, theta, target = "reflective")[1]


plt.scatter(x * 10 ** 6, dFz, label="LG_01 intensity",marker = "8")
plt.title('intensity plot with 300 numbers of grids',fontsize=20)

plt.legend(loc=1,fontsize=20)

plt.xlabel('X(um)',fontsize=20)
plt.ylabel('LG_01 intensity',fontsize=20)
plt.grid()
plt.show()
'''
'''
#17 plot Qx vs Rs/w0

a = 30 * 10 ** (-6)

w =  a * np.sqrt(2)

rho_0 = [a, a]

rho = np.linspace(15 * 10 ** (-6), a, 100)

radial_flisty = np.asarray(TQ.Fy_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y = radial_flisty  * c / ( n_0 * P )






plt.plot( rho * np.sqrt(2)/w, Q_y, lw=2, c="c", label="Qy")




radial_flistx = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x = radial_flistx  * c / ( n_0 * P )



plt.plot( rho * np.sqrt(2)/w, Q_x, lw=2, c="r", label="Qx")




Axial_flist_vs_Rs =  np.asarray(TQ.Fz_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))


Q_z = Axial_flist_vs_Rs * c / ( n_0 * P )




plt.plot(rho * np.sqrt(2)/w, Q_z , lw=2, c="g", label="Qz")


new_ticks1 = np.linspace(0.5, 1, 6) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.3, 0, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0.5))
ax.spines['bottom'].set_position(('data',0))

plt.legend(loc=4,fontsize=15)

plt.xlabel('rho/(w/sqrt(2))',fontsize=20)
plt.ylabel('Q',fontsize=20)
#plt.title('Qz vs radius/beam waist w0 at w = w0',fontsize=15)
plt.grid()
plt.show()

'''
#print R
'''
theta = np.linspace(0, np.pi/2, 100)
print (theta)

t = []
Ref = []

for i in range(len(theta)):
    theta_c = np.arcsin( ( n_s.real ) / n_0 )

    
    if theta[i] >= theta_c:
        
        
        theta_2 = np.pi/2
        
      
    else:
        
        theta_2 = np.arcsin(n_0*np.sin(theta[i])/n_s.real)
        
    print (theta[i])
    
    t.append(theta_2)
    
    print (np.around( TQ.rfcoe_sm(theta[i],theta_2, n_0, n_s), decimals= 15))
           
    R = np.around( TQ.rfcoe_sm(theta[i],theta_2, n_0, n_s), decimals= 15)
    
    Ref.append(np.around(TQ.rfcoe_sm(theta[i],theta_2, n_0, n_s), decimals= 15))
print (i)
print (t)    
print (Ref)     
'''

'''
#plot reproduction of Fy vs offset rho0y Roosen 1979
w = w_0
rho_0 = [0,0]
rho_0[1] = np.linspace(-25 * w, 25 *w,100)

radial_flisty = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w, z_R, P))

Q_y = radial_flisty * c / ( n_0 * P )


print(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w, z_R, P)[1])


plt.figure(2)
plt.plot( rho_0[1] /w, Q_y, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(-25, 25, 11) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1, 1, 5),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-25))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0y/w0',fontsize=20)
plt.ylabel('Qy',fontsize=20)
plt.grid()
plt.show()

'''

'''
#plot reproduction of Fx vs offset rho0x Roosen 1979
w = w_0
rho_0 = [0,0]
rho_0[0] = np.linspace(-25 * w, 25 *w,100)

radial_flistx = np.asarray(TQ.Fx_total_vs_rho0x_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, integration_method = integration_method, grid_size = grid_size))

Q_x = radial_flistx  * c / ( n_0 * P )





plt.figure(3)
plt.plot( rho_0[0]/w_0, Q_x, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(-25, 25, 11) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-1.5, 1.5, 7),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',-25))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho_0x/w0',fontsize=20)
plt.ylabel('Qx',fontsize=20)
plt.grid()
plt.show()
'''



'''
#plot vector field of optical forces in transverse xy plane

p = 0
l = 1

w = w_0
rho_0x,rho_0y = np.meshgrid(np.linspace(-30 * 10** (-6), 30 * 10** (-6), 10), np.linspace(-30 * 10** (-6), 30 * 10** (-6), 10))

u = 10 ** (10) * np.asarray(TQ.xy_plane_vectorplot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P)[0])
v = 10 ** (10) * np.asarray(TQ.xy_plane_vectorplot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P)[1])

theta = np.linspace(0, 2*np.pi, 100)


xc = rho * 10 ** 6 *np.cos(theta)
zc = rho * 10 ** 6 *np.sin(theta)

plt.quiver(rho_0x* 10** (6), rho_0y* 10** (6), u, v)

plt.plot(xc,zc)


new_ticks1 = np.linspace(-30, 30, 5) # plot axis
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-30, 30, 5),fontsize=20)

plt.xlabel('rho_0x(um)',fontsize=20)
plt.ylabel('rho_0y(um)',fontsize=20)

plt.title('Transverse Optical Force(10^(-10)N)', fontsize=12)
plt.grid()
plt.show()
'''


'''
#plot vector field of optical forces in xz plane
rho_0y = 0


rho_0x, w = np.meshgrid( np.linspace(-30 * 10 ** (-6), 30 * 10 ** (-6), 10 ), np.linspace(w_0, 50 * 10** (-6), 10))

Fx2 = 10 ** (10) * np.asarray( TQ.xz_plane_vectorplot(rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P)[0])
Fz = 10 ** (10) * (np.asarray( TQ.xz_plane_vectorplot(rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P)[1]))

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))


theta = np.linspace(0, 2*np.pi, 100)


xc = rho * 10 ** 6 *np.cos(theta)
zc = rho * 10 ** 6 *np.sin(theta)

plt.quiver(rho_0x* 10** (6), d* 10** (6), Fx2, -Fz)

plt.plot(xc,zc)


#circle1 = plt.Circle((0, 0), rho, color='r')

#plt.add_patch(circle1)

new_ticks1 = np.linspace(-30, 30, 5) # plot axis
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(0, 350, 8),fontsize=20)

plt.xlabel('rho_0x(um)',fontsize=20)
plt.ylabel('d(um)',fontsize=20)

plt.title('Axial Optical Force(10^(-10)N)', fontsize=15)



plt.grid()
plt.show()
'''




'''

#plot vector field of optical forces in xyz plane

rho_0x, rho_0y, w = np.meshgrid( np.linspace(-30 * 10 ** (-6), 30 * 10 ** (-6), 10 ), np.linspace(-30 * 10** (-6), 30 * 10** (-6), 10), np.linspace(w_0, 50 * 10** (-6), 5))

Fx = 5*10 ** (12) * np.asarray( TQ.xyz_plane_vectorplot(rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P)[0])
Fy = 5*10 ** (12) *np.asarray( TQ.xyz_plane_vectorplot(rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P)[1])
Fz = 5*10 ** (12) *np.asarray( TQ.xyz_plane_vectorplot(rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P)[2])

d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 ))

#Q_x = 3 * Fx  * c / ( n_0 * P ) #10*Q in definition
#Q_y = 3 * Fy  * c / ( n_0 * P )
#Q_z = 3 * Fz  * c / ( n_0 * P )

fig = plt.figure()
ax = fig.gca(projection='3d')

ax.quiver(rho_0x* 10** (6), rho_0y* 10** (6), d* 10** (6), Fx, Fy, -Fz)

ax.set_xlim([-30, 30])
ax.set_ylim([-30, 30])
ax.set_zlim([0, 350])

ax.set_xlabel('rho_0x(um)',fontsize=20)
ax.set_ylabel('rho_0y(um)',fontsize=20)
ax.set_zlabel('d(um)',fontsize=20)





plt.title('3D vector plot of optical forces', fontsize=20)



plt.grid()
plt.show()
'''



'''
#plot LG intensity


r = np.linspace(0, 4*10**(-6),50)
#I_LG01 = TQ.LG_01_Intensity(x,y,w_0)

#plt.plot(r, I_LG01)

P = TQ.power(w_0)

theta = np.linspace( 0, 2 * np.pi, 40)
r, theta = np.meshgrid(r, theta)

X = r * np.sin(theta)
Y = r * np.cos(theta)
Z = TQ.LG_01_Intensity(X, Y, w_0)

ax = plt.axes(projection='3d')
ax.plot_surface(X, Y,Z , rstride=1, cstride=1,
                cmap='viridis', edgecolor='none');
'''




'''
#plot Qy vs Rs/w
w = 5 * w_0

rho_0 = [ 0 , 10*w_0 ]

rho = np.linspace(10** (-6), 150 * w_0, 100)

radial_flisty = np.asarray(TQ.Fy_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, integration_method = integration_method, grid_size = grid_size))

Q_y = radial_flisty * c / ( n_0 * P )





plt.figure(2)
plt.plot(rho/w_0, Q_y, lw=2, c="c", label="My result")


new_ticks1 = np.linspace(0, 150, 4) # plot axis
print(new_ticks1)
plt.xticks(new_ticks1,fontsize=20)
plt.yticks(np.linspace(-0.5, 0.5, 3),fontsize=20)
ax = plt.gca()
ax.spines['top'].set_color('none')
ax.spines['right'].set_color('none')
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
ax.spines['left'].set_position(('data',0))
ax.spines['bottom'].set_position(('data',0))

#plt.legend(loc=1,fontsize=13)

plt.xlabel('rho/w_0',fontsize=20)
plt.ylabel('Qy',fontsize=20)

plt.title('Qy vs radius/beam waist w0 at w = 5w0 and rho_0y = 10*w0',fontsize=15)
plt.grid()
plt.show()
'''



'''
def f(x,y):
    
    return -(np.sin(x)*np.cos(y))
    
    #return -((np.sin(x)*np.cos(y)) **2 +(np.sin(x)*np.sin(y)-10)**2)* np.exp((np.sin(x)*np.cos(y)-10) **2 +(np.sin(x)*np.sin(y))**2) * np.sin(x) * np.cos(x) * np.sin(y) * np.sin(2*x)

Int = dblquad(f, 0, 2*np.pi, lambda y: 0, lambda phi: np.pi/2) #returns a tuple with first element the integral result and second element = upper bound error
print (Int)


def g(x,y):
    
    return -((np.sin(x)*np.cos(y)-10) **2 +(np.sin(x)*np.sin(y))**2)* np.exp((np.sin(x)*np.cos(y)-10) **2 +(np.sin(x)*np.sin(y))**2) * np.sin(x) * np.cos(x) * np.cos(y) * np.sin(2*x)

Intx = dblquad(g, 0, 2*np.pi, lambda y: 0, lambda phi: np.pi/2) #returns a tuple with first element the integral result and second element = upper bound error
print (Intx)
'''

'''

#plot Qy vs rho_0y/w0 for various w
w = [w_0, 2.5* w_0, 5 * w_0]

rho_0 = [ 0 , 0]

rho_0[1] = np.linspace(-25 * w_0, 25 *w_0,100)

#rho = np.linspace(10** (-6), 150 * w_0, 100)

radial_flisty0 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w[0], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y0 = radial_flisty0 * c / ( n_0 * P )


radial_flisty1 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w[1], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y1 = radial_flisty1* c / ( n_0 * P )

radial_flisty2 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w[2], z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y2 = radial_flisty2 * c / ( n_0 * P )


plt.figure(4)
plt.plot(rho_0[1]/w_0, Q_y0, lw=2, c="c", label="w = w0")
plt.plot(rho_0[1]/w_0, Q_y1, lw=2, c="r", label="w = 2.5w0")
plt.plot(rho_0[1]/w_0, Q_y2, lw=2, c="g", label="w = 5w0")

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

plt.xlabel('rho_0y/w_0',fontsize=20)
plt.ylabel('Qy',fontsize=20)

plt.title('Qy vs rho_0y/beam waist w0 rho_0x = 0',fontsize=15)
plt.grid()
plt.show()





'''




'''
#plot Qy vs rho_0y/w0 for various target radius rho
w = w_0

rho_0 = [ 0 , 0]

rho_0[1] = np.linspace(-25 * w_0, 25 *w_0,100)

rho = [0.75*w_0, w_0, 1.5*w_0, 11.2*w_0, 15*w_0]

radial_flisty0 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho[0], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y0 = radial_flisty0 * c / ( n_0 * P )


radial_flisty1 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho[1], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y1 = radial_flisty1* c / ( n_0 * P )

radial_flisty2 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho[2], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y2 = radial_flisty2 * c / ( n_0 * P )


radial_flisty3 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho[3], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y3 = radial_flisty3 * c / ( n_0 * P )

radial_flisty4 = np.asarray(TQ.Fy_total_vs_rho0y_plot(rho_0[0],rho_0[1],rho[4], n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y4 = radial_flisty4 * c / ( n_0 * P )

plt.figure(6)
plt.plot(rho_0[1]/w_0, Q_y0, lw=2, c="c", label="rho = 0.75w_0")
plt.plot(rho_0[1]/w_0, Q_y1, lw=2, c="r", label="rho = w_0")
plt.plot(rho_0[1]/w_0, Q_y2, lw=2, c="g", label="rho = 1.5w_0")
plt.plot(rho_0[1]/w_0, Q_y3, lw=2, c="y", label="rho = 11.2w_0")
plt.plot(rho_0[1]/w_0, Q_y4, lw=2, c="m", label="rho = 15w_0")

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

plt.xlabel('rho_0y/w_0',fontsize=20)
plt.ylabel('Qy',fontsize=20)

plt.title('Qy vs rho_0y/beam waist w0 rho_0x = 0',fontsize=15)
plt.grid()
plt.show()




'''



'''
#merging plots Qx, Qy, Qz vs rho/w0 plot
w = w_0

rho_0 = [ 0 , 0]

rho = np.linspace(10** (-6), 150 * w_0, 100)



Fz = np.asarray(TQ.Fz_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_z = Fz * c / ( n_0 * P )


Fy = np.asarray(TQ.Fy_total_vs_Rs_plot(rho_0[0],rho_0[1],rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_y = Fy* c / ( n_0 * P )

Fx = np.asarray(TQ.Fx_total_vs_Rs_plot(rho_0[0],rho_0[1], rho, n_0, n_s, w_0, w, z_R, P, target = "reflective", integration_method = integration_method, grid_size = grid_size))

Q_x = Fx * c / ( n_0 * P )




plt.figure(6)
plt.plot(rho/w_0, Q_z, lw=2, c="c", label="Qz")
plt.plot(rho/w_0, Q_y, lw=2, c="r", label="Qy")
plt.plot(rho/w_0, Q_x, lw=2, c="g", label="Qx")


new_ticks1 = np.linspace(0, 150, 4) # plot axis
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

plt.xlabel('rho/w_0',fontsize=20)
plt.ylabel('Q',fontsize=20)

plt.title('Qz, Qx, Qy merging plots vs rho/w0 at zero offsets',fontsize=15)
plt.grid()
plt.show()
'''


t_finish = time.time()

print(f"script executed in {t_finish-t_start:.2f} seconds")

#Trapping stability stiffness
