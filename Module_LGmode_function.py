# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 10:40:12 2020

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



def I_LG01(rou,z, w_0, z_R):
    '''LG01 intensity profile'''
    w = w_0 * np.sqrt(1 + (z/z_R)**2)
    return 2/(np.pi * w**2) * 2 * rou**2/w**2* np.exp(-2*rou**2/w**2)


#Axial integrands


def integrand(theta,phi,d, a, Rs, n_0, w_0, z_R):
    
    
    
    '''integrand of lower surface reflection'''
    
    c = 3 * 10**8
    
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta
    
    
    z = d
    #+ Rs * (1 - np.cos(theta))
    

    Reflection_coe = 1
    
    
    return (1/(2*c)) * n_0 * (1 + np.cos(2*theta)) * \
        I_LG01(rou,z, w_0, z_R) * Reflection_coe * Rs**2 * np.sin(2*theta)

    
    
#Calculation of axial force integral
    
def Axial_force_total(d,a,m, Rs, n_0, w_0, z_R):
    '''total axial forces calculated via integration'''
    
    g = 9.8
  
    forcenet = []
    for d_e in d:
        F_1tz = dblquad(integrand, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d_e, a, Rs, n_0, w_0, z_R))
        
    
        F_znet = F_1tz[0] - m * g #+ F_1tz[0] + F_2rz[0] + F_2tz[0]
        forcenet.append(F_znet * 10 ** 9)

    return forcenet



def disp_velo(d,a,m, Rs, n_0, w_0, z_R, t, y0, r0):
    
    g = 9.8
    
    d_list = []
    vz_list = []
    ya = [1,1]
    
    
    a_list = []
    vr_list = []
    ra = [1,1]
    
    for i in t: 
        if i == 0:
            F_1tz = dblquad(integrand, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (y0[0],r0[0], Rs, n_0, w_0, z_R)) #returns a tuple with first element the integral result and second element = upper bound error
            F_znet = F_1tz[0] - m * g
            b = F_znet
            disp_a = 0.5 * b/m * i**2 + y0[1] * i
            vz = y0[1] + b/m * i
            
            F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (y0[0],r0[0],Rs, n_0, w_0, z_R)) #returns a tuple with first element the integral result and second element = upper bound error #returns a tuple with first element the integral result and second element = upper bound error
            b = F_1rr[0]
            disp_r = 0.5 * b/m * i**2 + r0[1] * i
            vr = r0[1] + b/m * i
        else:
            F_1rz = dblquad(integrand, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (ya[0],r0[0], Rs, n_0, w_0, z_R)) #returns a tuple with first element the integral result and second element = upper bound error
            F_znet = F_1rz[0] - m * g
            b = F_znet
            disp_a = 0.5 * b/m * i**2 + ya[1] * i
            vz = ya[1] + b/m * i
            
            F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,ra[0],Rs, n_0, w_0, z_R)) #returns a tuple with first element the integral result and second element = upper bound error #returns a tuple with first element the integral result and second element = upper bound error
            b = F_1rr[0]
            disp_r = 0.5 * b/m * i**2 + ra[1] * i
            vr = ra[1] + b/m * i
    
        ya[0] = disp_a
        ya[1] = vz
        d_list.append(disp_a)
        vz_list.append(vz)
        
        ra[0] = disp_r
        ra[1] = vr
        a_list.append(disp_r)
        vr_list.append(vr)
    
    return d_list, vz_list, a_list, vr_list

def Radial_disp_velo(d,a,m, Rs, n_0, w_0, z_R, t, r0):
     
    
    a_list = []
    vr_list = []
    ra = [1,1]
    for i in t: 
        if i == 0:
            F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,r0[0],Rs, n_0, w_0, z_R)) #returns a tuple with first element the integral result and second element = upper bound error #returns a tuple with first element the integral result and second element = upper bound error
            b = F_1rr[0]
            disp = 0.5 * b/m * i**2 + r0[1] * i
            vr = r0[1] + b/m * i
        else:
            F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,ra[0],Rs, n_0, w_0, z_R)) #returns a tuple with first element the integral result and second element = upper bound error #returns a tuple with first element the integral result and second element = upper bound error
            b = F_1rr[0]
            disp = 0.5 * b/m * i**2 + ra[1] * i
            vr = ra[1] + b/m * i
    
        ra[0] = disp
        ra[1] = vr
        a_list.append(disp)
        vr_list.append(vr)
    return a_list, vr_list


#Radial force integrads

def integrand_1rr(theta,phi,d, a, Rs, n_0, w_0, z_R):
    
    '''integrand of lower surface reflection, note!!! the rfcoe not considered so it only works for reflection sphere now!!!'''
    c = 3 * 10**8
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta
    
    z = d #+ Rs * (1 - np.cos(theta))
    
    
    Reflection_coe = 1
    
    
    return - ((n_0/(2*c)) * np.sin(2*theta) * I_LG01(rou,z, w_0, z_R) * Reflection_coe * Rs**2 * np.cos(phi) * np.sin(2*theta))


        
#Calculation of radial force integral
def Radial_force_total(d,a,m, Rs, n_0, w_0, z_R):
    '''total radial forces calculated via integration'''
    
    
    forcenet_r = []
    for a_e in a:
        F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e,Rs, n_0, w_0, z_R)) #returns a tuple with first element the integral result and second element = upper bound error
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet * 10 ** 12)

    return forcenet_r