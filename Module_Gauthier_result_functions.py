# -*- coding: utf-8 -*-
"""
Created on Thu Jul 16 17:30:14 2020

@author: liuqi
"""


'''Gauther and Wallace 1995 paper reproduction functions
   NOTE: With current design, the function is only useful for reflecting sphere since parameters in rfcoe_sm is not added in the integrand functions'''

import scipy as sp
import numpy as np
import matplotlib.pylab as plt
import scipy.integrate as spi
from scipy.integrate import quad
import seaborn
from scipy.integrate import odeint
from scipy.integrate import dblquad

def r_s(theta, theta2, n_0, n_s):
    '''Fresnel coefficient of s-polarisation'''
    return (n_0*np.cos(theta) - n_s*np.cos(theta2))/(n_0*np.cos(theta) + n_s*np.cos(theta2))

def r_p(theta, theta2, n_0, n_s):
    '''Fresnel coefficient of p-polarisation'''
    return (n_0 * np.cos(theta2) - n_s * np.cos(theta))/ (n_s * np.cos(theta) + n_0 * np.cos(theta2))

def rfcoe_sm(theta,theta2):
    
    '''Average reflectance'''
    return (r_s(theta,theta2)**2 + r_p(theta,theta2)**2)/2


def I_GB(rou, z, w_0, z_R, P):
    '''GB intensity profile'''
    #Gaussian beam propagation parameters
    w = w_0 * np.sqrt(1 + (z/z_R)**2)
    #R = z + z_R**2/z
    return (2 * P/(np.pi*w**2)) * np.exp( - 2 * rou** 2 / w**2)

def I_LG01(rou,z, w_0, z_R):
    '''LG01 intensity profile'''
    w = w_0 * np.sqrt(1 + (z/z_R)**2)
    return 2/(np.pi * w**2) * 2 * rou**2/w**2* np.exp(-2*rou**2/w**2)


#Axial integrands


def integrand(theta,phi,d, a, Rs, n_0, n_s, w_0, z_R, P):
    
    
    
    '''integrand of lower surface reflection'''
    
    c = 3 * 10**8
    #rou = Rs * np.sin(theta) 
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta
    
    
    z = d
    #+ Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    

    
    #Reflection_coe = rfcoe_sm(theta,theta_2)
    Reflection_coe = 1
    
    
    return (1/(2*c)) * n_0 * (1 + np.cos(2*theta)) * \
        I_GB(rou, z, w_0, z_R, P) * Reflection_coe * Rs**2 * np.sin(2*theta)

def integrand_1tz(theta,phi,d, a, Rs, n_0, n_s, w_0, z_R, P):
    
    
    '''integrand of lower surface transmission'''
    
    c = 3 * 10**8
    
    #rou = Rs * np.sin(theta) #represent rou and z by theta
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi))
    
    z = d #+ Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
 
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    
    
    Transmission_coe = 1 - Reflection_coe
    
    return (1/(2*c)) * (n_0 - n_s*np.cos(theta - theta_2)) * \
        I_GB(rou, z, w_0, z_R, P) * Transmission_coe *Rs**2 * np.sin(2*theta)

def integrand_2rz(theta,phi,d, a,Rs, n_0, n_s, w_0, z_R, P):
    
    '''integrand of upper surface reflection'''
    
    #rou = Rs * np.sin(theta) #represent rou and z by theta
    c = 3 * 10**8
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi))
    
    z = d #+ Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    
    
    
    Transmission_coe = 1 - Reflection_coe
    
    return (1/(2*c)) * (n_s * (np.cos(theta - theta_2) + np.cos(3*theta_2 - theta)))* \
        I_GB(rou, z, w_0, z_R, P) * Reflection_coe * Transmission_coe * Rs**2 * np.sin(2*theta)

def integrand_2tz(theta,phi,d, a,Rs, n_0, n_s, w_0, z_R, P):
    
    '''integrand of upper surface transmission'''
    c = 3 * 10**8
    
    #rou = Rs * np.sin(theta) #represent rou and z by theta
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi))
    
    z = d #+ Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    Transmission_coe = 1 - Reflection_coe
    
    return (1/(2*c)) * (n_s * np.cos(theta - theta_2) - n_0 * np.cos(2*(theta - theta_2)))* \
    I_GB(rou, z, w_0, z_R, P) * Transmission_coe * Transmission_coe * Rs**2 * np.sin(2*theta)
    
    
#Calculation of axial force integral
    
def Axial_force_total(d,a,m, Rs, n_0, n_s, w_0, z_R, P):
    '''total axial forces calculated via integration'''
    
    g = 9.8
  
    forcenet = []
    for d_e in d:
        F_1rz = dblquad(integrand, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d_e, a, Rs, n_0, n_s, w_0, z_R, P))
        #F_1tz = dblquad(integrand_1tz, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d_e,a))
        #F_2rz = dblquad(integrand_2rz, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d_e,a))
        #F_2tz = dblquad(integrand_2tz, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d_e,a))
        #F_z = quad(integrand, 0, np.pi/2, args = (d_e)) #returns a tuple with first element the integral result and second element = upper bound error
        #F_1tz = quad(integrand_1tz, 0, np.pi/2, args = (d_e))
        #F_2rz = quad(integrand_2rz, 0, np.pi/2, args = (d_e))
        #F_2tz = quad(integrand_2tz, 0, np.pi/2, args = (d_e))
    
    
        F_znet = F_1rz[0] - m * g #+ F_1tz[0] + F_2rz[0] + F_2tz[0]
        forcenet.append(F_znet * 10 ** 12)

    return forcenet





#Radial force integrads

def integrand_1rr(theta,phi,d, a, Rs, n_0, n_s, w_0, z_R, P):
    
    '''integrand of lower surface reflection, note!!! the rfcoe not considered so it only works for reflection sphere now!!!'''
    c = 3 * 10**8
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta
    
    z = d #+ Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    #Reflection_coe = rfcoe_sm(theta,theta_2)
    
    Reflection_coe = 1
    
    
    return - ((n_0/(2*c)) * np.sin(2*theta) * I_GB(rou,z, w_0, z_R, P) * Reflection_coe * Rs**2 * np.cos(phi) * np.sin(2*theta))

def integrand_1tr(theta,phi,d, a, Rs, n_0, n_s, w_0, z_R, P):
    
    '''integrand of lower surface transmission'''
    c = 3 * 10**8
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta #represent rou and z by theta
    
    z = d #+ Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
 
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    
    Transmission_coe = 1 - Reflection_coe
    
    return n_s/(2*c) * np.sin(theta - theta_2) * I_GB(rou,z,w_0, z_R, P) * Transmission_coe * Rs**2 * np.cos(phi) * np.sin(2*theta)

def integrand_2rr(theta,phi,d, a, Rs, n_0, n_s, w_0, z_R, P):
    
    '''integrand of upper surface reflection'''
    c = 3 * 10**8
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta #represent rou and z by theta
    
    z = d #+ Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
   
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    
    Transmission_coe = 1 - Reflection_coe
    
    return n_s/(2*c) * (np.sin(3*theta_2 - theta) - np.sin(theta - theta_2)) * I_GB(rou,z,w_0, z_R, P) *\
        Reflection_coe * Transmission_coe * Rs**2 * np.cos(phi) * np.sin(2*theta)

def integrand_2tr(theta,phi,d, a, Rs, n_0, n_s, w_0, z_R, P):
    
    '''integrand of upper surface transmission'''
    c = 3 * 10**8
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta #represent rou and z by theta
    
    z = d + Rs * (1 - np.cos(theta))
    
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s)
    
    
    Reflection_coe = rfcoe_sm(theta,theta_2)
    
    
    Transmission_coe = 1 - Reflection_coe
    
    return 1/(2*c) * (n_0 * np.sin(2 * (theta - theta_2)) - n_s * np.sin(theta - theta_2)) * I_GB(rou,z,w_0, z_R, P) * \
        Transmission_coe * Transmission_coe * Rs**2 * np.cos(phi) * np.sin(2*theta)
        
#Calculation of radial force integral
def Radial_force_total(d,a,m, Rs, n_0, n_s, w_0, z_R, P):
    '''total radial forces calculated via integration'''
    
    
    forcenet_r = []
    for a_e in a:
        F_1rr = dblquad(integrand_1rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e,Rs, n_0, n_s, w_0, z_R, P)) #returns a tuple with first element the integral result and second element = upper bound error
        #F_1tr = dblquad(integrand_1tr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e))
       #F_2rr = dblquad(integrand_2rr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e))
        #F_2tr = dblquad(integrand_2tr, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (d,a_e))
    
    
        F_znet = F_1rr[0] #+ F_1tr[0] + F_2rr[0] + F_2tr[0]
        forcenet_r.append(F_znet * 10 ** 12)

    return forcenet_r