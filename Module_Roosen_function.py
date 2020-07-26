# -*- coding: utf-8 -*-
"""
Created on Fri Jul 17 21:25:45 2020

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
import cmath
#import Module_Gauthier_result_functions as GRF


def r_s(theta, theta2, n_0, n_s):
    '''Fresnel coefficient of s-polarisation'''
    return (n_0*np.cos(theta) - n_s*np.cos(theta2))/(n_0*np.cos(theta) + n_s*np.cos(theta2))

def r_p(theta, theta2, n_0, n_s):
    '''Fresnel coefficient of p-polarisation'''
    return (n_0 * np.cos(theta2) - n_s * np.cos(theta))/ (n_s * np.cos(theta) + n_0 * np.cos(theta2))

def rfcoe_sm(theta,theta2, n_0, n_s):
    
    '''Average reflectance'''
    return ( ( (r_s(theta,theta2, n_0, n_s).real)**2 + (r_s(theta,theta2, n_0, n_s).imag)**2) + ( (r_p(theta,theta2, n_0, n_s).real)**2 + (r_p(theta,theta2, n_0, n_s).imag)**2 ) ) / 2

def I_GB(rou, w, P, n_0):
    '''GB intensity profile'''
    #Gaussian beam propagation parameters
    
    c = 3 * 10**8
    Permittivity = 8.85 * 10**(-12)
    
    I = (2 * P/(np.pi*w**2)) * np.exp( - 2 * rou** 2 / w**2)
    E_modulussquare = 2 * I / (c * n_0 * Permittivity)
    
    return E_modulussquare

def TEM01_star(rou,w,P):
    '''TEM01* intensity profile'''
    mu_0 = 4*np.pi * 10**(-7)
    c = 3 * 10**8
    
    E_modulussquare2 = (8 * mu_0 * c * P / (np.pi * w**2)) * (rou**2 / w**2) * np.exp(- 2 * rou**2 / w**2)
    
    return E_modulussquare2

#Axial integrands


def axial_integrand_Roosen(theta,phi,a, Rs, n_0, n_s, w_0, w, z_R, P):
    
    
    
    '''integrand of total Fz'''
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 - 2*a*Rs*np.sin(theta) * np.sin(phi)) #represent rou and z by theta
    
    
    if n_s.real < n_0: #metallic reflective sphere or low index sphere
        theta_c = np.arcsin( ( n_s.real ) / n_0 )
    
    
        if theta > theta_c:
        
            theta_2 = np.pi/2
        
        else:
        
            theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) #for dielectric sphere
    
    else:
        theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) #for dielectric sphere

    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    if n_s.imag != 0:
        Tranmission = 0
        
        return -(n_0 *TEM01_star(rou,w,P)) * Rs**2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1)
        
        
    
    else:
        
        Tranmission = 1 - Reflection_coe
    
    
    
        return -(n_0 *TEM01_star(rou,w,P)) * Rs**2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1 - \
         Tranmission**2 * (np.cos(2*(theta - theta_2)) + Reflection_coe * np.cos(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe * np.cos(2*theta_2)) )
   
        
    


def radial_integrand_Roosen(theta,phi,a, Rs, n_0, n_s, w_0, w, z_R, P):
    
    
    
    '''integrand of total Fy'''
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 - 2*a*Rs*np.sin(theta) * np.sin(phi)) #represent rou and z by theta
    
                  
    if n_s.real < n_0: #metallic reflective sphere or low index sphere
        theta_c = np.arcsin( ( n_s.real ) / n_0 )
    
    
        if theta > theta_c:
        
            theta_2 = np.pi/2
        
        else:
        
            theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) #for dielectric sphere
    
    else:
        theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) #for dielectric sphere

    
    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    if n_s.imag != 0:
        Tranmission = 0
        
        return -(n_0 * TEM01_star(rou,w,P)) * Rs**2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta))
    
    
    
    else:
    
        Tranmission = 1 - Reflection_coe
    
    
    
        return - (n_0 *TEM01_star(rou,w,P)) * Rs**2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))


def axial_force_Ashkin(d, a, Rs, n_0, w_0, z_R, P):
    
    
    
    '''integrand of total Fz'''
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta
    
    
    z = z_R * np.sqrt( w**2 / w_0**2 - 1 )
    #+ Rs * (1 - np.cos(theta))
    

    Reflection_coe = 1
    
    
    return n_0 * P/c * (1 + Rs * np.cos(2*theta) - \
         Tranmission**2 * (np.cos(2*(theta - theta_2)) + Reflection_coe * np.cos(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe * np.cos(2*theta_2)) )


def radial_force_Ashkin(theta,phi, d, a, Rs, n_0, w_0, z_R, P):
    
    
    
    '''integrand of total Fy'''
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    
    
    rou = np.sqrt(a**2 + Rs**2 * (np.sin(theta))**2 + 2*a*Rs*np.sin(theta) * np.cos(phi)) #represent rou and z by theta
    
    
    z = d
    #+ Rs * (1 - np.cos(theta))
    

    Reflection_coe = 1
    
    
    return n_0 * P/c * (Rs * np.sin(2*theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))


def Axial_force_total(a,Rs, n_0, n_s, w_0, w, z_R, P):
    '''total axial forces calculated via integration'''
    
    g = 9.8
    F_z = dblquad(axial_integrand_Roosen, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = ( a, Rs, n_0, n_s, w_0, w, z_R, P))
    
    
    

    return F_z

def Radial_force_total(a, Rs, n_0, n_s, w_0, w, z_R, P):
    '''total radial forces calculated via integration'''
    
    
    
    F_y = dblquad(radial_integrand_Roosen, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = ( a ,Rs, n_0, n_s, w_0, w, z_R, P)) #returns a tuple with first element the integral result and second element = upper bound error
    

    return F_y

def Axial_force_total_plot(a,Rs, n_0, n_s, w_0, w, z_R, P):
    
    forcenet = []
    
    
    
    for Rs_e in Rs:
        F_1tz = Axial_force_total(a,Rs_e, n_0, n_s, w_0, w, z_R, P)
        
    
        F_znet = F_1tz[0]
        forcenet.append(abs(F_znet*10**9))

    return forcenet


def Radial_force_total_plot(a, Rs, n_0, n_s, w_0, w, z_R, P):
    forcenet_r = []
    for a_e in a:
        F_1rr = Radial_force_total(a_e, Rs, n_0, n_s, w_0, w, z_R, P)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet*10**9)

    return forcenet_r