# -*- coding: utf-8 -*-
"""
Created on Wed Aug  5 16:32:46 2020

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
import sympy as syp
from scipy.special import eval_genlaguerre as LGpoly
#import Module_Gauthier_result_functions as GRF
from scipy.misc import derivative
from sympy import *
from sympy import sin, cos, pi, exp, sqrt

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


def criticalangle(n_s, n_0):
    """determine critical angle for total internal reflection"""

    theta_c = np.arcsin( ( n_s.real ) / n_0 )
    
    return theta_c

def r_s(theta, theta2, n_0, n_s):
    """??? Fresnel coefficient of s-polarisation"""
    
    if n_s.real < n_0:
    
        if theta >= criticalangle(n_s, n_0): #total internal reflection, formula see LT6-7 Kenny notes
        
            delta_s = np.arctan(sqrt( (sin(theta)) ** 2 - n_s.real**2 / n_0 ** 2) / cos(theta)) #delta has to be real to ensure the modulus sqruare of rs is 1, if not -- attenuation 
        
            r_s1 = exp(2 * delta_s * 1j)
        
        else:                                #normal refraction formula:note the imaginary part is included                           
        
            r_s1 = (n_0*cos(theta) - n_s*cos(theta2))/(n_0*cos(theta) + n_s *cos(theta2))
    
    else:
        r_s1 = (n_0*cos(theta) - n_s*cos(theta2))/(n_0*cos(theta) + n_s *cos(theta2))
    
    return r_s1

def r_p(theta, theta2, n_0, n_s):
    """Fresnel coefficient of p-polarisation"""
    
    if n_s.real < n_0:
    
        if theta >= criticalangle(n_s, n_0):  #total internal reflection, formula see LT6-7 Kenny notes
        
            tan_delta_s = sqrt( (sin(theta)) ** 2 - n_s.real**2 / n_0 ** 2) / cos(theta) 
        
            delta_p = np.arctan( n_0 ** 2 / n_s.real ** 2) * tan_delta_s
        
            r_p1 = exp(2 * delta_p * 1j)
        
        else:                               #normal refraction formula:note the imaginary part is included                
            r_p1 = (n_s / cos(theta2) - n_0 /cos(theta))/ (n_0 / cos(theta) + n_s / cos(theta2))
    
    else: 
        r_p1 = (n_s / cos(theta2) - n_0 /cos(theta))/ (n_0 /cos(theta) + n_s /cos(theta2))
    
    return r_p1

def rfcoe_sm(theta,theta2, n_0, n_s):
    
    """Average reflectance"""
    
    r_smodsqr = abs(r_s(theta, theta2, n_0, n_s))**2
    
    r_pmodsqr = abs(r_p(theta, theta2, n_0, n_s))**2
    
    
    return (r_smodsqr + r_pmodsqr)/2


def LG_01_Intensity(r,w):
    
    """intensity of LG mode"""
    
    #r = np.sqrt(x**2 + y**2)
    
    return 2 / pi * (1 / w ** 2) * (2 *r**2 / (w ** 2)) *  exp(- 2 * r**2 / w**2)


def Fy_integrand_Roosen(theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    
    """integrand of Fy"""
    
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    r = sqrt( (rho * sin(theta) * cos(phi) - rho_0x) ** 2 + (rho * sin(theta) * sin(phi) - rho_0y) ** 2) #represent rou and z by theta
    
                  
    if n_s.real < n_0: 
        theta_c = arcsin( ( n_s.real ) / n_0 )
    
    
        if theta > theta_c:
        
            theta_2 = pi/2
        
        else:
        
            theta_2 = arcsin(n_0*sin(theta)/n_s.real) 
    
    else:
        theta_2 = arcsin(n_0*sin(theta)/n_s.real) 

    
    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    R = np.around(Reflection_coe, decimals= 15)
    
    if target == "reflective":
        Tranmission = 0
        
        return -(n_0 *LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * sin(theta) * cos(theta) * sin(phi) * \
        (Reflection_coe * sin(2 * theta))
    
    elif R == 1:
        Tranmission = 0
        
        return -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * sin(theta) * cos(theta) * \
        (Reflection_coe * cos(2*theta) + 1)
    
    else:
    
        Tranmission = 1 - Reflection_coe
    
    
    
        return - (n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * sin(theta) * cos(theta) * sin(phi) * \
        (Reflection_coe * sin(2 * theta) - \
         Tranmission**2 * (sin(2*(theta - theta_2)) + Reflection_coe * sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*cos(2*theta_2)))

def Fx_integrand_Roosen(theta, phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P,target):
    
    
    """integrand of Fx"""
   
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    r = np.sqrt( (rho * np.sin(theta) * np.cos(phi) - rho_0x) ** 2 + (rho * np.sin(theta) * np.sin(phi) - rho_0y) ** 2) #represent rou and z by theta
    
                  
    if n_s.real < n_0: 
        theta_c = np.arcsin( ( n_s.real ) / n_0 )
    
    
        if theta > theta_c:
        
            theta_2 = np.pi/2
        
        else:
        
            theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) 
    
    else:
        theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) 

    
    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    R = np.around(Reflection_coe, decimals= 15)
    
    if target == "reflective":
        Tranmission = 0
        
        return -(n_0 *LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.cos(phi) * \
        (Reflection_coe * np.sin(2 * theta))
    
    elif R == 1:
        Tranmission = 0
        
        return -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1)
    
    else:
    
        Tranmission = 1 - Reflection_coe
    
    
    
        return - (n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.cos(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))


def Fy_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """total radial forces calculated via integration"""
    
    theta, phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P = symbols("x y a b c d e f g h k")
    expr = Fy_integrand_Roosen(theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    F_y = expr.integrate(theta, phi, (theta, 0, pi/2), (phi, 0, 2*pi), manual = True) #returns a tuple with first element the integral result and second element = upper bound error
    

    return F_y

def Fx_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """total radial forces calculated via integration"""
    
    theta, phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P = symbols("x y a b c d e f g h k")
    
    expr = Fx_integrand_Roosen(theta, phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P,target)
    F_x = expr.integrate(theta, phi, (theta, 0, pi/2), (phi, 0, 2*pi), manual = True) #returns a tuple with first element the integral result and second element = upper bound error
    

    return F_x



