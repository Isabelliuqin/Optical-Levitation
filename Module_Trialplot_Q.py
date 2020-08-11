# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 20:22:34 2020

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
from sympy import *
from scipy.special import eval_genlaguerre as LGpoly
#import Module_Gauthier_result_functions as GRF
from scipy.misc import derivative


def criticalangle(n_s, n_0):
    """determine critical angle for total internal reflection"""

    theta_c = np.arcsin( ( n_s.real ) / n_0 )
    
    return theta_c

def r_s(theta, theta2, n_0, n_s):
    """??? Fresnel coefficient of s-polarisation"""
    
    if n_s.real < n_0:
    
        if theta >= criticalangle(n_s, n_0): #total internal reflection, formula see LT6-7 Kenny notes
        
            delta_s = np.arctan(np.sqrt( (np.sin(theta)) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(theta)) #delta has to be real to ensure the modulus sqruare of rs is 1, if not -- attenuation 
        
            r_s1 = np.exp(2 * delta_s * 1j)
        
        else:                                #normal refraction formula:note the imaginary part is included                           
        
            r_s1 = (n_0*np.cos(theta) - n_s*np.cos(theta2))/(n_0*np.cos(theta) + n_s *np.cos(theta2))
    
    else:
        r_s1 = (n_0*np.cos(theta) - n_s*np.cos(theta2))/(n_0*np.cos(theta) + n_s *np.cos(theta2))
    
    return r_s1

def r_p(theta, theta2, n_0, n_s):
    """Fresnel coefficient of p-polarisation"""
    
    if n_s.real < n_0:
    
        if theta >= criticalangle(n_s, n_0):  #total internal reflection, formula see LT6-7 Kenny notes
        
            tan_delta_s = np.sqrt( (np.sin(theta)) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(theta) 
        
            delta_p = np.arctan( n_0 ** 2 / n_s.real ** 2) * tan_delta_s
        
            r_p1 = np.exp(2 * delta_p * 1j)
        
        else:                               #normal refraction formula:note the imaginary part is included                
            r_p1 = (n_s / np.cos(theta2) - n_0 / np.cos(theta))/ (n_0 / np.cos(theta) + n_s / np.cos(theta2))
    
    else: 
        r_p1 = (n_s / np.cos(theta2) - n_0 / np.cos(theta))/ (n_0 / np.cos(theta) + n_s / np.cos(theta2))
    
    return r_p1

def rfcoe_sm(theta,theta2, n_0, n_s):
    
    """Average reflectance"""
    
    r_smodsqr = abs(r_s(theta, theta2, n_0, n_s))**2
    
    r_pmodsqr = abs(r_p(theta, theta2, n_0, n_s))**2
    
    
    return (r_smodsqr + r_pmodsqr)/2

def I_GB(r, w, P, n_0):
    """GB intensity profile"""
    #Gaussian beam propagation parameters
    
    c = 3 * 10**8
    Permittivity = 8.85 * 10**(-12)
    
    I = (2 * P/(np.pi*w**2)) * np.exp( - 2 * r** 2 / w**2)
    E_modulussquare = 2 * I / (c * n_0 * Permittivity)
    
    return E_modulussquare

def TEM01_star(r,w,P):  
    """TEM01* intensity profile"""
    mu_0 = 4*np.pi * 10**(-7)
    c = 3 * 10**8
    
    E_modulussquare2 = (8 * mu_0 * c * P / (np.pi * w**2)) * (r**2 / w**2) * np.exp(- 2 * r**2 / w**2)
    
    return E_modulussquare2


def LG_01_Intensity(r,w):
    
    """intensity of LG mode"""
    
    #r = np.sqrt(x**2 + y**2)
    
    return 2 / np.pi * (1 / w ** 2) * (2 *r**2 / (w ** 2)) *  np.exp(- 2 * r**2 / w**2)





def Fz_integrand_Roosen(theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    """integrand of Fz"""
    
    
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    r = np.sqrt( (rho * np.sin(theta) * np.cos(phi) - rho_0x) ** 2 + (rho * np.sin(theta) * np.sin(phi) - rho_0y) ** 2) #represent r(the parameter of laser beam) by theta
    
    
    if n_s.real < n_0:                              #for metallic reflective sphere or low index sphere
        theta_c = np.arcsin( ( n_s.real ) / n_0 )
    
    
        if theta > theta_c:                         #total internal reflection
        
            theta_2 = np.pi/2
        
        else:                                       #for non TIR
        
            theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) 
    
    else:                                           #for dielectric sphere
        theta_2 = np.arcsin(n_0*np.sin(theta)/n_s.real) 

    
    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    
    R = np.around(Reflection_coe, decimals= 15)     #turn 0.99999 to 1.0
    
    if target == "reflective":                      #reflective sphere has T = 0 by definition, R mignt not = 1
        Tranmission = 0
        
        return -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1)
        
        
    elif R == 1:
        Tranmission = 0
        
        return -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1)
        
    else:                                           #dielectric sphere ignore absorption
        
        Tranmission = 1 - Reflection_coe
    
    
    
        return -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1 - \
         Tranmission**2 * (np.cos(2*(theta - theta_2)) + Reflection_coe * np.cos(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe * np.cos(2*theta_2)) )
   
        
    


def Fy_integrand_Roosen(theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    
    """integrand of Fy"""
    
    
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
        
        return -(n_0 *LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta))
    
    elif R == 1:
        Tranmission = 0
        
        return -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) *\
        (Reflection_coe *np.sin(2 * theta))
    
    else:
    
        Tranmission = 1 - Reflection_coe
    
    
    
        return - (n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))

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
        
        return -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.cos(phi) *\
        (Reflection_coe *np.sin(2 * theta))
    
    else:
    
        Tranmission = 1 - Reflection_coe
    
    
    
        return - (n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.cos(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))

def Fz_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """total axial forces calculated via integration"""
    
    g = 9.8
    F_z = dblquad(Fz_integrand_Roosen, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target))
    
    
    

    return F_z

def Fy_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """total radial forces calculated via integration"""
    
    
    
    F_y = dblquad(Fy_integrand_Roosen, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)) #returns a tuple with first element the integral result and second element = upper bound error
    

    return F_y

def Fx_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """total radial forces calculated via integration"""
    
    
    
    F_x = dblquad(Fx_integrand_Roosen, 0, 2*np.pi, lambda phi: 0, lambda phi: np.pi/2, args = (rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)) #returns a tuple with first element the integral result and second element = upper bound error
    

    return F_x



def Fz_total_vs_d_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    forcenet = []
    
    
    for w_e in w:
        F_1tz = Fz_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
        
    
        F_znet = F_1tz[0]
        forcenet.append(F_znet)

    return forcenet

def Fy_total_vs_rho0y_plot(rho_0x,rho_0y,rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_r = []
    
    for rho_0_e in rho_0y:
        F_1rr = Fy_total(rho_0x, rho_0_e, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet)

    return forcenet_r


def Fx_total_vs_rho0x_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_r = []
    
    for rho_0_e in rho_0x:
        F_1rr = Fx_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet)

    return forcenet_r

def Fz_total_vs_rho0x_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_r = []
    
    for rho_0_e in rho_0x:
        F_1rr = Fz_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet)

    return forcenet_r


def Fx_total_vs_d_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    forcenet = []
    
    
    for w_e in w:
        F_1tz = Fx_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
        
    
        F_znet = F_1tz[0]
        forcenet.append(F_znet)

    return forcenet


def Fy_total_vs_d_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    forcenet = []
    
    
    for w_e in w:
        F_y = Fy_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
        
    
        F_znet = F_y[0]
        forcenet.append(F_znet)

    return forcenet


def xy_plane_vectorplot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_x = []
    
    
    forcenet_y = []
    
    
    for i in range(len(rho_0x)): #5
        forcenet_xe = []
        forcenet_ye = []
        
        for j in range(len(rho_0x[0])):
            
        
            F_xe = Fx_total(rho_0x[i][j],rho_0y[i][j], rho, n_0, n_s, w_0, w, z_R, P, target)
            F_x = F_xe[0]
            
            F_ye = Fy_total(rho_0x[i][j], rho_0y[i][j], rho, n_0, n_s, w_0, w, z_R, P, target)
            F_y = F_ye[0]
    
                   
            forcenet_xe.append(F_x)
            forcenet_ye.append(F_y)
            
        forcenet_x.append(forcenet_xe)
        forcenet_y.append(forcenet_ye)
    return forcenet_x, forcenet_y

def xz_plane_vectorplot(rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_x = []
    
    forcenet_d = []
    
    for i in range(len(rho_0x)): #5
        forcenet_xe = []
        forcenet_de = []
        
        for j in range(len(rho_0x[0])):
            F_xe = Fx_total(rho_0x[i][j], rho_0y, rho, n_0, n_s, w_0, w[i][j], z_R, P, target)
            F_x = F_xe[0]
    
            F_ze = Fz_total(rho_0x[i][j], rho_0y, rho, n_0, n_s, w_0, w[i][j], z_R, P, target)
            F_z = F_ze [0]
    
                   
            forcenet_xe.append(F_x)
            forcenet_de.append(F_z)
            
        forcenet_x.append(forcenet_xe)
        forcenet_d.append(forcenet_de)
    return forcenet_x, forcenet_d

def xyz_plane_vectorplot(rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_x = []
    forcenet_y = []
    forcenet_d = []
    
    for i in range(len(rho_0x)): #5
        forcenet_xe = []
        forcenet_ye = []
        forcenet_de = []
        
        for j in range(len(rho_0x[0])):
            forcenet_xee = []
            forcenet_yee = []
            forcenet_dee = []
            
            
            for k in range(len(rho_0x[0][0])):
                
                F_xe = Fx_total(rho_0x[i][j][k], rho_0y[i][j][k], rho, n_0, n_s, w_0, w[i][j][k], z_R, P, target)
    
                F_x = F_xe[0]
                
                F_ye = Fy_total(rho_0x[i][j][k], rho_0y[i][j][k], rho, n_0, n_s, w_0, w[i][j][k], z_R, P, target)
                F_y = F_ye[0]
    
                F_ze = Fz_total(rho_0x[i][j][k], rho_0y[i][j][k], rho, n_0, n_s, w_0, w[i][j][k], z_R, P, target)
                F_z = F_ze [0]
    
                   
                forcenet_xee.append(F_x)
                forcenet_yee.append(F_y)
                forcenet_dee.append(F_z)
            
            forcenet_xe.append(forcenet_xee)
            forcenet_ye.append(forcenet_yee)
            forcenet_de.append(forcenet_dee)
    
        forcenet_x.append(forcenet_xe)
        forcenet_y.append(forcenet_ye)
        forcenet_d.append(forcenet_de)
        
    return forcenet_x, forcenet_y, forcenet_d



def Fz_total_vs_Rs_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    forcenet = []
    
    
    for rho_e in rho:
        F_1tz = Fz_total(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target)
        
    
        F_znet = F_1tz[0]
        forcenet.append(F_znet)

    return forcenet


def Fy_total_vs_Rs_plot(rho_0x,rho_0y,rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_r = []
    
    for rho_e in rho:
        F_1rr = Fy_total(rho_0x, rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet)

    return forcenet_r


def Fx_total_vs_Rs_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    forcenet_r = []
    
    for rho_e in rho:
        F_1rr = Fx_total(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet)

    return forcenet_r  




def Fx_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """calculate x radial gradient"""
    
    
    
    forcenet_r = []
    
    for rho_0_e in rho_0x:
        F_x = Fx_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
        F_znet = F_x[0] 
        forcenet_r.append(F_znet)
    
    
    grad_Fx = []
    rho_0xlist = []
    
    for i in range( 1, len(forcenet_r) ):
        grad_Fxelement = ( forcenet_r[i] - forcenet_r[i-1] ) / (rho_0x[i] - rho_0x[i-1])
        
        rho_0xelement = (rho_0x[i] + rho_0x[i-1]) / (2 * (w / np.sqrt(2)) )
        
        grad_Fx.append(grad_Fxelement)
        
        rho_0xlist.append(rho_0xelement)
        
    return grad_Fx, rho_0xlist


def Fy_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """calculate x radial gradient"""
    
    
    
    forcenet_r = []
    
    for rho_0_e in rho_0y:
        F_y = Fy_total(rho_0x,rho_0_e, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
        F_znet = F_y[0] 
        forcenet_r.append(F_znet)
    
    
    grad_Fy = []
    rho_0ylist = []
    
    for i in range( 1, len(forcenet_r) ):
        grad_Fyelement = ( forcenet_r[i] - forcenet_r[i-1] ) / (rho_0y[i] - rho_0y[i-1])
        
        rho_0yelement = (rho_0y[i] + rho_0y[i-1]) / (2 * (w / np.sqrt(2)) )
        
        grad_Fy.append(grad_Fyelement)
        
        rho_0ylist.append(rho_0yelement)
        
    return grad_Fy, rho_0ylist

def Fz_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    """calculate x radial gradient"""
    
    
    
    forcenet = []
    
    
    for w_e in w:
        F_1tz = Fz_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
        
    
        F_znet = F_1tz[0]
        forcenet.append(F_znet)

    d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 )) 
    
    grad_Fz = []
    dlist = []
    
    for i in range( 1, len(forcenet) ):
        grad_Fzelement = ( forcenet[i] - forcenet[i-1] ) / (d[i] - d[i-1])
        
        delement = (d[i] + d[i-1]) /2
        
        grad_Fz.append(grad_Fzelement)
        
       
        dlist.append(delement)
        
    return grad_Fz, dlist


def Fx_total_gradient_vs_Rs_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    """stiffness stability plot"""
    forcenet_r = []
    
    
    for rho_e in rho:
        F_1rr = Fx_total_gradient(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target)
    
        
        forcenet_r.append(F_1rr)

    return forcenet_r  




