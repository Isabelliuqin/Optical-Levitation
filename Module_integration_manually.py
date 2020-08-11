# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 10:52:49 2020

@author: liuqi
"""


"""Do integration Mantually!"""



import scipy as sp
import numpy as np

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pylab as plt
import scipy.integrate as spi
from scipy.integrate import quad
import seaborn
from scipy.integrate import odeint
from scipy.integrate import dblquad



def criticalangle(n_s, n_0):
    """determine critical angle for total internal reflection"""

    theta_c = np.arcsin( ( n_s.real ) / n_0 )
    
    return theta_c



def Fz_integrand(Theta, Phi, theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    r = np.sqrt( (rho * np.sin(Theta) * np.cos(Phi) - rho_0x) ** 2 + (rho * np.sin(Theta) * np.sin(Phi) - rho_0y) ** 2) #represent r(the parameter of laser beam) by theta
    
    
    if n_s.real < n_0:                              #for metallic reflective sphere or low index sphere
        theta_c = np.arcsin( ( n_s.real ) / n_0 )
        
        
        theta_2 = [ [ None for p in range( len(theta) ) ]
                           for q in range( len(phi) ) ] 
        
        
        delta_s = [ [ None for p in range( len(theta) ) ]
                           for q in range( len(phi) ) ] 
        r_s1 = [ [ None for p in range( len(theta) ) ]
                        for q in range( len(phi) ) ] 
        
        
        tan_delta_s = [ [ None for p in range( len(theta) ) ]
                                   for q in range( len(phi) ) ] 
        delta_p = [ [ None for p in range( len(theta) ) ]
                                   for q in range( len(phi) ) ] 
                
        r_p1 = [ [ None for p in range( len(theta) ) ]
                                for q in range( len(phi) ) ] 
        
        
        r_smodsqr = [ [ None for p in range( len(theta) ) ]
                         for q in range( len(phi) ) ] 
        r_pmodsqr = [ [ None for p in range( len(theta) ) ]
                         for q in range( len(phi) ) ] 
    
        for i in range(len(theta)):
            
            for k in range(len(phi)):
                
                
             
                
                if Theta[k][i] > theta_c:
        
                    theta_2[k][i] = np.pi/2
                    
                    
                    
                    delta_s[k][i] = np.arctan(np.sqrt( (np.sin(Theta[k][i])) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(Theta[k][i])) #delta has to be real to ensure the modulus sqruare of rs is 1, if not -- attenuation 
        
                    r_s1[k][i] = np.exp(2 * delta_s[k][i] * 1j)
                    
                    
                    
                    tan_delta_s[k][i] = np.sqrt( (np.sin(Theta[k][i])) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(Theta[k][i]) 
        
                    delta_p[k][i] = np.arctan( n_0 ** 2 / n_s.real ** 2) * tan_delta_s[k][i]
        
                    r_p1[k][i] = np.exp(2 * delta_p[k][i] * 1j)
                    
                    r_smodsqr[k][i] = (abs( r_s1[k][i])) ** 2
            
                    r_pmodsqr[k][i] = (abs(r_p1[k][i])) ** 2
        
                else:
        
                    theta_2[k][i] = np.arcsin(n_0*np.sin(Theta[k][i])/n_s.real) 
                    
                    r_s1[k][i] = (n_0*np.cos(Theta[k][i]) - n_s*np.cos(theta_2[k][i]))/(n_0*np.cos(Theta[k][i]) + n_s *np.cos(theta_2[k][i]))
        
                    r_p1[k][i] = (n_s / np.cos(theta_2[k][i]) - n_0 / np.cos(Theta[k][i]))/ (n_0 / np.cos(Theta[k][i]) + n_s / np.cos(theta_2[k][i]))
                    
                    r_smodsqr[k][i] = (abs( r_s1[k][i])) ** 2
            
                    r_pmodsqr[k][i] = (abs(r_p1[k][i])) ** 2
    
    else:                                           #for dielectric sphere
        
        
        theta_2 = np.arcsin(n_0*np.sin(Theta)/n_s.real) 
        
        r_s1 = (n_0*np.cos(Theta) - n_s*np.cos(theta_2))/(n_0*np.cos(Theta) + n_s *np.cos(theta_2))
    
        r_p1 = (n_s / np.cos(theta_2) - n_0 / np.cos(Theta))/ (n_0 / np.cos(Theta) + n_s / np.cos(theta_2))
        
        r_smodsqr = (abs( r_s1)) ** 2
            
        r_pmodsqr = (abs(r_p1)) ** 2
        
    
    Reflection_coe =  (np.asarray(r_smodsqr) + np.asarray(r_pmodsqr)) / 2 
    
    R = np.around(Reflection_coe, decimals= 15) 
    
    if target == "reflective":                      #reflective sphere has T = 0 by definition, R mignt not = 1
        Tranmission = 0
        
        Fz_integrand = -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(Theta) * np.cos(Theta) * \
        (Reflection_coe * np.cos(2*Theta) + 1)
        
        return Fz_integrand
        
    
    else:   
        for i in range(len(theta)):
            
            for k in range(len(phi)):
                if R[k][i] == 1:
                    Tranmission = 0
        
                    Fz_integrand[k][i] =  -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(Theta) * np.cos(Theta) * \
                (Reflection_coe * np.cos(2*Theta) + 1)
        
                else:                                           #dielectric sphere ignore absorption
        
                    Tranmission = 1 - Reflection_coe
    
                    Fz_integrand[k][i] = -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(Theta) * np.cos(Theta) * \
        (Reflection_coe * np.cos(2*Theta) + 1 - \
         Tranmission**2 * (np.cos(2*(Theta - theta_2)) + Reflection_coe * np.cos(2*Theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe * np.cos(2*theta_2)) )
    
        return Fz_integrand
    
    
    

def Fy_integrand(Theta, Phi, theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):

    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    r = np.sqrt( (rho * np.sin(Theta) * np.cos(Phi) - rho_0x) ** 2 + (rho * np.sin(Theta) * np.sin(Phi) - rho_0y) ** 2) #represent r(the parameter of laser beam) by theta
    
    
    if n_s.real < n_0:                              #for metallic reflective sphere or low index sphere
        theta_c = np.arcsin( ( n_s.real ) / n_0 )
        
        
        theta_2 = [ [ None for p in range( len(theta) ) ]
                           for q in range( len(phi) ) ] 
        
        
        delta_s = [ [ None for p in range( len(theta) ) ]
                           for q in range( len(phi) ) ] 
        r_s1 = [ [ None for p in range( len(theta) ) ]
                        for q in range( len(phi) ) ] 
        
        
        tan_delta_s = [ [ None for p in range( len(theta) ) ]
                                   for q in range( len(phi) ) ] 
        delta_p = [ [ None for p in range( len(theta) ) ]
                                   for q in range( len(phi) ) ] 
                
        r_p1 = [ [ None for p in range( len(theta) ) ]
                                for q in range( len(phi) ) ] 
        
        
        r_smodsqr = [ [ None for p in range( len(theta) ) ]
                         for q in range( len(phi) ) ] 
        r_pmodsqr = [ [ None for p in range( len(theta) ) ]
                         for q in range( len(phi) ) ] 
    
        for i in range(len(theta)):
            
            for k in range(len(phi)):
                
                
             
                
                if Theta[k][i] > theta_c:
        
                    theta_2[k][i] = np.pi/2
                    
                    
                    
                    delta_s[k][i] = np.arctan(np.sqrt( (np.sin(Theta[k][i])) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(Theta[k][i])) #delta has to be real to ensure the modulus sqruare of rs is 1, if not -- attenuation 
        
                    r_s1[k][i] = np.exp(2 * delta_s[k][i] * 1j)
                    
                    
                    
                    tan_delta_s[k][i] = np.sqrt( (np.sin(Theta[k][i])) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(Theta[k][i]) 
        
                    delta_p[k][i] = np.arctan( n_0 ** 2 / n_s.real ** 2) * tan_delta_s[k][i]
        
                    r_p1[k][i] = np.exp(2 * delta_p[k][i] * 1j)
                    
                    r_smodsqr[k][i] = (abs( r_s1[k][i])) ** 2
            
                    r_pmodsqr[k][i] = (abs(r_p1[k][i])) ** 2
        
                else:
        
                    theta_2[k][i] = np.arcsin(n_0*np.sin(Theta[k][i])/n_s.real) 
                    
                    r_s1[k][i] = (n_0*np.cos(Theta[k][i]) - n_s*np.cos(theta_2[k][i]))/(n_0*np.cos(Theta[k][i]) + n_s *np.cos(theta_2[k][i]))
        
                    r_p1[k][i] = (n_s / np.cos(theta_2[k][i]) - n_0 / np.cos(Theta[k][i]))/ (n_0 / np.cos(Theta[k][i]) + n_s / np.cos(theta_2[k][i]))
                    
                    r_smodsqr[k][i] = (abs( r_s1[k][i])) ** 2
            
                    r_pmodsqr[k][i] = (abs(r_p1[k][i])) ** 2
    
    else:                                           #for dielectric sphere
        
        
        theta_2 = np.arcsin(n_0*np.sin(Theta)/n_s.real) 
        
        r_s1 = (n_0*np.cos(Theta) - n_s*np.cos(theta_2))/(n_0*np.cos(Theta) + n_s *np.cos(theta_2))
    
        r_p1 = (n_s / np.cos(theta_2) - n_0 / np.cos(Theta))/ (n_0 / np.cos(Theta) + n_s / np.cos(theta_2))
        
        r_smodsqr = (abs( r_s1)) ** 2
            
        r_pmodsqr = (abs(r_p1)) ** 2
        
    
    Reflection_coe =  (np.asarray(r_smodsqr) + np.asarray(r_pmodsqr)) / 2 
    
    R = np.around(Reflection_coe, decimals= 15)  

    if target == "reflective":
        Tranmission = 0
        
        Fy_integrand = -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta))
        
        return Fy_integrand
    
    else: 
        for i in range(len(theta)):
            
            for k in range(len(phi)):
        
                if R[k][i] == 1:
                    Tranmission = 0
                    
                    Fy_integrand[k][i] = -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) *\
        (Reflection_coe *np.sin(2 * theta))
        
        
    
                else:
    
                    Tranmission = 1 - Reflection_coe
                    
                    Fy_integrand[k][i] = - (n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))
    
    
        return Fy_integrand
        

def Fx_integrand(Theta, Phi, theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):


    c = 3 * 10**8
    mu_0 = 4*np.pi * 10**(-7)
    Permittivity = 8.85 * 10**(-12) 
    
    r = np.sqrt( (rho * np.sin(Theta) * np.cos(Phi) - rho_0x) ** 2 + (rho * np.sin(Theta) * np.sin(Phi) - rho_0y) ** 2) #represent r(the parameter of laser beam) by theta
    
    
    if n_s.real < n_0:                              #for metallic reflective sphere or low index sphere
        theta_c = np.arcsin( ( n_s.real ) / n_0 )
        
        
        theta_2 = [ [ None for p in range( len(theta) ) ]
                           for q in range( len(phi) ) ] 
        
        
        delta_s = [ [ None for p in range( len(theta) ) ]
                           for q in range( len(phi) ) ] 
        r_s1 = [ [ None for p in range( len(theta) ) ]
                        for q in range( len(phi) ) ] 
        
        
        tan_delta_s = [ [ None for p in range( len(theta) ) ]
                                   for q in range( len(phi) ) ] 
        delta_p = [ [ None for p in range( len(theta) ) ]
                                   for q in range( len(phi) ) ] 
                
        r_p1 = [ [ None for p in range( len(theta) ) ]
                                for q in range( len(phi) ) ] 
        
        
        r_smodsqr = [ [ None for p in range( len(theta) ) ]
                         for q in range( len(phi) ) ] 
        r_pmodsqr = [ [ None for p in range( len(theta) ) ]
                         for q in range( len(phi) ) ] 
    
        for i in range(len(theta)):
            
            for k in range(len(phi)):
                
                
             
                
                if Theta[k][i] > theta_c:
        
                    theta_2[k][i] = np.pi/2
                    
                    
                    
                    delta_s[k][i] = np.arctan(np.sqrt( (np.sin(Theta[k][i])) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(Theta[k][i])) #delta has to be real to ensure the modulus sqruare of rs is 1, if not -- attenuation 
        
                    r_s1[k][i] = np.exp(2 * delta_s[k][i] * 1j)
                    
                    
                    
                    tan_delta_s[k][i] = np.sqrt( (np.sin(Theta[k][i])) ** 2 - n_s.real**2 / n_0 ** 2) / np.cos(Theta[k][i]) 
        
                    delta_p[k][i] = np.arctan( n_0 ** 2 / n_s.real ** 2) * tan_delta_s[k][i]
        
                    r_p1[k][i] = np.exp(2 * delta_p[k][i] * 1j)
                    
                    r_smodsqr[k][i] = (abs( r_s1[k][i])) ** 2
            
                    r_pmodsqr[k][i] = (abs(r_p1[k][i])) ** 2
        
                else:
        
                    theta_2[k][i] = np.arcsin(n_0*np.sin(Theta[k][i])/n_s.real) 
                    
                    r_s1[k][i] = (n_0*np.cos(Theta[k][i]) - n_s*np.cos(theta_2[k][i]))/(n_0*np.cos(Theta[k][i]) + n_s *np.cos(theta_2[k][i]))
        
                    r_p1[k][i] = (n_s / np.cos(theta_2[k][i]) - n_0 / np.cos(Theta[k][i]))/ (n_0 / np.cos(Theta[k][i]) + n_s / np.cos(theta_2[k][i]))
                    
                    r_smodsqr[k][i] = (abs( r_s1[k][i])) ** 2
            
                    r_pmodsqr[k][i] = (abs(r_p1[k][i])) ** 2
    
    else:                                           #for dielectric sphere
        
        
        theta_2 = np.arcsin(n_0*np.sin(Theta)/n_s.real) 
        
        r_s1 = (n_0*np.cos(Theta) - n_s*np.cos(theta_2))/(n_0*np.cos(Theta) + n_s *np.cos(theta_2))
    
        r_p1 = (n_s / np.cos(theta_2) - n_0 / np.cos(Theta))/ (n_0 / np.cos(Theta) + n_s / np.cos(theta_2))
        
        r_smodsqr = (abs( r_s1)) ** 2
            
        r_pmodsqr = (abs(r_p1)) ** 2
        
    
    Reflection_coe =  (np.asarray(r_smodsqr) + np.asarray(r_pmodsqr)) / 2 
    
    R = np.around(Reflection_coe, decimals= 15) 


    if target == "reflective":
        Tranmission = 0
        
        Fx_integrand = -(n_0 *LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.cos(phi) * \
        (Reflection_coe * np.sin(2 * theta))
        
        return Fx_integrand
    
    else: 
        for i in range(len(theta)):
            
            for k in range(len(phi)):
    
                if R[k][i] == 1:
                    Tranmission = 0
                    
                    Fx_integrand[k][i] = -(n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.cos(phi) *\
        (Reflection_coe *np.sin(2 * theta))
        
        
    
                else:
    
                    Tranmission = 1 - Reflection_coe
                    
                    Fx_integrand[k][i] = - (n_0 * LG_01_Intensity(r,w)) * rho ** 2 / (2 * mu_0 * c**2) * np.sin(theta) * np.cos(theta) * np.cos(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))
    
    
    
        return Fx_integrand


def power_integrand(r,w, n_0):    
    
    c = 3 * 10**8
    Permittivity = 8.85 * 10**(-12)
    
    power_integrand = 0.5 * c * Permittivity * n_0 * 2 * np.pi * LG_01_Intensity(r,w) * r
    
    return power_integrand


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
    
    """modulus square of LG mode"""
    
    #r = np.sqrt(x**2 + y**2)
    
    return 2 / np.pi * (1 / w ** 2) * (2 *r**2 / (w ** 2)) *  np.exp(- 2 * r**2 / w**2)




#Integration

def dFz_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, x,y,z, target):
    
    Fz_slidelist = []
    
    x_list = []
    
    for i in range(1, len(phi)):
    
        thetad = np.linspace(theta, np.pi - theta, gridsize)
        phi = np.linspace(phi[i-1], phi[i], 100)
    
        Theta,Phi = np.meshgrid(thetad, phi)
    
        
        delta_phi = 2*np.pi * 1/(100 - 1)
        delta_theta = (np.pi - 2*theta) / (gridsize - 1)
    
    
    
        integrand = np.asarray(Fz_integrand(Theta, Phi, theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target))
    
        surf_element = delta_phi * delta_theta
    
        Fz_slide = sum(sum(integrand * surf_element))
        
        Fz_slidelist.append(Fz_slide)
        
        phi_ave = phi[i-1] + phi[1]/2
        
        x = rho * np.cos(phi_ave) * np.sin(theta)
        
        x_list.append(x)
    
    return Fz_slide, x_list





def Fz_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    theta = np.linspace(0, np.pi/2, gridsize)
    phi = np.linspace(0, np.pi * 2, gridsize)
    
    Theta,Phi = np.meshgrid(theta, phi)
    
    delta_theta = np.pi/2 * 1/(gridsize - 1)
    delta_phi = 2*np.pi * 1/(gridsize - 1)
    
    
    
    integrand = np.asarray(Fz_integrand(Theta, Phi, theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target))
    
    surf_element = delta_theta * delta_phi
    
    Fz_total = sum(sum(integrand * surf_element))
    
    return Fz_total
    
def Fy_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    theta = np.linspace(0, np.pi/2, gridsize)
    phi = np.linspace(0, np.pi * 2, gridsize)
    
    Theta,Phi = np.meshgrid(theta, phi)
    
    delta_theta = np.pi/2 * 1/(gridsize - 1)
    delta_phi = 2*np.pi * 1/(gridsize - 1)
    
    
    
    integrand = np.asarray(Fy_integrand(Theta, Phi, theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target))
    
    surf_element = delta_theta * delta_phi
    
    Fy_total = sum(sum(integrand * surf_element))
    
    return Fy_total

def Fx_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    theta = np.linspace(0, np.pi/2, gridsize)
    phi = np.linspace(0, np.pi * 2, gridsize)
    
    Theta,Phi = np.meshgrid(theta, phi)
    
    delta_theta = np.pi/2 * 1/(gridsize - 1)
    delta_phi = 2*np.pi * 1/(gridsize - 1)
    
    
    
    integrand = np.asarray(Fx_integrand(Theta, Phi, theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target))
    
    surf_element = delta_theta * delta_phi
    
    Fx_total = sum(sum(integrand * surf_element))
    
    return Fx_total

def Power_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, rrange, target):
    
    
    
    r = np.linspace(0, rrange, gridsize)
    
    
    
    delta_r = rrange * 1/(gridsize - 1)
    
    
    
    
    integrand = np.asarray(power_integrand(r,w, n_0))
    
    
    
    power_total = sum(integrand * delta_r )
    
    return power_total


#F vs own axis offset plot
def Fz_total_vs_w_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    """Fz total vs axial offset rho_0z plot"""
    
    forcenet = []
    
    
    for w_e in w:
        F_1tz = Fz_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, gridsize, target)
        print (w_e)
    
        F_znet = F_1tz
        forcenet.append(F_znet)

    return forcenet

def Fy_total_vs_rho0y_plot(rho_0x,rho_0y,rho, n_0, n_s, w_0, w, z_R, P,gridsize, target):
    
    """Fy total vs y axis offset plot"""
    forcenet_r = []
    
    for rho_0_e in rho_0y:
        F_1rr = Fy_Manual_integration(rho_0x,rho_0_e, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target)
    
        print (rho_0_e)
    
        F_znet = F_1rr 
        forcenet_r.append(F_znet)

    return forcenet_r


def Fx_total_vs_rho0x_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    
    """Fx total vs x axis offset plot"""
    forcenet_r = []
    
    for rho_0_e in rho_0x:
        F_1rr = Fx_Manual_integration(rho_0_e, rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target)
    
        print (rho_0_e)
    
        F_znet = F_1rr
        forcenet_r.append(F_znet)
        
    return forcenet_r
        
        
#crosstalk F vs other axis offset plot
        
def Fz_total_vs_rho0x_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
     """Fz total vs radial offset rho_0x plot"""
     forcenet_r = []
    
     for rho_0_e in rho_0x:
         F_1rr = Fz_Manual_integration(rho_0_e, rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target)
    
    
    
         F_znet = F_1rr 
         forcenet_r.append(F_znet)

     return forcenet_r

def Fy_total_vs_w_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    """Fx total vs axial offset plot"""
    
    forcenet = []
    
    
    for w_e in w:
        F_1tz = Fy_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, gridsize, target)
        
    
        F_znet = F_1tz
        forcenet.append(F_znet)

    return forcenet

def Fx_total_vs_w_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    """Fx total vs axial offset plot"""
    
    forcenet = []
    
    
    for w_e in w:
        F_1tz = Fx_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, gridsize, target)
        
    
        F_znet = F_1tz
        forcenet.append(F_znet)

    return forcenet



#F vs target radius
    
def Fz_total_vs_Rs_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    """Fz total vs radius rho plot"""
    
    forcenet = []
    
    
    for rho_e in rho:
        F_1tz = Fz_Manual_integration(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, gridsize, target)
        
    
        F_znet = F_1tz[0]
        forcenet.append(F_znet)

    return forcenet


def Fy_total_vs_Rs_plot(rho_0x,rho_0y,rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    """Fy total vs radius rho plot"""
    forcenet_r = []
    
    for rho_e in rho:
        F_1rr = Fy_Manual_integration(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, gridsize, target)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet)

    return forcenet_r


def Fx_total_vs_Rs_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    """Fx total vs radius rho plot"""
    forcenet_r = []
    
    for rho_e in rho:
        F_1rr = Fx_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target)
    
    
    
        F_znet = F_1rr[0] 
        forcenet_r.append(F_znet)

    return forcenet_r  

#Gradient plots
def Fx_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    
    """calculate x radial gradient"""
    
    
    
    forcenet_r = []
    
    for rho_0_e in rho_0x:
        F_x = Fx_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target)
    
    
    
        F_xnet = F_x[0] 
        forcenet_r.append(F_xnet)
    
    
    grad_Fx = []
    rho_0xlist = []
    
    for i in range( 1, len(forcenet_r) ):
        grad_Fxelement = ( forcenet_r[i] - forcenet_r[i-1] ) / (rho_0x[i] - rho_0x[i-1])
        
        rho_0xelement = (rho_0x[i] + rho_0x[i-1]) /2
        
        grad_Fx.append(grad_Fxelement)
        
        rho_0xlist.append(rho_0xelement)
        
    return grad_Fx, rho_0xlist


def Fy_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    """calculate y radial gradient"""
    
    
    
    forcenet = []
    
    
    for w_e in w:
        F_y = Fy_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target)
        
    
        F_ynet = F_y[0]
        forcenet.append(F_ynet)

    d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 )) 
    
    grad_Fy = []
    rho_0ylist = []
    
    for i in range( 1, len(forcenet) ):
        grad_Fyelement = ( forcenet[i] - forcenet[i-1] ) / (d[i] - d[i-1])
        
        rho_0yelement = (d[i] + d[i-1]) /2
        
        grad_Fy.append(grad_Fyelement)
        
       
        rho_0ylist.append(rho_0yelement)
        
    return grad_Fy, rho_0ylist



def Fz_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target):
    """calculate z axial gradient"""
    
    
    
    forcenet = []
    
    
    for w_e in w:
        F_z = Fz_Manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, gridsize, target)
        
    
        F_znet = F_z
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