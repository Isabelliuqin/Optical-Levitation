# -*- coding: utf-8 -*-
"""
Created on Sat Aug 29 16:49:26 2020

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
from scipy.misc import derivative


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


###################
##################

def LG_01_Intensity_sphere_coordinate(theta, phi, rho_0x, rho_0y, rho, w, n_0, P):
    
    """intensity of LG01 mode given sphere angular coordinates and offsets"""
    # theta - sphere theta, array
    # phi - sphere phi, array
    # rho_0x - beam displacement in x, scalar
    # rho_0y - beam displacement in y, scalar
    # rho - sphere radius, scalar
    # w - beam radius at target, scalar
    # n0 - refractive index of the surrounding
    # P - power of the laser beam
    
    
    r = beam_r(theta, phi, rho_0x, rho_0y, rho)
    
    return (P * 2 / np.pi) * (1 / w ** 2) * (2 *r**2 / (w ** 2)) *  np.exp(- 2 * r**2 / w**2)


def beam_r(theta, phi, rho_0x, rho_0y, rho):
    """ beam radial coordinate with respect to integration variables"""
    # theta - sphere theta, array
    # phi - sphere phi, array
    # rho_0x - beam displacement in x, scalar
    # rho_0y - beam displacement in y, scalar
    # rho - sphere radius, scalar

    return np.sqrt( (rho * np.sin(theta) * np.cos(phi) - rho_0x) ** 2 + (rho * np.sin(theta) * np.sin(phi) - rho_0y) ** 2)


def beam_w(theta, phi, rho, w, w_0, z_R):
    """ beam radial coordinate with respect to integration variables"""
    # theta - sphere theta, array
    # phi - sphere phi, array
    # w - beam radius at rho_0z: easier for the rest of the work
    # rho - sphere radius, scalar
    # w0 - beam waist
    # z_R - rayleigh range
    
    rho_0z = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 )) 
    
    axial_offset_at_artibrary_locationonsphere = rho_0z + rho*(1 - np.cos(theta))
    
    wnew = w_0 * np.sqrt(1 + axial_offset_at_artibrary_locationonsphere**2/z_R**2)

    return wnew

def spherical_to_cartesian_coordinates(r, theta, phi):
    """ convert spherical coordinates to cartesian"""
    # r - radial coordinate, array
    # theta - zenith angle, array
    # phi - azimuthal angle, array
    
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    
    return {'x' : x, 'y' : y, 'z' : z}


def r_s(theta, theta2, n_0, n_s):
    """Fresnel coefficient of s-polarisation"""
    # theta - incident angle, array
    # theta_2 - refracted angle, array
    # n_0 - incident refractive index, scalar
    # n_s - transmitted refractive index ,scalar
    
    
    theta = np.atleast_2d(theta)
    
    r_s1 = np.atleast_2d(np.array(theta, dtype='complex'))
    
    r_s1[:] = 0
        
    r_s1 = (n_0*np.cos(theta) - n_s*np.cos(theta2))/(n_0*np.cos(theta) + n_s *np.cos(theta2))
    
    return r_s1

def r_p(theta, theta2, n_0, n_s):
    """Fresnel coefficient of p-polarisation"""
    # theta - incident angle, array
    # theta_2 - refracted angle, array
    # n_0 - incident refractive index, scalar
    # n_s - transmitted refractive index ,scalar
    
        
    theta = np.atleast_2d(theta)
    
    r_p1 = np.atleast_2d(np.array(theta, dtype='complex'))
    
    r_p1[:] = np.NaN
        
    r_p1 = (n_0 * np.cos(theta2) - n_s * np.cos(theta))/ (n_0 * np.cos(theta2) + n_s * np.cos(theta))
    
    return r_p1

def rfcoe_sm(theta,theta2, n_0, n_s):#
    
    """Average reflectance"""
    
    r_smodsqr = abs(r_s(theta, theta2, n_0, n_s))**2
    
    r_pmodsqr = abs(r_p(theta, theta2, n_0, n_s))**2
    
    print (r_smodsqr)
    
    print (r_pmodsqr)
    
    return (r_smodsqr + r_pmodsqr)/2

def plot_intensity_sphere_overlap(theta, phi, rho_0x, rho_0y, rho, w, n_0, P):
    """ plot intensity profile and sphere outline"""
    # theta - sphere theta coordinate, array
    # phi - sphere phi coordinate, array
    # rho_0x - beam x displacement, scalar
    # rho_0y - beam y displacement, scalar
    # rho - sphere diameter, scalar
    
    intensity = LG_01_Intensity_sphere_coordinate(theta, phi, rho_0x, rho_0y, rho, w, n_0, P)
    
    cartesian_coordinates = spherical_to_cartesian_coordinates(5*rho, theta, phi)
    x = cartesian_coordinates['x']
    y = cartesian_coordinates['y']
    
    #x = np.linspace(-3 * rho, 3 * rho, 300)
    #y = np.linspace(-3 * rho, 3 * rho, 300)
    
    plt.contourf(x* 10 ** 6, y* 10 ** 6, intensity/np.max(intensity), 500, cmap='Reds')
    cbar = plt.colorbar()

    cbar.set_ticks([0,0.2, 0.4, 0.6, 0.8, 1])

    plt.xlabel('x(um)',fontsize=20)
    plt.ylabel('y(um)',fontsize=20)
    plt.xlim(-120, 120)
    plt.ylim(120, -120)
    '''
    plt.figure("beam_sphere_overlap_pcolormesh")
    plt.clf()
    plt.pcolormesh(x, y, intensity)
    '''
    # make sure coordiantes not distorted
    axes = plt.gca()
    axes.set_aspect('equal')
    
    sphere_outline = plt.Circle((0,0), radius = rho, facecolor = 'None', edgecolor = 'k')
    
    axes.add_artist(sphere_outline)
    
    #plt.figure("beam_sphere_overlap_surf")
    
    #plt.su
    return

#######################
def get_refracted_angle(theta, n_s, n_0):
    """ calculate angle of refracted ray"""
    # theta - angle of incidence (can be array)
    # n_s - refractive index of second medium, scalar, can be imaginary
    # n_0  - refractive index of first medium, scalar, must be real
    
    
    theta = np.atleast_2d(np.array(theta))  # make theta a numpy array
    theta_2 = np.copy(theta)   # initialise thea-2 array
        
        
    theta_2 = np.arcsin(n_0*np.sin(theta)/n_s) 
        
    return theta_2

def Fz_integrand_Roosen(theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    """integrand of Fz"""
    
    # theta - sphere theta coordinate, array
    # phi - sphere phi coordinate, array
    # rho_0x - beam x displacement, scalar
    # rho_0y - beam y displacement, scalar
    # rho - sphere diameter, scalar
    # n_0 - refractive index of medium
    # n_s - refractive index of sphere
    # w_0 - beam radius at waist
    # w - beam radius at target
    # z_R - beam rayleigh range
    # P - beam power (Watts)
    # target - string, set as 'reflective' if sphere is reflective, otherwise anything OK and will be ignored
    
    c = 3 * 10**8
    
    mu_0 = 4*np.pi * 10**(-7)
    
    Permittivity = 8.85 * 10**(-12) 
    
    P_norm = 0.5 * c * n_0 * Permittivity 
    
    r = beam_r(theta, phi, rho_0x, rho_0y, rho) #represent r(the parameter of laser beam) by theta
    
    wnew = beam_w(theta, phi, rho, w, w_0, z_R)    #represent w(parameter of laser beam a variable)
    

    theta_2 = get_refracted_angle(theta, n_s, n_0)
    
    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    
    I = LG_01_Intensity_sphere_coordinate(theta, phi, rho_0x, rho_0y, rho, wnew, n_0, P)   # beam intensity
    
    if target == "reflective":                      #reflective sphere has T = 0 by definition, R mignt not = 1
        Tranmission = 0
        
        return -n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1)
        
        
    elif R == 1:
        Tranmission = 0
        
        return -n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1)
        
    else:                                           #dielectric sphere ignore absorption
        
        Tranmission = 1 - Reflection_coe
    
    
    
        return -n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * \
        (Reflection_coe * np.cos(2*theta) + 1 - \
         Tranmission**2 * (np.cos(2*(theta - theta_2)) + Reflection_coe * np.cos(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe * np.cos(2*theta_2)) )
   
        
    


def Fy_integrand_Roosen(theta,phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target):
    
    
    """integrand of Fy"""
    
    # theta - sphere theta coordinate, array
    # phi - sphere phi coordinate, array
    # rho_0x - beam x displacement, scalar
    # rho_0y - beam y displacement, scalar
    # rho - sphere diameter, scalar
    # n_0 - refractive index of medium
    # n_s - refractive index of sphere
    # w_0 - beam radius at waist
    # w - beam radius at target
    # z_R - beam rayleigh range
    # P - beam power (Watts)
    # target - string, set as 'reflective' if sphere is reflective, otherwise anything OK and will be ignored
    
    c = 3 * 10**8
    
    mu_0 = 4*np.pi * 10**(-7)
    
    Permittivity = 8.85 * 10**(-12) 
    
    
    P_norm = 0.5 * c * n_0 * Permittivity 
    
    r = beam_r(theta, phi, rho_0x, rho_0y, rho) #represent rou and z by theta
    
    wnew = beam_w(theta, phi, rho, w, w_0, z_R)    #represent w(parameter of laser beam a variable)
    
    theta_2 = get_refracted_angle(theta, n_s, n_0)

    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    I = LG_01_Intensity_sphere_coordinate(theta, phi, rho_0x, rho_0y, rho, wnew, n_0, P)   # beam intensity
    
    if target == "reflective":
        Tranmission = 0
        
        return -n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta))
    
    elif R == 1:
        Tranmission = 0
        
        return -n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * np.sin(phi) *\
        (Reflection_coe *np.sin(2 * theta))
    
    else:
    
        Tranmission = 1 - Reflection_coe
    
    
    
        return - n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * np.sin(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))

def Fx_integrand_Roosen(theta, phi, rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P,target):
    
    
    """integrand of Fx"""
    # theta - sphere theta coordinate, array
    # phi - sphere phi coordinate, array
    # rho_0x - beam x displacement, scalar
    # rho_0y - beam y displacement, scalar
    # rho - sphere diameter, scalar
    # n_0 - refractive index of medium
    # n_s - refractive index of sphere
    # w_0 - beam radius at waist
    # w - beam radius at target
    # z_R - beam rayleigh range
    # P - beam power (Watts)
    # target - string, set as 'reflective' if sphere is reflective, otherwise anything OK and will be ignored
   
    
    c = 3 * 10**8
    
    mu_0 = 4*np.pi * 10**(-7)
    
    Permittivity = 8.85 * 10**(-12) 
    
    
    P_norm = 0.5 * c * n_0 * Permittivity 
    
    r = beam_r(theta, phi, rho_0x, rho_0y, rho) #represent rou and z by theta
    
    wnew = beam_w(theta, phi, rho, w, w_0, z_R)    #represent w(parameter of laser beam a variable)
        
                  
    theta_2 = get_refracted_angle(theta, n_s, n_0)

    
    Reflection_coe = rfcoe_sm(theta,theta_2, n_0, n_s)
    
    I = LG_01_Intensity_sphere_coordinate(theta, phi, rho_0x, rho_0y, rho, wnew, n_0, P)   # beam intensity
    
    if target == "reflective":
        Tranmission = 0
        
        return -n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * np.cos(phi) * \
        (Reflection_coe * np.sin(2 * theta))
    
    elif R == 1:
        Tranmission = 0
        
        return -n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * np.cos(phi) *\
        (Reflection_coe *np.sin(2 * theta))
    
    else:
    
        Tranmission = 1 - Reflection_coe
    
    
    
        return - n_0 * I * rho ** 2 / (2 * mu_0 * c**2 * P_norm) * np.sin(theta) * np.cos(theta) * np.cos(phi) * \
        (Reflection_coe * np.sin(2 * theta) - \
         Tranmission**2 * (np.sin(2*(theta - theta_2)) + Reflection_coe * np.sin(2*theta)) / (1 + Reflection_coe**2 + 2*Reflection_coe*np.cos(2*theta_2)))

###################################
     # my manual integrations
####################################
def F_total_manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, coordinate, grid_size = 400):
    """total forces calculated via manual integration"""
    # theta - sphere theta coordinate, array
    # phi - sphere phi coordinate, array
    # rho_0x - beam x displacement, scalar
    # rho_0y - beam y displacement, scalar
    # rho - sphere diameter, scalar
    # n_0 - refractive index of medium
    # n_s - refractive index of sphere
    # w_0 - beam radius at waist
    # w - beam radius at target
    # z_R - beam rayleigh range
    # P - beam power (Watts)
    # target - string, set as 'reflective' if sphere is reflective, otherwise anything OK and will be ignored
    # grid_size - integer, number of elements along square grid size, e.g. grid_size X grid_size
    # coordinate - str, 'x', 'y', or 'z', axis to return force for
    
    g = 9.8
    
    theta = np.linspace(0, np.pi/2, grid_size)
    
    theta_step = theta[1] - theta[0]
    
    phi = np.linspace(0, 2*np.pi, grid_size)
    phi_step = phi[1] - phi[0]
    
    theta, phi = np.meshgrid(theta, phi)
    
    if coordinate == 'z':
        
        force_function = Fz_integrand_Roosen
        
    elif coordinate == 'x':
        
        force_function = Fx_integrand_Roosen
        
    elif coordinate == 'y':
        
        force_function = Fy_integrand_Roosen
    
    force_grid = force_function(theta, phi, rho_0x, rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    force_total = force_grid.sum() * theta_step * phi_step
    
    
    I = LG_01_Intensity_sphere_coordinate(theta, phi, rho_0x, rho_0y, rho, w, n_0, P)
    
    cartesian_grid = spherical_to_cartesian_coordinates(rho, theta, phi)
    
    output = {'force_total' : force_total,
              'force_grid': force_grid,
              'intensity' : I,
              'theta' : theta,
              'phi' : phi,
              'x' : cartesian_grid['x'],
              'y' : cartesian_grid['y']
              }
    

    return output

####################################
    #dblquad integration
####################################
    
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

#############################################
    #plots
#############################################

def Fz_total_vs_d_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet = []
    for w_e in w: 
        if integration_method == 'integrated':
            
        
            F_1tz = Fz_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
        
            F_znet = F_1tz[0]
            forcenet.append(F_znet)
            
        elif integration_method == 'manual':
            
            print (w_e)
            integration_results = F_total_manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target, coordinate = 'z', grid_size = grid_size)
            
            forcenet.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
            
    return forcenet

def Fy_total_vs_rho0y_plot(rho_0x,rho_0y,rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet_r = []
    for rho_0_e in rho_0y:
        if integration_method == 'integrated':
            
    
            F_1rr = Fy_total(rho_0x, rho_0_e, rho, n_0, n_s, w_0, w, z_R, P, target)
    
            F_znet = F_1rr[0] 
            forcenet_r.append(F_znet)
            
        elif integration_method == 'manual':
            
            print (rho_0_e)
            integration_results = F_total_manual_integration(rho_0x,rho_0_e, rho, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'y', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
        

    return forcenet_r


def Fx_total_vs_rho0x_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    forcenet_r = []
    
    for rho_0_e in rho_0x:
        if integration_method == 'integrated':
            
    
            F_1rr = Fx_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
        
            F_znet = F_1rr[0] 
            forcenet_r.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'x', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    

    return forcenet_r

def Fz_total_vs_rho0x_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet_r = []
    for rho_0_e in rho_0x:
        if integration_method == 'integrated':
            
    
            F_1rr = Fz_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
            F_znet = F_1rr[0] 
            forcenet_r.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'z', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    

    return forcenet_r


def Fx_total_vs_d_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet = []
    
    for w_e in w:
        if integration_method == 'integrated':
            
    
            F_1tz = Fx_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
        
    
            F_znet = F_1tz[0]
            forcenet.append(F_znet)
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target, coordinate = 'x', grid_size = grid_size)
            
            forcenet.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    
    
    return forcenet


def Fy_total_vs_d_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet = []
    
    for w_e in w:
        if integration_method == 'integrated':
            
    
            F_y = Fy_total(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
        
    
            F_znet = F_y[0]
            forcenet.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target, coordinate = 'y', grid_size = grid_size)
            
            forcenet.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    
    return forcenet

###################################################################################################
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



def Fz_total_vs_Rs_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet = []
    
    for rho_e in rho:
        if integration_method == 'integrated':
            
    
            F_1tz = Fz_total(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target)
        
            F_znet = F_1tz[0]
            forcenet.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'z', grid_size = grid_size)
            
            forcenet.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    
    return forcenet


def Fy_total_vs_Rs_plot(rho_0x,rho_0y,rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet_r = []
    
    for rho_e in rho:
        if integration_method == 'integrated':
            F_1rr = Fy_total(rho_0x, rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target)
        
            F_znet = F_1rr[0] 
            forcenet_r.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'y', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
            forcenet_r = []
    
    return forcenet_r


def Fx_total_vs_Rs_plot(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    
    # integration_method - str, 'integrated' or 'manual', use integrated or manual integration method
    # grid_size - integer, grid side length
    
    forcenet_r = []
    for rho_e in rho:
        if integration_method == 'integrated':
            
            F_1rr = Fx_total(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target)
    
            F_znet = F_1rr[0] 
            forcenet_r.append(F_znet)
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'x', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
        
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    
    return forcenet_r  
###################################################################################################################################

##################################################
    #Force gradients plots
##################################################




def Fx_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    """calculate x radial gradient"""
    
    
    
    forcenet_r = []
    
    for rho_0_e in rho_0x:
        
        if integration_method == 'integrated':
            F_x = Fx_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
            F_znet = F_x[0] 
            forcenet_r.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'x', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
            
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    
    grad_Fx = []
    rho_0xlist = []
    
    for i in range( 1, len(forcenet_r) ):
        grad_Fxelement = ( forcenet_r[i] - forcenet_r[i-1] ) / (rho_0x[i] - rho_0x[i-1])
        
        rho_0xelement = (rho_0x[i] + rho_0x[i-1]) / (2 * (w / np.sqrt(2)) )
        
        grad_Fx.append(grad_Fxelement)
        
        rho_0xlist.append(rho_0xelement)
        
    output = {'Fx_grad' : grad_Fx,
              'rho_0x': rho_0xlist,
              }
        
    return output
        
    


def Fy_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    """calculate x radial gradient"""
    
    
    
    forcenet_r = []
    
    for rho_0_e in rho_0y:
        
        if integration_method == 'integrated':
            F_y = Fy_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
            F_znet = F_y[0] 
            forcenet_r.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0x,rho_0_e, rho, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'y', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
            
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
        
    
    
    grad_Fy = []
    rho_0ylist = []
    
    for i in range( 1, len(forcenet_r) ):
        grad_Fyelement = ( forcenet_r[i] - forcenet_r[i-1] ) / (rho_0y[i] - rho_0y[i-1])
        
        rho_0yelement = (rho_0y[i] + rho_0y[i-1]) / (2 * (w / np.sqrt(2)) )
        
        grad_Fy.append(grad_Fyelement)
        
        rho_0ylist.append(rho_0yelement)
        
    output = {'Fy_grad' : grad_Fy,
              'rho_0y': rho_0ylist,
              }
        
    return output

def Fz_total_gradient(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    """calculate x radial gradient"""
    
    
    
    forcenet = []
    
    
    for w_e in w:
        
        if integration_method == 'integrated':
            F_z = Fz_total(rho_0_e,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target)
    
    
    
            F_znet = F_z[0] 
            forcenet.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0x,rho_0y, rho, n_0, n_s, w_0, w_e, z_R, P, target, coordinate = 'z', grid_size = grid_size)
            
            forcenet.append(integration_results['force_total'])
            
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
        

    d = np.sqrt( (z_R ** 2 / w_0 ** 2) * ( w**2 - w_0**2 )) 
    
    grad_Fz = []
    dlist = []
    
    for i in range( 1, len(forcenet) ):
        grad_Fzelement = ( forcenet[i] - forcenet[i-1] ) / (d[i] - d[i-1])
        
        delement = (d[i] + d[i-1]) /2
        
        grad_Fz.append(grad_Fzelement)
        
       
        dlist.append(delement)
        
        
    output = {'Fz_grad' : grad_Fz,
              'rho_0z': dlist,
              }
        
    return output


def Fz_total_gradient_vs_radialoffset(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, integration_method = "integrated", grid_size = 500):
    """calculate x radial gradient"""
    
    
    
    forcenet_r = []
    
    for rho_0xe in rho_0x:
        
        if integration_method == 'integrated':
            F_z = Fz_total(rho_0xe,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target)
    
    
    
            F_znet = F_z[0] 
            forcenet_r.append(F_znet)
            
        elif integration_method == 'manual':
            integration_results = F_total_manual_integration(rho_0xe,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, target, coordinate = 'z', grid_size = grid_size)
            
            forcenet_r.append(integration_results['force_total'])
            
        else:
            print(F"WRONG INTEGRATION METHOD SPECIFIED OF: {integration_method}")
    
    grad_Fz = []
    rho_0xlist = []
    
    for i in range( 1, len(forcenet_r) ):
        grad_Fzelement = ( forcenet_r[i] - forcenet_r[i-1] ) / (rho_0x[i] - rho_0x[i-1])
        
        rho_0xelement = (rho_0x[i] + rho_0x[i-1]) / 2
        
        grad_Fz.append(grad_Fzelement)
        
        rho_0xlist.append(rho_0xelement)
        
    output = {'Fz_grad' : grad_Fz,
              'rho_0x': rho_0xlist,
              }
        
    return output






