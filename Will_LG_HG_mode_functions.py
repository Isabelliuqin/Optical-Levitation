# -*- coding: utf-8 -*-
"""
Created on Wed Jul  8 17:10:45 2020

@author: Will_work

functions to create Laguerre-Gaussian and Hermite-Gaussian modes

"""

import numpy as np
import warnings
from scipy.special import eval_genlaguerre as LGpl
from scipy.special import eval_hermite as HGn

def Laguerre_Gaussian(p, l, x, y, wavelength = None, q = None):
    # returns Laguerre Gaussian mode for radial index p, azimuthal index l at array of points x, y
    # if wavelength and complex q parameter not specified x and y assumed normalised to beam waist radius
      
    # if selecting wavelength or q the other must be set (^ is XOR)
    if (wavelength is None) ^ (q is None):
        warnings.warn('both wavelength and q must be set or unset')
        return []
    
    if wavelength is None:
        w = 1.
        psi = 0
        wavelength = np.pi
        q = get_q(w,wavelength)
    else:
        w = np.sqrt(-1/np.imag(1/q) * wavelength /np.pi)      # from q definition
        psi = np.arctan(np.real(q)/np.imag(q))       # from q(z) = z + zR*i
           
    k = 2*np.pi/wavelength    
        
    phi = np.arctan2(y,x)       #must use arctan2, arctan wraps around prematurely (doesn't preserve quadrant)
    rsq = x**2 + y**2

    factorial = np.math.factorial
    
    return (-1)**p * np.sqrt( 2*factorial(p) / ( np.pi*w**2 *factorial(p+abs(l)) ) ) * (np.sqrt(2*rsq)/w)**abs(l) * LGpl( p, abs(l), 2*rsq/(w**2)) * np.exp( -1j*k*rsq/(2*q) - 1j*(2*p+abs(l)+1)*psi - 1j*l*phi )

def Hermite_Gaussian2D(n, m, x, y, wavelength = None , qx = None, qy = None):
    # returns Hermite Gaussian mode for horizontal index n, vertical index m at array of points x, y
    # if wavelength and complex q parameter not specified x and y assumed normalised to beam waist radius

    if qy is None:
        qy = qx
        
    return Hermite_Gaussian1D(n, x, wavelength = wavelength, q = qx) * Hermite_Gaussian1D(m, y, wavelength = wavelength, q = qy)
    
def Hermite_Gaussian1D(n, x, wavelength = None, q = None):
    
    # if selecting wavelength or q the other must be set (^ is XOR)
    if (wavelength is None) ^ (q is None):
        warnings.warn('both wavelength and q must be set or unset')
        return []    
    
    if wavelength is None:
        k = 0
        w = 1
        psi = 0
    else:        
        k = 2*np.pi/wavelength
        w = np.sqrt(-1/np.imag(1/q) * wavelength /np.pi)      # from q definition
        psi = np.arctan(np.real(q)/np.imag(q))       # from q(z) = z + zR*i
    
    factorial = np.math.factorial
    
    return (2/np.pi)**0.25 * np.sqrt( np.exp(-1j*(2*n+1)*psi) / (2**n *factorial(n)*w) ) * np.exp(-1j*k* x**2 / (2*q) ) * HGn(n, np.sqrt(2)*x/w)

def get_q(w, wavelength, R = np.inf, M2 = 1, n = 1):
    return 1/( 1/R - 1j*wavelength*M2/(np.pi*n*w**2) )

def get_zR(q):
    return np.imag(q)