# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 15:31:40 2020

@author: liuqi
"""


#################################################
#axial stiffness at particular parameter functions
#################################################

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
import Will_Module_Trialplot_Q as WMTQ

def Fz_stiffness_at_parameter(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, resolution, target, integration_method, grid_size = 500):
    
    chosen_resolution = resolution
    w_minus_2 = w - 2 * resolution
    w_minus_1 = w - 1 * resolution
    
    w_plus_1 = w + 1 * resolution
    w_plus_2 = w + 2 * resolution
    
    w_input = np.asarray([w_minus_2, w_minus_1, w, w_plus_1, w_plus_2])
    
    Fz_grad = WMTQ.Fz_total_gradient(rho_0x, rho_0y, rho, n_0, n_s, w_0, w_input, z_R, P, target, integration_method, grid_size)
    Fz_grad_output = np.asarray(Fz_grad['Fz_grad'])
    
    
    return Fz_grad_output


def Fz_stiffness_vs_rho_plots(rho_0x,rho_0y, rho, n_0, n_s, w_0, w, z_R, P, resolution, target, integration_method, grid_size = 500):
    
    Fz_grad = []
    
    for rho_e in rho:
        Fz_grad_element = np.mean(Fz_stiffness_at_parameter(rho_0x,rho_0y, rho_e, n_0, n_s, w_0, w, z_R, P, resolution, target, integration_method, grid_size))
        
        Fz_grad.append(Fz_grad_element)
        
    return Fz_grad