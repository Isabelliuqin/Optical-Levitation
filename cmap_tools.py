# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 12:06:00 2020

@author: wrk10_home

functions to manipulate colormap objects and data

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import warnings

def cmap_add_alpha(cmap_in, alpha_upper, alpha_lower = 0):
    # adds a linear alpha fade to the lower end of a colormap
    # cmap_in  = string name or array of values to modify
    # alpha_upper = upper start point of fade
    # alpha_lower = lower end point of fade
    # returns - (listed)colormap object, can be used in imsave etc.
        
    if type(cmap_in) is str:
        cmap = plt.get_cmap(cmap_in)
        cmap = cmap(np.linspace(0,1,num=256))
        cmap_name = cmap_in
    else:
        cmap_name = ''


    alpha_upper_idx = int(round(cmap.shape[0]*alpha_upper))
    alpha_lower_idx = int(round(cmap.shape[0]*alpha_lower))

    cmap[0:alpha_lower_idx,3] = 0
    cmap[alpha_lower_idx:alpha_upper_idx,3] = np.linspace(0,1,num=alpha_upper_idx-alpha_lower_idx)
    
    cmap = ListedColormap(cmap,name=[cmap_name+'_alpha_modified'])
    return cmap

def cmap_alpha_to_RGB(cmap, background=(1,1,1)):
    
    try:
        cmap_name = cmap.name[0]
    except:
        cmap_name = cmap
        pass
    
    cmap = plt.get_cmap(cmap)
    
    # get RGBa values of cmap
    cmap_array = cmap(np.linspace(0,1,256))
    
    # extract and tile alpha into RGB columns
    alpha_array = np.tile(np.atleast_2d(cmap_array[:,3]).T,(1,3))
    
    # store background colours
    background_array = np.tile(np.array([background]),(256,1))
    
    cmap_array[:,0:3] = cmap_array[:,0:3]*alpha_array + background_array*(1-alpha_array)
    cmap_array[:,3] = 1
    
    cmap = ListedColormap(cmap_array,name=[cmap_name+'alpha_modified_RGB'])
    return cmap

def convert_to_colormap_object(cmap):
    # cmap is RGB(A) array of colors, string name, or already listed_colormap object
    #
    
    if type(cmap) is ListedColormap:
        return cmap
    if type(cmap) is str:
        try:
            cmap = get_cmap_from_string(cmap)
            return ListedColormap(cmap)
        except ValueError:
            print(f'colormap name {cmap} not found')
            return
    else:
        try:
            return ListedColormap(cmap)
        except:
            print('colormap could not be converted')
    
    return
    
    

def export_colorbar(cmap, export_name = None, save_extension = 'png', save_folder = '', npixels_long = 256, npixels_wide = 10):
    # export an image of the cmap given
    # cmap = str or Colormap object to export
    # export_name = name of export file, if not given defaults to colorbar_(cmap name if available)
    # save_extension = file type to export to
    # save_folder = path of save folder (include / at end)
    # npixels_long = number of pixels along colormap gradient
    # npixels_wide = number of pixels wide
    
    if export_name is None:
        save_name = 'colorbar'
        try:
            save_name = save_name + '_' + cmap.name[0]
        except:
            pass
    else:
        save_name = export_name
    
    save_name = str(save_folder) + save_name + '.' + save_extension
    
    gradient = np.arange(0,npixels_long,1)
    gradient = np.tile(gradient,(npixels_wide,1))
    
    # pixel color is at mid point value
    gradient = gradient + 0.5
    
    cmap = convert_to_colormap_object(cmap)
    
    plt.imsave(save_name, gradient, cmap = cmap)
    return

def export_cmap_LUT(cmap,export_name = None,normalised = False,samples=256,number_format='%.18e',delimiter='\t'):
    # export LUT file of given colormap
    # cmap = str or colormap object ro export
    # export_name = name of exported LUT file
    # normalised = if True, values 0 -> 1, else 0 -> 255 integers
    
    if export_name is None:
        try:
            export_name = cmap.name[0] + '.LUT'
        except:
            export_name = 'cmap.LUT'
            pass
    else:
        export_name = export_name + '.LUT'
        
    array = cmap(np.linspace(0,1,num=samples))
    
    if not normalised:
        array = array*255
        array = array.round()
        array = array.astype('int')
        
        number_format = '%1d'
    
    np.savetxt(export_name,array,fmt=number_format,delimiter=delimiter)
    return


def convert_rgba_to_intensity(array, rgba_cmap):
    # convert rgba values to normalised intensity corresponding to given colormap
    # array - np array of rgba values, in range 0 -> 255
    # rgba_cmap - [N,4] numpy array of RGBa intensity values with value range 0 -> 255
    # https://stackoverflow.com/questions/3720840/how-to-reverse-a-color-map-image-to-scalar-values
    import scipy.cluster.vq as scv
    
    array = array / 255
    
    gradient = rgba_cmap

    # Reshape arr to something like (240*240, 4), all the 4-tuples in a long list...
    arr2 = array.reshape((array.shape[0]*array.shape[1], array.shape[2]))

    # Use vector quantization to shift the values in arr2 to the nearest point in
    # the code book (gradient).
    code, dist = scv.vq(arr2, gradient)

    # code is an array of length arr2 (240*240), holding the code book index for
    # each observation. (arr2 are the "observations".)
    # Scale the values so they are from 0 to 1.
    values = code.astype('float') / (gradient.shape[0] - 1)  # maximum index is (shape[0] - 1) due to 0 indexing in python

    # Reshape values back to (240,240)
    values = values.reshape(array.shape[0], array.shape[1])
    
    return values

def get_cmap_from_string(name):
    # search colormap databases for colormap name
    # name - string of cmap name
    
    cmap_found = False
    
    try:
        cmap = get_matplotlib_cmap(name)
        cmap_found = True
    except ValueError:
        pass
    
    try:
        cmap = get_cinogy_cmap(name)
        cmap_found = True
    except ValueError:
        pass
    
    try:
        cmap = get_matlab_cmap(name)
        cmap_found = True
    except ValueError:
        pass
    
    try:
        cmap = get_colorcet_cmap(name)
        cmap_found = True
    except ImportError:
        warnings.warn(f"colorcet module not installed, was not searched for colormap '{name}'")
        pass
    except ValueError:
        pass
    
    if cmap_found:
        return cmap
    else:
        raise ValueError(f'colormap of name {name} not found')
    
def get_cinogy_cmap(name):
    # returns LUT from Cinogy Ray-Ci in 0 -> 1 scale
    # if cmap not found returns False
    
    filename = 'LUT_cinogy/' + name + '.LUT'
    
    try:
        cmap = np.loadtxt(filename)
    except IOError:
        raise ValueError(f"cmap '{name}' not found")
        return
        
    
    # Cinogy LUT in 0 -> 255 scale
    
    cmap_out = np.ones((cmap.shape[0], 4))
    cmap_out[:, 0:3] = cmap[:, 0:3] / 255
    
    return cmap_out

def get_matlab_cmap(name):
    # returns LUT from MATLAB in 0 -> 1 scale
    
    filename = 'LUT_MATLAB/' + name + '.LUT'
    
    try:
        cmap = np.loadtxt(filename)
    except IOError:
        raise ValueError(f"cmap '{name}' not found")
        return
    # MATLAB LUT in 0 -> 1 scale
    
    cmap_out = cmap
    
    return cmap_out

def get_matplotlib_cmap(name):
    # returns Matplotlib cmap array in 0 -> 1 scale
    # if not found returns False
    
    import matplotlib.pyplot as plt
    
    try:
        cmap = getattr(plt.cm, name)
    except AttributeError:
        # cmap name not found
        raise ValueError(f"cmap '{name}' not found")
        return
    
    return cmap(np.arange(cmap.N))

def get_colorcet_cmap(name):
    # returns colorcet_cmap array in 0 -> 1 scale
    # if not found returns False
    
    try:
        import colorcet as cc
    except ImportError:
        raise ImportError(f"colorcet module not installed")
        return
    
    try:
        cmap = getattr(cc.cm, name)
    except KeyError:
        # cmap name not found
        raise ValueError(f"cmap '{name}' not found")
        return
    
    return cmap(np.arange(cmap.N))

def get_linear_cmap(color_1, color_2, n_steps = 255):
    # returns linear gradient colormap between color_1 and color_2 RGBa values
    # color_1 - starting RGBa value, list or array
    # color_2 - stopping RGBa value, list or array
    # n_steps - number of values in gradient
    
    color_array = np.ones([n_steps, 4])
    
    for i in range(4):
        color_array[:,i] = np.linspace(color_1[i], color_2[i], num = n_steps)
    
    return ListedColormap(color_array)

    

if __name__ == '__main__':
    name = 'fire'
    
    
    print(convert_to_colormap_object(name))
        