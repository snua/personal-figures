# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 17:14:36 2018

@author: Anuar
"""
import numpy as np
import matplotlib.pyplot as plt

def save_image(data=np.arange(100).reshape((10, 10)), cm='jet', fn='out1.png'):
   # from https://fengl.org/2014/07/09/matplotlib-savefig-without-borderframe/
    """
    Saves figures without border. Change the height and width if required.
    """
    sizes = np.shape(data)
    height = float(sizes[0])*10
    width = float(sizes[1])*10
     
    fig = plt.figure()
    fig.set_size_inches(width/height, 1, forward=False)
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
     
    ax.imshow(data, cmap=cm)
    plt.savefig(fn, dpi = height) 
    plt.close()
    
save_image()