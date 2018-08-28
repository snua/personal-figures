#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 10 16:23:54 2018

@author: anuar
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def generate_perm(ngridx, ngridy, set_seed=587329, dist_type='uniform', mu=0, sigma=1):

    np.random.seed(seed=set_seed)   
    
    if dist_type=='uniform':
        perm = np.random.rand(ngridx,ngridy)
    elif dist_type=='normal':   
        perm = sigma * np.random.randn(ngridx,ngridy) + mu
    
    return perm

def save_perm(perm, prm_file='test.prm'):
    
    ngridx, ngridy = np.shape(perm)
    total_grid = ngridx*ngridy
    perm_mistress_format = np.zeros((total_grid,6))
    
    #generate the grid numbering in Mistress format:
    # KBLX, KBLY, KTRX, KTRY, AKXVAL, AKYVAL
    for j in range(ngridy):
        for i in range(ngridx):
    
            n = j*ngridx + i
            
            perm_mistress_format[n,0] = i+1
            perm_mistress_format[n,1] = j+1
            perm_mistress_format[n,2] = i+1
            perm_mistress_format[n,3] = j+1
            
    # now write the permeability in x-direction
    perm_mistress_format[:,4] = np.reshape(perm,total_grid)
    # assume same permeability in y-direction
    perm_mistress_format[:,5] = perm_mistress_format[:,4] 
        
    np.savetxt(prm_file, perm_mistress_format,\
               fmt=['%-4d', '%-4d', '%-4d', '%-4d', '%-1.3f', '%-1.3f'],\
               delimiter=' ')

def scale_perm(perm, refine_level=2):
    
    ngridx, ngridy = np.shape(perm)   
    scaled_perm = np.full((ngridx*refine_level, ngridy*refine_level), np.nan)
    
    for i in range(ngridx):        
        for j in range(ngridy):
            for m in range(refine_level):
                scaled_perm[i*refine_level,j*refine_level+m]=perm[i,j]
    
    for i in range(ngridx):       
        for m in range(refine_level):
            scaled_perm[i*refine_level+m,:]=scaled_perm[i*refine_level,:] 
            
    return scaled_perm

def save_conc(conc, conc_file='test.conc'):
    
    ngridx, ngridy = np.shape(conc)
    total_grid = ngridx*ngridy

    conc_mistress_format  = pd.DataFrame(np.nan, index=np.arange(total_grid), columns=np.arange(6))
    conc = conc.flatten()
    for j in range(ngridy):
        for i in range(ngridx):
    
            n = j*ngridx + i
            
            conc_mistress_format.iloc[n,0] = 'MODCINIT'
            conc_mistress_format.iloc[n,1] = str(int(i+1))
            conc_mistress_format.iloc[n,2] = str(int(i+1))
            conc_mistress_format.iloc[n,3] = str(int(j+1))
            conc_mistress_format.iloc[n,4] = str(int(j+1))
            conc_mistress_format.iloc[n,5] = np.round(conc[n],3)
    
    conc_mistress_format.to_csv(conc_file,sep=' ', index=False, header=False,float_format="%.3f")


def diagonal_grid(n_grid, prm_file, base_perm=1, low_perm=0.00001):
    # create diagonal grid quarter by quarter
    quarter_grid_1 = np.ones((n_grid , n_grid))*base_perm
    quarter_grid_1 *= np.tri(*quarter_grid_1.shape)
    quarter_grid_2 = np.flipud(quarter_grid_1)
    quarter_grid_1_2 = np.concatenate((quarter_grid_1, quarter_grid_2), axis=0)   
    quarter_grid_3_4 = np.fliplr(quarter_grid_1_2)
    # merge all the quarters
    merged_quarters = np.concatenate((quarter_grid_3_4, quarter_grid_1_2), axis=1)
    # delete one duplicate row/column in the middle 
    diag_grid = np.delete(merged_quarters,n_grid,0)
    diag_grid = np.delete(diag_grid,n_grid,1)
    # take note the active cells
    active_grid = int(np.sum(diag_grid))
    #substitute zeros with low perm value
    diag_grid[diag_grid==0] = low_perm
        
    # write in Mistress format
    ngridx = len(diag_grid)
    ngridy = len(diag_grid)
    total_grid = ngridx*ngridy
    perm_mistress_format = np.zeros((total_grid,6))
    
    #generate the grid numbering in Mistress format:
    # KBLX, KBLY, KTRX, KTRY, AKXVAL, AKYVAL
    for j in range(ngridy):
        for i in range(ngridx):
    
            n = j*ngridx + i
            
            perm_mistress_format[n,0] = i+1
            perm_mistress_format[n,1] = j+1
            perm_mistress_format[n,2] = i+1
            perm_mistress_format[n,3] = j+1
            
    # now write the permeability in x-direction
    perm_mistress_format[:,4] = np.reshape(diag_grid,total_grid)
    # assume same permeability in y-direction
    perm_mistress_format[:,5] = perm_mistress_format[:,4] 

    np.savetxt(prm_file, perm_mistress_format,\
           fmt=['%-4d', '%-4d', '%-4d', '%-4d', '%-1.5f', '%-1.5f'],\
           delimiter=' ')

    return diag_grid, ngridx, ngridy, active_grid
#Input-----------------------------------
n_grid = 6 # Diagonal grid blocks at the border
prm_file = 'test_diag.prm' # perm file name
#----------------------------------------
diag_grid, ngridx, ngridy, active_grid = diagonal_grid(n_grid, prm_file)

plt.imshow(diag_grid)
print('Diagonal grid blocks at the border='+str(n_grid))
print('Grid blocks to be specified in Mistress='+str(ngridx)+'x'+str(ngridy))
print('Active grid blocks='+str(active_grid))
print('Active pore volume='+str(active_grid)+'/('+str(ngridx)+'x'+str(ngridy)+')='+str(active_grid/(ngridx*ngridy)))
plt.title('Grid')
plt.show()

#--------------------------------------------------------------------------
#---------------------------------
# specify the following parameters
random_box_size = 3 # define how big is the concentration plume at the injection well
set_seed=38271 # specify different seed number to get different realisation
dist_type='uniform' #specify either 'uniform' or 'normal'
conc_file = 'test' # perm file base name
#----------------------------------

#initialise all conc to be zero
conc_1 = np.zeros((ngridx,ngridy)) 
# create random data from which the random concentration is drawn from
temp_conc_1 = generate_perm(ngridx, ngridy, set_seed=set_seed, dist_type=dist_type)  

for i in range(random_box_size):
    for j in range(random_box_size):
        conc_1[i+n_grid-j-1,j+i] = temp_conc_1[i+n_grid-j-1,j+i] 
        conc_1[i+n_grid-j-1,j+i+1] = temp_conc_1[i+n_grid-j-1,j+i+1] 

# scale by 0.1 and shift by +0.9, so that so that distribution now is between 0.9 to 1.0
#conc_1 = conc_1*0.1+0.9
## manually set some blocks at the circumference to be low concentration 
#conc_1[0,4]=0.1
#conc_1[4,0]=0.1
#conc_1[2,3]=0.1
#conc_1[3,2]=0.1
#
##manually create the quarter circle
#conc_1[1:4,4]=0
#conc_1[4,1:4]=0
#conc_1[3,3]=0
#conc_1[4,4]=0

# view it
plt.figure()
plt.imshow(np.flipud(conc_1))
plt.imshow((conc_1))
plt.title('Base concentration')

# save it
save_conc(conc_1, conc_file=conc_file+'1.conc')
#-----------------------------------------------------------------------