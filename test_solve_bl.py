#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 18 11:11:53 2018

@author: anuar
"""

import thesis_tools as tt
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------
# Plot Mistress results 

#file_name = 'BL' # vtk file base name <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#vtk_no = 3 # no. of vtlk files <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#vtk_time_step = 0.1 # time between the vtk files <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
#
#x_grid, y_grid, z_grid = vtk_dim(file_name+'.1.vtk')
#x = np.linspace(0,1,x_grid)
#
#plt.figure(figsize=(10,5))
#plt.subplot(1,2,2)
#
#for vtk in range(vtk_no):  
#    td = (vtk+1)*vtk_time_step
#    sat = extract_2d_vtk(file_name+'.'+str(vtk+1)+'.vtk')   
#    avg_sat = np.mean(sat, axis=0)
#    plt.plot(x/td,avg_sat,'--r', label='td='+str(round(td,2)))
#--------------------------------------------------------------------------------    
    
#------------------------------------------------------------------------------------
# Analytical solution   
    
nw=5 # Corey's exponent (water) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
no=1 # Corey's exponent (oil) <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
mw=1 # water viscosity <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
mo=100 # oil viscosity <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Swc=0 # Connate water <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
Sor=0 # Irreducible oil <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

S = np.linspace(Swc+1/100,1-Sor,100)
Swf, vwf, S_solution, v_solution = tt.solve_bl(nw=nw,no=no,mw=mw,mo=mo,Swc=Swc,Sor=Sor)
#------------------------------------------------------------------------------------

plt.subplot(1,2,1)
plt.plot(S,tt.frac_flow(S,nw=nw,no=no,mw=mw,mo=mo,Swc=Swc,Sor=Sor),'-k')
plt.plot(S,vwf*(S-Swc),'--k')
plt.scatter(Swf,tt.frac_flow(Swf,nw=nw,no=no,mw=mw,mo=mo,Swc=Swc,Sor=Sor),c='k')
plt.scatter(Swc,0,c='k')
plt.xlim(0,1)
plt.ylim(0,1)
plt.xlabel('Sw')
plt.ylabel('Fw')

plt.subplot(1,2,2)
plt.plot(v_solution,S_solution,'-k',label='Analytical')
plt.plot([vwf,vwf], [Swc,Swf], '-k')
plt.scatter(vwf,Swf,c='k')
plt.scatter(vwf,Swc,c='k')
plt.xlim(0,vwf*1.5)
plt.ylim(0,1)
plt.xlabel('v=x/t')
plt.ylabel('Sw')
plt.legend()

plt.savefig('BL.png')