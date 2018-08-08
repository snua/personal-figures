#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  6 11:21:47 2018

@author: anuar
"""

import thesis_tools as tt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

C = tt.extract_2d_vtk('test.1.vtk')
plt.imshow(C,cmap='jet')
plt.show()
from skimage import measure


level = 0.5
# Find contours at a constant value of 0.8
contours_scipy = measure.find_contours(C, level)
#plt.plot(contours[0][:,1],contours[0][:,0])
for n, contour_scipy in enumerate(contours_scipy):
    plt.plot(contour_scipy[:,1],contour_scipy[:,0],'k')


C_test = np.full((len(C), len(C)), np.nan)
C_vertex = np.full((len(C)+1, len(C)+1), np.nan)

points = []

for i in range(len(C_test)):
    for j in range(len(C_test)):               
        if C[i,j]<level:
            C_test[i,j] = 0
        else:
            C_test[i,j] = 1

for i in range(len(C_test)-1):
    for j in range(len(C_test)-1): 

        if C_test[i,j]==0 and C_test[i+1,j]==0 and C_test[i,j+1]==0 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 0
        elif C_test[i,j]==0 and C_test[i+1,j]==1 and C_test[i,j+1]==0 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 1
            points.append([i+0.5,j,i+1,j+0.5])
        elif C_test[i,j]==0 and C_test[i+1,j]==0 and C_test[i,j+1]==0 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 2
            points.append([i+1,j+0.5,i+0.5,j+1])
        elif C_test[i,j]==0 and C_test[i+1,j]==1 and C_test[i,j+1]==0 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 3
            points.append([i+0.5,j,i+0.5,j+1])

        elif C_test[i,j]==0 and C_test[i+1,j]==0 and C_test[i,j+1]==1 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 4
            points.append([i,j+0.5,i+0.5,j+1])
        elif C_test[i,j]==0 and C_test[i+1,j]==1 and C_test[i,j+1]==1 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 5
            points.append([i,j+0.5,i+0.5,j])
            points.append([i+1,j+0.5,i+0.5,j+1])
        elif C_test[i,j]==0 and C_test[i+1,j]==0 and C_test[i,j+1]==1 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 6
            points.append([i,j+0.5,i+1,j+0.5])
        elif C_test[i,j]==0 and C_test[i+1,j]==1 and C_test[i,j+1]==1 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 7
            points.append([i,j+0.5,i+0.5,j])
            
        if C_test[i,j]==1 and C_test[i+1,j]==0 and C_test[i,j+1]==0 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 8
            #points.append([i,j+0.5,i+0.5,j])
            points.append([i+0.5,j,i,j+0.5])
        elif C_test[i,j]==1 and C_test[i+1,j]==1 and C_test[i,j+1]==0 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 9
            points.append([i,j+0.5,i+1,j+0.5])
        elif C_test[i,j]==1 and C_test[i+1,j]==0 and C_test[i,j+1]==0 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 10
            points.append([i+0.5,j,i+1,j+0.5])
            points.append([i,j+0.5,i+0.5,j+1])
        elif C_test[i,j]==1 and C_test[i+1,j]==1 and C_test[i,j+1]==0 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 11
            #points.append([i,j+0.5,i+0.5,j+1])
            points.append([i+0.5,j+1,i,j+0.5])

        elif C_test[i,j]==1 and C_test[i+1,j]==0 and C_test[i,j+1]==1 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 12
            points.append([i+0.5,j,i+0.5,j+1])
            #points.append([i+0.5,j+1,i+0.5,j])
        elif C_test[i,j]==1 and C_test[i+1,j]==1 and C_test[i,j+1]==1 and C_test[i+1,j+1]==0:
            C_vertex[i+1,j+1] = 13
            points.append([i+1,j+0.5,i+0.5,j+1])
        elif C_test[i,j]==1 and C_test[i+1,j]==0 and C_test[i,j+1]==1 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 14
            points.append([i+0.5,j,i+1,j+0.5])
        elif C_test[i,j]==1 and C_test[i+1,j]==1 and C_test[i,j+1]==1 and C_test[i+1,j+1]==1:
            C_vertex[i+1,j+1] = 15
            #points.append([i+,j+,i+,j+])



#        C_vertex[i+1,j+1] = (C[i,j] + C[i+1,j] + C[i,j+1] + C[i+1,j+1])/4            
points = np.asarray(points)
plt.scatter(points[:,1], points[:,0],c='r')
plt.scatter(points[:,3], points[:,2],marker='x')
plt.imshow(C_test)


contour = np.copy(points)

def swap_points(points):
    
    temp_points = np.copy(points)
    points[0:2] = temp_points[2:4]
    points[2:4] = temp_points[0:2]
    
    return points

sorted_points = np.full((len(contour),4), np.nan)

# first the first point
min_point_a = contour[np.argmin(contour[:,0]),0] 
min_point_b = contour[np.argmin(contour[:,2]),2] 

if min_point_a<min_point_b:
    min_point_loc = np.argmin(contour[:,0]) 
else:
    min_point_loc = np.argmin(contour[:,2])

min_point = np.copy(contour[min_point_loc,:])
if min_point[0]>min_point[2]:
    min_point = swap_points(min_point)

sorted_points[0,:] = min_point

ints = np.arange(len(contour))

# remove the entry already used
#ints = np.delete(ints, min_point_loc)
ints = np.delete(ints,np.argwhere(ints==min_point_loc))

for k in range(1,len(contour)):    
    
    for i,j in enumerate(ints):
           if all(contour[j,0:2]==sorted_points[k-1,2:4]) or all(contour[j,2:4]==sorted_points[k-1,2:4]):
               next_point_loc = j
           #else:
               #print('Cant find match')
    
    next_point = np.copy(contour[next_point_loc,:])
    
    if any(next_point[0:2]!=sorted_points[k-1,2:4]):
        next_point = swap_points(next_point)
    
    sorted_points[k,:] = next_point
    
    # remove the entry already used
    ints = np.delete(ints,np.argwhere(ints==next_point_loc))

plt.plot(sorted_points[:,1],sorted_points[:,0],'g')