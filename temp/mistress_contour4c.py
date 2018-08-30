#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  8 08:40:02 2018

@author: anuar
"""

#import thesis_tools as tt
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def swap_points(points):
    
    temp_points = np.copy(points)
    points[0:2] = temp_points[2:4]
    points[2:4] = temp_points[0:2]
    
    return points

def import_contour(file_name):
    
    df = pd.read_csv(file_name, delim_whitespace=True)
    time_array = df['Time'].unique()
    contour_level = df['ContourValue'].unique()
    
    return df, time_array, contour_level

def extract_contour_points(df, time, level):

    points = df.loc[df['Time'] == time]
    points = points.loc[points['ContourValue'] == level]

    # the contour format from mistress: Time ContourLevel X1 Y1 X2 Y2
    contour = np.copy(points.iloc[:,2:6].values)
    
    return contour
    
def sort_contour_points(contour):
    
    multi_line = [0]
    
    sorted_points = np.full((len(contour),4), np.nan)
    
    # find the first point i.e. the point with minimum y
    min_point_a = contour[np.argmin(contour[:,1]),1] 
    min_point_b = contour[np.argmin(contour[:,3]),3] 
    
    if min_point_a<min_point_b:
        min_point_loc = np.argmin(contour[:,1]) 
    else:
        min_point_loc = np.argmin(contour[:,3])
    
    min_point = np.copy(contour[min_point_loc,:])
    if min_point[1]>min_point[3]:
        min_point = swap_points(min_point)
    
    sorted_points[0,:] = min_point
    
    # the index of all the points
    ints = np.arange(len(contour))
    
    # remove the entry already used
    ints = np.delete(ints,np.argwhere(ints==min_point_loc))
    
    # now process the rest of the points
    for k in range(1,len(contour)):    
        
        for i,j in enumerate(ints):
               if all(contour[j,0:2]==sorted_points[k-1,2:4]) or \
               all(contour[j,2:4]==sorted_points[k-1,2:4]):
                   break           

        next_point_loc = j
        next_point = np.copy(contour[next_point_loc,:])
                   
        # make sure the order of the points are correct
        if any(next_point[0:2]!=sorted_points[k-1,2:4]):
            next_point = swap_points(next_point)
    
        sorted_points[k,:] = next_point
        
        # remove the entry already used
        ints = np.delete(ints,np.argwhere(ints==next_point_loc))
        
        #check where the lines are disconnected, by marking where the line breaks
        if any(sorted_points[k,0:2]!=sorted_points[k-1,2:4]):
            multi_line = np.append(multi_line,k)
    #also add the last point to the array
    multi_line = np.append(multi_line,len(contour))
    
    # find the location with the most points
    longest_line_loc  = np.argmax(np.diff(multi_line)) + 1
    longest_line = sorted_points[multi_line[longest_line_loc-1]:multi_line[longest_line_loc],:]
    
    return sorted_points, multi_line, longest_line

file_name = 'earlytime_5'
df, time_array, contour_level = import_contour(file_name+'.cntr')

store = np.full((1000,2,3), np.nan)
store_fingers = np.full((len(time_array),2), np.nan)

multiplier = 5
for Time in range(int(len(time_array)/multiplier)):
    for c_lvl in range(1,2):
        
        x_dist = np.full((500,2), np.nan)
        x_dist2 = np.full((500,2), np.nan)

        time = time_array[Time*multiplier]
        contour = extract_contour_points(df, time, contour_level[c_lvl])
                   
        sorted_points, multi_line, longest_line = sort_contour_points(contour)
        
#        plt.scatter(contour[:,0], contour[:,1],marker='.',c='g')
        plt.plot(longest_line[:,0],longest_line[:,1],'r')
        
        #---------------------------------------------------------------------
        # find where the change of sign occur
        a = longest_line[:,0] - np.roll(longest_line[:,0], 1) 
        asign = np.sign(a)
        
        for i in range(len(asign)):
            if asign[i]==0:
                asign[i]=asign[i-1]
        
        signchange = ((np.roll(asign, 1) - asign) != 0).astype(int)
        
        arg_signchange = np.argwhere(signchange==1) 
        nfingers = np.sum(signchange)/2 #ignore first point
        
        for s, enum in enumerate(arg_signchange-1):
#            plt.scatter(longest_line[enum,0], longest_line[enum,1])
            distance = [longest_line[enum,0], longest_line[enum,1]]
            x_dist[s,:] = distance
        
        for s, enum in enumerate(arg_signchange-2):
            if abs(x_dist[s+1,0] - x_dist[s,0])>2:
                x_dist2[s,:] = x_dist[s+1,:]
                plt.scatter(longest_line[enum,0], longest_line[enum,1])
        
        nfingers = np.count_nonzero(~np.isnan(x_dist2[:,0]))/2
#        plt.show()
    #-------------------------------------------------------------------
         #only count the fingers if contour's peak-to-peak is significant enough
#        if np.max(contour[:,0])-np.min(contour[:,0])<3:
#            store[Time,0,c_lvl], store[Time,1,c_lvl] = round(time,5), nfingers
#        else:
#            store[Time,0,c_lvl], store[Time,1,c_lvl] = round(time,5), nfingers
        
        store[Time,0,c_lvl], store[Time,1,c_lvl] = round(time,5), nfingers
        store_fingers[Time,:] = [round(time,5), nfingers]

plt.figure()
plt.plot(store[:,0,0], store[:,1,0],'-^r',label='C=0.1')
plt.plot(store[:,0,1], store[:,1,1],'-^b',label='C=0.5')
plt.plot(store[:,0,2], store[:,1,2],'-^g',label='C=0.9')
plt.legend()
plt.xlabel('Time')
plt.ylabel('No. of fingers')
plt.savefig(file_name+'.png')

np.savetxt(file_name+'_fingers.txt',store_fingers,fmt="%2.4f %2.2f")