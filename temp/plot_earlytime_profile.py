#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 24 08:26:14 2018

@author: anuar
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

M=20
DT = [0.002, 0.004, 0.006, 0.008, 0.01]
color = ['r','g','b','c','m']
color_dotted = ['--r','--g','--b','--c','--m']
#file_name = earlytime
case_no = 5
#plt.subplots(3,2)
plt.figure(figsize=(10,12))
plt.subplot(321)
for i in range(case_no):
    analytical = 1/((2**(5/2))*np.pi*(M+1)/(M-1)*DT[i])
    fingers = np.loadtxt('earlytime_'+str(i+1)+'_fingers.txt')
    plt.plot(fingers[:,0], fingers[:,1], color[i], label='$D_L$=0, $D_T$='+str(DT[i]))
    plt.plot(fingers[:,0], np.full(len(fingers),analytical), color_dotted[i])
plt.xlabel('T')
plt.ylabel('No. of fingers')
plt.legend()
#plt.show()

plt.subplot(322)
for i in range(case_no):
    analytical = 1/((2**(5/2))*np.pi*(M+1)/(M-1)*DT[i])
    df = pd.read_csv('earlytime_'+str(i+1)+'.conc', header=None, skiprows=1, delim_whitespace=True)
    c = df.values
    mix_len = np.full((len(df)), np.nan)
    for j in range(len(df)):
        mix_len[j] = np.min(np.where(c[j,1:]<0.1))-np.min(np.where(c[j,1:]<0.9))
    plt.plot(c[:,0],mix_len[:]/(np.size(c,axis=1)-1), color[i], label='$D_L$=0, $D_T$='+str(DT[i]))
plt.xlabel('T')
plt.ylabel('Mixing length')
#plt.legend()
#plt.show()

plt.subplot(323)
for i in range(case_no):
    analytical = 1/((2**(5/2))*np.pi*(M+1)/(M-1)*DT[i])
    fingers = np.loadtxt('earlytime_'+str(i+1)+'_fingers.txt')
    plt.plot(fingers[:,0]/DT[i], fingers[:,1]/analytical, color[i], label='$D_L$=0, $D_T$='+str(DT[i]))
    plt.plot(fingers[:,0]/DT[i], np.full(len(fingers),analytical)/analytical, color_dotted[i])
    
plt.xlabel('T/$D_T$')
plt.ylabel('No. of fingers/Analytical fingers')
#plt.legend()
#plt.show()

plt.subplot(324)
for i in range(case_no):
    analytical = 1/((2**(5/2))*np.pi*(M+1)/(M-1)*DT[i])
    df = pd.read_csv('earlytime_'+str(i+1)+'.conc', header=None, skiprows=1, delim_whitespace=True)
    c = df.values
    mix_len = np.full((len(df)), np.nan)
    for j in range(len(df)):
        mix_len[j] = np.min(np.where(c[j,1:]<0.1))-np.min(np.where(c[j,1:]<0.9))
    plt.plot(c[:,0]/DT[i],mix_len[:]/(np.size(c,axis=1)-1)/analytical, color[i], label='$D_L$=0, $D_T$='+str(DT[i]))
plt.xlabel('T/$D_T$')
plt.ylabel('Mixing length/Analytical fingers')
#plt.legend()
#plt.show()

plt.subplot(325)
for i in range(case_no):
    analytical = 1/((2**(5/2))*np.pi*(M+1)/(M-1)*DT[i])
    fingers = np.loadtxt('earlytime_'+str(i+1)+'_fingers.txt')
    plt.scatter(DT[i], analytical, c='r', marker='^')
    plt.scatter(DT[i], np.nanmax(fingers[:,1]), c='k', marker='x')
#    plt.plot(fingers[:,0]/DT[i], np.full(len(fingers),analytical)/analytical, color_dotted[i])
    
plt.xlabel('$D_T$')
plt.ylabel('No. of fingers')
plt.xlim(0,0.012)
plt.legend(['Analytical','Mistress'])
#plt.show()
plt.suptitle('M=20, $D_L=0$, grid=400x400', fontsize=15)
plt.savefig('early_time.png')
