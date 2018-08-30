# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 14:07:57 2018

@author: Anuar
"""
import thesis_tools as tt
import matplotlib.pyplot as plt

test = tt.extract_conc('unit_test.conc')
plt.plot(test[:,0], test[:,1])