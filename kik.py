# -*- coding: utf-8 -*-
"""
Created on Wed Oct  7 20:49:34 2020

@author: erlin
"""
import hyperspy.api as hs
import matplotlib.pyplot as plt
import numpy as np
import kikuchipy as kp 

filepath = 'filename' #Filepath to raw data
s = kp.load(filepath, lazy=False) #Loading raw EBSD patterns

s.remove_static_background(operation="subtract", relative=True) #Removing static background from patterns

s.remove_dynamic_background() #Removing dynamic background from patterns

w_gauss = kp.filters.Window(window="gaussian",std=1) #Creating gaussian filter 
s.average_neighbour_patterns(window=w_gauss) #Neighbour pattern averaging the patterns, using the gaussian filter above

s.save("output_filepath") #Saving the averaged patterns

   