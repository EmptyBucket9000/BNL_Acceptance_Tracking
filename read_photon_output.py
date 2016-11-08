# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 11:21:46 2016

@author: Eric Schmidt
"""

import numpy as np
import csv
import os
import glob
    
def main(): 
#==============================================================================
# Particles
#==============================================================================
    
    particle_file = glob.glob("%s/../Output/particle_matrix.csv"%(os.getcwd()))
    particle_file = particle_file[0]
    
    with open(particle_file, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
        stuff = list(reader)
        
        '''
        Photon #                        0
        Steps                           1
        Kill Event                      2
        Starting Global x-Position      3
        Starting Global y-Position      4
        Starting Global z-Position      5
        Ending Calorimeter x (mm)       6
        Ending Calorimeter y (mm)       7
        Energy (GeV)                    8
        Steps Inside Matter             9
        Distance Inside Matter (cm)     10
        dt                              11
        Kill Timestamp                  12
        '''
      
