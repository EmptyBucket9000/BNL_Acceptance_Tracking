# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:42:14 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import glob


def main():
    
    file_name = "output_particles"
    files = glob.glob("../%s/Output/%s.csv"%(os.getcwd(),file_name))
    file = files[0]
    i = 0

    with open(file, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
    
        stuff = list(reader)
        N = len(stuff)
        particle = np.chararray((N,2))  # () [Particle #, Kill event]
        x = np.zeros((N,3))         # (m) [x,y,z] Local muon position at decay
        sigma = np.zeros((N,2))     # (m) [x,y,z] Beam width
        xbar = np.zeros((N,2))      # (m) [x,y,z] Beam centroid
        p = np.zeros((N,3))         # (eV/c) [Starting, Ending, Difference]
        inside = np.zeros((N,2))    # Inside matter [steps, distance (cm)]
        photon_count = np.zeros((N,2)) # [Total, # above minimum energy]
        
        for row in stuff:
            
            particle[i,0] = row[0]
            particle[i,1] = row[1]
            
            x[i,0] = row[2]
            xbar[i,0] = row[3]
            sigma[i,0] = row[4]
            
            x[i,1] = row[5]
            xbar[i,1] = row[6]
            sigma[i,1] = row[7]
            
            x[i,2] = row[8]
            
            p[i,0] = row[9]
            p[i,1] = row[10]
            p[i,2] = row[11]
            
            inside[i,0] = row[12]
            inside[i,1] = row[13]
            
            photon_count[i,0] = row[14]
            photon_count[i,1] = row[15]
            
            i = i + 1
            
    print('Photons released per centimeter in matter: %0.3f'
        %(sum(photon_count[:,0])/sum(inside[:,1])))
    print('High-energy photons released per centimeter in matter: %0.3f'
        %(sum(photon_count[:,1])/sum(inside[:,1])))
            
#==============================================================================
#     Plotting
#==============================================================================
        
    n = 0
    
    plt.figure(n)
    n = n + 1
    
    ax = plt.subplot(111)
    ax.scatter(x[:,0], x[:,1], color='g', label='Contact')
    ax.grid(True)
    lgd = ax.legend(bbox_to_anchor=(1.73,1.11))
    ax.set_title("Starting Position")
    ax.set_xlabel('x-Position (m)')
    ax.set_ylabel('y-position (m)')
    
if __name__ == '__main__':

    main()