# -*- coding: utf-8 -*-
"""
Created on Thu Nov  3 11:21:46 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import glob
    
def main():
    
    ts = 12
#    extra = "_angle" # Note the underscore that should be added
    extra = ""
    
#==============================================================================
# Particles
#==============================================================================
    
    particle_file = glob.glob("%s/../Output/photon_matrix%s_%d.csv"%(
                                os.getcwd(),extra,ts))
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
        x Calorimeter Angle             13
        y Calorimeter Angle             14
        Total Calorimeter Angle         15
        '''
        
        i = 0
        
        # of photons
        N_photons = len(stuff)
        
        # [Kill Event, Energy, dt, Kill Timestamp]
        # Possible kill events:
        # Pair-Production in Short Quad
        # Pair-Production in Long Quad
        # Pair-Production in HV Standoff
        # Pair-Production in HV Standoff Screw
        # Pair-Production in Standoff Plate
        photon = np.zeros((N_photons,4),dtype=object)
        
        # [x,y,z] Starting position (global coords)
        x = np.zeros((N_photons,3))
        
        # [x,y] Calorimeter contact position (local calorimeter coords)
        x_cal = np.zeros((N_photons,3))
        
        # [Steps inside matter, Distance inside matter]
        in_matter = np.zeros((N_photons,2))
        
        # Counters
        
        pp_sqel_counter = 0
        pp_dqel_counter = 0
        pp_sp_counter = 0
        pp_so_counter = 0
        pp_sos_counter = 0
        cal_counter = 0
        
        for row in stuff:
            
            photon[i,0] = row[2]
            photon[i,1] = row[8]
            photon[i,2] = row[11]
            photon[i,3] = row[12]
            
            if photon[i,0] == "Pair-Production in Short Quad":
                pp_sqel_counter = pp_sqel_counter + 1
            
            if photon[i,0] == "Pair-Production in Long Quad":
                pp_dqel_counter = pp_dqel_counter + 1
            
            if photon[i,0] == "Pair-Production in Standoff Plate":
                pp_sp_counter = pp_sp_counter + 1
            
            if photon[i,0] == "Pair-Production in HV Standoff":
                pp_so_counter = pp_so_counter + 1
            
            if photon[i,0] == "Pair-Production in HV Standoff Screw":
                pp_sos_counter = pp_sos_counter + 1
            
            if photon[i,0] == "Calorimeter Contact":
                cal_counter = cal_counter + 1
            
            x[i,0] = row[3]
            x[i,1] = row[4]
            x[i,2] = row[5]
            
            x_cal[i,0] = row[6]
            x_cal[i,1] = row[7]
            
            in_matter[i,0] = row[9]
            in_matter[i,1] = row[10]
            
            i = i + 1
            
#==============================================================================
# Data Processing and Display
#==============================================================================
            
    pp_tot = pp_sqel_counter + pp_dqel_counter + pp_sp_counter + \
            pp_so_counter + pp_sos_counter
    in_matter = np.array(in_matter,dtype=float)
    d_tot = sum(in_matter[:,1])
            
    print('Total # of photons created: %d'%N_photons)
    print('Total distance inside matter: %0.2f cm'%(d_tot*100))
    print('Total # of pair-production events: %d'%pp_tot)
    print('# of pair-production events in sqel: %d'%pp_sqel_counter)
    print('# of pair-production events in dqel: %d'%pp_dqel_counter)
    print('# of pair-production events in sp: %d'%pp_sp_counter)
    print('# of pair-production events in so: %d'%pp_so_counter)
    print('# of pair-production events in sos: %d'%pp_sos_counter)
    print('# of calorimeter contacts: %d'%cal_counter)
            
#==============================================================================
#     Plotting
#==============================================================================
        
    n = 0
    
    plt.figure(n)
    n = n + 1
    
    x_cal = np.array(x_cal, dtype = float)
    x_cal = x_cal[np.any(x_cal != 0, axis = 1)]
    ax = plt.subplot(1,1,1)
    ax.scatter(x_cal[:,0]*100, x_cal[:,1]*100,
               color='g',
               s = 0.7,
               label='Contact Points')
    ax.plot([-11.25,11.25],[7,7],'k-',label='Perimeter')
    ax.plot([-11.25,11.25],[-7,-7],'k-')
    ax.plot([-11.25,-11.25],[-7,7],'k-')
    ax.plot([11.25,11.25],[-7,7],'k-')
    plt.xlim(-24,24)
    plt.ylim(-15.5,15.5)
    ax.grid(True)
    ax.legend(bbox_to_anchor=(1.33,1.11))
    ax.set_title("Calorimeter Contact Position")
    ax.set_xlabel('x-Position (cm)')
    ax.set_ylabel('y-position (cm)')
    plt.axis('equal') # Prevents a skewed look
    
if __name__ == '__main__':

    main()