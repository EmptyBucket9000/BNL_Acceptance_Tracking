# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:42:14 2016

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
        # Particle #                        0
        Steps                               1
        Kill Event                          2
        Charge                              3
        Starting Global x-Position (mm)	 4
        Starting Global y-Position (mm)	 5
        Starting Global z-Position (mm)     6
        Starting Momentum (GeV/c)           7
        Ending Momentum (GeV/c)             8
        Delta Momentum (GeV/c)              9
        Steps Inside Short Quad             10
        Distance Inside Short Quad (cm)	 11
        Total # of Photons Released         12
        # of Detectable Photons Released	 13
        Steps Inside Long Quad              14
        Distance Inside Long Quad (cm)	 15
        Total # of Photons Released         16
        # of Detectable Photons Released	 17
        Steps Inside Standoff Plate         18
        Distance Inside Standoff Plate (cm) 19
        Total # of Photons Released         20
        # of Detectable Photons Released	 21
        Steps Inside HV Standoff            22
        Distance Inside HV Standoff (cm)	 23
        Total # of Photons Released         24
        # of Detectable Photons Released	 25
        dt                                  26
        Kill Timestamp                      27
        '''
        N_particles = len(stuff)
        i = 0
                
        particle = np.zeros((N_particles,3),dtype=object) # () [Kill event,dt,charge]
        x = np.zeros((N_particles,3))         # (mm) [x,y,z] Global muon position at decay
        p = np.zeros((N_particles,3))         # (GeV/c) [Starting, Ending, Difference]
        
        # Inside matter [steps, distance (cm),total photons, HE photons]
        
        in_sqel = np.zeros((N_particles,4))
        in_dqel = np.zeros((N_particles,4))
        in_sp = np.zeros((N_particles,4))
        in_so = np.zeros((N_particles,4))
                
        # Counters
        
        cal_con_particle = 0
        cal_con_photon = 0
        so_contact = 0
        sqel_contact = 0
        dqel_contact = 0
        so_contact = 0
        sp_contact = 0
        cal_con_particle_so = 0     # Cal con and touch standoff
        for row in stuff:
            
            particle[i,0] = row[2]
            particle[i,1] = row[26]
            particle[i,2] = row[3]
            
            if row[2] == "Calorimeter Contact":
                cal_con_particle = cal_con_particle + 1
            
            x[i,0] = row[4]         
            x[i,1] = row[5]         
            x[i,2] = row[6]
            
            p[i,0] = row[7]
            p[i,1] = row[8]
            p[i,2] = row[9]
            
            in_sqel[i,0] = row[10]
            in_sqel[i,1] = row[11]
            in_sqel[i,2] = row[12]
            in_sqel[i,3] = row[13]
            
            if float(row[10]) > 0:
                sqel_contact = sqel_contact + 1
            
            in_dqel[i,0] = row[14]
            in_dqel[i,1] = row[15]
            in_dqel[i,2] = row[16]
            in_dqel[i,3] = row[17]
            
            if float(row[14]) > 0:
                dqel_contact = dqel_contact + 1
            
            in_sp[i,0] = row[18]
            in_sp[i,1] = row[19]
            in_sp[i,2] = row[20]
            in_sp[i,3] = row[21]
            
            if float(row[18]) > 0:
                sp_contact = sp_contact + 1
            
            in_so[i,0] = row[22]
            in_so[i,1] = row[23]
            in_so[i,2] = row[24]
            in_so[i,3] = row[25]
            
            if float(row[22]) > 0:
                so_contact = so_contact + 1
            
            if row[2] == "Calorimeter Contact" and float(row[22]) > 0:
                cal_con_particle_so = cal_con_particle_so + 1
            
            i = i + 1
            
#==============================================================================
# Data Processing
#==============================================================================
        
    ## Remove 'zero' rows
        
    x = x[0:i:1]
    p = p[0:i:1]
    in_sqel = in_sqel[0:i:1]
    in_dqel = in_dqel[0:i:1]
    in_sp = in_sp[0:i:1]
    in_so = in_so[0:i:1]
        
    total_particles = len(x)
    total_photons = sum(in_sqel[:,2]) + sum(in_dqel[:,2]) + \
                    sum(in_sp[:,2]) + sum(in_so[:,2])
    total_in_matter = sum(in_sqel[:,1]) + sum(in_dqel[:,1]) + \
                      sum(in_sp[:,1]) + sum(in_so[:,1])
    print('Total muon decays: %d'%total_particles)
    print('Total particles: %d'%total_particles)
    print('Total distance in matter: %0.3f cm'%total_in_matter)
    print('Total photons: %d'%total_photons)
    print('Total single quad contacts: %d'%sqel_contact)
    print('Total double quad contacts: %d'%dqel_contact)
    print('Total standoff plate contacts: %d'%sp_contact)
    print('Total HV standoff contacts: %d'%so_contact)
    print('Total particle calorimeter contacts: %d'%cal_con_particle)
    print('Total HV standoff contacts that hit the calorimeter: %d'%cal_con_particle_so)
        
#==============================================================================
# Photons
#==============================================================================
            
#==============================================================================
# Data Processing
#==============================================================================
        
    ## Remove 'zero' rows
            
#==============================================================================
#     Plotting
#==============================================================================
        
#    n = 0
#    
#    plt.figure(n)
#    n = n + 1
#    
#    ax = plt.subplot(111)
#    ax.scatter(x[:,0]/1000, x[:,1]/1000, color='g')
#    ax.grid(True)
##    lgd = ax.legend(bbox_to_anchor=(1.33,1.11))
#    ax.set_title("Starting Position")
#    ax.set_xlabel('x-Position (m)')
#    ax.set_ylabel('y-position (m)')
#    plt.axis('equal') # Prevents a skewed look
    
if __name__ == '__main__':

    main()