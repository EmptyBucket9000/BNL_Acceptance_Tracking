# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:42:14 2016

@author: Eric Schmidt
"""

"""
See README.md for information.

"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import glob
import geometry
from datetime import datetime
    
def main():
    
    startTime = datetime.now()
    
    save_plots = 1                      # Set to 1 to save plots, 0 otherwise
    save_dir = "../Output/Images"       # Set save directory
    image_dpi = 300                     # Set saved image dpi
    
    ts = 13
    extra = "_group_1" # E.g. "_group_2" Note the beginning underscore
    
#    particle_file = glob.glob("%s/../Output/particle_matrix%s_%d.csv"%(
#                                os.getcwd(),extra,ts))
    particle_file = glob.glob("%s/../Output/combined_particle_matrix_%d.csv"%(
                                os.getcwd(),ts))
    
#    photon_file = glob.glob("%s/../Output/photon_matrix%s_%d.csv"%(
#                                os.getcwd(),extra,ts))
    photon_file = glob.glob("%s/../Output/combined_photon_matrix_%d.csv"%(
                                os.getcwd(),ts))
    
    # Bin width in GeV for the energy histogram
    binwidth = 0.1
    
    show_angle_histograms = 0    
    
    # Number of bins in the calorimeter contact angle histogram
    calorimeter_angle_hist = 60
    
    geo_pack = geometry.geo()               # All the permanent geometries

    # Unpack 'geo_pack'
    
    cal_rad = geo_pack[2]
    so_rad = geo_pack[6]
    sp_rad = geo_pack[10]
    sqel_rad = geo_pack[15]
    dqel_rad = geo_pack[17]
    R = geo_pack[19]*100
    rail_rad = geo_pack[25]
    
#==============================================================================
# Photons
#==============================================================================
    
    photon_file = photon_file[0]
    
    with open(photon_file, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
        stuff = list(reader)
        
        i = 0
        
        # of photons
        N_photons = len(stuff)
        # [kill event,energy,ending x,ending y,muon #,muon #,muon set #]
        photon = np.zeros((N_photons,8),dtype=object)
        
        photon_cal_con = np.zeros((N_photons))
        
        for row in stuff:
            photon[i,0] = row[2]
            photon[i,1] = float(row[8])
            photon[i,2] = float(row[6])
            photon[i,3] = float(row[7])
            photon[i,4] = int(row[16])
            photon[i,5] = int(row[17])
            photon[i,6] = float(row[3])
            photon[i,7] = float(row[4])
            
            if photon[i,0] == "Calorimeter Contact" or \
                photon[i,0] == "Calorimeter Edge Contact":
                    
                photon_cal_con[i] = photon[i,1]
                
            i = i + 1
    
    photon_dist = np.copy(photon[photon[:,2] != 0,:])
#    photon_dist = np.copy(photon[0:i:1])
    photon = photon[0:i:1]
#==============================================================================
# Particles
#==============================================================================
    
    particle_file = particle_file[0]
    
    with open(particle_file, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
        stuff = list(reader)
        
        '''
        Particle #                          0
        Steps                               1
        Kill Event                          2
        Charge                              3
        Starting Global x-Position (mm)	 4
        Starting Global y-Position (mm)	 5
        Starting Global z-Position (mm)     6
        Ending Calorimeter x (mm)           7
        Ending Calorimeter y (mm)           8
        Starting Momentum (GeV/c)           9
        Ending Momentum (GeV/c)             10
        Delta Momentum (GeV/c)              11
        Steps Inside Short Quad             12
        Distance Inside Short Quad (cm)	 13
        Total # of Photons Released         14
        # of Detectable Photons Released	 15
        Steps Inside Long Quad              16
        Distance Inside Long Quad (cm)	 17
        Total # of Photons Released         18
        # of Detectable Photons Released	 19
        Steps Inside Standoff Plate         20
        Distance Inside Standoff Plate (cm) 21
        Total # of Photons Released         22
        # of Detectable Photons Released	 23
        Steps Inside HV Standoff            24
        Distance Inside HV Standoff (cm)	 25
        Total # of Photons Released         26
        # of Detectable Photons Released	 27
        Steps Inside HV Standoff Screws     28
        Distance Inside HV Standoff Screws	 29
        Total # of Photons Released         30
        # of Detectable Photons Released	 31
        dt                                  32
        Pair Produced (0 or 1)              33
        Kill Timestamp                      34
        x Calorimeter Angle                 35
        y Calorimeter Angle                 36
        Total Calorimeter Angle             37
        Starting Local x (mm)               38
        Startine Local y (mm)               39
        Starting Local x-prime (mrad)       40
        Starting Local y-prime (mrad)       41
        Muon #                              42
        Muon set #                          43
        '''
        
        N_particles = len(stuff)
        i = 0
        
        # [Kill event,dt,charge,pair-produced]
        particle = np.zeros((N_particles,4),dtype=object)
        
        # [part end x,part end y,photon end x,photon end y,energy,brem. loc.]
        part_phot = np.zeros((N_particles,8))
        
        # Momentum of particles that contact the calorimeter
        p_cal_con = np.zeros((N_particles))
        
        # [x,y,z] Global muon position at decay
        x = np.zeros((N_particles,3))       # (mm)
        
        # [Starting, Ending, Difference]
        p = np.zeros((N_particles,3))       # (GeV/c)
                
        # [Starting, Ending, Difference] Not from pp
        p_orig = np.zeros((N_particles,3))
        
        # [Starting, Ending, Difference] pp only particles
        p_pp = np.zeros((N_particles,3))    # (GeV/c)
        
        # Calorimeter contact position [x,y]
        x_cal = np.zeros((N_particles,2))
        
        # Inside matter [steps, distance (cm),total photons, HE photons]
        
        in_sqel = np.zeros((N_particles,4))
        in_dqel = np.zeros((N_particles,4))
        in_sp = np.zeros((N_particles,4))
        in_so = np.zeros((N_particles,4))
        in_sos = np.zeros((N_particles,4))
        
        # Calorimeter contact angles [x,y,total]
        angles = np.zeros((N_particles,3))
                
        # Counters
        
        rail_contact = 0            # Trolly rail contact
        cal_con = 0                 # Calorimeter contact
        so_contact = 0              # HV standoff contact
        sqel_contact = 0            # Single quad contact
        dqel_contact = 0            # Double quad contact
        sos_contact = 0             # HV standoff screw contact
        sp_contact = 0              # Standoff plate contact
        cal_con_so = 0              # Calorimeter and HV standoff contact
        cal_edge_con = 0        # Calorimeter edge contact (detected)
        cal_end_edge_contact = 0    # Cal end edge contact (not detected)
        pp_con = 0                  # Cal contact and pair-produced
        connectedCount = 0          # Counter for connecting particles and phot
        
        starting_pos = np.zeros((1,10),dtype=int)
        data_qel_contact = np.zeros((10))
        data_no_qel_contact = np.zeros((10))
        yerr_qel_contact = np.zeros((10))
        yerr_no_qel_contact = np.zeros((10))
        
        # sqel or dqel contact
        through_quad = np.zeros((1,10),dtype=int)
        no_through_quad = np.zeros((1,10),dtype=int)
        
        # sqel or dqel contact then cal contact
        through_quad_contact = np.zeros((1,10),dtype=int)
        
        # sqel or dqel contact then cal edge contact
        through_quad_edge_contact = np.zeros((1,10),dtype=int)
        
        # No sqel or dqel contact then cal contact
        no_through_quad_contact = np.zeros((1,10),dtype=int)
        
        # Color array for gammas vs. distance plot
        phot_color = np.zeros((len(part_phot)),dtype=object)
        
        # Color array for positron vs. distance plot
        pos_color = np.zeros((len(part_phot)),dtype=object)
        
        pos_energy_distance_alternate = 0
        
        for row in stuff:
            
            particle[i,0] = row[2]
            particle[i,1] = row[32]
            particle[i,2] = row[3]
            particle[i,3] = row[33]
            
            if (row[2] == "Calorimeter Contact" or \
                row[2] == "Calorimeter Edge Contact"):
                p_cal_con[i] = float(row[10])
            
            if (row[2] == "Calorimeter Contact"):
                    
                temp = [k for k, x in enumerate(photon_dist[:,4])
                        if x == int(row[42])]
                temp2 = [k for k in temp if photon_dist[k,5] == int(row[43])]
                
                if temp2 != []:
                    for t in temp2:
                        part_phot[connectedCount,0] = float(row[7])
                        part_phot[connectedCount,1] = float(row[8])
                        part_phot[connectedCount,2] = float(photon_dist[t,2])
                        part_phot[connectedCount,3] = float(photon_dist[t,3])
                        part_phot[connectedCount,4] = float(photon_dist[t,1])
                        part_phot[connectedCount,5] = float(row[10])
            
                        # sqel
                        if ((float(row[12]) > 0 or float(row[16]) > 0) and
                            float(row[20]) == 0 and float(row[24]) == 0 and
                            float(row[28]) == 0):
                            pos_color[connectedCount] = "Red"
            
                        # sp
                        elif (float(row[12]) == 0 and float(row[16]) == 0 and
                            float(row[20]) > 0 and float(row[24]) == 0 and
                            float(row[28]) == 0):
                            pos_color[connectedCount] = "Green"
            
                        # so
                        elif (float(row[12]) == 0 and float(row[16]) == 0 and
                            float(row[20]) == 0 and float(row[24]) > 0 and
                            float(row[28]) == 0):
                            pos_color[connectedCount] = "Purple"
                            
                        # sos
                        elif (float(row[12]) == 0 and float(row[16]) == 0 and
                            float(row[20]) == 0 and float(row[24]) > 0 and
                            float(row[28]) > 0):
                            pos_color[connectedCount] = "Blue"
                        
                        # None (Shouldn't exist)
                        elif (float(row[12]) == 0 and float(row[16]) == 0 and
                            float(row[20]) == 0 and float(row[24]) == 0 and
                            float(row[28]) == 0):
                                pos_color[connectedCount] == "Yellow"
                        
                        # Multiple
                        else:
                            pos_color[connectedCount] = "Orange"
                                
                        # Uncomment for alternate positron energy v distance
#                        if float(row[20]) > 0:
#                            pos_color[connectedCount] = "Green"
#                        pos_energy_distance_alternate = 1
                        
                        r_phot = np.sqrt((float(photon_dist[t,6]))**2 + 
                                    (float(photon_dist[t,7]))**2)
                        
                        if ((r_phot > sqel_rad[0] and r_phot < sqel_rad[1]) or
                            (r_phot > dqel_rad[0] and r_phot < dqel_rad[1])):
                            phot_color[connectedCount] = "Red"
                        
                        elif r_phot > sp_rad[0] and r_phot < sp_rad[1]:
                            phot_color[connectedCount] = "Green"
                        
                        elif r_phot > so_rad[0] and r_phot < so_rad[1]:
                            phot_color[connectedCount] = "Purple"
                                                
                        elif ((so_rad[0] < r_phot and 
                            (so_rad[0] + 0.00635) > r_phot) or
                            ((so_rad[1] - 0.00635) < r_phot and 
                            so_rad[1] > r_phot)):
                            phot_color[connectedCount] = "Blue"
                            
                        else:
                            phot_color[connectedCount] = "Orange"
                        
                        connectedCount = connectedCount + 1
            if row[33] == 1 and (row[2] == "Calorimeter Contact" or \
                row[2] == "Calorimeter Edge Contact"):
                pp_con = pp_con + 1
            
            if row[2] == "Calorimeter Contact":
                cal_con = cal_con + 1
            
            if row[2] == "Trolley Rail Contact":
                rail_contact = rail_contact + 1
                
            if row[2] == "Calorimeter Edge Contact":
                cal_edge_con = cal_edge_con + 1
                
            if row[2] == "Calorimeter End Edge Contact":
                cal_end_edge_contact = cal_end_edge_contact + 1
            
            x[i,0] = float(row[4])
            x[i,1] = float(row[5])
            x[i,2] = float(row[6])
            
            r = np.sqrt(x[i,0]**2 + x[i,1]**2)/10   # (cm) Starting radius
            d = (r - R)                               # (cm) Starting local x
#            d = float(row[38])/10
#            print(d)
            
            x_cal[i,0] = row[7]
            x_cal[i,1] = row[8]
            
            p[i,0] = row[9]
            p[i,1] = row[10]
            p[i,2] = row[11]
            
            if float(row[33]) == 0:
            
                p_orig[i,0] = row[9]
                p_orig[i,1] = row[10]
                p_orig[i,2] = row[11]
            
            if float(row[33]) == 1:
            
                p_pp[i,0] = row[9]
                p_pp[i,1] = row[10]
                p_pp[i,2] = row[11]
            
            in_sqel[i,0] = row[12]
            in_sqel[i,1] = row[13]
            in_sqel[i,2] = row[14]
            in_sqel[i,3] = row[15]
            
            if float(row[12]) > 0:
                sqel_contact = sqel_contact + 1
            
            in_dqel[i,0] = row[16]
            in_dqel[i,1] = row[17]
            in_dqel[i,2] = row[18]
            in_dqel[i,3] = row[19]
            
            if float(row[16]) > 0:
                dqel_contact = dqel_contact + 1
            
            in_sp[i,0] = row[20]
            in_sp[i,1] = row[21]
            in_sp[i,2] = row[22]
            in_sp[i,3] = row[23]
            
            if float(row[20]) > 0:
                sp_contact = sp_contact + 1
            
            in_so[i,0] = row[24]
            in_so[i,1] = row[25]
            in_so[i,2] = row[26]
            in_so[i,3] = row[27]
            
            if float(row[24]) > 0:
                so_contact = so_contact + 1
            
            in_sos[i,0] = row[28]
            in_sos[i,1] = row[29]
            in_sos[i,2] = row[30]
            in_sos[i,3] = row[31]
            
            if float(row[28]) > 0:
                sos_contact = sos_contact + 1
            
            if row[2] == "Calorimeter Contact" and \
                (float(row[24]) > 0 or float(row[28]) > 0):
                cal_con_so = cal_con_so + 1
                
            angles[i,0] = float(row[35])*180/np.pi
            angles[i,1] = float(row[36])*180/np.pi
            angles[i,2] = float(row[37])*180/np.pi
            
            ## For different starting x-position limits, check if calorimeter
            ## contact was made and if so, if the positron passed through any
            ## quads on the way
            
            if float(particle[i,3]) == 0:
                
                # If starting x-position (cm) is greater 3
                if d >= 3:
                    
                    # Add 1 to the counter for the total # of particles
                    starting_pos[0,8] = starting_pos[0,8] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,8] = through_quad[0,8] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,8] = \
                                through_quad_contact[0,8] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,8] = \
                                through_quad_edge_contact[0,8] + 1
                    
                    # Check if the particle hit the calorimter and missed quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,8] = no_through_quad[0,8] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,8] = \
                                no_through_quad_contact[0,8] + 1
                
                if d < 3 and d >= 2:
                    
                    # Add 1 to the counter for the total # of particles
                    starting_pos[0,7] = starting_pos[0,7] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,7] = through_quad[0,7] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,7] = \
                                through_quad_contact[0,7] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,7] = \
                                through_quad_edge_contact[0,7] + 1
                    
                    # Check if the particle hit the calorimter and missed quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,7] = no_through_quad[0,7] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,7] = \
                                no_through_quad_contact[0,7] + 1
                
                if d < 2 and d >= 1:
                    
                    starting_pos[0,6] = starting_pos[0,6] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,6] = through_quad[0,6] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,6] = \
                                through_quad_contact[0,6] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,6] = \
                                through_quad_edge_contact[0,6] + 1
                    
                    # Check if the particle hit the calorimter and missed quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,6] = no_through_quad[0,6] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,6] = \
                                no_through_quad_contact[0,6] + 1
                        
                if d < 1 and d >= 0:
                    
                    starting_pos[0,5] = starting_pos[0,5] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,5] = through_quad[0,5] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,5] = \
                                through_quad_contact[0,5] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,5] = \
                                through_quad_edge_contact[0,5] + 1
                    
                    # Check if the particle hit the calorimter and missed all quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,5] = no_through_quad[0,5] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,5] = \
                                no_through_quad_contact[0,5] + 1
                
                if d < 0 and d >= -1:
                    
                    starting_pos[0,4] = starting_pos[0,4] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,4] = through_quad[0,4] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,4] = \
                                through_quad_contact[0,4] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,4] = \
                                through_quad_edge_contact[0,4] + 1
                    
                    # Check if the particle hit the calorimter and missed quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,4] = no_through_quad[0,4] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,4] = \
                                no_through_quad_contact[0,4] + 1
                
                if d < -1 and d >= -2:
                    
                    starting_pos[0,3] = starting_pos[0,3] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,3] = through_quad[0,3] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,3] = \
                                through_quad_contact[0,3] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,3] = \
                                through_quad_edge_contact[0,3] + 1
                    
                    # Check if the particle hit the calorimter and missed quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,3] = no_through_quad[0,3] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,3] = \
                                no_through_quad_contact[0,3] + 1
                
                if d < -2 and d >= -3:
                    
                    starting_pos[0,2] = starting_pos[0,2] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,2] = through_quad[0,2] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,2] = \
                                through_quad_contact[0,2] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,2] = \
                                through_quad_edge_contact[0,2] + 1
                    
                    # Check if the particle hit the calorimter and missed quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,2] = no_through_quad[0,2] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,2] = \
                                no_through_quad_contact[0,2] + 1
                
                if d < -3:
                    
                    starting_pos[0,1] = starting_pos[0,1] + 1
                    
                    # Check if the particle hit the calorimeter and a quad
                    
                    if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        through_quad[0,1] = through_quad[0,1] + 1
                        
                        if particle[i,0] == "Calorimeter Contact":
                            through_quad_contact[0,1] = \
                                through_quad_contact[0,1] + 1
                        
                        if particle[i,0] == "Calorimeter Edge Contact":
                            through_quad_edge_contact[0,1] = \
                                through_quad_edge_contact[0,1] + 1
                    
                    # Check if the particle hit the calorimter and missed quads
                    
                    if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        no_through_quad[0,1] = no_through_quad[0,1] + 1
                        
                        if (particle[i,0] == "Calorimeter Contact" or 
                            particle[i,0] == "Calorimeter Edge Contact"):
                                
                            # Count the particle
                            no_through_quad_contact[0,1] = \
                                no_through_quad_contact[0,1] + 1
                
            i = i + 1
            
#==============================================================================
# Data Processing
#==============================================================================
    
    # Remove rows of all zeros
    part_phot = part_phot[np.any(part_phot != 0,axis=1)]
    phot_color = phot_color[phot_color != 0]
    phot_color = phot_color.astype('str')
    pos_color = pos_color[pos_color != 0]
    pos_color = pos_color.astype('str')
    angles = angles[np.any(angles != 0, axis = 1)]
    p_pp = p_pp[np.any(p_pp != 0,axis=1)]
    p_orig = p_orig[np.any(p_orig != 0,axis=1)]
    p_cal_con = p_cal_con[p_cal_con != 0]
    photon_cal_con = photon_cal_con[photon_cal_con != 0]
         
    dist = np.zeros((len(part_phot))) 
    k = 0
    for row in part_phot:
        distVec = np.array([row[0]-row[2],row[1]-row[3]])
        dist[k] = np.sqrt(np.dot(distVec,distVec))*100
        k = k + 1
    
    angles_mean = np.mean(angles[:,2])
    
    total_particles = len(x)
    total_photons = sum(in_sqel[:,2]) + sum(in_dqel[:,2]) + \
                    sum(in_sp[:,2]) + sum(in_so[:,2]) + sum(in_sos[:,2])
                      
    print('Total particles: %d'%total_particles)
#    print('Total distance in matter: %0.3f cm'%(total_in_matter))
    print('Total photons: %d'%total_photons)
    print('Total single quad contacts: %d'%sqel_contact)
    print('Total # of sqel photons released: %d'%sum(in_sqel[:,2]))
    print('Total double quad contacts: %d'%dqel_contact)
    print('Total # of dqel photons released: %d'%sum(in_dqel[:,2]))
    print('Total standoff plate contacts: %d'%sp_contact)
    print('Total # of sp photons released: %d'%sum(in_sp[:,2]))
    print('Total HV standoff contacts: %d'%so_contact)
    print('Total # of so photons released: %d'%sum(in_so[:,2]))
    print('Total HV standoff screw contacts: %d'%sos_contact)
    print('Total # of sos photons released: %d'%sum(in_sos[:,2]))
    print('Total # of trolley rail contacts: %d'%rail_contact)
    print('Total # of calorimeter edge contacts: %d'%cal_edge_con)
    print('Average calorimeter contact angle: %0.3f'%angles_mean)
    print('Total particle calorimeter contacts: %d'%(cal_con + cal_edge_con))
    print('# of calorimeter contacts after qel contact: %d'%(np.sum(
            through_quad_contact) + np.sum(through_quad_edge_contact)))
    print('# of calorimeter contacts without qel contact: %d'%np.sum(
            no_through_quad_contact))
    print('Total SO/SO screw contacts that hit the calorimeter: %d'\
            %cal_con_so)
    print('Total calorimeter contact from p.p. particles: %d'%pp_con)
    print('Positron/gamma pairs on cal contact: %d'%(len(part_phot)))
#==============================================================================
#     Plotting
#==============================================================================
            
    i = 0
    
    while i < np.shape(starting_pos)[1]:
        
        # Get the fraction that passed through a quad and hit a calorimeter to
        # the total # in that starting x-position range
        data_qel_contact[i] = (through_quad_contact[0,i] + 
                    through_quad_edge_contact[0,i]) / through_quad[0,i]
        
        # Get Poisson uncertainties
        yerr_qel_contact[i] = np.sqrt(through_quad_contact[0,i] + 
                                    through_quad_edge_contact[0,i]) / \
                                    through_quad[0,i]
        
        # Get the fraction that did not pass through a quad and hit a
        # calorimeter to the total # in that starting x-position range
        data_no_qel_contact[i] = (no_through_quad_contact[0,i]) / \
                                    no_through_quad[0,i]
        
        # Poisson uncertainties
        yerr_no_qel_contact[i] = np.sqrt(no_through_quad_contact[0,i]) / \
                                    no_through_quad[0,i]
        print('%d-Through, front contact: %d'%(i,through_quad_contact[0,i]))
        print('%d-Through, edge contact: %d'%(
            i,through_quad_edge_contact[0,i]))
        print('%d-No through, contact: %d'%(i,no_through_quad_contact[0,i]))
        print('%d-Total in starting range: %d'%(i,starting_pos[0,i]))
        i = i + 1
    
    print(data_qel_contact)
    print(yerr_qel_contact)
    print(data_no_qel_contact)
    print(yerr_no_qel_contact)
    
    # Convert string to float
    x_cal = np.array(x_cal, dtype = float)
    
    # Remove 'zero' rows
    x_cal = x_cal[np.any(x_cal != 0, axis = 1)]
        
    n = 0
    
    # Create a scatter plot of the fraction of particles that passed through
    # a quad and hit a calorimeter to the total # in that starting x-position
    # range, and create a scatter plot of the fraction that did not pass
    # through a quad to the total # in that starting x-position range.
    # The fit parameters were found using least-squares fitting.
    
    # Independent axis values for the fit functions
    xxfit = np.linspace(-4.2,4.2,500)
    
    # Quadratic fit parameters [a,b,c] in a + b*x + c*x**2
#    qfit_through = np.array([0.2559,-0.00543762,0.000746657])
#    qfit_no_through = np.array([0.508965,-0.00339719,-0.00186841])
    
    # Linear fit parameters [a,b] in a + b*x
    lfit_through = np.array([0.761623,-0.00577427])
    
    lfit_no_through = np.array([0.793202,-0.0123672])
    
    plt.figure(n)
    n = n + 1
    
    xx = np.arange(-4.5,5)
    ax = plt.subplot(1,1,1)
#    ax.errorbar(xx,data_no_qel_contact,yerr=yerr_no_qel_contact,fmt='.',
#                color='g')
    ax.scatter(xx,data_no_qel_contact,
                color='g')
    ax.errorbar(xx,data_qel_contact,yerr=yerr_qel_contact,fmt='.',
                color='b')
#    ax.plot(xxfit,qfit[0] + qfit[1]*xxfit + qfit[2] * xxfit**2,
#            label='Quadratic fit')
    ax.plot(xxfit,lfit_through[0] + lfit_through[1]*xxfit,
            label='Through electrode')
    ax.plot(xxfit,lfit_no_through[0] + lfit_no_through[1]*xxfit,
            label='No through')
    ax.legend(bbox_to_anchor=(1.4,1.1))
    ax.set_title("Position vs. Acceptance")
    ax.set_xlabel("x-Position (cm)")
    ax.set_ylabel("Particles Detected / Total")
    
    if save_plots == 1:
            plt.savefig('%s/acceptance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
    plt.figure(n)
    n = n + 1
    
    xx = np.arange(-4.5,5)
    ax = plt.subplot(1,1,1)
    ax.errorbar(xx,data_no_qel_contact,yerr=yerr_no_qel_contact,fmt='.',
                color='g')
    ax.plot(xxfit,lfit_no_through[0] + lfit_no_through[1]*xxfit)
#    ax.legend(bbox_to_anchor=(1,1))
    ax.set_title("Position vs. No through-quad acceptance")
    ax.set_xlabel("x-Position (cm)")
    ax.set_ylabel("Particles Detected / Total")
    
    if save_plots == 1:
            plt.savefig('%s/no_through_acceptance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
    plt.figure(n)
    n = n + 1
    
    xx = np.arange(-4.5,5)
    ax = plt.subplot(1,1,1)
    ax.errorbar(xx,data_qel_contact,yerr=yerr_qel_contact,fmt='.',
                color='b')
    ax.plot(xxfit,lfit_through[0] + lfit_through[1]*xxfit)
#    ax.legend(bbox_to_anchor=(1,1))
    ax.set_title("Position vs. Through-quad acceptance")
    ax.set_xlabel("x-Position (cm)")
    ax.set_ylabel("Particles Detected / Total")
    
    if save_plots == 1:
            plt.savefig('%s/through_acceptance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
#    plt.figure(n)
#    n = n + 1
#    
#    ax = plt.subplot(1,1,1)
#    ax.errorbar(xx,data_no_qel_contact,yerr=yerr_no_qel_contact,fmt='.',
#                color='g')
##    ax.plot(xxfit,qfit[0] + qfit[1]*xxfit + qfit[2] * xxfit**2,
##            label='Quadratic fit')
#    ax.plot(xxfit,lfit[0] + lfit[1]*xxfit,
#            label='Linear fit\n $\chi^2_R$ = ??')
#    ax.legend(bbox_to_anchor=(0.7,0.31))
#    ax.set_title("No through-quad acceptance vs. position")
#    ax.set_xlabel("x-Position (cm)")
#    ax.set_ylabel("Particles Detected / Total")
#    
#    if save_plots == 1:
#            plt.savefig('%s/acceptance_no_through.png'%save_dir,
#                        bbox_inches='tight',dpi=image_dpi)
    
    # Plot the contact points on the calorimeter in the calorimeter local
    # coordinate system
    
    plt.figure(n)
    n = n + 1
    
    ax = plt.subplot(1,1,1)
    ax.scatter(x_cal[:,0]*100, x_cal[:,1]*100,
               color='g',
               s = 0.7,
               label='Contact Points')
    ax.plot([-11.25,11.25],[7,7],'b-',label='Perimeter')
    ax.plot([-11.25,11.25],[-7,-7],'b-')
    ax.plot([-11.25,-11.25],[-7,7],'b-')
    ax.plot([11.25,11.25],[-7,7],'b-')
    ax.grid(True)
    ax.legend(bbox_to_anchor=(1.33,1.11))
    ax.set_title("Calorimeter Contact Position (Particles)")
    ax.set_xlabel('x-Position (cm)')
    ax.set_ylabel('y-position (cm)')
    plt.axis('equal')
    
    if save_plots == 1:
            plt.savefig('%s/particle_calorimeter_contact.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
      
    plt.figure(n)
    n = n + 1
    ax = plt.subplot(1,1,1)
    ax.hist(photon[:,1],
            bins=np.arange(min(photon[:,1]), max(photon[:,1]) + binwidth, 
                           binwidth), color="Red",
            label='Gammas at birth > 0.04 GeV')
    ax.hist(photon_cal_con,
            bins=np.arange(min(photon_cal_con), max(photon_cal_con) + binwidth, 
                           binwidth), color="Turquoise",
            label='Gammas at contact > 0.2 GeV')
#    ax.hist(p_pp[:,0],
#            bins=np.arange(min(p_pp[:,0]), max(p_pp[:,0]) + binwidth, 
#                           binwidth),
#            label='Pair-producded')
    ax.legend(bbox_to_anchor=(1,1))
    ax.set_title('Gamma Energy Histograms')
    ax.set_xlabel('Energy (GeV)')
    ax.set_ylabel('Count')
    
    if save_plots == 1:
            plt.savefig('%s/gamma_energy_distribution_hist.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
      
    plt.figure(n)
    n = n + 1
    ax = plt.subplot(1,1,1)
    ax.hist(p_orig[:,0],
            bins=np.arange(min(p_orig[:,0]), max(p_orig[:,0]) + binwidth, 
                           binwidth),
            label='Positron at birth')
    ax.hist(p_cal_con,
            bins=np.arange(min(p_cal_con), max(p_cal_con) + binwidth, 
                           binwidth),
            label='Positron at contact')
    ax.legend(bbox_to_anchor=(0.5,1))
    ax.set_title('Positron Momemtum')
    ax.set_xlabel('Energy (GeV/c)')
    ax.set_ylabel('Count')
    
    if save_plots == 1:
            plt.savefig('%s/positron_energy_distribution_hist.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
      
    plt.figure(n)
    n = n + 1
    binwidth = 0.25
    ax = plt.subplot(1,1,1)
    ax.hist(dist,
            bins=np.arange(min(dist), max(dist) + binwidth,binwidth))
    ax.set_title('Distance Between Positron and Gamma Contact Points')
    ax.set_xlabel('Distance (cm)')
    ax.set_ylabel('Count')
    if save_plots == 1:
            plt.savefig('%s/cal_distance_hist.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
    plt.figure(n)
    n = n + 1
    ax = plt.subplot(1,1,1)
    ax.scatter(part_phot[:,4],dist,s=6,c=phot_color,edgecolors="None")
    ax.set_title('Distance vs. Gamma Energy')
    ax.set_ylabel('Distance (cm)')
    ax.set_xlabel('Gamma Energy (GeV)')
    ax.text(1.9, 15, 'Red: Electrode \nGreen: Cage plate\nPurple: HV standoff\nBlue: HV standoff screw\nOrange: Multiple',
            bbox={'facecolor':'white', 'alpha':1, 'pad':3})
    ax.grid(True)
    if save_plots == 1:
            plt.savefig('%s/gamma_energy_distance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
    plt.figure(n)
    n = n + 1
    ax = plt.subplot(1,1,1)
    ax.scatter(part_phot[:,5],dist,s=6,c=pos_color,edgecolors="None")
    ax.set_ylabel('Distance (cm)')
    ax.set_xlabel('Positron Energy (GeV)')
    ax.grid(True)
    
    if pos_energy_distance_alternate == 0:
        ax.text(2.2, 15, 'Red: Electrode \nGreen: Cage plate\nPurple: HV standoff\nBlue: HV standoff screw\nOrange: Multiple',
                bbox={'facecolor':'white', 'alpha':1, 'pad':3})
        ax.set_title('Distance vs. Positron Energy')
        if save_plots == 1:
            plt.savefig('%s/positron_energy_distance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
        
    else:
        ax.text(2.2, 15, 'Red: Electrode \nGreen: At least cage plate\nPurple: HV standoff\nBlue: HV standoff screw\nOrange: Multiple (not CP)',
            bbox={'facecolor':'white', 'alpha':1, 'pad':3})
        ax.set_title('Distance vs. Positron Energy (Alternate)')
        if save_plots == 1:
            plt.savefig('%s/positron_energy_distance_alt.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
                        
    plt.figure(n)
    n = n + 1
    i = 40
    ax = plt.subplot(1,1,1)
    while i < 50:
        x = np.array([part_phot[i,0],part_phot[i,2]])*100
        y = np.array([part_phot[i,1],part_phot[i,3]])*100
        ax.plot(x,y,'b-')
        ax.scatter(x[1],y[1],color='r',s = 4)
        ax.scatter(x[0],y[0],color='b',s = 4)
        i = i + 1
    ax.plot([-11.25,11.25],[7,7],'b-',label='Perimeter')
    ax.plot([-11.25,11.25],[-7,-7],'b-')
    ax.plot([-11.25,-11.25],[-7,7],'b-')
    ax.plot([11.25,11.25],[-7,7],'b-')
    ax.grid(True)
    ax.set_title('Sample Set of Positron/Gamma Contacts')
    ax.set_xlabel('Distance (cm)')
    ax.set_ylabel('Distance (cm)')
    plt.axis('equal')
    
    if save_plots == 1:
            plt.savefig('%s/cal_separation.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
    # Plot a histogram of calorimeter contact angles where the angle is from
    # the positive x-axis
                        
    if show_angle_histograms == 1:
    
        plt.figure(n)
        n = n + 1
        
        ax = plt.subplot(1,1,1)
        ax.hist(angles[:,0],calorimeter_angle_hist)
        ax.set_title("x Calorimter Contact Angles (Particles)")
        ax.set_xlabel('Angle (deg)')
        ax.set_ylabel('Count')
        
        if save_plots == 1:
                plt.savefig('%s/particle_calorimeter_contact_x_angle.png'%
                            save_dir,bbox_inches='tight',dpi=image_dpi)
    
        # Plot a histogram of calorimeter contact angles where the angle is
        # from the positive y-axis
        
        plt.figure(n)
        n = n + 1
        
        ax = plt.subplot(1,1,1)
        ax.hist(angles[:,1],calorimeter_angle_hist)
        ax.set_title("y Calorimter Contact Angles (Particles)")
        ax.set_xlabel('Angle (deg)')
        ax.set_ylabel('Count')
        
        if save_plots == 1:
                plt.savefig('%s/particle_calorimeter_contact_y_angle.png'%
                            save_dir,bbox_inches='tight',dpi=image_dpi)
    
        # Plot a histogram of calorimeter contact angles where the angle is
        # from the projection of the final momentum vector onto the calorimeter
        # plane to the the final momentum vector (in local calorimeter
        # coordinates)
        
        plt.figure(n)
        n = n + 1
        
        ax = plt.subplot(1,1,1)
        ax.hist(angles[:,2],calorimeter_angle_hist)
        ax.set_title("Total Calorimter Contact Angles (Particles)")
        ax.set_xlabel('Angle (deg)')
        ax.set_ylabel('Count')
        
        if save_plots == 1:
                plt.savefig('%s/particle_calorimeter_contact_angle.png'%
                            save_dir,bbox_inches='tight',dpi=image_dpi)

    plt.show()
    print(datetime.now() - startTime)
    
if __name__ == '__main__':

    main()