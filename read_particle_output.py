# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:42:14 2016

@author: Eric Schmidt
"""

"""
Reads the particle and photon csv files and outputs the data to plots.

See README.md for more information.

"""

import numpy as np
import matplotlib.pyplot as plt
from find_fit_chisquare import fit as ffch
import csv
import os
import glob
import geometry
from datetime import datetime
    
def main():
    
    startTime = datetime.now()
    
    save_plots = 0                      # Set to 1 to save plots, 0 otherwise
    save_dir = "../Output/Images"       # Set save directory
    image_dpi = 300                     # Set saved image dpi
    num_slots = 20                      # Number of points for acceptance plots
    show_angle_histograms = 0           # Show a set of histograms
    show_old_plots = 0                  # Show other plots
    
    ts = 13
    
    # E.g. "_group_2" Note the beginning underscore, only used if a combined
    # csv file is not used.
    extra = "_group_1" 
    
    # Uncomment either the line with 'extra' or the line with the combined
    # file in each pair below.
    
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
    
    # Number of bins in the calorimeter contact angle histogram
    calorimeter_angle_hist = 60
    
    geo_pack = geometry.geo()               # All the permanent geometries

    # Unpack 'geo_pack'
    
    so_rad = geo_pack[6]
    sp_rad = geo_pack[10]
    sqel_rad = geo_pack[15]
    dqel_rad = geo_pack[17]
    R = geo_pack[19]*100
    
#==============================================================================
# Photons
#==============================================================================
    
    # See 'read_photon_output.py' for information on the photon csv file.
    
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
        
        # Create the photon array
        
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
    
    ## Create specific version of the 'photon' array depending on what is being
    ## searched for. This is done to drastically reduce search times later in
    ## the code.
    
    # Create an array of only those photons that hit the calorimeter front
    photon_dist = np.copy(photon[photon[:,2] != 0,:])
    
    # Photon array of photons that hit the side of the calorimeter
    photon_edge = np.copy(photon[photon[:,0] == "Calorimeter Edge Contact"])
    
    # Remove unused parts of the 'photon' array.
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
        
        # [x-prime,y-prime] for particles that contact the calorimeter
        prime_cal_con = np.zeros((N_particles,2))
        
        #[x-prime,y-prime]
        prime = np.zeros((N_particles,2))
        
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
        
        # Used in creating the sections for plotting acceptance vs x-prime

        max_prime = 0.004365
        min_prime = -0.004434        # [x-prime]
        prime_step = (max_prime - min_prime) / num_slots
        
        
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
        cal_edge_con = 0            # Calorimeter edge contact (detected)
        cal_end_edge_contact = 0    # Cal end edge contact (not detected)
        pp_con = 0                  # Cal contact and pair-produced
        connectedCount = 0          # Counter for connecting particles and phot
        
        # Counters for specific circumstances
        
        total_t = 0
        total_no_rail_t = 0
        rail_contact_t = 0            # Trolly rail contact
        cal_con_t = 0                 # Calorimeter contact
        cal_edge_con_t = 0
        so_contact_t = 0              # HV standoff contact
        sqel_contact_t = 0            # Single quad contact
        dqel_contact_t = 0            # Double quad contact
        sos_contact_t = 0             # HV standoff screw contact
        sp_contact_t = 0              # Standoff plate contact
        rail_contact_t = 0
        starting_energy_t = 0
        ending_energy_t = 0
        
        max_pos = 4.183048
        min_pos = -4.101644
#        min_pos = -max_pos

        pos_step = (max_pos - min_pos) / num_slots
        starting_pos = np.zeros((1,num_slots),dtype=int)
        
        # sqel or dqel contact
        through_quad = np.zeros((1,num_slots),dtype=int)
        through_matter = np.zeros((1,num_slots),dtype=int)
        no_through_quad = np.zeros((1,num_slots),dtype=int)
        no_through_matter = np.zeros((1,num_slots),dtype=int)
        
        # sqel or dqel contact then cal contact
        through_quad_contact = np.zeros((1,num_slots),dtype=int)
        through_matter_contact = np.zeros((1,num_slots),dtype=int)
        
        # sqel or dqel contact then cal edge contact
        through_quad_edge_contact = np.zeros((1,num_slots),dtype=int)
        through_matter_edge_contact = np.zeros((1,num_slots),dtype=int)
        
        # No sqel or dqel contact then cal contact
        no_through_quad_contact = np.zeros((1,num_slots),dtype=int)
        no_through_matter_contact = np.zeros((1,num_slots),dtype=int)
        
        # Color array for gammas vs. distance plot
        phot_color = np.zeros((len(part_phot)),dtype=object)
        
        # Color array for positron vs. distance plot
        pos_color = np.zeros((len(part_phot)),dtype=object)
        
        starting_pos_prime = np.zeros((1,num_slots),dtype=int)
        data_matter_contact_prime = np.zeros((num_slots))
        data_no_matter_contact_prime = np.zeros((num_slots))
        yerr_matter_contact_prime = np.zeros((num_slots))
        yerr_no_matter_contact_prime = np.zeros((num_slots))
        
        # sqel or dqel contact
        through_matter_prime = np.zeros((1,num_slots),dtype=int)
        no_through_matter_prime = np.zeros((1,num_slots),dtype=int)
        
        # sqel or dqel contact then cal contact
        through_matter_contact_prime = np.zeros((1,num_slots),dtype=int)
        
        # sqel or dqel contact then cal edge contact
        through_matter_edge_contact_prime = np.zeros((1,num_slots),dtype=int)
        
        # No sqel or dqel contact then cal contact
        no_through_matter_contact_prime = np.zeros((1,num_slots),dtype=int)
        
        pos_energy_distance_alternate = 0
        
        for row in stuff:
            
            particle[i,0] = row[2]
            particle[i,1] = row[32]
            particle[i,2] = row[3]
            particle[i,3] = row[33]
            
            if (particle[i,0] == "Calorimeter Contact" or \
                particle[i,0] == "Calorimeter Edge Contact"):
                p_cal_con[i] = float(row[10])
                prime_cal_con[i,0] = float(row[40])
                prime_cal_con[i,1] = float(row[41])
                
            ## Section for plotting calorimeter position-related data
            
            if (particle[i,0] == "Calorimeter Contact"):
                    
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
            
            x_init = float(row[38])*100          # (cm) Initial x-position
            
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
            
            prime[i,0] = float(row[40])
            prime[i,1] = float(row[41])
            
            # Look at statistics for starting x-position range
            
            if x_init > -4.5 and x_init < 4.5:
                
                total_t = total_t + 1
            
                if float(row[12]) > 0:
                    sqel_contact_t = sqel_contact_t + 1
                
                if float(row[16]) > 0:
                    dqel_contact_t = dqel_contact_t + 1
                
                if float(row[20]) > 0:
                    sp_contact_t = sp_contact_t + 1
                
                if float(row[24]) > 0:
                    so_contact_t = so_contact_t + 1
                
                if float(row[28]) > 0:
                    sos_contact_t = sos_contact_t + 1
            
                if particle[i,0] == "Trolley Rail Contact":
                    rail_contact_t = rail_contact_t + 1
                    
                if particle[i,0] != "Trolley Rail Contact":
                    starting_energy_t = starting_energy_t + float(row[9])
                    ending_energy_t = ending_energy_t + float(row[10])
                    total_no_rail_t = total_no_rail_t + 1
                
                if particle[i,0] == "Calorimeter Contact":
                    cal_con_t = cal_con_t + 1
                    
                if particle[i,0] == "Calorimeter Edge Contact":
                    cal_edge_con_t = cal_edge_con_t + 1
                
            
            if float(particle[i,3]) == 0:
                
                k = 0
                while k < num_slots :
            
                    ## For different starting x-prime limits, check if
                    ## calorimeter contact was made and if so, if the positron
                    ## passed through any quads on the way
                
                    if ((prime[i,0] >= \
                            min_prime + (k)*prime_step) and \
                        (prime[i,0] < \
                            min_prime + (k+1)*prime_step)):
                    
                        # Add 1 to the counter for the total # of particles
                        starting_pos_prime[0,k] = starting_pos_prime[0,k] + 1
                        
                        # Check if the particle hit the calorimeter and matter
                        
                        if (float(in_sqel[i,0]) > 0 or \
                            float(in_dqel[i,0]) > 0) or \
                            float(in_sp[i,0]) > 0 or \
                            float(in_so[i,0]) > 0 or \
                            float(in_sos[i,0]) > 0 or \
                            particle[i,0] == "Trolley Rail Contact":
                                
                            through_matter_prime[0,k] = through_matter_prime[0,k] + 1
                            
                            if particle[i,0] == "Calorimeter Contact":
                
                                if distancePartPhot(photon_dist,photon_edge,row):
                                    through_matter_contact_prime[0,k] = \
                                        through_matter_contact_prime[0,k] + 1
                            
                            if particle[i,0] == "Calorimeter Edge Contact":
                
                                if edgePartPhot(photon_edge,row):
                                    through_matter_edge_contact_prime[0,k] = \
                                        through_matter_edge_contact_prime[0,k] + 1
                        
                        # Check if the particle hit the calorimter and missed
                        # all matter
                        
                        if (float(in_sqel[i,0]) == 0 and \
                            float(in_dqel[i,0]) == 0 and \
                            float(in_sp[i,0]) == 0 and \
                            float(in_so[i,0]) == 0 and \
                            float(in_sos[i,0]) == 0) and \
                            particle[i,0] != "Trolley Rail Contact":
                                
                            no_through_matter_prime[0,k] = no_through_matter_prime[0,k] + 1
                            
                            if (particle[i,0] == "Calorimeter Contact" or 
                                particle[i,0] == "Calorimeter Edge Contact"):
                                    
                                # Count the particle
                                no_through_matter_contact_prime[0,k] = \
                                    no_through_matter_contact_prime[0,k] + 1
            
                    ## For different starting x-position limits

                    if x_init >= min_pos + k*pos_step and \
                        x_init < min_pos + (k+1)*pos_step:
                        
                        # Add 1 to the counter for the total # of particles
                        starting_pos[0,k] = starting_pos[0,k] + 1
                        
                        # Check if the particle hit the calorimeter and an elec
                        
                        if (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                            through_quad[0,k] = through_quad[0,k] + 1
                            
                            if particle[i,0] == "Calorimeter Contact":
                
                                if distancePartPhot(photon_dist,photon_edge,row):
                                    through_quad_contact[0,k] = \
                                        through_quad_contact[0,k] + 1
                            
                            if particle[i,0] == "Calorimeter Edge Contact":
                
                                if edgePartPhot(photon_edge,row):
                                    through_quad_edge_contact[0,k] = \
                                        through_quad_edge_contact[0,k] + 1
                        
                        # Check if the particle hit the calorimeter and matter
                        
                        if (float(in_sqel[i,0]) > 0 or \
                            float(in_dqel[i,0]) > 0) or \
                            float(in_sp[i,0]) > 0 or \
                            float(in_so[i,0]) > 0 or \
                            float(in_sos[i,0]) > 0 or \
                            particle[i,0] == "Trolley Rail Contact":
                                
                            through_matter[0,k] = through_matter[0,k] + 1
                            
                            if particle[i,0] == "Calorimeter Contact":
                
                                if distancePartPhot(photon_dist,photon_edge,row):
                                    through_matter_contact[0,k] = \
                                        through_matter_contact[0,k] + 1
                            
                            if particle[i,0] == "Calorimeter Edge Contact":
                
                                if edgePartPhot(photon_edge,row):
                                    through_matter_edge_contact[0,k] = \
                                        through_matter_edge_contact[0,k] + 1
                        
                        # Check if the particle hit the calorimter and missed quads
                        
                        if (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                            no_through_quad[0,k] = no_through_quad[0,k] + 1
                            
                            if (particle[i,0] == "Calorimeter Contact" or 
                                particle[i,0] == "Calorimeter Edge Contact"):
                                    
                                # Count the particle
                                no_through_quad_contact[0,k] = \
                                    no_through_quad_contact[0,k] + 1
                        
                        # Check if the particle hit the calorimter and missed
                        # all matter
                        
                        if (float(in_sqel[i,0]) == 0 and \
                            float(in_dqel[i,0]) == 0 and \
                            float(in_sp[i,0]) == 0 and \
                            float(in_so[i,0]) == 0 and \
                            float(in_sos[i,0]) == 0) and \
                            particle[i,0] != "Trolley Rail Contact":
                                
                            no_through_matter[0,k] = no_through_matter[0,k] + 1
                            
                            if (particle[i,0] == "Calorimeter Contact" or 
                                particle[i,0] == "Calorimeter Edge Contact"):
                                    
                                # Count the particle
                                no_through_matter_contact[0,k] = \
                                    no_through_matter_contact[0,k] + 1
                                    
                    k = k + 1
                
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
    prime_cal_con = prime_cal_con[np.any(prime_cal_con != 0,axis=1)]
         
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
    print('Minimum x-prime: %0.6f'%(min(prime_cal_con[:,0])))
    print('Maximum x-prime: %0.6f'%(max(prime_cal_con[:,0])))
    print('%% of %d positrons that passed through the following:'%total_t)
    print('El: %0.3f'%((sqel_contact_t + dqel_contact_t)*100/total_t))
    print('Sp: %0.3f'%(sp_contact_t*100/total_t))
    print('So: %0.3f'%(so_contact_t*100/total_t))
    print('Sos: %0.3f'%(sos_contact_t*100/total_t))
    print('Rail: %0.3f'%(rail_contact_t*100/total_t))
    print('Ave. starting momentum: %0.3f'%(starting_energy_t/total_no_rail_t))
    print('Ave. ending momentum: %0.3f'%(ending_energy_t/total_no_rail_t))
    print('# of calorimeter contacts: %d'%cal_con_t)
    print('# of edge contacts: %d'%cal_edge_con_t)
#==============================================================================
#     Plotting
#==============================================================================
        
    data_matter_contact_prime,yerr_matter_contact_prime, \
    data_no_matter_contact_prime,yerr_no_matter_contact_prime = \
        getDataVsAcceptance(through_matter_contact_prime,
                            through_matter_edge_contact_prime,
                            through_matter_prime,
                            no_through_matter_contact_prime,
                            no_through_matter_prime,
                            num_slots,
                            starting_pos_prime)
    
    print(data_matter_contact_prime)
    print(yerr_matter_contact_prime)
    print(data_no_matter_contact_prime)
    print(yerr_no_matter_contact_prime)
        
    
    # Convert string to float
    x_cal = np.array(x_cal, dtype = float)
    
    # Remove 'zero' rows
    x_cal = x_cal[np.any(x_cal != 0, axis = 1)]
    
    # Plot counter
    n = 0
    
    ## Position vs. Acceptance
    
    # Create a scatter plot of the fraction of particles that passed through
    # a quad and hit a calorimeter to the total # in that starting x-position
    # range, and create a scatter plot of the fraction that did not pass
    # through a quad to the total # in that starting x-position range.
    # The fit parameters were found using least-squares fitting.
    
    # Set x-axis for scatter plots
    
    prime_step = (max_prime - min_prime) / (num_slots)
    xx_prime = np.arange(min_prime + prime_step/2,max_prime,prime_step)
    
    init_guess_lin = np.array([0,0])
    init_guess_quad = np.array([0,0,0])
    fit_type = "poly"
    lfit_t = ffch(xx_prime,data_matter_contact_prime,yerr_matter_contact_prime,
                        init_guess_lin,fit_type)
    lfit_not = ffch(xx_prime,data_no_matter_contact_prime,
                    yerr_no_matter_contact_prime,init_guess_lin,fit_type)
    qfit_t = ffch(xx_prime,data_matter_contact_prime,
                  yerr_matter_contact_prime,init_guess_quad,fit_type)
    qfit_not = ffch(xx_prime,data_no_matter_contact_prime,
                    yerr_no_matter_contact_prime,init_guess_quad,fit_type)
                  
    print('lfit_t: %s'%lfit_t.message)  
    print('qfit_t: %s'%qfit_t.message)  
    print('lfit_not: %s'%lfit_not.message)  
    print('qfit_not: %s'%qfit_not.message)
    
    # Independent axis values for the fit functions
    xxfit_p = np.linspace(min_prime,max_prime,500)
    
    # Position vs. Through-matter acceptance
    
    fig = plt.figure(n)
    n = n + 1
    
    x2lin = lfit_t.fun/(len(xx_prime)-len(init_guess_lin)-1)
    dx2lin = np.sqrt(2*(len(xx_prime)-len(init_guess_lin)-1)) / \
            (len(xx_prime)-len(init_guess_lin)-1)
    x2quad = qfit_t.fun/(len(xx_prime)-len(init_guess_quad)-1)
    dx2quad = np.sqrt(2*(len(xx_prime)-len(init_guess_quad)-1)) / \
                (len(xx_prime)-len(init_guess_quad)-1)
    
    ax = plt.subplot(1,1,1)
    ax.errorbar(xx_prime,data_matter_contact_prime,
                yerr=yerr_matter_contact_prime,fmt='.',color='b')
    ax.plot(xxfit_p,lfit_t.x[0] + lfit_t.x[1]*xxfit_p)
    ax.plot(xxfit_p,qfit_t.x[0] + qfit_t.x[1]*xxfit_p + qfit_t.x[2]*xxfit_p**2)
    fig.text(0.24, -0.25,
             "Lin fit: y = %0.4f + %0.4f$x$\n"
             "$\chi^2_R$ = %0.2f$\pm$%0.2f\n"
             "Quad fit: y = %0.4f + %0.4f$x$ + %0.4f$x^2$\n"
             "$\chi^2_R$ = %0.2f$\pm$%0.2f"%(
                lfit_t.x[0],
                lfit_t.x[1],
                x2lin,dx2lin,
                qfit_t.x[0],
                qfit_t.x[1],
                qfit_t.x[2],
                x2quad,dx2quad
             ),
                bbox={'facecolor':'white', 'alpha':1, 'pad':3})
#    ax.legend(bbox_to_anchor=(1,1))
    ax.set_title("x-Prime vs. Through-matter acceptance")
    ax.set_xlabel("x-Prime (trans. mom. / long. mom.)")
    ax.set_ylabel("Particles Detected / Total")
    
    if save_plots == 1:
            plt.savefig('%s/through_acceptance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
    # Position vs. No through-matter acceptance
    
    fig = plt.figure(n)
    n = n + 1
    
    x2lin = lfit_not.fun/(len(xx_prime)-len(init_guess_lin)-1)
    dx2lin = np.sqrt(2*(len(xx_prime)-len(init_guess_lin)-1)) / \
            (len(xx_prime)-len(init_guess_lin)-1)
    x2quad = qfit_not.fun/(len(xx_prime)-len(init_guess_quad)-1)
    dx2quad = np.sqrt(2*(len(xx_prime)-len(init_guess_quad)-1)) / \
                (len(xx_prime)-len(init_guess_quad)-1)
    
    ax = plt.subplot(1,1,1)
    ax.errorbar(xx_prime,data_no_matter_contact_prime,
                yerr=yerr_no_matter_contact_prime,fmt='.',color='g')
    ax.plot(xxfit_p,lfit_not.x[0] + lfit_not.x[1]*xxfit_p)
    ax.plot(xxfit_p,qfit_not.x[0] + qfit_not.x[1]*xxfit_p + \
            qfit_not.x[2]*xxfit_p**2)
    fig.text(0.24, -0.25,
             "Lin fit: y = %0.4f + %0.4f$x$\n"
             "$\chi^2_R$ = %0.2f$\pm$%0.2f\n"
             "Quad fit: y = %0.4f + %0.4f$x$ + %0.4f$x^2$\n"
             "$\chi^2_R$ = %0.2f$\pm$%0.2f"%(
                lfit_not.x[0],
                lfit_not.x[1],
                x2lin,dx2lin,
                qfit_not.x[0],
                qfit_not.x[1],
                qfit_not.x[2],
                x2quad,dx2quad
             ),
          bbox={'facecolor':'white', 'alpha':1, 'pad':3})
#    ax.legend(bbox_to_anchor=(1,1))
    ax.set_title("x-Prime vs. No through-matter acceptance")
    ax.set_xlabel("x-Prime (trans. mom. / long. mom.)")
    ax.set_ylabel("Particles Detected / Total")
    
    if save_plots == 1:
            plt.savefig('%s/no_through_acceptance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
                        
#==============================================================================
# Old plots
#==============================================================================
    
    if show_old_plots == 1:
        
        data_qel_contact,yerr_qel_contact, \
        data_no_qel_contact,yerr_no_qel_contact = \
            getDataVsAcceptance(through_quad_contact,
                                through_quad_edge_contact,
                                through_quad,
                                no_through_quad_contact,
                                no_through_quad,
                                num_slots,
                                starting_pos)
            
        data_matter_contact,yerr_matter_contact, \
        data_no_matter_contact,yerr_no_matter_contact = \
            getDataVsAcceptance(through_matter_contact,
                                through_matter_edge_contact,
                                through_matter,
                                no_through_matter_contact,
                                no_through_matter,
                                num_slots,
                                starting_pos)
    
        # Position vs. Through-matter acceptance
    
        # Set x-axis for scatter plots
        
        t_step = (max_pos - min_pos) / (num_slots)
        xx = np.arange(min_pos+t_step/2,max_pos,t_step)/100
    
        # Independent axis values for the fit functions
        xxfit = np.linspace(min_pos,max_pos,500)/100
        
        fig = plt.figure(n)
        n = n + 1
        
        x2lin = lfit_t.fun/(len(xx)-len(init_guess_lin)-1)
        dx2lin = np.sqrt(2*(len(xx)-len(init_guess_lin)-1)) / \
                (len(xx)-len(init_guess_lin)-1)
        x2quad = qfit_t.fun/(len(xx)-len(init_guess_quad)-1)
        dx2quad = np.sqrt(2*(len(xx)-len(init_guess_quad)-1)) / \
                    (len(xx)-len(init_guess_quad)-1)
        
        ax = plt.subplot(1,1,1)
        ax.errorbar(xx,data_matter_contact,yerr=yerr_matter_contact,fmt='.',
                    color='b')
        ax.plot(xxfit,lfit_t.x[0] + lfit_t.x[1]*xxfit)
        ax.plot(xxfit,qfit_t.x[0] + qfit_t.x[1]*xxfit + qfit_t.x[2]*xxfit**2)
        fig.text(0.24, -0.25,
                 "Lin fit: y = %0.4f + %0.4f$x$\n"
                 "$\chi^2_R$ = %0.2f$\pm$%0.2f\n"
                 "Quad fit: y = %0.4f + %0.4f$x$ + %0.4f$x^2$\n"
                 "$\chi^2_R$ = %0.2f$\pm$%0.2f"%(
                    lfit_t.x[0],
                    lfit_t.x[1],
                    x2lin,dx2lin,
                    qfit_t.x[0],
                    qfit_t.x[1],
                    qfit_t.x[2],
                    x2quad,dx2quad
                 ),
                    bbox={'facecolor':'white', 'alpha':1, 'pad':3})
    #    ax.legend(bbox_to_anchor=(1,1))
        ax.set_title("Position vs. Through-matter acceptance")
        ax.set_xlabel("x-Position (m)")
        ax.set_ylabel("Particles Detected / Total")
        
        if save_plots == 1:
                plt.savefig('%s/through_acceptance_x_prime.png'%save_dir,
                            bbox_inches='tight',dpi=image_dpi)
        
        # Position vs. No through-matter acceptance
        
        fig = plt.figure(n)
        n = n + 1
        
        x2lin = lfit_not.fun/(len(xx)-len(init_guess_lin)-1)
        dx2lin = np.sqrt(2*(len(xx)-len(init_guess_lin)-1)) / \
                (len(xx)-len(init_guess_lin)-1)
        x2quad = qfit_not.fun/(len(xx)-len(init_guess_quad)-1)
        dx2quad = np.sqrt(2*(len(xx)-len(init_guess_quad)-1)) / \
                    (len(xx)-len(init_guess_quad)-1)
        
        ax = plt.subplot(1,1,1)
        ax.errorbar(xx,data_no_matter_contact,yerr=yerr_no_matter_contact,fmt='.',
                    color='g')
        ax.plot(xxfit,lfit_not.x[0] + lfit_not.x[1]*xxfit)
        ax.plot(xxfit,qfit_not.x[0] + qfit_not.x[1]*xxfit + qfit_not.x[2]*xxfit**2)
        fig.text(0.24, -0.25,
                 "Lin fit: y = %0.4f + %0.4f$x$\n"
                 "$\chi^2_R$ = %0.2f$\pm$%0.2f\n"
                 "Quad fit: y = %0.4f + %0.4f$x$ + %0.4f$x^2$\n"
                 "$\chi^2_R$ = %0.2f$\pm$%0.2f"%(
                    lfit_not.x[0],
                    lfit_not.x[1],
                    x2lin,dx2lin,
                    qfit_not.x[0],
                    qfit_not.x[1],
                    qfit_not.x[2],
                    x2quad,dx2quad
                 ),
              bbox={'facecolor':'white', 'alpha':1, 'pad':3})
    #    ax.legend(bbox_to_anchor=(1,1))
        ax.set_title("Position vs. No through-matter acceptance")
        ax.set_xlabel("x-Position (m)")
        ax.set_ylabel("Particles Detected / Total")
        
        if save_plots == 1:
                plt.savefig('%s/no_through_acceptance_x_prime.png'%save_dir,
                            bbox_inches='tight',dpi=image_dpi)
    
        # Position vs. Through-quad acceptance
        
        fig = plt.figure(n)
        n = n + 1
        
        ax = plt.subplot(1,1,1)
        ax.errorbar(xx,data_qel_contact,yerr=yerr_qel_contact,fmt='.',
                    color='b')
        ax.plot(xxfit,lfit_t[0] + lfit_t[1]*xxfit)
        fig.text(0.45, 0.85, 'Fit: y = %0.4f - %0.4fx'%(lfit_t[0],
                                                   np.abs(lfit_t[1])),
                    bbox={'facecolor':'white', 'alpha':1, 'pad':3})
    #    ax.legend(bbox_to_anchor=(1,1))
        ax.set_title("Position vs. Through-quad acceptance")
        ax.set_xlabel("x-Position (m)")
        ax.set_ylabel("Particles Detected / Total")
        
        if save_plots == 1:
                plt.savefig('%s/through_acceptance.png'%save_dir,
                            bbox_inches='tight',dpi=image_dpi)
        
        # Position vs. No through-quad acceptance
        
        fig = plt.figure(n)
        n = n + 1
        
        ax = plt.subplot(1,1,1)
        ax.errorbar(xx,data_no_qel_contact,yerr=yerr_no_qel_contact,fmt='.',
                    color='g')
        ax.plot(xxfit,lfit_not[0] + lfit_not[1]*xxfit)
        fig.text(0.45, 0.85,
                'Fit: y = %0.4f - %0.4fx'%(lfit_not[0],
                                       np.abs(lfit_not[1])),
              bbox={'facecolor':'white', 'alpha':1, 'pad':3})
    #    ax.legend(bbox_to_anchor=(1,1))
        ax.set_title("Position vs. No through-quad acceptance")
        ax.set_xlabel("x-Position (m)")
        ax.set_ylabel("Particles Detected / Total")
        
        if save_plots == 1:
                plt.savefig('%s/no_through_acceptance.png'%save_dir,
                            bbox_inches='tight',dpi=image_dpi)
        
        ## Calorimeter Contact Position (Particles)
        
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
          
        ## Gamma Energy Histograms
        
        plt.figure(n)
        n = n + 1
        ax = plt.subplot(1,1,1)
        ax.hist(photon[:,1],
                bins=np.arange(min(photon[:,1]), max(photon[:,1]) + binwidth, 
                               binwidth), color="Red",
                label='Gammas at birth > 0.04 GeV')
        ax.hist(photon_cal_con,
                bins=np.arange(min(photon_cal_con),
                               max(photon_cal_con) + binwidth, 
                               binwidth), color="Turquoise",
                label='Gammas at contact > 0.2 GeV')
        ax.legend(bbox_to_anchor=(1,1))
        ax.set_title('Gamma Energy Histograms')
        ax.set_xlabel('Energy (GeV)')
        ax.set_ylabel('Count')
        
        if save_plots == 1:
                plt.savefig('%s/gamma_energy_distribution_hist.png'%save_dir,
                            bbox_inches='tight',dpi=image_dpi)
          
        ## Positron Momemtum Histograms
        
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
          
        ## Distance Between Positron and Gamma Contact Points
        
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
        
        ## Distance vs. Gamma Energy
        
        plt.figure(n)
        n = n + 1
        ax = plt.subplot(1,1,1)
        ax.scatter(part_phot[:,4],dist,s=6,c=phot_color,edgecolors="None")
        ax.set_title('Distance vs. Gamma Energy')
        ax.set_ylabel('Distance (cm)')
        ax.set_xlabel('Gamma Energy (GeV)')
        ax.text(1.9, 15,
                'Red: Electrode \nGreen: Cage plate\nPurple: HV standoff\nBlue: HV standoff screw\nOrange: Multiple',
                bbox={'facecolor':'white', 'alpha':1, 'pad':3})
        ax.grid(True)
        if save_plots == 1:
                plt.savefig('%s/gamma_energy_distance.png'%save_dir,
                            bbox_inches='tight',dpi=image_dpi)
                            
        ## Distance vs. Positron Energy
        
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
    
        ## Sample Set of Positron/Gamma Contacts             
                            
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
    
def getDataVsAcceptance(through_matter_contact,through_matter_edge_contact,
                                through_matter,no_through_matter_contact,
                                no_through_matter,num_slots,starting):
                                    
    ## Set up data and uncertainties for "Position vs. Acceptance" plots                     
    
    data_matter_contact = np.zeros((num_slots))
    data_no_matter_contact = np.zeros((num_slots))
    yerr_matter_contact = np.zeros((num_slots))
    yerr_no_matter_contact = np.zeros((num_slots))
        
    i = 0
    
    while i < num_slots:
                    
        # Set variables used in determining uncertainties
        
        a = (through_matter_contact[0,i] + 
            through_matter_edge_contact[0,i])
            
        # Poisson uncertainty of 'a'
        delta_a = np.sqrt(a)
        b = np.copy(through_matter[0,i])
        
        # Poisson uncertainty of 'b'
        delta_b = np.sqrt(b)
        
        # Get the fraction that passed through matter and hit a calorimeter
        # to the total # in that starting x-position range that passed through
        data_matter_contact[i] = a / b
        
        # Get uncertainties                                    
        yerr_matter_contact[i] = np.sqrt(
                ((1/b)*delta_a)**2 + 
                ((a/b**2)*delta_b)**2 + 2*(1/b)*(-a/b**2)*delta_a**2
            )
        
        # Get the fraction that did not pass through matter and hit a
        # calorimeter to the total # in that starting x-position range
                                    
        a = np.copy(no_through_matter_contact[0,i])
        delta_a = np.sqrt(a)
        b = np.copy(no_through_matter[0,i])
        delta_b = np.sqrt(b)
        data_no_matter_contact[i] = a / b
        
        # Uncertainties                                    
        yerr_no_matter_contact[i] = np.sqrt(
                ((1/b)*delta_a)**2 + 
                ((a/b**2)*delta_b)**2 + 2*(1/b)*(-a/b**2)*delta_a**2
            )
            
        k = np.where(np.isnan(yerr_no_matter_contact))[0]
        
        if k.size:
            for el in k:
                yerr_no_matter_contact[el] = \
                    1/np.sqrt(no_through_matter_contact[0,el])
                print(no_through_matter_contact[0,el])
            
        print('%d-Through, front contact: %d'%(i,
                                               through_matter_contact[0,i]))
        print('%d-Through, edge contact: %d'%(
            i,through_matter_edge_contact[0,i]))
        print('%d-No through, contact: %d'%(i,
                                            no_through_matter_contact[0,i]))
        print('%d-Total in starting range: %d'%(i,starting[0,i]))
        
        i = i + 1
        
    return data_matter_contact,yerr_matter_contact, \
        data_no_matter_contact,yerr_no_matter_contact
        
def distancePartPhot(photon_dist,photon_edge,row):
    
    k = 0
    
    # Maximum distance to count as within one cell where actual cell width
    # is 2.5 cm
    md = 0.02       # (m)
    
    # Total energy of positrons and gammas if contact is within max distance
    total_energy = momentum2Energy(float(row[10]))
                    
    temp = [kk for kk, x in enumerate(photon_dist[:,4])
            if x == int(row[42])]
    temp2 = [kk for kk in temp if photon_dist[kk,5] == int(row[43])]
    
    if temp2 != []:
        
        # [part end x,part end y,photon end x,photon end y,photon energy]
        ppt = np.zeros((len(temp2),5))
        
        for t in temp2:
            ppt[k,0] = float(row[7])
            ppt[k,1] = float(row[8])
            ppt[k,2] = float(photon_dist[t,2])
            ppt[k,3] = float(photon_dist[t,3])
            ppt[k,4] = float(photon_dist[t,1])
            
            vec = np.array([ppt[k,0]-ppt[k,2],ppt[k,1]-ppt[k,3]])
            d = np.sqrt(np.dot(vec,vec))
            
            if d < md:
                total_energy = total_energy + ppt[k,4]
            
            k = k + 1
            
    k = 0
    
    temp = None
    temp2 = None
                    
    temp = [kk for kk, x in enumerate(photon_edge[:,4])
            if x == int(row[42])]
    temp2 = [kk for kk in temp if photon_edge[kk,5] == int(row[43])]
    
    if temp2 != []:
        for t in temp2:
            if float(row[7]) <= -0.0875:
                total_energy = total_energy + float(photon_edge[t,1])
            
            k = k + 1
            
    if total_energy > 1.8:
        return True
    else:
        return False
        
def edgePartPhot(photon_edge,row):
    
    k = 0
    # Total energy of positrons and gammas if contact is within max distance
    total_energy = momentum2Energy(float(row[10]))
                    
    temp = [k for k, x in enumerate(photon_edge[:,4])
            if x == int(row[42])]
    temp2 = [k for k in temp if photon_edge[k,5] == int(row[43])]
    
    if temp2 != []:        
        for t in temp2:            
            total_energy = total_energy + float(photon_edge[t,1])
            k = k + 1
    if total_energy > 1.8:
        return True
    else:
        return False

## Get total energy from momentum and mass

def momentum2Energy(p):
    
    m = 0.510999                      # (GeV/c**2) Particle mass
    energy = np.sqrt(np.dot(p,p) + m**2)
 
    return energy
    
if __name__ == '__main__':

    main()