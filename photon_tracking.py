# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

'''
See README.md for information.
'''

import numpy as np
import callable_functions as cf

def track(particle_pos,particle_matrix,particle_proc,photon_pos,photon_proc,
          photon_matrix,dt,steps,m,B,k_min,k_max,geo_pack,particle_count,
          photon_steps,photon_dt,photon_row_index,min_tracking,muon_number):
    
    c = 2.99792458*10**8                    # (m/s) Speed of light
    photon_kill_event_text = "Unknown Failure"
    photon_energy = photon_proc[photon_row_index,6]
    
    # Counter for number of steps a particle is inside matter
    steps_inside = 0
    
    step_counter = photon_proc[photon_row_index,10]
    
    x = np.zeros((photon_steps,3))         # Initialize photon position array
    p = np.zeros((3))                      # Initialize particle momentum array
    
    # Get the position of the photon and the momentum of the particle that
    # created it at the time of creation. The momentum vector is used for
    # directionality only.
    
    x[0,0] = photon_proc[photon_row_index,0]
    x[0,1] = photon_proc[photon_row_index,1]
    x[0,2] = photon_proc[photon_row_index,2]
    
    p[0] = photon_proc[photon_row_index,3]
    p[1] = photon_proc[photon_row_index,4]
    p[2] = photon_proc[photon_row_index,5]

    # Unpack 'geo_pack'
    
    cal_theta = geo_pack[0]
    cal_rad = geo_pack[2]
    cal_box_theta = geo_pack[3]
    so_z_max = geo_pack[5]
    so_rad = geo_pack[6]
    so_theta = geo_pack[7]
    sp_rad = geo_pack[10]
    sp_theta = geo_pack[11]
    qel_z_max = geo_pack[14]
    sqel_rad = geo_pack[15]
    sqel_theta = geo_pack[16]
    dqel_rad = geo_pack[17]
    dqel_theta = geo_pack[18]
    R = geo_pack[19]
    R_i = geo_pack[20]
    cal_width = geo_pack[21]
    cal_height = geo_pack[22]
    cal_theta_glob = geo_pack[23]
    rail_height = geo_pack[24]
    rail_rad = geo_pack[25]
    cal_det_theta = geo_pack[26]
    
    # Get the normalized momentum vector to create the photon velocity vector
    
    p_norm = p / cf.mag(p)
    v_photon = p_norm*c
            
    # Initialize the calorimeter contact position array
    cal_con_x = np.zeros((2))

    i = 0
    
    while i < photon_steps - 1:
        
        x[i+1] = x[i] + v_photon*photon_dt
        
        i = i + 1
        step_counter = step_counter + 1
        
        # Kill if energy below tracking energy        
        
        if photon_energy <= min_tracking:
            photon_kill_event_text = "Energy Below Minimum"
            break

        # Calorimeter front contact

        if np.abs(x[i,2]) < cal_height:
        
            if cf.noPassthroughElementContact(x[i],cal_rad,cal_theta):
                photon_kill_event_text = "Calorimeter Contact"
                cal_con_x = cf.getPositionOnCalorimeter(x[i],cal_width,R_i)
                break

        # Calorimeter top
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_det_theta):
            photon_kill_event_text = "Calorimeter Edge Contact"
            break
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_box_theta):
            photon_kill_event_text = "Calorimeter End Edge Contact"
            break
        
        # Outer radius limit
        
        if cf.outerLimit(x[i],R):
            photon_kill_event_text = "Heading Out"
            break
        
        # Trolly rail contact
        
        if cf.railContact(x[i],R,rail_height,rail_rad):
            photon_kill_event_text = "Trolley Rail Contact"
            break
        
        ## Check for and create pair-production events

        # Check if the photon's z-position is within the maximum height of the
        # HV standoff
        if x[i,2] < so_z_max:
                    
            # Check if the photon is inside an HV standoff
            if cf.passthroughHVStandoff(x[i],so_rad,so_theta):
                steps_inside = steps_inside + 1
                
                # Check if a pair-production event occurs
                if cf.ifPairProduction(photon_energy,photon_dt,"Ma"):
                    
                    # Add the new particles to the particle_proc array and
                    # increase the particle count
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = "Pair-Production in HV Standoff"
                    break
                    
            # Check if inside an HV standoff screw
            if cf.passthroughHVStandoffScrews(x[i],so_rad,so_theta):
                steps_inside = steps_inside + 1
                if cf.ifPairProduction(photon_energy,photon_dt,"SiBr"):
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = \
                        "Pair-Production in HV Standoff Screw"
                    break

        # Check if the photon's z-position is within the maximum height of the
        # electrodes
        if x[i,2] < qel_z_max:

            # Check if inside a double quad
            if cf.passthroughElementContact(x[i],dqel_rad,dqel_theta):
                steps_inside = steps_inside + 1
                if cf.ifPairProduction(photon_energy,photon_dt,"Al"):
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = "Pair-Production in Long Quad"
                    break
                
            # Check if inside a single quad    
            if cf.passthroughElementContact(x[i],sqel_rad,sqel_theta):
                steps_inside = steps_inside + 1
                if cf.ifPairProduction(photon_energy,photon_dt,"Al"):
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = "Pair-Production in Short Quad"
                    break

        # Check if inside a standoff plate
        if cf.passthroughElementContact(x[i],sp_rad,sp_theta):
            steps_inside = steps_inside + 1
            if cf.ifPairProduction(photon_energy,photon_dt,"Al"):
                particle_proc,particle_count = \
                    cf.doPairProduction(photon_energy,particle_count,
                                        particle_proc,m,p_norm,x[i],
                                        step_counter)
                photon_kill_event_text = "Pair-Production in Standoff Plate"
                break
            
    # Update the photon_proc array to indicate how many steps this photon took
    # and that it has been tracked
    photon_proc[photon_row_index,8] = i
    photon_proc[photon_row_index,9] = 1
    
    ang_x = 0
    ang_y = 0
    ang_tot = 0
    
    # Get the angle at which the incident particle/x-ray hits the calorimeter.
    
    if photon_kill_event_text == "Calorimeter Contact":
        
        # Get the projeced position on the calorimeter from the 2nd to last
        # photon position, this does not take into account A3
        cal_con_pre_x = cf.getPositionOnCalorimeter(x[i-1],cal_width,R_i)
        
        # Set this to pi/2 for now as sin(ang_tot) is needed in the second
        # iteration below but not in the first
        ang_tot = np.pi/2
        
        # See similar section in particle_tracking.py for comments on the loop
        k = 0
        while k < 2:
            
            # Add a z-component to the projection array from above, the 2nd
            # iteration of the code improves this value
            cal_con_pre_x[2] = np.sin(ang_tot) * c * dt
            
            ang_x,ang_y,ang_tot = \
                cf.getAnglesFromCalorimeter(cal_con_pre_x,cal_con_x,
                                            cal_theta_glob)
            
            k = k + 1
    
    # Add the new photon to the photon matrix for output to a single file
    photon_matrix[photon_row_index + 1] = np.array(
                         [photon_row_index,i,photon_kill_event_text,
                          x[0,0],
                          x[0,1],
                          x[0,2],
                          cal_con_x[0],
                          cal_con_x[1],
                          photon_energy/10**9,
                          steps_inside,steps_inside*c*photon_dt,photon_dt,
                          step_counter*photon_dt,
                          ang_x,
                          ang_y,
                          ang_tot,
                          muon_number])
    
    return particle_pos,particle_matrix,particle_proc,photon_proc, \
           particle_count,photon_matrix