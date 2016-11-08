# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

import numpy as np
import callable_functions as cf

def track(particle_pos,particle_matrix,particle_proc,photon_pos,photon_proc,
          photon_matrix,dt,steps,m,B,k_min,k_max,geo_pack,particle_count,
          photon_steps,photon_dt,photon_row_index):
    
    c = 2.99792458*10**8                    # (m/s) Speed of light
    photon_kill_event_text = "Unknown Failure"
    photon_energy = photon_proc[photon_row_index,6]
#    print('photon energy: %0.3f'%(photon_energy/10**9))
    
    # Counter for number of steps a particle is inside matter
    steps_inside = 0
    
    step_counter = photon_proc[photon_row_index,10]
    
    x = np.zeros((photon_steps,3))         # Initialize photon position array
    p = np.zeros((3))                      # Initialize particle momentum array
    
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
    
    # Get the normalized velocity vector           
    p_norm = p / cf.mag(p)
    v_photon = p_norm*c

    i = 0
    
    while i < photon_steps - 1:
            
        cal_con_x = np.zeros((2))
        
        x[i+1] = x[i] + v_photon*photon_dt
        
        i = i + 1
        step_counter = step_counter + 1

        # Calorimeter front contact

        if np.abs(x[i,2]) < cal_height:
        
            if cf.noPassthroughElementContact(x[i],cal_rad,cal_theta):
                photon_kill_event_text = "Calorimeter Contact"
                r = cf.getParticleRadialPosition(x[i])
                cal_con_x = np.array([R_i + cal_width/2 - r,x[i,2]])
                break

        # Calorimeter top
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_box_theta):
            photon_kill_event_text = "Calorimeter Top Contact"
            break
        
        # Outer radius limit
        
        if cf.outerLimit(x[i],R + 0.2):
            photon_kill_event_text = "Heading Out"
            break

        if x[i,2] < so_z_max:
                    
            if cf.passthroughHVStandoff(x[i],so_rad,so_theta):
                steps_inside = steps_inside + 1
                if cf.ifPairProduction(photon_energy,photon_dt,"Ma"):
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = "Pair-Production in HV Standoff"
                    break
                    
            if cf.passthroughHVStandoffScrews(x[i],so_rad,so_theta):
                steps_inside = steps_inside + 1
                if cf.ifPairProduction(photon_energy,photon_dt,"SiBr"):
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = "Pair-Production in HV Standoff Screw"
                    break

        if x[i,2] < qel_z_max:

            if cf.passthroughElementContact(x[i],dqel_rad,dqel_theta):
                steps_inside = steps_inside + 1
                if cf.ifPairProduction(photon_energy,photon_dt,"Al"):
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = "Pair-Production in Long Quad"
                    break
    
            if cf.passthroughElementContact(x[i],sqel_rad,sqel_theta):
                steps_inside = steps_inside + 1
                if cf.ifPairProduction(photon_energy,photon_dt,"Al"):
                    particle_proc,particle_count = \
                        cf.doPairProduction(photon_energy,particle_count,
                                            particle_proc,m,p_norm,x[i],
                                            step_counter)
                    photon_kill_event_text = "Pair-Production Short Quad"
                    break

        if cf.passthroughElementContact(x[i],sp_rad,sp_theta):
            steps_inside = steps_inside + 1
            if cf.ifPairProduction(photon_energy,photon_dt,"Al"):
                particle_proc,particle_count = \
                    cf.doPairProduction(photon_energy,particle_count,
                                        particle_proc,m,p_norm,x[i],
                                        step_counter)
                photon_kill_event_text = "Pair-Production in Standoff Plate"
                break
            
    
    photon_proc[photon_row_index,8] = i
    photon_proc[photon_row_index,9] = 1
    
    photon_matrix[photon_row_index + 1] = np.array(
                         [photon_row_index,i,photon_kill_event_text,
                          x[0,0],
                          x[0,1],
                          x[0,2],
                          cal_con_x[0],
                          cal_con_x[1],
                          photon_energy/10**9,
                          steps_inside,steps_inside*c,photon_dt,
                          step_counter*photon_dt])
    
    return particle_pos,particle_matrix,particle_proc,photon_proc, \
           particle_count,photon_matrix