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

def track(particle_pos,particle_matrix,particle_proc,photon_pos,
                 photon_proc,dt,steps,m,B,k_min,energy,geo_pack,
                 particle_count,photon_count,particle_row_index,muon_number,
                 min_tracking):

    c = 2.99792458*10**8                    # (m/s) Speed of light
    
    x = np.zeros((steps,3))                 # Initialize position array
    p = np.zeros((steps,3))                 # Initialize velocity array
    
    # Particle charge
    q = particle_proc[particle_row_index,6]
    
    # The number of steps into the photon's tracking when the pair-production
    # occured
    step_counter = particle_proc[particle_row_index,8]
    
    # Get the starting position and momentum of the particle
    
    x[0,0] = particle_proc[particle_row_index,0]
    x[0,1] = particle_proc[particle_row_index,1]
    x[0,2] = particle_proc[particle_row_index,2]
    
    p[0,0] = particle_proc[particle_row_index,3]
    p[0,1] = particle_proc[particle_row_index,4]
    p[0,2] = particle_proc[particle_row_index,5]
    
    # Get the energy of the particle
    energy = np.sqrt(np.dot(p[0],p[0]) + m**2)
    
    # (eV) Minimum energy detectable by the calorimeter
    min_detectable_energy = 0.2*10**9
    
    # Counter for number of steps a particle is inside matter
    # [sqel,dqel,sp,so,sos]
    steps_inside = np.zeros((5))
    
    # How far the particle has traveled in matter
    # [sqel,dqel,sp,so,sos]
    d_matter = np.zeros((5))                             # (m)
    
    # Radiation lengths
    X0_al = 0.08897                  # (m) Radiation length of aluminum
    X0_ma = 0.05198                  # (m) Radiation length of macor
    X0_sibr = 0.01468                # (m) Radiation length of silicon bronze
#    X0_al = 10
#    X0_ma = 10
#    X0_sibr = 10  

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
    
    n = .142                                # () Used in E-field
    
    # Electric field
    E = np.zeros((3))                       # (V/m) Initialize E-field
        
    if cf.isInSQuad(x[0],sqel_theta,R) or \
        cf.isInDQuad(x[0],dqel_theta,R):
            inEField = 1

    else:
        inEField = 0
        
    E = cf.getElectricField(x[0],B,R,n,inEField)       # Set initial E-field values
    
    # Event text
    kill_event_text = "Unknown failure" # In case nothing happens
    
    # of photons released [total, above min. energy]
    sqel_photon_count = np.zeros([2], dtype=int)
    dqel_photon_count = np.zeros([2], dtype=int)
    sp_photon_count = np.zeros([2], dtype=int)
    so_photon_count = np.zeros([2], dtype=int)
    sos_photon_count = np.zeros([2], dtype=int)     # HV standoff screws
    total_photon_count = 0
        
    # Initialize the calorimeter contact position array
    cal_con_x = np.zeros((2))
    
#==============================================================================
#   Tracking by Runga-Kutta 4th
#==============================================================================
    
    # Loop counter
    i = 0
    
    while i < steps - 1:
            
        a = q*c**2*(E + np.cross(p[i]/energy,B))
        dp1 = a*dt
        
        a = q*c**2*(E + np.cross((p[i] + dp1/2)/energy,B))
        dp2 = a*dt
        
        a = q*c**2*(E + np.cross((p[i] + dp2/2)/energy,B))
        dp3 = a*dt
        
        a = q*c**2*(E + np.cross((p[i] + dp3)/energy,B))
        dp4 = a*dt
        
        dp = (dp1 + 2*dp2 + 2*dp3 + dp4) / 6
    
        p[i+1] = p[i] + dp
        
        x[i+1] = x[i] + ((p[i+1]) / energy)*c*dt
        
        i = i + 1        
        energy = cf.momentum2Energy(p[i],m)        
        step_counter = step_counter + 1
        
        ''' Check for contact with permanent geometries '''
        
        # Quad electrodes, first check if z-position within range        
        if np.abs(x[i,2]) < qel_z_max:
        
            # Single-quad electrode            
            if cf.passthroughElementContact(x[i],sqel_rad,sqel_theta):
                
                # Update the variables to track time spent inside this matter
                steps_inside[0],d_matter[0] = \
                    cf.updateInsideMatter(p[i],energy,dt,steps_inside[0],
                                          d_matter[0])
                
                # Check if any photons are released due Bremsstrahlung
                photons_released = \
                    cf.isPhotonReleased(k_min,energy,X0_al,p[i],dt,m)
                    
                # If any photons are released
                if photons_released > 0:
                    
                    k = 0
                    
                    # Check each photon
                    while k < photons_released:
            
                        # Make sure the current particle energy is greater than
                        # the minimal energy we care about
                        if energy > k_min: # Nonsense if k_min > energy
                            p_norm = p[i]/cf.mag(p[i])
                            
                            # Do the Bremsstrahlung event, adjusting the 
                            # particle's momentum accordingly
                            p[i],sqel_photon_count, \
                            photon_energy,total_photon_count = \
                                cf.bremsstrahlung(p[i],m,k_min,energy,
                                                  i,sqel_photon_count,
                                                  min_detectable_energy,
                                                  total_photon_count)
                                                  
                            # Create the new line in the photon tracking array
                            photon_proc[total_photon_count-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          p_norm[0],p_norm[1],p_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                        k = k + 1
            
            # Double-quad electrode
            
            if cf.passthroughElementContact(x[i],dqel_rad,dqel_theta):
                
                steps_inside[1],d_matter[1] = \
                    cf.updateInsideMatter(p[i],energy,dt,steps_inside[1],
                                          d_matter[1])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,energy,X0_al,p[i],dt,m)
                    
                if photons_released > 0:
                    
                    k = 0                    
                    while k < photons_released:
            
                        if energy > k_min: # Nonsense if k_min > energy
                            p_norm = p[i]/cf.mag(p[i])
                            p[i], dqel_photon_count, \
                            photon_energy,total_photon_count = \
                                cf.bremsstrahlung(p[i],m,k_min,energy,
                                                  i,dqel_photon_count,
                                                  min_detectable_energy,
                                                  total_photon_count)
                            photon_proc[total_photon_count-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          p_norm[0],p_norm[1],p_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                        k = k + 1
        
        # Side support plate
        
        if cf.passthroughElementContact(x[i],sp_rad,sp_theta):
                
                steps_inside[2],d_matter[2] = \
                    cf.updateInsideMatter(p[i],energy,dt,steps_inside[2],
                                          d_matter[2])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,energy,X0_al,p[i],dt,m)
                    
                if photons_released > 0:
                    
                    k = 0                    
                    while k < photons_released:
            
                        if energy > k_min: # Nonsense if k_min > energy
                            p_norm = p[i]/cf.mag(p[i])
                            p[i], sp_photon_count, \
                            photon_energy,total_photon_count = \
                                cf.bremsstrahlung(p[i],m,k_min,energy,
                                                  i,sp_photon_count,
                                                  min_detectable_energy,
                                                  total_photon_count)
                            photon_proc[total_photon_count-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          p_norm[0],p_norm[1],p_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                        k = k + 1
        
        # High-voltage standoff
        
        if x[i,2] < so_z_max:
            if cf.passthroughHVStandoff(x[i],so_rad,so_theta):
                
                steps_inside[3],d_matter[3] = \
                    cf.updateInsideMatter(p[i],energy,dt,steps_inside[3],
                                          d_matter[3])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,energy,X0_ma,p[i],dt,m)
                    
                if photons_released > 0:
                    
                    k = 0                    
                    while k < photons_released:
            
                        if energy > k_min: # Nonsense if k_min > energy
                            p_norm = p[i]/cf.mag(p[i])
                            p[i], so_photon_count, \
                            photon_energy,total_photon_count = \
                                cf.bremsstrahlung(p[i],m,k_min,energy,
                                                  i,so_photon_count,
                                                  min_detectable_energy,
                                                  total_photon_count)
                            photon_proc[total_photon_count-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          p_norm[0],p_norm[1],p_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                        k = k + 1
        
            # High-voltage standoff screws
        
            if cf.passthroughHVStandoffScrews(x[i],so_rad,so_theta):
                
                steps_inside[4],d_matter[4] = \
                    cf.updateInsideMatter(p[i],energy,dt,steps_inside[4],
                                          d_matter[4])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,energy,X0_sibr,p[i],dt,m)
                    
                if photons_released > 0:
                    
                    k = 0                    
                    while k < photons_released:
            
                        if energy > k_min: # Nonsense if k_min > energy
                            p_norm = p[i]/cf.mag(p[i])
                            p[i], sos_photon_count, \
                             photon_energy,total_photon_count = \
                                cf.bremsstrahlung(p[i],m,k_min,energy,
                                                  i,sos_photon_count,
                                                  min_detectable_energy,
                                                  total_photon_count)
                            photon_proc[total_photon_count-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          p_norm[0],p_norm[1],p_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                        k = k + 1
            
        # Break if particle energy below detectability/10
        
        energy = cf.momentum2Energy(p[i],m)
        
        if energy <= min_tracking:
            kill_event_text = "Energy Below Minimum"
            break
                
        photons_released = 0 # Reset for next step in particle motion
        
        if cf.isInSQuad(x[i],sqel_theta,R) or \
            cf.isInDQuad(x[i],dqel_theta,R):
                inEField = 1
                
        else:
            inEField = 0
        
        # Get the electric field based on position
        E = cf.getElectricField(x[i],B,R,n,inEField)
        
        # Check if the particle has come into contact with the calorimeter and
        # kill the loop if it has
        
        if np.abs(x[i,2]) < cal_height/2:
        
            if cf.noPassthroughElementContact(x[i],cal_rad,cal_theta):
                kill_event_text = "Calorimeter Contact"
                cal_con_x = cf.getPositionOnCalorimeter(x[i],cal_width,R_i)
                break
            
        # Calorimeter top
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_det_theta):
            kill_event_text = "Calorimeter Edge Contact"
            break
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_box_theta):
            kill_event_text = "Calorimeter End Edge Contact"
            break
                        
        # Inner radius limit
        
        if cf.innerLimit(x[i],R_i):
            kill_event_text = "Inner Limit Reached"
            break
        
        # Outer radius limit
        
        if cf.outerLimit(x[i],R):
            kill_event_text = "Into the Iron"
            break
        
        # Trolly rail contact
        
        if cf.railContact(x[i],R,rail_height,rail_rad):
            kill_event_text = "Trolly Rail Contact"
            break
        
    # Prepare the data to be saved into the particle matrix for later saving
    # to a single file.
    
    p_init_mag = cf.mag(p[0])
    p_end_mag = cf.mag(p[i])
    particle_pos[particle_row_index] = np.copy(x)
    charge = particle_proc[particle_row_index,6]
    pp = particle_proc[particle_row_index,9]
    
    ang_x = 0
    ang_y = 0
    ang_tot = 0
    
    # Get the angle at which the incident particle/x-ray hits the calorimeter
    
    if kill_event_text == "Calorimeter Contact":
        
        # Get the projeced position on the calorimeter from the 2nd to last
        # particle position, this does not take into account A3
        cal_con_pre_x = cf.getPositionOnCalorimeter(x[i-1],cal_width,R_i)
        
        # Set this to pi/2 for now as sin(ang_tot) is needed in the second
        # iteration below but not in the first
        ang_tot = np.pi/2
        
        k = 0
        
        # Loop twice to provide better accuracy as the first loop assumes the
        # z-component of the particle vector is equal to the magnitude of the
        # velocity of the particle * dt, when in fact it is v*dt*sin(angle),
        # where angle is the angle of the particle vector with respect to the
        # calorimeter plane. The second iteration used the angle found in the
        # first iteration to provide a better z-component. The further the
        # angle between the vector and the plane is from pi/2, the less
        # accurate the first iteration is.
        while k < 2:
            
            # Add a z-component to the projection array from above, the 2nd
            # iteration of the code improves this value
            cal_con_pre_x[2] = p_end_mag/cf.momentum2Energy(p_end_mag,m) * \
                                np.sin(ang_tot) * c * dt
            
            ang_x,ang_y,ang_tot = \
                cf.getAnglesFromCalorimeter(cal_con_pre_x,cal_con_x,
                                            cal_theta_glob)
            
            k = k + 1
        
    # Add the new line to the particle matrix array that will be saved to file
    particle_matrix[particle_row_index + 1] = np.array(
                         [particle_row_index,                       # 0
                          i,                                        # 1
                          '%s'%kill_event_text,                     # 2
                          charge,                                   # 3
                          '%0.3f'%(x[0,0]*10**3),                   # 4
                          '%0.3f'%(x[0,1]*10**3),                   # 5
                          '%0.3f'%(x[0,2]*10**3),                   # 6
                          cal_con_x[0],                             # 7
                          cal_con_x[1],                             # 8
                          '%0.7f'%(p_init_mag/10**9),               # 9
                          '%0.7f'%(p_end_mag/10**9),                # 10
                          '%0.7f'%((p_end_mag - p_init_mag)/10**9), # 11
                          steps_inside[0],                          # 12
                          d_matter[0]*100,                          # 13
                          sqel_photon_count[0],                     # 14
                          sqel_photon_count[1],                     # 15
                          steps_inside[1],                          # 16
                          d_matter[1]*100,                          # 17          
                          dqel_photon_count[0],                     # 18
                          dqel_photon_count[1],                     # 19
                          steps_inside[2],                          # 20
                          d_matter[2]*100,                          # 21
                          sp_photon_count[0],                       # 22
                          sp_photon_count[1],                       # 23
                          steps_inside[3],                          # 24
                          d_matter[3]*100,                          # 25
                          so_photon_count[0],                       # 26
                          so_photon_count[1],                       # 27
                          steps_inside[4],                          # 28
                          d_matter[4]*100,                          # 29
                          sos_photon_count[0],                      # 30
                          sos_photon_count[1],                      # 31
                          '%2e'%dt,                                 # 32
                          pp,                                       # 33
                          '%5e'%(step_counter*dt),                  # 34
                          ang_x,                                    # 35
                          ang_y,                                    # 36
                          ang_tot])                                 # 37
                          
    # Update the particle tracking array to account for having fully tracked
    # this particle
    particle_count = particle_count + 1
    particle_proc[particle_row_index,7] = 1
                            
    return particle_pos,particle_matrix,particle_proc,photon_count, \
            photon_proc