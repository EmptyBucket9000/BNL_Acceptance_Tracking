# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

import numpy as np
import callable_functions as cf

def track(particle_pos,particle_matrix,particle_proc,photon_pos,
                 photon_proc,dt,steps,m,B,k_min,k_max,geo_pack,
                 particle_count,photon_count,particle_row_index,muon_number):

    c = 2.99792458*10**8                    # (m/s) Speed of light
    
    x = np.zeros((steps,3))                 # Initialize position array
    v = np.zeros((steps,3))                 # Initialize velocity array
    
    q = particle_proc[particle_row_index,6]
    
    step_counter = particle_proc[particle_row_index,8]
    
    x[0,0] = particle_proc[particle_row_index,0]
    x[0,1] = particle_proc[particle_row_index,1]
    x[0,2] = particle_proc[particle_row_index,2]
    
    v[0,0] = particle_proc[particle_row_index,3]
    v[0,1] = particle_proc[particle_row_index,4]
    v[0,2] = particle_proc[particle_row_index,5]
    
    min_detectable_energy = 0.2*10**9       # (eV) Minimum energy detectable

    p_init = cf.beta2Momentum(v[0]/c,m)
    
    # Counter for number of steps a particle is inside matter
    # [sqel,dqel,sp,so]
    steps_inside = np.zeros((4))
    
    # How far the particle has traveled in matter
    # [sqel,dqel,sp,so]
    d_matter = np.zeros((4))                             # (m)
    
    # Radiation lengths
#    X0_al = 0.08897                         # (m) Radiation length of aluminum
#    X0_ma = 0.131                           # (m) Radiation length of macor
    X0_al = 0.005
    X0_ma = 0.005

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
    
    n = .142                                # () Used in E-field
    
    # Electric field
    E = np.zeros((3))                       # (V/m) Initialize E-field
        
    if cf.isInSQuad(x[0],sqel_theta,R) or \
        cf.isInDQuad(x[0],dqel_theta,R):
            loc = "In"
            
    else:
        loc = "Out"
        
    E = cf.getElectricField(x[0],B,R,n,loc)           # Set initial E-field values
    
    # Force vector 
    F = np.zeros((3))                       # (N) Initialze force array
    F = cf.forceDueFields(v[0],B,E,q)       # Set initial force due fields

    # Magnitude of relativistic beta vector
    beta = v[0]/c
    beta_mag = np.zeros(1)                  # ()
    beta_mag[0] = cf.mag(beta)
    
    # Relativistic gamma        
    gamma = np.zeros(1)                     # ()
    gamma[0] = cf.beta2Gamma(beta_mag[0])
    
    # Event text
    kill_event_text = "Unknown failure" # In case nothing happens
    
    photon_count = np.zeros([2], dtype='int') # of photons released
    
#==============================================================================
#   Tracking by Runga-Kutta 4th
#==============================================================================
    
    # Loop counter
    i = 0
    
    ''' RK4 work '''
    
    while i < steps - 1:
        
        # Relativistic mass to find the acceleration from the force
        rmass = m*gamma/c**2
        
        a = F/(rmass)
        dv1 = dt * a
        
        a = cf.forceDueFields(v[i] + dv1/2,B,E,q)/(rmass)
        dv2 = dt * a
        
        a = cf.forceDueFields(v[i] + dv2/2,B,E,q)/(rmass)
        dv3 = dt * a
        
        a = cf.forceDueFields(v[i] + dv3,B,E,q)/(rmass)
        dv4 = dt * a
        
        dv = (dv1 + 2*dv2 + 2*dv3 + dv4) / 6
        
        v[i+1] = v[i] + dv
        
        # New position vector
        x[i+1] = x[i] + dt * v[i+1]
        
        k_max = cf.velocity2Energy(v[i],m)
        
        i = i + 1
        step_counter = step_counter + 1
        
        ''' Check for contact with permanent geometries '''
        
        # Quad electrodes, first check if z-position within range
        
        if np.abs(x[i,2]) < qel_z_max:
        
            # Single-quad electrode
            
            if cf.passthroughElementContact(x[i],sqel_rad,sqel_theta):
                
                steps_inside[0],d_matter[0] = \
                    cf.updateInsideMatter(v[i],dt,steps_inside[0],d_matter[0])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,k_max,X0_al,v[i],dt,m)
                    
                if photons_released > 0:
                    
#                    particle_passthrough_text = "Single-Quad Electrode"
                    k = 0
#                    print('particle %d energy SQ: %0.3f'%(particle_row_index,k_max/10**9))
                    
                    while k < photons_released:
            
                        if k_max > k_min: # Nonsense if k_min > k_max
                            v_norm = v[i]/cf.mag(v[i])
                            p, k_max, v[i], photon_count, \
                            gamma, photon_energy = \
                                cf.bremsstrahlung(v[i],m,k_min,k_max,
                                                  i,photon_count,
                                                  min_detectable_energy)
                            photon_proc[photon_count[0]-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          v_norm[0],v_norm[1],v_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                            k = k + 1
            
            # Double-quad electrode
            
            if cf.passthroughElementContact(x[i],dqel_rad,dqel_theta):
                
                steps_inside[1],d_matter[1] = \
                    cf.updateInsideMatter(v[i],dt,steps_inside[1],d_matter[1])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,k_max,X0_al,v[i],dt,m)
                    
                if photons_released > 0:
                    
#                    particle_passthrough_text = "Double-Quad Electrode"
                    k = 0
#                    print('particle %d energy DQ: %0.3f'%(particle_row_index,k_max/10**9))
                    
                    while k < photons_released:
            
                        if k_max > k_min: # Nonsense if k_min > k_max
                            v_norm = v[i]/cf.mag(v[i])
                            p, k_max, v[i], photon_count, \
                            gamma, photon_energy = \
                                cf.bremsstrahlung(v[i],m,k_min,k_max,
                                                  i,photon_count,
                                                  min_detectable_energy)
                            photon_proc[photon_count[0]-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          v_norm[0],v_norm[1],v_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                            k = k + 1
        
        # Side support plate
        
        if cf.passthroughElementContact(x[i],sp_rad,sp_theta):
                
                steps_inside[2],d_matter[2] = \
                    cf.updateInsideMatter(v[i],dt,steps_inside[2],d_matter[2])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,k_max,X0_al,v[i],dt,m)
                    
                if photons_released > 0:
                    
#                    particle_passthrough_text = "Standoff Plate"
                    k = 0
#                    print('particle %d energy SP: %0.3f'%(particle_row_index,k_max/10**9))
                    
                    while k < photons_released:
            
                        if k_max > k_min: # Nonsense if k_min > k_max
                            v_norm = v[i]/cf.mag(v[i])
                            p, k_max, v[i], photon_count, \
                            gamma, photon_energy = \
                                cf.bremsstrahlung(v[i],m,k_min,k_max,
                                                  i,photon_count,
                                                  min_detectable_energy)
                            photon_proc[photon_count[0]-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          v_norm[0],v_norm[1],v_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                            k = k + 1
        
        # High-voltage standoff
        
        if x[i,2] < so_z_max:
        
            if cf.passthroughElementContact(x[i],so_rad,so_theta):
                
                steps_inside[3],d_matter[3] = \
                    cf.updateInsideMatter(v[i],dt,steps_inside[3],d_matter[3])
                
                photons_released = \
                    cf.isPhotonReleased(k_min,k_max,X0_ma,v[i],dt,m)
                    
                if photons_released > 0:
                    
#                    particle_passthrough_text = "HV Standoff"
                    k = 0
#                    print('particle %d energy HVS: %0.3f'%(particle_row_index,k_max/10**9))
                    
                    while k < photons_released:
            
                        if k_max > k_min: # Nonsense if k_min > k_max
                            v_norm = v[i]/cf.mag(v[i])
                            p, k_max, v[i], photon_count, \
                            gamma, photon_energy = \
                                cf.bremsstrahlung(v[i],m,k_min,k_max,
                                                  i,photon_count,
                                                  min_detectable_energy)
                            photon_proc[photon_count[0]-1] = \
                                np.array([x[i,0],x[i,1],x[i,2],
                                          v_norm[0],v_norm[1],v_norm[2],
                                          photon_energy,particle_row_index,
                                          0,0,step_counter])
                            k = k + 1
            
        # Break if particle energy below detectability/10
        
        cur_energy = cf.velocity2Energy(v[i],m)
        
        if cur_energy <= k_min:
            kill_event_text = "Energy Very Low"
            break
                
        photons_released = 0 # Reset for next step in particle motion
        
        if cf.isInSQuad(x[i],sqel_theta,R) or \
            cf.isInDQuad(x[i],dqel_theta,R):
                loc = "In"
                
        else:
            loc = "Out"
        
        # Get the electric field based on position
#        E = cf.getElectricField(x[i],B,R,n,loc)
        E = np.array([0,0,0])
        
        # New force vector
        F = cf.forceDueFields(v[i],B,E,q)\
        
        # Calorimeter front
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_theta):
            kill_event_text = "Calorimeter Contact"
            break
        
        # Calorimeter top
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_box_theta):
            kill_event_text = "Calorimeter Top Contact"
            break
                        
        # Inner radius limit
        
        if cf.innerLimit(x[i],R_i):
            kill_event_text = "Inner Limit Reached"
            break
        
        # Outer radius limit
        
        if cf.outerLimit(x[i],R + 0.2):
            kill_event_text = "Into the Iron"
            break
        
    print(kill_event_text)
        
    p_init_mag = cf.mag(p_init)
    p_end_mag = cf.mag(cf.beta2Momentum(v[i]/c,m))
    particle_pos[particle_row_index] = np.copy(x)
    charge = particle_proc[particle_row_index,6]
        
    particle_matrix[particle_row_index + 1] = np.array(
                         [particle_row_index,i,
                          '%s'%kill_event_text,charge,
                          '%0.3f'%(x[0,0]*10**3),
                          '%0.3f'%(x[0,1]*10**3),
                          '%0.3f'%(x[0,2]*10**3),
                          '%0.7f'%(p_init_mag/10**9),
                          '%0.7f'%(p_end_mag/10**9),
                          '%0.7f'%((p_end_mag - p_init_mag)/10**9),
                          steps_inside[0],d_matter[0]*100,
                          steps_inside[1],d_matter[1]*100,
                          steps_inside[2],d_matter[2]*100,
                          steps_inside[3],d_matter[3]*100,
                          photon_count[0],
                          photon_count[1],'%5e'%(step_counter*dt)])
    particle_count = particle_count + 1
    particle_proc[particle_row_index,7] = 1
                            
    return particle_pos,particle_matrix,particle_proc,photon_count, \
            photon_proc