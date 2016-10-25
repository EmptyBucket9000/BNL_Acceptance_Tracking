# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
import callable_functions as cf
import random
import csv

#==============================================================================
# Two coordinate systems are used, global and local:
#
#   Global is the entire ring in 3-dimensions cartesian coordinates where the
#    origin is the center of the ring and the z-direction indicates 'up' and
#   'down', from the perspective of a person standing in the ring, and is
#   parallel to the y-direction in the local coordinate system.
#
#   Local is a 2-dimensional x,y cartesian coordinate where the origin is the
#   beam centroid, the negative x-direction points towards the center of the
#   ring and the y-direction is parallel to the z-direction in the global
#   system. Local is primarily used in determining positron initial conditions,
#   i.e. the conditions describing the muon at decay
#
# Variable names formatted 'm_*' are for muons, all others are for positrons
#
# The muon density function is Gaussian, following:
#   (np.sqrt(2*np.pi*sigma**2)**(-1))*np.exp(-((x-xbar)**2)/(2*sigma**2))
#
# Magic momentum is set to 3.09435 GeV/c giving a magnetic field of 1.4513 T
# for a radius of 7.112 m.
#==============================================================================

#==============================================================================
#==============================================================================
# #   'main' function
#==============================================================================
#==============================================================================

def main():

    ''' Begin editable variables '''

    # Output
    
    print_text = 1              # Set to 1 to print output text, 0 for none
    make_plots = 1              # Set to 1 to display plots
    save_plots = 0              # Set to 1 to save plots as images
    save_output = 0             # Set to 1 to save data output to csv
    
    # Idealizations in local coordinates. Set to 1 in appropriate array index
    # if an idealization is desired, i.e. set first element in array to 1 for
    # ideal x, set the second element to 1 for ideal y. Ideal gives no
    # oscillations for that variable.

    m_xbar_ideal = np.array([1,1])      # Ideal mean position
    m_sigma_ideal = np.array([1,1])     # Ideal beam width
        
    m_theta_set = 1                     # 1, use m_theta below, 0 random
    
    N = 1                              # Number of muons in beam
    steps = 10**6                       # Nnumber of steps for integration
    dt = 10**-12                        # Timestep for integration
    p_range = np.array([.2,3.09435])*10**9 # (eV/c) Possible positron momentums
#    p_range = np.array([0.6,0.6])*10**9
    p_n = 1000                  # Number of possilbe values in range    
    
    if m_theta_set == 1:
        m_theta = 1.5*np.pi/8   # (rad) Muon azimuth. position in global coords
    elif m_theta_set == 0:
        m_theta_min = 0
        m_theta_max = 2*np.pi
    
    # Amplitude of mean position oscillation
    m_xbar_amp = np.array([2,0])*10**-3 # (m)
    
    # Position around which xbar oscillates
    m_xbar_0 = np.array([0,0])*10**-3
    
    # Position phase offset
    m_xphi = np.array([0,0])  # (rad)
    
    # Beam width phase offset
    m_sigmaphi = np.array([0,0]) # (rad)
    
    # Amplitude of beam width oscillation
    m_sigma_amp = np.array([11.156/2,0])*10**-3
    
    # Value around which x_sigma oscillates
    m_sigma_0 = np.array([15,0])*10**-3
    
    ''' Permanent constants '''
    
    q = 1                                   # (e) Positron charge
#    X0_al = 0.08897                         # (m) Radiation length of aluminum
#    X0_ma = 0.0965                          # (m) Radiation length of macor
    X0_al = 0.001
    X0_ma = 0.001
    c = 2.99792458*10**8                    # (m/s) Speed of light
    m = 0.510999*10**6                      # (eV/c**2) Positron mass
#    m_m = 105.65837*10**6                  # (eV/c**2) Muon mass
    B = np.array([0,0,1.4513])              # (T) Magnetic field
    min_detectable_energy = 0.2*10**9       # (eV) Minimum energy detectable
    
    ''' Other variables '''
    
    particle_number = 1         # Used as a counter for multiple particles
    
    # Output for each bremsstrahlung event
    out_brem = np.array([["Particle #","Bremsstrahlung Location",
                         "Photon Energy (eV)","Photon Kill Event",
                         "New Particle Momentum (GeV/c)"]])
                         
    # Output for each particle
    out_particle = np.array([["Particle #","Kill Event",
                              "Starting Local x-Position (mm)",
                              "Mean x (mm)","x Beam Width (mm)",
                              "Starting Local y-Position (mm)",
                              "Mean y (mm)","y Beam Width (mm)",
                              "Starting Global Theta (deg)",
                              "Starting Momentum (GeV/c)",
                              "Ending Momentum (GeV/c)",
                              "Delta Momentum (GeV/c)",
                              "Steps Inside Matter",
                              "Distance Inside matter (cm)",
                              "Total # of Photons Released",
                              "# of Detectable Photons Released"]])

#==============================================================================
# Setting up permanent geometries
#==============================================================================

    R = 7.112                   # (m) Radius of the ring
    R_i = R - 0.5525            # (m) Inner radius limit for tracking
    
    ''' Calorimeters '''
    
    cal_len = 0.36576           # (m) Calorimeter length
    #cal_phi = 0                # (rad) Angle from radial line (Unused)
    cal_depth = 0.4572          # (m) Depth of calorimeter
    
    # (rad) Location as a function of theta. 0.001 subtracted from
    # 'cal_theta_end' to prevent particles from being counted as contact with
    # the calorimeter by contacting the 'top' of the calorimeter.
    
    cal_theta_start = np.linspace(0,2*np.pi-np.pi/12,24) + np.pi/12
#    cal_theta_start = np.array([0])
    cal_box_theta_end = cal_theta_start - cal_depth/(R_i+cal_len)
    cal_theta_end = cal_theta_start - .001
    cal_theta = np.column_stack((cal_theta_start,cal_theta_end))
    cal_box_theta = np.column_stack((cal_theta_start -0.001,cal_box_theta_end))
    
    # Radial position [r_min,r_max]
    cal_rad = np.array([R_i,R_i+cal_len]) # (m) 
    
    ''' HV Standoff '''
    
    so_z_max = 3.8*10**-3       # (m) Top of standoff
    so_length = 7.62*10**-3     # (m) Width of standoff
    so_depth = 27.5*10**-3      # (m) Standoff depth
    so_rad_start = 50.5*10**-3  # (m) Starting distance in from R
    so_rad = np.array([R - (so_rad_start + so_depth), R - so_rad_start])
    so_theta_start_base = np.array([2.48,8.48,13.48,15.48,21.48,27.23])
    so_theta_start = np.array([2.48,8.48,13.48,15.48,21.48,27.23])
        
    i = 1
    while i < 12:
        so_theta_start = np.concatenate((
            so_theta_start, so_theta_start_base + i*30
        ))
        i = i + 1
    so_theta_start = so_theta_start*np.pi/180
    so_theta_end = so_theta_start - so_length/(R-so_rad_start)
    
    # Theta position [theta_min, theta_max]
    so_theta = np.column_stack((so_theta_start,so_theta_end))
    
    ''' Support plates '''
    
    sp_length = 33*10**-3        # (m) Support plate width
    sp_depth = 1.7*10**-3        # (m) Support plate thickness
    sp_rad_start = 78*10**-3     # (m) Distance in from R
    sp_rad = np.array([R - (sp_rad_start + sp_depth), R - sp_rad_start]) # (m)
    sp_theta_start_base = np.array([2.48,8.48,13.48,15.48,21.48,27.23])+0.1
    sp_theta_start = np.array([2.48,8.48,13.48,15.48,21.48,27.23])+0.1
        
    i = 1
    while i < 12:
        sp_theta_start = np.concatenate((
            sp_theta_start, sp_theta_start_base + i*30
        ))
        i = i + 1
    sp_theta_start = sp_theta_start*np.pi/180
    sp_theta_end = sp_theta_start - sp_length/(R-sp_rad_start)
    
    # Theta position [theta_min, theta_max]
    sp_theta = np.column_stack((sp_theta_start,sp_theta_end))
    
    ''' Variables for both quads '''
    
    qel_z_max = 23.5*10**-3     # (m) Top of electrodes
    qel_depth = 0.5*10**-3      # (m) Electrode thickness
    qel_rad_start = 50*10**-3   # (m) Starting distance in from R
    
    ''' Single-quad electrodes (without edge curls) '''
    
    sqel_rad = np.array([R - (qel_rad_start + qel_depth),R - qel_rad_start])
    sqel_theta_base = 90 - np.array([31.89,44.89])-4
    sqel_theta = 90 - np.array([31.89,44.89])-4
    
    i = 1
    while i < 4:
        sqel_theta = np.row_stack((
            sqel_theta, sqel_theta_base + i*90
        ))
        i = i + 1
        
    sqel_theta = sqel_theta*np.pi/180
    
    ''' Double-quad electrodes (without edge curls) '''
    
    dqel_rad = np.array([R - (qel_rad_start + qel_depth),R - qel_rad_start])
    dqel_theta_base = 90 - np.array([48.89,74.89])-4
    dqel_theta = 90 - np.array([48.89,74.89])-4
    
    i = 1
    while i < 4:
        dqel_theta = np.row_stack((
            dqel_theta, dqel_theta_base + i*90
        ))
        i = i + 1
        
    dqel_theta = dqel_theta*np.pi/180  
    
    # Pack up the geometry variables for easier transfer to the run() function
    geo_pack = np.array([cal_theta,cal_theta_start,cal_rad,
                         cal_box_theta,cal_box_theta_end,
                         so_z_max,so_rad,so_theta,so_theta_start,so_theta_end,
                         sp_rad,sp_theta,sp_theta_start,sp_theta_end,
                         qel_z_max,
                         sqel_rad,sqel_theta,
                         dqel_rad,dqel_theta,
                         R,R_i])
    
#==============================================================================
#   Run the code!!
#==============================================================================
                             
    while particle_number <= N:
    
#==============================================================================
#   Set all initial local muon variables
#==============================================================================
    
        '''Unless otherwise stated, 3-element arrays are of the form x,y,z.'''
        
        # 'fit' is the type of function that will fit the muon distribution
        # density. Curreantly 'Gaussian' is the only option. ([x,y,z])
        m_xfit = np.array(["Gaussian","",""])
        
        # Number of possible muon positions within range, used for building the
        # muon distribution arrays
        m_xnum = np.array([1000,1000])
        m_sigmanum = np.array([1000,1000])
    
        if m_theta_set == 0:
            m_theta = (m_theta_max - m_theta_min)*random.random() + m_theta_min
        
        # (m) Mean muon position (centroid position)
        m_xbar_x = cf.getParticleXBar(m_xbar_amp[0],m_xbar_ideal[0],m_xnum[0],
                                      m_xbar_0[0],'x')
        m_xbar_y = cf.getParticleXBar(m_xbar_amp[1],m_xbar_ideal[1],m_xnum[1],
                                      m_xbar_0[1],'y')
        
        m_xbar = np.array([m_xbar_x,m_xbar_y])
    
        ''' Minimum/maximum muon position in beam '''
    
        m_xmin = np.zeros((2))                # (m)
        m_xmax = np.zeros((2))                # (m)
        
        # Beam distribution width
        
        m_sigma = cf.getParticleSigma(m_sigma_amp,m_sigma_ideal,
                                         m_sigma_0,m_sigmanum)
        
        if m_sigma_ideal[0] == 1:
            m_xmin[0] = m_xbar[0]             # (m)
            m_xmax[0] = m_xbar[0]             # (m)
        else:
            m_xmin[0] = m_xbar[0] - 41*10**-3 # (m)
            m_xmax[0] = m_xbar[0] + 41*10**-3 # (m)
        
        # Muon x, y, z prime values at decay
        m_xprime = np.array([0,0,0])
    
        ''' Set the muon position array '''
    
        # x-position
        m_pos_x = cf.getXParticlePositions(
            m_xbar[0],m_sigma[0],m_xmin[0],m_xmax[0],m_xnum[0],m_xfit[0]
        )
        
        # y-position
        m_pos_y = cf.getYParticlePositions(
            m_xbar[1],m_sigma[1],m_xmin[1],m_xmax[1],m_xnum[1],m_xfit[1]
        )
        
        # Full position vector
        m_x = np.array([m_pos_x,m_pos_y])
    
        # Run the code where everything happens
        out_brem,out_particle = \
            run(geo_pack,m_x,m,c,out_brem,min_detectable_energy,out_particle,
                print_text,make_plots,save_plots,save_output,N,steps,dt,
                m_theta,p_range,p_n,particle_number,q,B,X0_al,X0_ma,m_xbar,
                m_sigma,'original')
                
        print("Percent complete: %0.1f%%"%(particle_number*(100/N)), end="\r")
                
        particle_number = particle_number + 1
        
    # Save data to file
        
    if save_output == 1:
        path = "Output/output_particles.csv"%()
        with open(path, "w", newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            for row in out_particle:
                writer.writerow(row)
                
        path = "Output/output_brem.csv"%()
        with open(path, "w", newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            for row in out_brem:
                writer.writerow(row)

#==============================================================================
#==============================================================================
# #   Begin 'run' function, where everything happens
#==============================================================================
#==============================================================================

def run(geo_pack,m_x,m,c,out_brem,min_detectable_energy,out_particle,
        print_text,make_plots,save_plots,save_output,N,steps,dt,m_theta,
        p_range,p_n,particle_number,q,B,X0_al,X0_ma,m_xbar,m_sigma,part_type):
    
#==============================================================================
#   Initialization and setting variables
#==============================================================================

    # Unpack 'geo_pack'    
    
    cal_theta = geo_pack[0]
    cal_theta_start = geo_pack[1]
    cal_rad = geo_pack[2]
    cal_box_theta = geo_pack[3]
    cal_box_theta_end = geo_pack[4]
    so_z_max = geo_pack[5]    
    so_rad = geo_pack[6]
    so_theta = geo_pack[7]
    so_theta_start = geo_pack[8]    
    so_theta_end = geo_pack[9]
    sp_rad = geo_pack[10]
    sp_theta = geo_pack[11]
    sp_theta_start = geo_pack[12]
    sp_theta_end = geo_pack[13]
    qel_z_max = geo_pack[14]    
    sqel_rad = geo_pack[15]
    sqel_theta = geo_pack[16]
    dqel_rad = geo_pack[17]
    dqel_theta = geo_pack[18]
    R = geo_pack[19]
    R_i = geo_pack[20]
    
#==============================================================================
#   If particle is directly from decay, skip if created by pair production
#==============================================================================
    
    if part_type == 'original':
    
        # Momentum array for a single muon away from the magic momentum
        m_p = cf.getMuonMomentumAtDecay() # (eV/c)
        
        '''Use local muon variables to get inital global positron variables'''
        
        # Radial distance from center of ring
        r = R + m_x[0]                          # (m)
        
        # Set x,y,z positions
        
        x = np.zeros((steps,3))                 # Initialize positon array
        x[0,0] = r*np.cos(m_theta)              # (m) x-position
        x[0,1] = r*np.sin(m_theta)              # (m) y-position
        x[0,2] = m_x[1]                         # (m) z-position
        
        # Positron momentum at decay
        v = np.zeros((steps,3))                 # (m/s) Initialize velocity
        
        # (eV/c) Momentum
        p = p_initial = cf.getPositronMomentumAtDecay(x,m_theta,m_p,m,p_range,
                                                      p_n)
        
        beta = cf.momentum2Beta(p,m)            # () Relativistic beta
        v[0] = beta*c              # (m/s) Initial velcoty  
    
#==============================================================================
#   Initialize misc. variables
#==============================================================================
    
    # Magnitude of relativistic beta vector
    
    beta_mag = np.zeros(1)                  # ()
    beta_mag[0] = cf.mag(beta)
    
    # Relativistic gamma
    
    gamma = np.zeros(1)                     # ()
    gamma[0] = cf.beta2Gamma(beta_mag[0])
    
    # Electric field
    
    E = np.zeros((steps,3))                 # (V/m) Initialize electric field
    E[0] = cf.getElectricField(x[0])        # Set initial electric field values
    
    # Force vector
    
    F = np.zeros((steps,3))                 # (N) Initialze force array
    F[0] = cf.forceDueFields(v[0],B,E[0],q) # Set initial force due fields
    
    # Counter for photons released by Bremsstrahlung
    photons_count = 0
    
    # Counter for detectable photons released by Bremsstrahlung   
    photons_count_min = 0
    
    # Possible energies of Bremsstrahlung photons
    k_min = 10**6                           # (eV)
    k_max = cf.momentum2Energy(p_initial,m) # (eV)
    
    # Array of indices for Bremsstrahlung points
    brem_array = np.zeros(0, dtype=np.int)
    photon_energies = np.zeros((0,2))
    
    # For checking if bremsstrahlung happened or not
    photons_released = 0
            
    photon_steps = 10**3
    photon_dt = 10**-11
    
    # Photon position array
    photon_x_base = np.zeros((1,photon_steps,3))
    photon_x = np.zeros((0,photon_steps,3))
    
    # Counter for number of steps a positron is inside matter
    steps_inside = 0
    
    # How far the positron has traveled in matter
    d_matter = 0                            # (m)
    
#==============================================================================
#   Tracking by Runga-Kutta 4th
#==============================================================================
    
    # Event text
    kill_event_text = "Unknown failure" # In case nothing happens
    
    # Loop counter
    i = 0
    
    ''' RK4 work '''
    
    while i < steps - 1:
        
        # Relativistic mass to find the acceleration from the force
        rmass = m*gamma/c**2   
        
        a = F[i]/(rmass)
        dv1 = dt * a
        
        a = cf.forceDueFields(v[i] + dv1/2,B,E[i],q)/(rmass)
        dv2 = dt * a
        
        a = cf.forceDueFields(v[i] + dv2/2,B,E[i],q)/(rmass)
        dv3 = dt * a
        
        a = cf.forceDueFields(v[i] + dv3,B,E[i],q)/(rmass)
        dv4 = dt * a
        
        # New velocity vector        
        v[i+1] = v[i] + (dv1 + 2*dv2 + 2*dv3 + dv4) / 6
        
        # New position vector
        x[i+1] = x[i] + dt * v[i+1]
        
        # Get the electric field based on position
        E[i+1] = cf.getElectricField(x[i+1])
        
        # New force vector
        F[i+1] = cf.forceDueFields(v[i+1],B,E[i+1],q)
        
        i = i + 1
        
        ''' Check for contact with permanent geometries '''
        
        # Quad electrodes, first check if z-position within range
        
        if np.abs(x[i,2]) < qel_z_max:
        
            # Single-quad electrode
            
            if cf.passthroughElementContact(x[i],sqel_rad,sqel_theta):
                
                steps_inside,d_matter = \
                    cf.updateInsideMatter(v[i],dt,steps_inside,d_matter)
                
                if k_max > k_min: # Nonsense if k_min > k_max
                
                    photons_released = \
                        cf.isPhotonReleased(k_min,k_max,X0_al,v[i],dt,m)
                    
                    if photons_released > 0:
                        
                        text = "Single-Quad Electrode"
                        k = 0
                        
                        while k < photons_released:
                        
                            p, k_max, v[i], brem_array, photons_count, \
                            gamma, photons_count_min, \
                            photon_energies = \
                                cf.bremsstrahlung(v[i],m,k_min,k_max,
                                                  brem_array,i,photons_count,
                                                  min_detectable_energy,
                                                  photons_count_min,
                                                  photon_energies)
                            k = k + 1
            
            # Double-quad electrode
            
            if cf.passthroughElementContact(x[i],dqel_rad,dqel_theta):
                
                steps_inside,d_matter = \
                    cf.updateInsideMatter(v[i],dt,steps_inside,d_matter)
                
                if k_max > k_min: # Nonsense if k_min > k_max
                
                    photons_released = \
                        cf.isPhotonReleased(k_min,k_max,X0_al,v[i],dt,m)
                    
                    if photons_released > 0:
                        
                        text = "Double-Quad Electrode"
                        k = 0

                        while k < photons_released:
                        
                            p, k_max, v[i], brem_array, photons_count, \
                            gamma, photons_count_min, \
                            photon_energies = \
                                cf.bremsstrahlung(v[i],m,k_min,k_max,
                                                  brem_array,i,photons_count,
                                                  min_detectable_energy,
                                                  photons_count_min,
                                                  photon_energies)
                            k = k + 1
        
        # Side support plate
        
        if cf.passthroughElementContact(x[i],sp_rad,sp_theta):
                
            steps_inside,d_matter = \
                cf.updateInsideMatter(v[i],dt,steps_inside,d_matter)
            
            if k_max > k_min: # Nonsense if k_min > k_max
                
                photons_released = \
                    cf.isPhotonReleased(k_min,k_max,X0_al,v[i],dt,m)
                
                if photons_released > 0:
                    
                    text = "Side Support Plate"
                    k = 0

                    while k < photons_released:
                    
                        p, k_max, v[i], brem_array, photons_count, \
                        gamma, photons_count_min, \
                            photon_energies = \
                            cf.bremsstrahlung(v[i],m,k_min,k_max,
                                              brem_array,i,photons_count,
                                                  min_detectable_energy,
                                                  photons_count_min,
                                                  photon_energies)
                        k = k + 1
        
        # High-voltage standoff
        
        if x[i,2] < so_z_max:
        
            if cf.passthroughElementContact(x[i],so_rad,so_theta):
                
                steps_inside,d_matter = \
                    cf.updateInsideMatter(v[i],dt,steps_inside,d_matter)
                
                if k_max > k_min: # Nonsense if k_min > k_max
                
                    photons_released = \
                        cf.isPhotonReleased(k_min,k_max,X0_al,v[i],dt,m)
                    
                    if photons_released > 0:
                        
                        text = "HV Standoff"
                        k = 0
    
                        while k < photons_released:
                        
                            p, k_max, v[i], brem_array, photons_count, \
                            gamma, photons_count_min, \
                            photon_energies = \
                                cf.bremsstrahlung(v[i],m,k_min,k_max,
                                                  brem_array,i,photons_count,
                                                  min_detectable_energy,
                                                  photons_count_min,
                                                  photon_energies)
                            k = k + 1
                        
        # Adds the bremsstrahlung event to the output array
                        
        if photons_released > 0:
                     
            photon_energies = photon_energies[~np.all(photon_energies == 0,
                                                      axis=1)]
            
            photon_kill_event_text = "Low Energy"
            
            # Get the normalized velocity vector           
            v_norm = v[i-1] / cf.mag(v[i-1])
            v_photon = v_norm*c
            
            k = 0
            
            while k < photons_released:
            
                le = len(brem_array)
                
#                print(photon_energies)
                
                if photon_energies[le-k-1,0] > min_detectable_energy:
                    photon_x = np.append(photon_x,photon_x_base,axis=0)
                
                    # Set index of newest photon_x page
                    cur_el = np.shape(photon_x)[0]-1
                    
                    # Get position where the photon was created
                    photon_x[cur_el,0] = np.array([x[brem_array[le-k-1],0],
                                              x[brem_array[le-k-1],1],
                                              x[brem_array[le-k-1],2]])
                
                    l = 0
                    
                    while l < photon_steps-1:
                        
                        photon_x[cur_el,l+1] = \
                            photon_x[cur_el,l] + v_photon*photon_dt
                        
                        l = l + 1
                
                        # Calorimeter front contact
                        
                        if cf.noPassthroughElementContact(photon_x[cur_el,l],
                                                          cal_rad,cal_theta):
                            photon_kill_event_text = \
                                "Photon Calorimeter Contact"
                            break
                
                        # Calorimeter top
                        
                        if cf.noPassthroughElementContact(photon_x[cur_el,l],
                                                          cal_rad,
                                                          cal_box_theta):
                            photon_kill_event_text = \
                                "Photon Calorimeter Top Contact"
                            break
                        
                        # Outer radius limit
                        
                        if cf.outerLimit(photon_x[cur_el,l],R + 0.1):
                            photon_kill_event_text = "Photon Heading Out"
                            break
        
                        if photon_x[cur_el,l,2] < so_z_max:
                        
                            if cf.passthroughElementContact(photon_x[cur_el,l],
                                                            so_rad,so_theta):
                                print('Photon Standoff')
            
                        if cf.passthroughElementContact(photon_x[cur_el,l],
                                                        dqel_rad,dqel_theta):
                            print('Photon Double Quad')
#                            cf.pairProduction()
                
                out_brem = np.append(out_brem,
                                     [[particle_number,'%s'%text,
                                       '%e'%photon_energies[le-k-1,0],
                                       '%s'%photon_kill_event_text,
                                       '%e'%(cf.mag(p)/10**9)]], 
                                     axis=0
                                    )
                                    
                k = k + 1
                
            photons_released = 0 # Reset for next step in positron motion
        
        # Calorimeter front
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_theta):
            kill_event_text = "Calorimeter Contact"
            break
        
        # Calorimeter top
        
        if cf.noPassthroughElementContact(x[i],cal_rad,cal_box_theta):
            kill_event_text = "Calorimeter Top Contact"
            break
            
        # Break if particle energy below detectability
        
        cur_energy = cf.velocity2Energy(v[i],m)
        
        if cur_energy < min_detectable_energy:
            kill_event_text = "Energy Below Minimum"
            break
                        
        # Outer radius limit
        
        if cf.innerLimit(x[i],R_i):
            kill_event_text = "Inner Limit Reached"
            break
        
#==============================================================================
#   Cleanup and some useful variable creation
#==============================================================================
        
    # Removes unused parts of the matricies (the unchanged zeros)
        
    E = E[0:i+1:1]
    x = x[0:i+1:1]
    v = v[0:i+1:1]
    
    # Get final momentum and beta
    
    beta_final = v[i]/c
    p = cf.beta2Momentum(beta_final,m)
    
    if print_text == 1:
        
        # Print what happened to the positron
        np.set_printoptions(formatter={'float': lambda x: format(x, '6.3E')})
        print('Starting Momentum: %0.6f GeV/c'%(cf.mag(p_initial)/10**9))
        print('Kill Event: %s'%kill_event_text)
        print('Final Momentum: %0.6f GeV/c'%(cf.mag(p)/(10**9)))
        print('Photons Released: %d'%photons_count)
        print('Photon Energies:')
        print(photon_energies[:,0])
        
#==============================================================================
#   Plotting and other output
#==============================================================================

    if make_plots == 1:
    
        ''' Setting useful variables '''
        
        # Counter used for setting multiple figures
        n = 0
        
        # Used in adding a text box to the plots
        props = dict(boxstyle='square', facecolor='wheat', alpha=0.8)
        
        ''' Plotting '''
        
        # Figure
        plt.figure(n)
        n = n + 1
        
        ax = plt.subplot(1,1,1)
        ax.plot(x[:,0],x[:,1], label='Particle Track')
        plt.xlabel('x-position (m)')
        plt.ylabel('y-position (m)')  
        plt.title('Positron Position')
        
        ''' Add geometries to the figure '''
        
        # Add inner radius
        
        xt = np.linspace(-R_i,R_i,steps)    # (m)
        ax.plot(xt,np.sqrt((R_i)**2 - xt**2),'-.r')
        ax.plot(xt,-np.sqrt(R_i**2 - xt**2),'-.r')
        
        # Add support plates
        
        count = len(sp_theta_start)
        
        k = 0
        while k < count:
            ax.plot(
                [sp_rad[0]*np.cos(sp_theta_start[k]),
                 sp_rad[0]*np.cos(sp_theta_end[k])],
                [sp_rad[0]*np.sin(sp_theta_start[k]),
                 sp_rad[0]*np.sin(sp_theta_end[k])],
                'k'
            )
            ax.plot(
                [sp_rad[1]*np.cos(sp_theta_start[k]),
                 sp_rad[1]*np.cos(sp_theta_end[k])],
                [sp_rad[1]*np.sin(sp_theta_start[k]),
                 sp_rad[1]*np.sin(sp_theta_end[k])],
                'k'
            )
            
            k = k + 1
        
        # Add HV Standoff
        
        count = len(so_theta_start)
        k = 0
        
        while k < count:
            ax.plot(
                [so_rad[0]*np.cos(so_theta_start[k]),
                 so_rad[1]*np.cos(so_theta_start[k])],
                [so_rad[0]*np.sin(so_theta_start[k]),
                 so_rad[1]*np.sin(so_theta_start[k])],
                'k'
            )
            ax.plot(
                [so_rad[0]*np.cos(so_theta_end[k]),
                 so_rad[1]*np.cos(so_theta_end[k])],
                [so_rad[0]*np.sin(so_theta_end[k]),
                 so_rad[1]*np.sin(so_theta_end[k])],
                'k'
            )
            
            k = k + 1
        
        # Add single-quad electrodes      
        
        count = len(sqel_theta)
        k = 0
        M = 100
        
        # Plot those for y > 0
        
        while k < count/2:
            
            xt = np.linspace(
                sqel_rad[0]*np.cos(sqel_theta[k,0]),
                sqel_rad[0]*np.cos(sqel_theta[k,1]),
                M
            )
            
            ax.plot(xt,np.sqrt(sqel_rad[0]**2 - xt**2),'k')
            
            xt = np.linspace(
                sqel_rad[1]*np.cos(sqel_theta[k,0]),
                sqel_rad[1]*np.cos(sqel_theta[k,1]),
                M
            )
            ax.plot(xt,np.sqrt(sqel_rad[1]**2 - xt**2),'k')
            
            k = k + 1
            
        # Plot those for y < 0
            
        while  k < count:
            
            xt = np.linspace(
                sqel_rad[0]*np.cos(sqel_theta[k,0]),
                sqel_rad[0]*np.cos(sqel_theta[k,1]),
                M
            )
            
            ax.plot(xt,-np.sqrt(sqel_rad[0]**2 - xt**2),'k')
            
            xt = np.linspace(
                sqel_rad[1]*np.cos(sqel_theta[k,0]),
                sqel_rad[1]*np.cos(sqel_theta[k,1]),
                M
            )
            ax.plot(xt,-np.sqrt(sqel_rad[1]**2 - xt**2),'k')
            
            k = k + 1
        
        # Add double-quad electrodes      
        
        count = len(dqel_theta)
        k = 0
        M = 100
        
        # Plot those for y > 0
        
        while k < count/2:
            
            xt = np.linspace(
                dqel_rad[0]*np.cos(dqel_theta[k,0]),
                dqel_rad[0]*np.cos(dqel_theta[k,1]),
                M
            )            
            ax.plot(xt,np.sqrt(float(dqel_rad[0])**2 - xt**2),'k')
            
            xt = np.linspace(
                dqel_rad[1]*np.cos(dqel_theta[k,0]),
                dqel_rad[1]*np.cos(dqel_theta[k,1]),
                M
            )
            ax.plot(xt,np.sqrt(dqel_rad[1]**2 - xt**2),'k')
            
            k = k + 1
            
        # Plot those for y < 0
            
        while  k < count:
            
            xt = np.linspace(
                dqel_rad[0]*np.cos(dqel_theta[k,0]),
                dqel_rad[0]*np.cos(dqel_theta[k,1]),
                M
            )            
            ax.plot(xt,-np.sqrt(dqel_rad[0]**2 - xt**2),'k')
            
            xt = np.linspace(
                dqel_rad[1]*np.cos(dqel_theta[k,0]),
                dqel_rad[1]*np.cos(dqel_theta[k,1]),
                M
            )
            ax.plot(xt,-np.sqrt(dqel_rad[1]**2 - xt**2),'k')
            
            k = k + 1
        
        # Add calorimeters
        
        count = len(cal_theta_start)        
        
        k = 0
        while k < count:
            ax.plot(
                [cal_rad[0]*np.cos(cal_theta_start[k]),
                cal_rad[1]*np.cos(cal_theta_start[k])],
                [cal_rad[0]*np.sin(cal_theta_start[k]),
                cal_rad[1]*np.sin(cal_theta_start[k])],
                'k-'
            )
                
            ax.plot(
                [cal_rad[0]*np.cos(cal_box_theta_end[k]),
                cal_rad[1]*np.cos(cal_box_theta_end[k])],
                [cal_rad[0]*np.sin(cal_box_theta_end[k]),
                cal_rad[1]*np.sin(cal_box_theta_end[k])],
                'k-'
            )
             
            k = k + 1
            
        # Add Bremsstrahlung photons if they were created
            
        if photons_count_min > 0:
            
            k = 0
            
            while k < photons_count_min:
                
                temp_x = photon_x[k]
                temp_x = temp_x[~np.all(temp_x == 0, axis=1)]
#                print(temp_x)
#                L = len(temp_x)
#                temp_x = temp_x[0:L:1]
                ax.plot(temp_x[:,0],temp_x[:,1],'r')
                
                k = k + 1
            
        # Add text box

        textstr = ''
        ax.text(0.75, 0.95, textstr, transform=ax.transAxes,
                fontsize=12, verticalalignment='top', bbox=props)

        # Save the plot(s) if save_plots == 1

        if save_plots == 1:
            plt.savefig('Output/Images/plot.png', bbox_inches='tight', dpi=300)

        # Add legend
#        ax.legend()

        # Set axes limits based on min/max particle positions

        plt.axis('equal') # Prevents a skewed look
        plt.xlim(min(x[:,0]-0.2),max(x[:,0]+0.2))
        plt.ylim(min(x[:,1]-1),max(x[:,1]+0.2))
#        plt.xlim(4.5,5)
#        plt.ylim(4.5,5)
        plt.grid()
    
        # Show plot
        plt.show()
        
    out_p_initial = cf.mag(p_initial)/10**9
    out_p_end = cf.mag(p)/10**9
        
    out_particle = np.append(out_particle,
                         [[particle_number,'%s'%kill_event_text,
                           '%0.3f'%(m_x[0]*10**3),'%0.3f'%(m_xbar[0]*10**3),
                            '%0.3f'%(m_sigma[0]*10**3),
                           '%0.3f'%(m_x[1]*10**3),'%0.3f'%(m_xbar[1]*10**3),
                            '%0.3f'%(m_sigma[1]*10**3),
                            m_theta*180/np.pi,
                            '%0.7f'%(out_p_initial),
                            '%0.7f'%(out_p_end),
                            '%0.7f'%(out_p_end - out_p_initial),
                            steps_inside,
                            d_matter*100,
                            photons_count,
                            photons_count_min]], 
                            axis=0
                            )
    
    return out_brem,out_particle 
        
#==============================================================================
#   End 'run' function
#==============================================================================
    
if __name__ == '__main__':
    
    main()