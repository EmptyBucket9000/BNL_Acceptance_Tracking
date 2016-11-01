# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
import callable_functions as cf
import geometry
import particle_tracking as pt
import photon_tracking as pht
import muon_position as mp
import plot_geometries as pg
import csv

"""
See README.md for information.

"""

#==============================================================================
#==============================================================================
# #   'main' function
#==============================================================================
#==============================================================================

def main():

    ''' Begin editable variables '''

    # Output
    
#    print_text = 1              # Set to 1 to print output text, 0 for none
    make_plots = 1              # Set to 1 to display plots
    save_plots = 0              # Set to 1 to save plots as images
    save_output = 1             # Set to 1 to save data output to csv
    
    # Idealizations in local coordinates. Set to 1 in appropriate array index
    # if an idealization is desired, i.e. set first element in array to 1 for
    # ideal x, set the second element to 1 for ideal y. Ideal gives no
    # oscillations for that variable.

    m_xbar_ideal = np.array([0,1])      # Ideal mean position
    m_sigma_ideal = np.array([1,1])     # Ideal beam width
        
    m_theta_set = 1                     # 1, use m_theta below, 0 random
    
    N = 1                              # Number of muons in beam
    steps = 2*10**5                       # Nnumber of steps for integration
    dt = 10**-12                        # Timestep for integration
    
    p_range = np.array([.2,3.09435])*10**9 # (eV/c) Possible particle momentums
#    p_range = np.array([3.09435,3.09435])*10**9
    p_n = 1000                  # Number of possilbe values in range
    
    if m_theta_set == 1:
        m_theta = 2.3*np.pi/8   # (rad) Muon azimuth. position in global coords
    else:
        m_theta = -99
        
    # Min/max muon theta values, only used if m_theta_set == 0
        
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
    
    # Maximum prime values
    xprimemax = np.array([4.5,2.5])
    
    # How far from 0 does the beam reach (m)
    m_xlimit = 41*10**-3
    
    ''' Permanent constants '''
    
    q = 1                                   # (e) Particle charge
    c = 2.99792458*10**8                    # (m/s) Speed of light
    m = 0.510999*10**6                      # (eV/c**2) Particle mass
#    m_m = 105.65837*10**6                  # (eV/c**2) Muon mass
    B = np.array([0,0,1.4513])              # (T) Magnetic field
    
    geo_pack = geometry.geo()               # All the permanent geometries
    
    ''' Other variables '''
    
    muon_number = 0         # Used as a counter for multiple particles
                         
    # Output for each particle
    particle_matrix = np.array([["Particle #","Steps","Kill Event","Charge",
                              "Starting Global x-Position (mm)",
                              "Starting Global y-Position (mm)",
                              "Starting Global z-Position (mm)",
                              "Starting Momentum (GeV/c)",
                              "Ending Momentum (GeV/c)",
                              "Delta Momentum (GeV/c)",
                              "Steps Inside Short Quad",
                              "Distance Inside Short Quad (cm)",
                              "Steps Inside Long Quad",
                              "Distance Inside Long Quad (cm)",
                              "Steps Inside Standoff Plate",
                              "Distance Inside Standoff Plate (cm)",
                              "Steps Inside HV Standoff",
                              "Distance Inside HV Standoff (cm)",
                              "Total # of Photons Released",
                              "# of Detectable Photons Released",
                              "Kill Timestamp"]])
          
    # Output for each photon
    photon_matrix = np.array([["Photon #","Steps","Kill Event",
                               "Starting Global x-Position",
                               "Starting Global y-Position",
                               "Starting Global z-Position",
                               "Energy (GeV)","Steps Inside Matter",
                               "Distance Inside Matter (cm)",
                               "Kill Timestamp"]])
    
#==============================================================================
#   Run the code!!
#==============================================================================
                             
    while muon_number < N:
    
        # Set all initial local muon variables
    
        part_type = 1
    
        m_x,m_sigma,m_theta,m_xbar,m_xprime = \
            mp.muon(m_theta_set,m_theta_min,m_theta_max,m_xbar_amp,m_xbar_0,
                    m_xphi,m_sigmaphi,m_sigma_amp,m_sigma_0,m_xbar_ideal,
                    m_sigma_ideal,m_theta,xprimemax,m_xlimit)
        
    
        # Run the code where everything happens
        particle_matrix = \
            run(geo_pack,m_x,m_sigma,m_theta,m_xbar,m_xprime,m,c,photon_matrix,
                particle_matrix,make_plots,save_plots,save_output,
                N,steps,dt,p_range,p_n,q,B,part_type)
                
        print("Percent complete: %0.1f%%"%(muon_number*(100/N)),
              end="\r")
                
        muon_number = muon_number + 1

#==============================================================================
#==============================================================================
# #   Begin 'run' function, where everything happens
#==============================================================================
#==============================================================================

def run(geo_pack,m_x,m_sigma,m_theta,m_xbar,m_xprime,m,c,photon_matrix,
        particle_matrix,make_plots,save_plots,save_output,N,steps,dt,
        p_range,p_n,q,B,part_type):
    
#==============================================================================
#   Initialization and setting variables
#==============================================================================

    R = geo_pack[19]
    
#==============================================================================
#   If particle is directly from decay; skip if created by pair production
#==============================================================================
    
    if part_type == 1:
            
        photon_steps = 3*10**4
        photon_dt = 10**-12
        
        particle_count = 0
    
        # Momentum array for a single muon away from the magic momentum
        m_p = cf.getMuonMomentumAtDecay()       # (eV/c)
        
        '''Use local muon variables to get inital global particle variables'''
        
        # Radial distance from center of ring
        r = R + m_x[0]                          # (m)
        
        # Set x,y,z positions
        
        x = np.zeros((1,3))                   # Initialize position array
        x[0,0] = r*np.cos(m_theta)            # (m) x-position
        x[0,1] = r*np.sin(m_theta)            # (m) y-position
        x[0,2] = m_x[1]                       # (m) z-position
        
        # Particle velocity at decay
        v = np.zeros((steps,3))                 # (m/s) Initialize velocity
        
        # (eV/c) Momentum
        p = p_initial = \
            cf.getParticleMomentumAtDecay(x,m_theta,m_p,m,p_range,p_n)
        
        beta = cf.momentum2Beta(p,m)            # () Relativistic beta
        v[0] = beta*c                           # (m/s) Initial velocity
    
        # Convert dx/ds to dx/dt
        m_xprime[0] = m_xprime[0] * cf.mag(v[0])/1000
        
        xprime = np.zeros((3))
        xprime[0] = m_xprime[0]*np.cos(m_theta)
        xprime[1] = m_xprime[0]*np.sin(m_theta)
        
        v[0] = cf.addVelocities(v[0],xprime)
        
        particle_pos = np.zeros((500,steps,3))
        photon_pos = np.zeros((photon_steps,3))
        particle_proc = np.zeros((500,9))
        particle_proc_old = np.zeros((500,9))
        particle_proc[0] = np.array([x[0,0],x[0,1],x[0,2],
                                    v[0,0],v[0,1],v[0,2],1,0,0])
        photon_proc = np.zeros((500,11))
        photon_proc_old = np.zeros((500,11))
    
#==============================================================================
#   Initialize misc. variables
#==============================================================================
    
    # [0] Counter for photons released by Bremsstrahlung
    # [1] Counter for detectable photons released by Bremsstrahlung 
    photon_count = np.zeros([2],dtype='int')
    
    # Possible energies of Bremsstrahlung photons
    k_min = 0.04*10**9                          # (eV)
    k_max = cf.momentum2Energy(p_initial,m)     # (eV)
    
    # Photon position array
    
    particle_row_index = 0
    photon_row_index = 0
    kk = 0
#    while particle_proc_change == 1 or photon_proc_change == 1:
    while 1 < 2:
        
        particle_proc_old = np.copy(particle_proc)
            
        for row in particle_proc:
        
            if (np.any(row != 0)) and row[7] == 0:
                print('particle #: %d'%particle_row_index)
                particle_pos,particle_matrix,particle_proc,photon_count,\
                    photon_proc = \
                    pt.track(particle_pos,particle_matrix,particle_proc,
                             photon_pos,photon_proc,dt,steps,m,B,k_min,
                             k_max,geo_pack,particle_count,photon_count,
                             particle_row_index)
                             
            particle_row_index = particle_row_index + 1
        
#        if np.array_equal(particle_proc,particle_proc_old):
#            particle_proc_change = 0
            
        photon_proc_old = np.copy(photon_proc)
            
        for row in photon_proc:
            
            if (np.any(row != 0)) and row[9] == 0:
                print('photon #: %d'%photon_row_index)
                particle_pos,particle_matrix,particle_proc,photon_proc, \
                    particle_count,photon_matrix = \
                    pht.track(particle_pos,particle_matrix,particle_proc,
                              photon_pos,photon_proc,photon_matrix,dt,steps,
                              m,B,k_min,k_max,geo_pack,particle_count,
                              photon_steps,photon_dt,photon_row_index)
                             
            photon_row_index = photon_row_index + 1
            
        if np.array_equal(photon_proc,photon_proc_old) and \
            np.array_equal(particle_proc,particle_proc_old):
            break
        
        
        kk = kk + 1
        particle_row_index = 0
        photon_row_index = 0
                
        
#==============================================================================
#   Cleanup and some useful variable creation
#==============================================================================
        
    # Removes unused parts of the matricies (the unchanged zeros)
        
    particle_proc = particle_proc[~np.all(particle_proc == 0,axis=1)]
    particle_pos = particle_pos[0:len(particle_proc):1]
    photon_proc = photon_proc[~np.all(photon_proc == 0, axis=1)]
            
#==============================================================================
#   Plotting and other output
#==============================================================================
    
    if make_plots == 1:
        
        print('Plotting...')
            
        # Counter used for setting multiple figures
        n = 0
            
        # Figure
        plt.figure(n)
        n = n + 1
        ax = plt.subplot(1,1,1)
        
        pg.plot(geo_pack,steps,ax) # Adds the geometries to the figure
        
        part_index = 1 # As first row in particle_matrix is headers
    
        for particle in particle_pos:
            
            final_index = int(particle_matrix[part_index,1])
            
            x = particle[0:final_index:1]
            
#            if print_text == 1:
#                
#                # Print what happened to the particle
#                print('Particle')
##                np.set_printoptions(formatter={'float': lambda x: format(x, '6.3E')})
#                print('Starting Momentum: %0.6f GeV/c'%float(particle_matrix[part_index,7]))
#                print('Kill Event: %s'%particle_matrix[part_index,2])
#                print('Final Momentum: %0.6f GeV/c'%float(particle_matrix[part_index,8]))
#                print('Photons Released: %d'%int(particle_matrix[part_index,12]))
#    #            print('Photon Energies:')
#    #            print(photon_energies[:,0])
                
            part_index = part_index + 1
                
            ax.plot(x[:,0],x[:,1], label='Particle Track',lw = 0.5)
            
        # Add Bremsstrahlung photons if they were created
            
        for photon in photon_proc:
            x_end = photon[0]+photon[3]*c*photon[8]*photon_dt
            y_end = photon[1]+photon[4]*c*photon[8]*photon_dt
            ax.plot([photon[0],x_end],[photon[1],y_end],'r-')
            
#        if photon_count_min > 0:
#            
#            k = 0
#            
#            while k < photon_count_min:
#                
#                temp_x = photon_x[k]
#                temp_x = temp_x[~np.all(temp_x == 0, axis=1)]
#                ax.plot(temp_x[:,0],temp_x[:,1],'r')
#                
#                k = k + 1
        
        # Used in adding a text box to the plots
        props = dict(boxstyle='square', facecolor='wheat', alpha=0.8)
        plt.xlabel('x-position (m)')
        plt.ylabel('y-position (m)')  
        plt.title('Particle Position')
            
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
        plt.xlim(min(x[:,0]-0.1),max(x[:,0]+0.2))
        plt.ylim(min(x[:,1]-1),max(x[:,1]+1))
#        plt.xlim(5,5.1)
#        plt.ylim(5,5.12)
        plt.grid()
    
        # Show plot
        plt.show()
        
        if save_output == 1:
    
            output_dir = "../Output/"
            path = output_dir + "particle_matrix.csv"
            with open(path, "w", newline='') as csv_file:
                writer = csv.writer(csv_file, delimiter=',')
                for row in particle_matrix:
                    writer.writerow(row)
    
            output_dir = "../Output/"
            path = output_dir + "photon_matrix.csv"
            with open(path, "w", newline='') as csv_file:
                writer = csv.writer(csv_file, delimiter=',')
                for row in photon_matrix:
                    writer.writerow(row)
    
    return particle_matrix
        
#==============================================================================
#   End 'run' function
#==============================================================================
    
if __name__ == '__main__':
    
    main()