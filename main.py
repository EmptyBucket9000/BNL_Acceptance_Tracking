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
import muon_data as md
import plot_geometries as pg
import os
import glob
import csv
from process_single_files import process

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
    make_plots = 0              # Set to 1 to display plots
    save_plots = 0              # Set to 1 to save plots as images
    save_output = 1             # Set to 1 to save data output to csv

    # Name of csv containing muon data    
    file_name = "EndOfTracking_phase_space.csv"
        
    m_theta_set = 0                     # 1, use m_theta below, 0 random
    
    N = 5171                              # Number of muons in beam
#    N = 1
    steps = 5*10**6                     # Nnumber of steps for integration
#    steps = 3
    ts = 13
    dt = 10**-ts                        # Timestep for integration

    # Used in naming spcecial output files.
    # Folder within 'Output' must exist with this name and a subfolder with the
    # name of the value of 'ts' must also exist.
    extra = ""
    
    p_magic = 3.09435*10**9             # (eV/c) Muon magic momentum
    
#    if m_theta_set == 1:
#        m_theta = 3*np.pi/8   # (rad) Muon azimuth. position in global coords
#    else:
#        m_theta = -99
        
    m_theta_array = np.array([0,2*np.pi])
    
    ''' Permanent constants '''
    
    q = 1                                   # (e) Particle charge
    c = 2.99792458*10**8                    # (m/s) Speed of light
    m = 0.510999*10**6                      # (eV/c**2) Particle mass
    B = np.array([0,0,1.4513])              # (T) Magnetic field
    
    geo_pack = geometry.geo()               # All the permanent geometries
    
    ''' Other variables '''
                         
    # Output for each particle
    particle_matrix_header = np.array(["Particle #","Steps","Kill Event",
                                 "Charge",
                                 "Starting Global x-Position (mm)",
                                 "Starting Global y-Position (mm)",
                                 "Starting Global z-Position (mm)",
                                 "Ending Calorimeter x-Position (mm)",
                                 "Ending Calorimeter y-Position (mm)",
                                 "Starting Momentum (GeV/c)",
                                 "Ending Momentum (GeV/c)",
                                 "Delta Momentum (GeV/c)",
                                 "Steps Inside Short Quad",
                                 "Distance Inside Short Quad (cm)",
                                 "Total # of Photons Released",
                                 "# of Detectable Photons Released",
                                 "Steps Inside Long Quad",
                                 "Distance Inside Long Quad (cm)",
                                 "Total # of Photons Released",
                                 "# of Detectable Photons Released",
                                 "Steps Inside Standoff Plate",
                                 "Distance Inside Standoff Plate (cm)",
                                 "Total # of Photons Released",
                                 "# of Detectable Photons Released",
                                 "Steps Inside HV Standoff",
                                 "Distance Inside HV Standoff (cm)",
                                 "Total # of Photons Released",
                                 "# of Detectable Photons Released",
                                 "Steps Inside HV Standoff Screws",
                                 "Distance Inside HV Standoff Screws (cm)",
                                 "Total # of Photons Released",
                                 "# of Detectable Photons Released",
                                 "dt","Pair Produced",
                                 "Kill Timestamp"])
          
    # Output for each photon
    photon_matrix_header = np.array(["Photon #","Steps","Kill Event",
                               "Starting Global x-Position",
                               "Starting Global y-Position",
                               "Starting Global z-Position",
                               "Ending Calorimeter x-Position (mm)",
                               "Ending Calorimeter y-Position (mm)",
                               "Energy (GeV)","Steps Inside Matter",
                               "Distance Inside Matter (cm)",
                               "dt",
                               "Kill Timestamp"])
                               
    N_part_mat = len(particle_matrix_header)                               
    N_phot_mat = len(photon_matrix_header)
    
    muon_number = 0         # Used as a counter for multiple particles
    
    # Estimated maximum # of particles and photons that could be created from
    # each muon decay (smaller -> less memory usage but too small and an
    # error could be thrown)
    N_particles = 50
    N_photons = 50
    
    ## Delete all the current single files
    
    single_files = glob.glob("%s/../Output/Single_Files/*.csv"%(os.getcwd()))
    if len(single_files) > 0:
        for f in single_files:
            os.remove(f)
#==============================================================================
#   Run the code!!
#==============================================================================
                             
    while muon_number < N:
        
        m_x_list,m_p_list,m_theta = \
            md.muon(N,file_name,p_magic,m_theta_array,m_theta_set)
        
        m_x = m_x_list[muon_number]
        m_p = m_p_list[muon_number]
        
        particle_matrix = np.zeros((N_particles,N_part_mat),dtype=object)
        photon_matrix = np.zeros((N_photons,N_phot_mat),dtype=object)
    
        particle_matrix[0] = particle_matrix_header
        photon_matrix[0] = photon_matrix_header
        
        # Set all initial local muon variables
                
        particle_matrix,photon_matrix = \
            run(geo_pack,m_x,m_p,m_theta,m,c,photon_matrix,
                particle_matrix,make_plots,save_plots,save_output,
                N,steps,dt,q,B,muon_number,N_particles,N_photons,ts)
                
        muon_number = muon_number + 1
                
        print("Percent complete: %0.3f%%"%(muon_number*(100/N)),end="\r")
              
    # Finally, convert all the single particle_matrix and photon_matrix files
    # into single files.
              
    if save_output == 1:
        process(particle_matrix_header,photon_matrix_header,N_part_mat,
                N_phot_mat,ts,extra)

#==============================================================================
#==============================================================================
# #   Begin 'run' function, where everything happens
#==============================================================================
#==============================================================================

def run(geo_pack,m_x,m_p,m_theta,m,c,photon_matrix,
        particle_matrix,make_plots,save_plots,save_output,
        N,steps,dt,q,B,muon_number,N_particles,N_photons,ts):

#==============================================================================
#   Initialization and setting variables
#==============================================================================

    R = geo_pack[19]
        
    photon_steps = 3*10**5
    photon_dt = 10**-11
    
    particle_count = 0
    
    '''Use local muon variables to get inital global particle variables'''
    
    m_m = 105.65837*10**6                  # (eV/c**2) Muon mass
    
    # Radial distance from center of ring
    r = R + m_x[0]                          # (m)
    
    # Set x,y,z positions
    
    x = np.zeros((1,3))                   # Initialize position array
    x[0,0] = r*np.cos(m_theta)            # (m) x-position
    x[0,1] = r*np.sin(m_theta)            # (m) y-position
    x[0,2] = m_x[1]                       # (m) z-position
    p_range = np.array([1.8*10**9,cf.mag(m_p)])
#        p_range = np.array([1.2,1.2])*10**9
    
    p_s = cf.getParticleMomentumAtDecay(p_range,m_p,m_theta,m_m)
    
    p = np.zeros((3))
    p[0] = p_s[0] + m_p[0]*np.cos(m_theta)
    p[1] = p_s[1] + m_p[0]*np.sin(m_theta)
    p[2] = np.copy(m_p[1])
    
    particle_pos = np.zeros((N_particles,steps,3))
    photon_pos = np.zeros((photon_steps,3))
    particle_proc = np.zeros((N_particles,10))
    particle_proc_old = np.zeros((N_particles,10))
    particle_proc[0] = np.array([x[0,0],x[0,1],x[0,2],
                                p[0],p[1],p[2],1,0,0,0])
    photon_proc = np.zeros((N_photons,11))
    photon_proc_old = np.zeros((N_photons,11))
    
    # [0] Counter for photons released by Bremsstrahlung
    # [1] Counter for detectable photons released by Bremsstrahlung 
    photon_count = np.zeros([2],dtype='int')
    
    # Possible energies of Bremsstrahlung photons
    k_min = 0.04*10**9                          # (eV)
    k_max = cf.momentum2Energy(p,m)     # (eV)
    
    # Photon position array
    
    particle_row_index = 0
    photon_row_index = 0
    
    while True:
        
        particle_proc_old = np.copy(particle_proc)
            
        for row in particle_proc:
        
            if (np.any(row != 0)) and row[7] == 0:
                
                particle_pos,particle_matrix,particle_proc,photon_count,\
                    photon_proc = \
                    pt.track(particle_pos,particle_matrix,particle_proc,
                             photon_pos,photon_proc,dt,steps,m,B,k_min,
                             k_max,geo_pack,particle_count,photon_count,
                             particle_row_index,muon_number)
                             
            particle_row_index = particle_row_index + 1
            
        photon_proc_old = np.copy(photon_proc)
            
        for row in photon_proc:
            
            if (np.any(row != 0)) and row[9] == 0:
                
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
            
        # Counter used for setting multiple figures
        n = 0
        
        plt.figure(n)
        ax = plt.subplot(1,1,1)
        
        index = 1 # As first row in particle_matrix is headers
    
        for particle in particle_pos:
                     
            x = float(particle_matrix[index,7])
            y = float(particle_matrix[index,8])
            index = index + 1
            if x != 0.0 or y != 0.0:
                ax.scatter(x,y,c='b',label='Particle on Calorimeter')
        
        index = 1 # As first row in particle_matrix is headers
    
        for photon in photon_proc:
                       
            x = float(photon_matrix[index,6])
            y = float(photon_matrix[index,7])
            index = index + 1
            if x != 0.0 or y != 0.0:
                ax.scatter(x,y,c='r',label='Photon on Calorimeter')
            
        plt.axis('equal')
        plt.grid()
        plt.xlim(-0.2,0.2)
        plt.ylim(-0.15,0.15)
        plt.title('Particle and Photon Positions on Calorimeter')
        
        n = n + 1
            
        # Figure
        plt.figure(n)
        n = n + 1
        ax = plt.subplot(1,1,1)
        
        pg.plot(geo_pack,steps,ax) # Adds the geometries to the figure
        
        part_index = 1 # As first row in particle_matrix is headers
    
        for particle in particle_pos:
            
            final_index = int(particle_matrix[part_index,1])
            x = particle[0:final_index:1]
            part_index = part_index + 1
            ax.plot(x[:,0],x[:,1], label='Particle Track',lw = 0.5)
            
        # Add Bremsstrahlung photons if they were created
            
        for photon in photon_proc:
            x_end = photon[0]+photon[3]*c*photon[8]*photon_dt
            y_end = photon[1]+photon[4]*c*photon[8]*photon_dt
            ax.plot([photon[0],x_end],[photon[1],y_end],'r-')
        
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

        # Add legend
#        ax.legend()

        # Set axes limits based on min/max particle positions

        plt.axis('equal') # Prevents a skewed look
        plt.xlim(min(x[:,0]-2),max(x[:,0]+2))
        plt.ylim(min(x[:,1]-2),max(x[:,1]+2))
#        plt.xlim(min(x[:,0]-0.1),max(x[:,0]+0.1))
#        plt.ylim(min(x[:,1]-0.1),max(x[:,1]+0.1))
#        plt.xlim(3.5,3.7)
#        plt.ylim(-6.1,-6.0)
        plt.grid()

        if save_plots == 1:
            plt.savefig('../Output/Images/plot.png', bbox_inches='tight', dpi=500)
    
        # Show plot
        plt.show()
        
    if save_output == 1:
        
        particle_matrix_print = \
            particle_matrix[np.any(particle_matrix != 0,axis=1)]
        photon_matrix_print = \
            photon_matrix[np.any(photon_matrix != 0,axis=1)]

        output_dir = "../Output/Single_Files/%d/"%ts
        path = output_dir + "particle_matrix_%d.csv"%muon_number
        with open(path, "w", newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            for row in particle_matrix_print:
                writer.writerow(row)

        output_dir = "../Output/Single_Files/%d/"%ts
        path = output_dir + "photon_matrix_%d.csv"%muon_number
        with open(path, "w", newline='') as csv_file:
            writer = csv.writer(csv_file, delimiter=',')
            for row in photon_matrix_print:
                writer.writerow(row)
    
    return particle_matrix,photon_matrix
        
#==============================================================================
#   End 'run' function
#==============================================================================
    
if __name__ == '__main__':
    
    main()