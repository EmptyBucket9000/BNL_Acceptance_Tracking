# -*- coding: utf-8 -*-
"""
Created on Tue Oct 11 14:42:14 2016

@author: Eric Schmidt
"""

import numpy as np
import matplotlib.pyplot as plt
import csv
import os
import glob
    
def main():
    
    save_plots = 0                      # Set to 1 to save plots, 0 otherwise
    save_dir = "../Output/Images"       # Set save directory
    image_dpi = 500                     # Set saved image dpi
    
    ts = 12
#    extra = "_angle" # Note the underscore that should be added
    extra = "_group_2"
    
    R = 7.112*100
    
#==============================================================================
# Particles
#==============================================================================
    
    particle_file = glob.glob("%s/../Output/particle_matrix%s_%d.csv"%(
                                os.getcwd(),extra,ts))
#    particle_file = glob.glob("%s/particle_matrix%s_%d.csv"%(
#                                os.getcwd(),extra,ts))
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
        '''
        
        N_particles = len(stuff)
        i = 0
        
        # [Kill event,dt,charge]
        particle = np.zeros((N_particles,3),dtype=object)
        
        # [x,y,z] Global muon position at decay
        x = np.zeros((N_particles,3))       # (mm) 
        
        # [Starting, Ending, Difference]
        p = np.zeros((N_particles,3))       # (GeV/c)
        
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
        
        through_quad = np.zeros((1,8),dtype=int)
        data_qel_contact = np.zeros((8))
        data_no_qel_contact = np.zeros((8))
        yerr_qel_contact = np.zeros((8))
        yerr_no_qel_contact = np.zeros((8))
        
        # sqel or dqel contact then cal contact
        through_quad_contact = np.zeros((1,8),dtype=int)
        
        # No sqel or dqel contact then cal contact
        no_through_quad_contact = np.zeros((1,8),dtype=int)
        
        for row in stuff:
            
            particle[i,0] = row[2]
            particle[i,1] = row[32]
            particle[i,2] = row[3]
            
            if row[2] == "Calorimeter Contact":
                cal_con = cal_con + 1
            
            if row[2] == "Trolly Rail Contact":
                rail_contact = rail_contact + 1
            
            x[i,0] = float(row[4])
            x[i,1] = float(row[5])
            x[i,2] = float(row[6])
            
            r = np.sqrt(x[i,0]**2 + x[i,1]**2)/10   # (cm) Starting radius
            d = R - r                               # (cm) Starting local x
            
            x_cal[i,0] = row[7]
            x_cal[i,1] = row[8]
            
            p[i,0] = row[9]
            p[i,1] = row[10]
            p[i,2] = row[11]
            
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
            
            # If starting x-position (cm) is greater 3
            if d >= 3:
                
                # Add 1 to the counter for the total # of particles
                through_quad[0,7] = through_quad[0,7] + 1
                
                # Check if the particle hit the calorimeter and a quad
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    # Count the particle
                    through_quad_contact[0,7] = through_quad_contact[0,7] + 1
                
                # Check if the particle hit the calorimter and missed all quads
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    # Count the particle
                    no_through_quad_contact[0,7] = no_through_quad_contact[0,7] + 1
            
            if d < 3 and d >= 2:
                
                through_quad[0,6] = through_quad[0,6] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    through_quad_contact[0,6] = through_quad_contact[0,6] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    no_through_quad_contact[0,6] = no_through_quad_contact[0,6] + 1
            
            if d < 2 and d >= 1:
                
                through_quad[0,5] = through_quad[0,5] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    through_quad_contact[0,5] = through_quad_contact[0,5] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    no_through_quad_contact[0,5] = no_through_quad_contact[0,5] + 1
                    
            if d < 1 and d >= 0:
                
                through_quad[0,4] = through_quad[0,4] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    through_quad_contact[0,4] = through_quad_contact[0,4] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    no_through_quad_contact[0,4] = no_through_quad_contact[0,4] + 1
            
            if d < 0 and d >= -1:
                
                through_quad[0,3] = through_quad[0,3] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    through_quad_contact[0,3] = through_quad_contact[0,3] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    no_through_quad_contact[0,3] = no_through_quad_contact[0,3] + 1
            
            if d < -1 and d >= -2:
                
                through_quad[0,2] = through_quad[0,2] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    through_quad_contact[0,2] = through_quad_contact[0,2] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    no_through_quad_contact[0,2] = no_through_quad_contact[0,2] + 1
            
            if d < -2 and d >= -3:
                
                through_quad[0,1] = through_quad[0,1] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    through_quad_contact[0,1] = through_quad_contact[0,1] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    no_through_quad_contact[0,1] = no_through_quad_contact[0,1] + 1
            
            if d < -3:
                
                through_quad[0,0] = through_quad[0,0] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) > 0 or float(in_dqel[i,0]) > 0):
                        
                    through_quad_contact[0,0] = through_quad_contact[0,0] + 1
                
                if particle[i,0] == "Calorimeter Contact" and \
                    (float(in_sqel[i,0]) == 0 and float(in_dqel[i,0]) == 0):
                        
                    no_through_quad_contact[0,0] = no_through_quad_contact[0,0] + 1   
                
            i = i + 1
            
#==============================================================================
# Data Processing
#==============================================================================
    
    # Remove rows of all zeros
    angles = angles[np.any(angles != 0, axis = 1)]
    angles_mean = np.mean(angles[:,2])
    
    total_particles = len(x)
#    total_photons = sum(in_sqel[:,2]) + sum(in_dqel[:,2]) + \
#                    sum(in_sp[:,2]) + sum(in_so[:,2]) + sum(in_sos[:,2])
#    total_in_matter = sum(in_sqel[:,1]) + sum(in_dqel[:,1]) + \
#                      sum(in_sp[:,1]) + sum(in_so[:,1] + sum(in_sos[:,1]))
                      
    print('Total particles: %d'%total_particles)
#    print('Total distance in matter: %0.3f cm'%(total_in_matter))
#    print('Total photons: %d'%total_photons)
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
    print('Average calorimeter contact angle: %0.3f'%angles_mean)
    print('Total particle calorimeter contacts: %d'%cal_con)
    print('# of calorimeter contacts after qel contact: %d'%np.sum(
            through_quad_contact))
    print('# of calorimeter contacts without qel contact: %d'%np.sum(
            no_through_quad_contact))
    print('Total SO/SO screw contacts that hit the calorimeter: %d'\
            %cal_con_so)
#==============================================================================
#     Plotting
#==============================================================================
            
    i = 0
    
    while i < np.shape(through_quad)[1]:
        
        # Get the fraction that passed through a quad and hit a calorimeter to
        # the total # in that starting x-position range
        data_qel_contact[i] = through_quad_contact[0,i] / through_quad[0,i]
        
        # Get Poisson uncertainties
        yerr_qel_contact[i] = np.sqrt(through_quad_contact[0,i]) / \
                                through_quad[0,i]
        
        # Get the fraction that did not pass through a quad and hit a
        # calorimeter to the total # in that starting x-position range
        data_no_qel_contact[i] = no_through_quad_contact[0,i] / through_quad[0,i]
        
        # Poisson uncertainties
        yerr_no_qel_contact[i] = np.sqrt(no_through_quad_contact[0,i]) / \
                                    through_quad[0,i]
        print('%d: %d'%(i,through_quad_contact[0,i]))
        print('%d: %d'%(i,no_through_quad_contact[0,i]))
        print('%d: %d'%(i,through_quad[0,i]))
        i = i + 1
            
    # Number of bins in the calorimeter contact angle histogram
    calorimeter_angle_hist = 40
    
    # Convert string to float
    x_cal = np.array(x_cal, dtype = float)
    
    # Remove 'zero' rows
    x_cal = x_cal[np.any(x_cal != 0, axis = 1)]
        
    n = 0
    
    # Create a scatter plot of the fraction of particles that passed through
    # a quad and hit a calorimeter to the total # in that starting x-position
    # range, and create a scatter plot of the fraction that did not pass
    # through a quad to the total # in that starting x-position range
    
    plt.figure(n)
    n = n + 1
    
    xx = np.arange(-3.5,4)
    ax = plt.subplot(1,1,1)
    ax.errorbar(xx,data_no_qel_contact,yerr=yerr_no_qel_contact,fmt='.',
                color='g',
                label='No through-quad')
    ax.errorbar(xx,data_qel_contact,yerr=yerr_qel_contact,fmt='.',
                color='b',
                label='Through-quad')
    ax.legend(bbox_to_anchor=(0.48,0.26))
    ax.set_title("Admittance as a function of muon position")
    ax.set_xlabel("x-Position (cm)")
    ax.set_ylabel("Particles Detected / Total")
    
    if save_plots == 1:
            plt.savefig('%s/admittance.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
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
    plt.xlim(-24,24)
    plt.ylim(-15.5,15.5)
    ax.grid(True)
    ax.legend(bbox_to_anchor=(1.33,1.11))
    ax.set_title("Calorimeter Contact Position (Particles)")
    ax.set_xlabel('x-Position (cm)')
    ax.set_ylabel('y-position (cm)')
    plt.axis('equal')
    
    if save_plots == 1:
            plt.savefig('%s/particle_calorimeter_contact.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
                        
    # Plot a histogram of calorimeter contact angles where the angle is from
    # the positive x-axis
    
    plt.figure(n)
    n = n + 1
    
    ax = plt.subplot(1,1,1)
    ax.hist(angles[:,0],calorimeter_angle_hist)
    ax.set_title("x Calorimter Contact Angles (Particles)")
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Count')
    
    if save_plots == 1:
            plt.savefig('%s/particle_calorimeter_contact_x_angle.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)

    # Plot a histogram of calorimeter contact angles where the angle is from
    # the positive y-axis
    
    plt.figure(n)
    n = n + 1
    
    ax = plt.subplot(1,1,1)
    ax.hist(angles[:,1],calorimeter_angle_hist)
    ax.set_title("y Calorimter Contact Angles (Particles)")
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Count')
    
    if save_plots == 1:
            plt.savefig('%s/particle_calorimeter_contact_y_angle.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)

    # Plot a histogram of calorimeter contact angles where the angle is from
    # the projection of the final momentum vector onto the calorimeter plane to
    # the the final momentum vector (in local calorimeter coordinates)
    
    plt.figure(n)
    n = n + 1
    
    ax = plt.subplot(1,1,1)
    ax.hist(angles[:,2],calorimeter_angle_hist)
    ax.set_title("Total Calorimter Contact Angles (Particles)")
    ax.set_xlabel('Angle (deg)')
    ax.set_ylabel('Count')
    
    if save_plots == 1:
            plt.savefig('%s/particle_calorimeter_contact_angle.png'%save_dir,
                        bbox_inches='tight',dpi=image_dpi)
    
    
if __name__ == '__main__':

    main()