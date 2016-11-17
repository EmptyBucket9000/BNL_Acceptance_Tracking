# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 15:26:56 2016

@author: Eric Schmidt
"""

"""
See README.md for information.

"""

import numpy as np
import csv
import os
import glob


def process(particle_matrix_header,photon_matrix_header,N_part_mat,N_phot_mat,
            ts,extra):
    output_dir = "../Output/"
                
    if extra == "":
        extra_out = ""
    else:
        extra_out = "_" + extra[:-1]
    
#==============================================================================
# Particle Files
#==============================================================================
        
    particle_files = \
        glob.glob("%s/../Output/Single_Files/%s%d/particle_*.csv"%(
                os.getcwd(),extra,ts))
    path = output_dir + "particle_matrix%s_%d.csv"%(extra_out,ts)
                                                                   
    total_muons = len(particle_files)
    N_particles = total_muons*5
    i = 1
        
    particle_matrix_full = np.zeros((N_particles,N_part_mat),dtype=object)
                         
    # Output for each particle
    particle_matrix_full[0] = particle_matrix_header
    
    for file in particle_files:
    
        with open(file, "rt") as inf:
            reader = csv.reader(inf, delimiter=',')
            next(reader, None)  # skip the headers
            stuff = list(reader)
            for row in stuff:
                particle_matrix_full[i] = row
                i = i + 1
                
    particle_matrix_full = particle_matrix_full[0:i:1]

    with open(path, "w", newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for row in particle_matrix_full:
            writer.writerow(row)
    
#==============================================================================
# Photon Files
#==============================================================================
    
    photon_files = \
        glob.glob("%s/../Output/Single_Files/%s%d/photon_*.csv"%(
                    os.getcwd(),extra,ts))
    path = output_dir + "photon_matrix%s_%d.csv"%(extra_out,ts)
            
    total_files = len(photon_files)
    N_photons = total_files*10
    i = 1
    
    photon_matrix_full = np.zeros((N_photons,N_phot_mat),dtype=object)
              
    # Output for each photon
    photon_matrix_full[0] = photon_matrix_header
    
    for file in photon_files:
    
        with open(file, "rt") as inf:
            reader = csv.reader(inf, delimiter=',')
            next(reader, None)  # skip the headers
            stuff = list(reader)
            for row in stuff:
                photon_matrix_full[i] = row
                i = i + 1
                
    photon_matrix_full = photon_matrix_full[0:i:1]

    with open(path, "w", newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for row in photon_matrix_full:
            writer.writerow(row)