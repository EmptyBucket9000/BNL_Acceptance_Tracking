# -*- coding: utf-8 -*-
"""
Created on Wed Nov  2 15:26:56 2016

@author: Eric Schmidt
"""

import numpy as np
import csv
import os
import glob


def process(particle_matrix_header,photon_matrix_header):
    
#==============================================================================
# Particle Files
#==============================================================================
    
    particle_files = \
        glob.glob("%s/../Output/Single_Files/particle_*.csv"%(os.getcwd()))
    total_muons = len(particle_files)
    N_particles = total_muons*5
    i = 1
        
    particle_matrix_full = np.zeros((N_particles,33),dtype=object)
                         
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

    output_dir = "../Output/"
    path = output_dir + "particle_matrix.csv"
    with open(path, "w", newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for row in particle_matrix_full:
            writer.writerow(row)
    
#==============================================================================
# Photon Files
#==============================================================================
    
    photon_files = \
        glob.glob("%s/../Output/Single_Files/photon_*.csv"%(os.getcwd()))
    total_files = len(photon_files)
    N_photons = total_files*5
    i = 1
        
    photon_matrix_full = np.zeros((N_photons,11),dtype=object)
              
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

    output_dir = "../Output/"
    path = output_dir + "photon_matrix.csv"
    with open(path, "w", newline='') as csv_file:
        writer = csv.writer(csv_file, delimiter=',')
        for row in photon_matrix_full:
            writer.writerow(row)