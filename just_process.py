# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 06:21:29 2016

@author: schmidte
"""

import numpy as np
from process_single_files import process

ts = 13
#extra = "angle/"
extra = ""
                         
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

process(particle_matrix_header,photon_matrix_header,N_part_mat,N_phot_mat,ts,extra)