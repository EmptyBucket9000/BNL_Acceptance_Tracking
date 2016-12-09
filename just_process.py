# -*- coding: utf-8 -*-
"""
Created on Mon Nov  7 06:21:29 2016

@author: Eric Schmidt
"""

'''
'just_process.py' is used only when the 'main.py' code is still running (so you
can see the current progress) or the file processing failed for some reason.

To use, simply set 'ts' and 'extra' then run.

*** Important ***
If you updated the header files in 'main.py', you must also update them here.

See README.md for more information.
'''

import numpy as np
import process_single_files as psf

ts = 13
#extra = "angle/" # Note the forward slash that must be added
extra = "group_8/"
                         
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
                             "Kill Timestamp",
                             "x Calorimeter Angle",
                             "y Calorimeter Angle",
                             "Total Calorimeter Angle",
                             "Starting Local x (m)",
                             "Starting Local y (m)",
                             "Starting Local x-prime (rad)",
                             "Starting Local y-prime (rad)",
                             "Muon #"])
      
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
                           "Kill Timestamp",
                           "x Calorimeter Angle",
                           "y Calorimeter Angle",
                           "Total Calorimeter Angle",
                           "Muon #"])
                           
N_part_mat = len(particle_matrix_header)                               
N_phot_mat = len(photon_matrix_header)

# Calls the function that would have been called normally after all muon decays
# had been tracked.
psf.process(particle_matrix_header,photon_matrix_header,N_part_mat,N_phot_mat,
            ts,extra)