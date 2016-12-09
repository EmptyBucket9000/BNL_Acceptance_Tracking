# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 07:12:50 2016

@author: Eric Schmidt
"""

"""
This script is used if the same muon data set is used for multiple runs in
order to reduce statistical uncertainties.

Combines all the files created by the process_single_files.py script to be
read as usual.

**All files to be combined must be copied to a single location, by default,
it's: ' ../Output/Combined/ts/' where 'ts' is set below.**

"""

import numpy as np
import csv
import os
import glob

ts = 13
out_extra = ""          # E.g. "_group_2", be sure to begin with "_" (for output)
output_dir = "../Output/"
N_max = 100000      # Maximum # of tracked particles for pre-allocation

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
                                 "Muon #",
                                 "Muon Set #"])
                             
                             
      
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
                               "Muon #",
                               "Muon Set #"])
                           
N_part_mat = len(particle_matrix_header)                               
N_phot_mat = len(photon_matrix_header)

#==============================================================================
# Particle Files
#==============================================================================

# Get list of files to be combined                                                               
particle_files = sorted(glob.glob("%s/../Output/Combined/%d/particle_*.csv"%(
                            os.getcwd(),ts)), key=os.path.getmtime)

# Set output file name and location
path = output_dir + "combined_particle_matrix%s_%d.csv"%(out_extra,ts)

# Row counter, starts at 1 as the output file has a header
i = 1
setNum = 0

# Array to be written to file

particle_matrix_full = np.zeros((N_max,N_part_mat),dtype=object)
particle_matrix_full[0] = particle_matrix_header

# Loop through each file
for file in particle_files:

    with open(file, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
        stuff = list(reader)
        N = len(stuff)
        for row in stuff:
            row = np.append(row,setNum)
            # Write each row in file to the new output file
            particle_matrix_full[i] = row
            i = i + 1

    setNum = setNum + 1
    
# Remove unused rows
particle_matrix_full = particle_matrix_full[0:i:1]

# Write the output array to file
with open(path, "w", newline='') as csv_file:
    writer = csv.writer(csv_file, delimiter=',')
    for row in particle_matrix_full:
        writer.writerow(row)

#==============================================================================
# Photon Files
#==============================================================================

# See the above particle section for comments

photon_files = sorted(glob.glob("%s/../Output/Combined/%d/photon_*.csv"%(
                        os.getcwd(),ts)), key=os.path.getmtime)
path = output_dir + "combined_photon_matrix%s_%d.csv"%(out_extra,ts)
i = 1
setNum = 0
photon_matrix_full = np.zeros((N_max,N_phot_mat),dtype=object)
photon_matrix_full[0] = photon_matrix_header

for file in photon_files:
    with open(file, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
        stuff = list(reader)
        for row in stuff:
            row = np.append(row,setNum)
            photon_matrix_full[i] = row
            i = i + 1

    setNum = setNum + 1
            
photon_matrix_full = photon_matrix_full[0:i:1]

with open(path, "w", newline='') as csv_file:
    writer = csv.writer(csv_file, delimiter=',')
    for row in photon_matrix_full:
        writer.writerow(row)