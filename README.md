# BNL_Admittance_Tracking
## Eric Schmidt

######################################
# Initial Notes
######################################

First, I should say that if I were to do this code again, I would make it more object-oriented and less procedural. When creating this code, I had only just learned about classes and therefore did not use them. In particular, a variable called 'geo_pack' is passed around between multiple functions. This file contains all the permanent geometry information about the ring. A geometry class would have been much more user-friendly.

If any changes are made to the data stored that require the editing of either 'particle_matrix_header' or 'photon_matrix_header' variables, they must be updated is all of the following files:

- combine.py
- just_process.py
- main.py

In addition, the following files must be updated to read the new changes:

- read_particle_output.py
- read_photon_output.py

If you do make these changes, it's advised to add the new entries to the end of the variables to make updating the 'read' files much easier.

For future updates, see the 'Future Updates' section at the bottom of this readme.

######################################
# Known Bugs/Issues
######################################

### Unique Output Files (Priority: low)
Adding a clarifying string to the output files by making 'extra' non-empty allows the code to be run multiple times without overwriting previous output. Improve this process by having the required folders automatically created.

### Pair-Produced Particle Contact (Priority: medium)
Pair-produced positrons are tracked but they are not counted if they contact the calorimeter.

######################################
# Assumptions
######################################

This section covers assumptions made and whether or not I plan to remove them.

### Bremsstrahlung (Priority: low)
The x-ray momentum vector released via Bremsstrahlung was assumed to be parallel to the particle momentum vector \cite{bib:brem_angle}.

### Compton Scattering  (Priority: medium)
Due to the high-energy of the x-rays produced from Bremsstrahlung, Compton scattering was ignored \cite{bib:scattering_rates}.

### Ring Geometries (Priority: low)
Two ring structures, high-voltage standoffs and standoff plates, were geometrically approximated by allowing them to 'bend' radially around the ring (to greatly simplify their coding). However, due to their small sizes, this approximation causes deviations of only small fractions of a millimeter from their true shapes.

### Electrode Curls (Priority: low)
The edges of the electrodes curl around. This is more a complicated geometry to code plus it only has a very small effect on tracking. Therefore it is common practice in G-2 to ignore them, which I have done.

### Trolley Rails (Priority: medium)
Due to the odd shape (being very difficult to represent digitally), it was assumed that if any particle or x-ray made contact with them, the particle or x-ray would no longer be tracked, i.e. it was 'killed'. This is a valid assumption because 1) the rails were large enough to cause a very large reduction in energy of the particle passing through them (preventing them from contacting the calorimeter), and 2) only about 0.2\% of particles make contact with them.

######################################
# Running the Code
######################################

### Editable variables ###

The variables that need editing are at the top of the 'main.py' file and are commented, but they are:

'make_plots' - Set to 0 or 1. Determines if you want to display the plots from a single muon run. Plots are 1) global position of initial positron and any possible subsequest x-rays and pair-produced particles, and 2) location on calorimeter of particle(s)/x-ray(s) if they make contact.

'save_plots' - Set to 0 or 1. Determines whether or not to save the plots created, 'make_plots' must be set to 1 for this to work.

'save_output' - Set to 0 or 1. Determines if output files are created for the particles and x-rays. See below in the section 'Variables and Functions' for the output.

'file_name' - The name of the .csv file (including extension) containing the muon data (see the 'Muons' section below for details).

'N' - Number of muons you want to create, must be <= the number of muons in the muon data file.

'ts' - Timestep power. For example, to have a timestep of 10**-13, set ts = 13.

'extra' - String to add uniqueness to the names of a certain set of output files if two different runs are being made at the same time. A folder must be created with this text inside the 'Single_Files' folder and the folder of the 'ts' value must exist within it. For example, if you set extra = "example/" (note the forward slash) and ts = 13, you must have the folder "../Output/Single_Files/example/13" where the origin of this path is the location of the code. The output files will also have 'example' such as 'particle_matrix_example_13.csv' will be created. 'extra' will generally be left as extra = "".

### IMPORTANT: if you have 'extra' set to non-empty or you change 'ts', the following files must be edited:

- just_process.py
- read_particle_output.py
- read_photon_output.py
- combine.py (optional, only if you will use)

The edits are very obvious at the top of the files. This will hopefully get updated at some later time to be more automated.

### How the code works

The code takes in one .csv file of the muon data (see the section 'Muons' below) then for each muon, creates a single positron that it tracks through the ring. The positron can create photons from Bremsstrahlung and thos photons can create electron-positron pairs from pair-production. Each photon and particle is tracked until it is killed from losing too much energy so it can be ignored, by coming into contact with the calorimeter or inner limit of the ring, or by leaving the outer limit of the ring.

For each particle or x-ray created, a row in the applicable array (for photon or particle) is added with the starting position and momentum vectors plus an entry called 'processed' (initially set to 0). These rows are added at the time of creation. Since multiple particles and photons can be created from each muon decay (after Bremsstrahlung and pair-production) these array are used to keep track of the particles and photons that have or have not been tracked. Once a particle/photon has been fully tracked, it's 'processed' entry is set to 1 and is thus ignored in the future when checks of those arrays are conducted to see if untracked particles/photons exist.

The output from each particle and photon is saved to its own file in the directory '../Output/Single_Files/ts' (where ts is the timestamp power as described above). It is does this way (instead of creating one large output array then writing that to file at the end) to keep memory usage down and to prevent having to restart from the beginning if the code crashes for some reason (e.g. power outage).

Finally, after all the particles and photons are finished being tracked, all the single files are then read and two final files are created, 'particle_matrix_ts.csv' and 'photon_matrix_ts.csv'. These two files can then be read by using the 'read_particle_output.py' and 'read_photon_output.py' files (or your own).

### Combining multiple runs into one readable file

If you have multiple runs complete and you want to combine them into a single particle matrix file and a single photon matrix file, use combine.py. To use this script, a copy of all the 'particle_matrix_*.csv' and 'photon_matrix_*.csv' files that you want combined must be added to a single directory, by default this is located (relative to main.py) at '../Output/Combined/ts' where 'ts' is the timestep set and can be changed inside combine.py. You should not combine different sets of files with different timesteps. If the directory does not exist, be sure to create it. The new files will be in the Output directory and will be called 'combined_particle_matrix*.csv' and 'combined_photon_matrix*.csv' where * represents a string you can set in the combine.py file.

######################################
# Files Used in the Project
######################################

The files used are the following:

-	main.py - The primary file that is run and where the editable variables live.

-	callable_functions.py - A file of functions used by the other files.

-	combine.py - If using mulitple runs of the same set of muons, combine.py will take all the files created by process_single_files.py and combine them into one final output to be read as normal. See the script for how to prepare your files to be used.

-	find_fit_chisquare.py - Note: Not fully commented as it's still being updated/changed regularly. Finds polynomial fits and chi-square values for the data. Used by particle_tracking.py.

-	geometry.py - Builds the geometry of the ring to be used in determining if a particle/x-ray is inside it, and used for plotting.

-	muon_data.py - Returns the muon position and momentum vectors in the local coordinate system.

-	just_process.py - Used only for checking the single output files of completed muons while the code is running.

-	particle_tracking.py - Note: Not fully commented as it's still being updated/changed regularly. Tracks the movement of the positrons and electrons.

-	photon_tracking.py - Note: Not fully commented as it's still being updated/changed regularly. Tracks the movement of the x-rays.

-	plot_geometries.py - Plots the geometries if make_plots == 1.

-	process_single_files.py - Converts all the single output files (one for each muon) to a single files for particles (particle_matrix_ts.csv) and a single file for x-rays (photon_matrix_ts.csv) where 'ts' in both files is the timestamp set in the variables section. 

-	read_particle_output.py - Reads the particle file created by 'process_single_files.py.

-	read_photon_output.py - Reads the photon file created by 'process_single_files.py.

-	README.md - This file.

######################################
# Coordinate Systems
######################################

Two primary coordinate systems are used, global and local, plus a third, temporary system:
 
## Global
Global is the entire ring in 3-dimensions cartesian coordinates where the origin is the center of the ring and the z-direction indicates 'up' and 'down', from the perspective of a person standing in the ring, and is parallel to the y-direction in the local coordinate system.

## Local  
Local is a 2-dimensional x,y cartesian coordinate where the origin is the beam centroid (assume point of view of the muons), the positive x-direction points towards the center of the ring and the y-direction is parallel to the z-direction in the global system. Local is used in determining particle initial conditions, i.e. the conditions describing the muon at decay

## Other Local
In addition, a different type of local coordinate system is used when determining if a particle is inside matter. In some cases, a useful local coordinate system is used for the comparison, which is some x,y cartesian system used only for that particular matter, such as the high-voltage standoff. Another is used when determining contact information with the front face of the calorimeters.

######################################
# Muons
######################################

The muon data, i.e. momentum and position vectors (in local coordinates), come from a data file from the simulation group that contains muon information after 100 microseconds in the ring. The first 5 columns should be (including sample data from 2 muons):

|x			|xp			|y			|yp			|pz			|		
|-----------|-----------|-----------|-----------|-----------|
|0.0033485	|0.0025989	|-0.0091195	|-0.00037342|0.0018116 	|
|-0.013607	|-0.00085284|-0.011613	|-0.0014243	|0.00037121	|

where x,y are positions in the local coordinate system, xp and yp are xprime and yprime (from x,y phase-spaces), and pz is (momentum - magic momentum)/magic momentum.

######################################
# Basic Variables and Functions
######################################
 
Variable names formatted 'm_*' are for muons, all others are for particles
 
Magic momentum is set to 3.09435 GeV/c giving a magnetic field of 1.4513 T
for a ring radius of 7.112 m.

The output variable that is saved to a file contains:

## For particles
(the text is the column header and the number is the index in the array used to generate the output .csv):

|Column Header							|Index	|
|---------------------------------------|-------|
|Particle #                          	|0		|
|Steps                               	|1		|
|Kill Event                          	|2		|
|Charge                              	|3		|
|Starting Global x-Position (mm)	 	|4		|  
|Starting Global y-Position (mm)	 	|5		|  
|Starting Global z-Position (mm)     	|6		|  
|Ending Calorimeter x (mm)           	|7		|  
|Ending Calorimeter y (mm)           	|8		|  
|Starting Momentum (GeV/c)           	|9		|  
|Ending Momentum (GeV/c)             	|10		|  
|Delta Momentum (GeV/c)              	|11		|
|Steps Inside Short Quad             	|12		|  
|Distance Inside Short Quad (cm)	 	|13		|  
|Total # of Photons Released         	|14		|  
|# of Detectable Photons Released		|15		|  
|Steps Inside Long Quad              	|16		|
|Distance Inside Long Quad (cm)	 		|17		|  
|Total # of Photons Released         	|18		|  
|# of Detectable Photons Released		|19		|  
|Steps Inside Standoff Plate         	|20		|  
|Distance Inside Standoff Plate (cm) 	|21		|  
|Total # of Photons Released         	|22		|  
|# of Detectable Photons Released		|23		|  
|Steps Inside HV Standoff            	|24		|  
|Distance Inside HV Standoff (cm)		|25		|  
|Total # of Photons Released         	|26		|  
|# of Detectable Photons Released		|27		|  
|Steps Inside HV Standoff Screws     	|28		|  
|Distance Inside HV Standoff Screws		|29		|  
|Total # of Photons Released         	|30		|  
|# of Detectable Photons Released		|31		|  
|dt (particle timestep)		        	|32		|  
|Pair Produced (0 or 1)              	|33		|  
|Kill Timestamp                      	|34		|
|x Calorimeter Angle					|35		|
|y Calorimeter Angle					|36		|
|Total Calorimeter Angle				|37		|
|Starting Local x (mm)               	|38		|
|Startine Local y (mm)               	|39		|
|Starting Local x-prime (mrad)       	|40		|
|Starting Local y-prime (mrad)       	|41		|

## For photons

|Column Header							|Index	|
|---------------------------------------|-------|
|Photon #                        		|0		|
|Steps                           		|1		|
|Kill Event                      		|2		|  
|Starting Global x-Position      		|3		|  
|Starting Global y-Position      		|4		|  
|Starting Global z-Position      		|5		|  
|Ending Calorimeter x (mm)       		|6		|  
|Ending Calorimeter y (mm)       		|7		|  
|Energy (GeV)                    		|8		|  
|Steps Inside Matter             		|9		|  
|Distance Inside Matter (cm)     		|10		|  
|dt (photon timestep)            		|11		|
|Kill Timestamp                  		|12		|  
|x Calorimeter Angle					|13		|
|y Calorimeter Angle					|14		|
|Total Calorimeter Angle				|15		|

######################################
# Future Updates
######################################

These updates are not set in stone and may or may not happen depending on need and time.

### Calorimeters (Likelihood: Medium)
Differentiate between contact with different calorimeters from the same decay. It's rare that contact is made with more than one calorimeter from the same muon decay.

### x-prime (Likelihood: High)
Determine acceptance as a function of x-prime (transverse momentum / longitudinal momentum).

### Spin (S_x) (Likelihood: Medium)
Determine acceptace as a function of muon spin.

### y, y-prime, S_y (Likelihood: Medium)
Repeat the above but in the y phase-space and with the y-component of spin.
