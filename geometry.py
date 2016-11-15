# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

'''
A file 'geo_pack' is created here, to unpack, use:

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
    cal_width = geo_pack[21]
    cal_height = geo_pack[22]
    cal_theta_glob = geo_pack[23]
    
See README.md for information.

'''

import numpy as np

def geo():

    R = 7.112                   # (m) Radius of the ring
    R_i = 6.805
    
    ''' Calorimeters '''
    
    cal_width = 0.225           # (m) Calorimeter length (width)
    cal_depth = 0.4572          # (m) Depth of calorimeter
    cal_height = 0.14           # (m) Height of the calorimeter
    cal_theta_glob = 0.070      # (rad) Angle of cal w.r.t. radial
    
    # (rad) Location as a function of theta. 0.001 subtracted from
    # 'cal_theta_end' to prevent particles from being counted as contact with
    # the calorimeter by contacting the 'top' of the calorimeter.
    
    cal_theta_start = np.linspace(0,2*np.pi-np.pi/12,24) + np.pi/12
    cal_box_theta_end = cal_theta_start - cal_depth/(R_i+cal_width)
    cal_theta_end = cal_theta_start - .001
    cal_theta = np.column_stack((cal_theta_start,cal_theta_end))
    cal_box_theta = np.column_stack((cal_theta_start -0.001,cal_box_theta_end))
    
    # Radial position [r_min,r_max]
    cal_rad = np.array([R_i,R_i+cal_width]) # (m)
    
    ''' HV Standoff '''
    
    rot = 9.23 + 0.48
    so_z_max = 3.8*10**-3       # (m) Top of standoff
    so_length = 7.62*10**-3     # (m) Width of standoff
    so_depth = 27.5*10**-3      # (m) Standoff depth
    so_rad_start = 50.5*10**-3  # (m) Starting distance in from R
    so_rad = np.array([R - (so_rad_start + so_depth), R - so_rad_start])
    so_theta_start_base = np.array([2,8,13,15,21,26.75,32,38,43])+rot
    so_theta_start = np.array([2,8,13,15,21,26.75,32,38,43])+rot
        
    i = 1
    while i < 5:
        so_theta_start = np.concatenate((
            so_theta_start, so_theta_start_base + i*90
        ))
        i = i + 1
    so_theta_start = so_theta_start*np.pi/180 + so_length/(2*(R-so_rad_start))
    so_theta_end = so_theta_start - so_length/(R-so_rad_start)
    
    # Theta position [theta_min, theta_max]
    so_theta = np.column_stack((so_theta_start,so_theta_end))
    
    ''' Support plates '''
    
    rot = 9.23
    sp_length = 33*10**-3        # (m) Support plate width
    sp_depth = 1.7*10**-3        # (m) Support plate thickness
    sp_rad_start = 78*10**-3     # (m) Distance in from R
    sp_rad = np.array([R - (sp_rad_start + sp_depth), R - sp_rad_start]) # (m)
    sp_theta_start_base = np.array([2.48,8.48,13.48,15.48,21.48,27.23])+rot
    sp_theta_start = np.array([2.48,8.48,13.48,15.48,21.48,27.23])+rot
    
    i = 1
    while i < 12:
        sp_theta_start = np.concatenate((
            sp_theta_start, sp_theta_start_base + i*30
        ))
        i = i + 1
    sp_theta_start = sp_theta_start*np.pi/180 + sp_length/(2*(R-sp_rad_start))
    sp_theta_end = sp_theta_start - sp_length/(R-sp_rad_start)
    
    # Theta position [theta_min, theta_max]
    sp_theta = np.column_stack((sp_theta_start,sp_theta_end))
    
    ''' Variables for both quads '''
    
    qel_z_max = 23.5*10**-3     # (m) Top of electrodes
    qel_depth = 0.5*10**-3      # (m) Electrode thickness
    qel_rad_start = 50*10**-3   # (m) Starting distance in from R
    
    ''' Single-quad electrodes (without edge curls) '''
    rot = 3
    sqel_rad = np.array([R - (qel_rad_start + qel_depth),R - qel_rad_start])
    sqel_theta_base = 90 - np.array([33.39,46.39])-rot
    sqel_theta = 90 - np.array([33.39,46.39])-rot
    
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
    
    ''' Trolly Rail '''
    
    rail_height = 40*10**-3                  # (m) Where the rail starts in y
    rail_rad = R - np.array([56,46])*10**-3  # (m) Distance from R 
    
    
#     Pack up the geometry variables for easier transfer
    geo_pack = np.array([cal_theta,cal_theta_start,cal_rad,
                         cal_box_theta,cal_box_theta_end,
                         so_z_max,so_rad,so_theta,so_theta_start,so_theta_end,
                         sp_rad,sp_theta,sp_theta_start,sp_theta_end,
                         qel_z_max,
                         sqel_rad,sqel_theta,
                         dqel_rad,dqel_theta,
                         R,R_i,cal_width,cal_height,cal_theta_glob,
                         rail_height,rail_rad])
                         
    return geo_pack