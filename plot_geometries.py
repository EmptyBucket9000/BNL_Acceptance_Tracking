# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

import numpy as np

def plot(geo_pack,steps,ax):
    
    lw = 0.2

    # Unpack 'geo_pack'    
    
    cal_theta_start = geo_pack[1]
    cal_rad = geo_pack[2]
    cal_box_theta_end = geo_pack[4]
    so_rad = geo_pack[6]
    so_theta_start = geo_pack[8]
    so_theta_end = geo_pack[9]
    sp_rad = geo_pack[10]
    sp_theta_start = geo_pack[12]
    sp_theta_end = geo_pack[13]
    sqel_rad = geo_pack[15]
    sqel_theta = geo_pack[16]
    dqel_rad = geo_pack[17]
    dqel_theta = geo_pack[18]
    R_i = geo_pack[20]
    
    # Add inner radius
                
    xt = np.linspace(-R_i,R_i,steps)    # (m)
    ax.plot(xt,np.sqrt((R_i)**2 - xt**2),'-.r',lw=lw)
    ax.plot(xt,-np.sqrt(R_i**2 - xt**2),'-.r',lw=lw)
    
    # Add support plates
    
    count = len(sp_theta_start)
    
    k = 0
    while k < count:
        ax.plot(
            [sp_rad[0]*np.cos(sp_theta_start[k]),
             sp_rad[0]*np.cos(sp_theta_end[k])],
            [sp_rad[0]*np.sin(sp_theta_start[k]),
             sp_rad[0]*np.sin(sp_theta_end[k])],
            'k',lw=lw
        )
        ax.plot(
            [sp_rad[1]*np.cos(sp_theta_start[k]),
             sp_rad[1]*np.cos(sp_theta_end[k])],
            [sp_rad[1]*np.sin(sp_theta_start[k]),
             sp_rad[1]*np.sin(sp_theta_end[k])],
            'k',lw=lw
        )
        
        k = k + 1
    
    # Add HV Standoff
    
    count = len(so_theta_start)
    k = 0
    
    while k < count:
        ax.plot(
            [so_rad[0]*np.cos(so_theta_start[k]),
             so_rad[1]*np.cos(so_theta_start[k])],
            [so_rad[0]*np.sin(so_theta_start[k]),
             so_rad[1]*np.sin(so_theta_start[k])],
            'k',lw=lw
        )
        ax.plot(
            [so_rad[0]*np.cos(so_theta_end[k]),
             so_rad[1]*np.cos(so_theta_end[k])],
            [so_rad[0]*np.sin(so_theta_end[k]),
             so_rad[1]*np.sin(so_theta_end[k])],
            'k',lw=lw
        )
        
        k = k + 1
    
    # Add single-quad electrodes      
    
    count = len(sqel_theta)
    k = 0
    M = 100
    
    # Plot those for y > 0
    
    while k < count/2:
        
        xt = np.linspace(
            sqel_rad[0]*np.cos(sqel_theta[k,0]),
            sqel_rad[0]*np.cos(sqel_theta[k,1]),
            M
        )
        
        ax.plot(xt,np.sqrt(sqel_rad[0]**2 - xt**2),'k',lw=lw)
        
        xt = np.linspace(
            sqel_rad[1]*np.cos(sqel_theta[k,0]),
            sqel_rad[1]*np.cos(sqel_theta[k,1]),
            M
        )
        ax.plot(xt,np.sqrt(sqel_rad[1]**2 - xt**2),'k',lw=lw)
        
        k = k + 1
        
    # Plot those for y < 0
        
    while  k < count:
        
        xt = np.linspace(
            sqel_rad[0]*np.cos(sqel_theta[k,0]),
            sqel_rad[0]*np.cos(sqel_theta[k,1]),
            M
        )
        
        ax.plot(xt,-np.sqrt(sqel_rad[0]**2 - xt**2),'k',lw=lw)
        
        xt = np.linspace(
            sqel_rad[1]*np.cos(sqel_theta[k,0]),
            sqel_rad[1]*np.cos(sqel_theta[k,1]),
            M
        )
        ax.plot(xt,-np.sqrt(sqel_rad[1]**2 - xt**2),'k',lw=lw)
        
        k = k + 1
    
    # Add double-quad electrodes      
    
    count = len(dqel_theta)
    k = 0
    M = 100
    
    # Plot those for y > 0
    
    while k < count/2:
        
        xt = np.linspace(
            dqel_rad[0]*np.cos(dqel_theta[k,0]),
            dqel_rad[0]*np.cos(dqel_theta[k,1]),
            M
        )            
        ax.plot(xt,np.sqrt(float(dqel_rad[0])**2 - xt**2),'k',lw=lw)
        
        xt = np.linspace(
            dqel_rad[1]*np.cos(dqel_theta[k,0]),
            dqel_rad[1]*np.cos(dqel_theta[k,1]),
            M
        )
        ax.plot(xt,np.sqrt(dqel_rad[1]**2 - xt**2),'k',lw=lw)
        
        k = k + 1
        
    # Plot those for y < 0
        
    while  k < count:
        
        xt = np.linspace(
            dqel_rad[0]*np.cos(dqel_theta[k,0]),
            dqel_rad[0]*np.cos(dqel_theta[k,1]),
            M
        )            
        ax.plot(xt,-np.sqrt(dqel_rad[0]**2 - xt**2),'k',lw=lw)
        
        xt = np.linspace(
            dqel_rad[1]*np.cos(dqel_theta[k,0]),
            dqel_rad[1]*np.cos(dqel_theta[k,1]),
            M
        )
        ax.plot(xt,-np.sqrt(dqel_rad[1]**2 - xt**2),'k',lw=lw)
        
        k = k + 1
    
    # Add calorimeters
    
    count = len(cal_theta_start)        
    
    k = 0
    while k < count:
        ax.plot(
            [cal_rad[0]*np.cos(cal_theta_start[k]),
            cal_rad[1]*np.cos(cal_theta_start[k])],
            [cal_rad[0]*np.sin(cal_theta_start[k]),
            cal_rad[1]*np.sin(cal_theta_start[k])],
            'k-',lw=lw
        )
            
        ax.plot(
            [cal_rad[0]*np.cos(cal_box_theta_end[k]),
            cal_rad[1]*np.cos(cal_box_theta_end[k])],
            [cal_rad[0]*np.sin(cal_box_theta_end[k]),
            cal_rad[1]*np.sin(cal_box_theta_end[k])],
            'k-',lw=lw
        )
         
        k = k + 1
        
    
