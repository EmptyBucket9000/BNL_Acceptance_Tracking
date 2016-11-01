# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:51 2016

@author: Eric Schmidt
"""

import callable_functions as cf
import numpy as np
import random
    
def muon(theta_set,theta_min,theta_max,xbar_amp,xbar_0,xphi,
             sigmaphi,sigma_amp,sigma_0,xbar_ideal,sigma_ideal,theta,
             xprimemax,m_xlimit):

    '''Unless otherwise stated, 3-element arrays are of the form x,y,z.'''
    
    # 'fit' is the type of function (in local coordinates) that will fit the
    # muon distribution density. Curreantly 'Gaussian' is the only option.
    # ([x,y])
    xfit = np.array(["Gaussian","Gaussian"])
    
    # Number of possible muon positions within range, used for building the
    # muon distribution arrays
    xnum = np.array([1000,1000])
    sigmanum = np.array([1000,1000])
    
    # Prime array
    xprime = np.zeros((2))

    ''' Minimum/maximum muon position in beam '''

    xmin = np.zeros((2))                # (m)
    xmax = np.zeros((2))                # (m)
    
    # (m) Mean muon position (centroid position)
    xbar_x = cf.getParticleXBar(xbar_amp[0],xbar_ideal[0],
                                     xnum[0],xbar_0[0],'x')
    xbar_y = cf.getParticleXBar(xbar_amp[1],xbar_ideal[1],
                                     xnum[1],xbar_0[1],'y')
    
    xbar = np.array([xbar_x,xbar_y])
    
    if sigma_ideal[0] == 1:
        xmin[0] = xbar[0]             # (m)
        xmax[0] = xbar[0]             # (m)
    else:
        xmin[0] = xbar[0] - m_xlimit # (m)
        xmax[0] = xbar[0] + m_xlimit # (m)
    
    if theta_set == 0:
        theta = (theta_max - theta_min)*random.random() + \
                theta_min

# Beam distribution width
    
    sigma = cf.getParticleSigma(sigma_amp,sigma_ideal,
                                     sigma_0,sigmanum)
    
    ''' Set the muon position array '''

    # x-position
    pos_x = cf.getXParticlePositions(
        xbar[0],sigma[0],xmin[0],xmax[0],xnum[0],
        xfit[0]
    )
    
    # y-position
    pos_y = cf.getYParticlePositions(
        xbar[1],
        sigma[1],
        xmin[1],xmax[1],xnum[1],
        xfit[1]
    )
    # Full position vector
    x_pos = np.array([pos_x,pos_y])
    
    # Muon x, y, z prime values at decay (dx/ds)
    xprime[0] = cf.getParticleXPrime(xbar[0],pos_x,m_xlimit,xprimemax[0])
    
    return np.array([x_pos,sigma,theta,xbar,xprime])