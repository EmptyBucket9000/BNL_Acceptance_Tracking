# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 08:00:03 2016

@author: schmidte
"""

import numpy as np
import matplotlib.pyplot as plt
import csv

def main():

    file_name = "EndOfTracking_phase_space.csv"
    i = 0

    with open(file_name, "rt") as inf:
        reader = csv.reader(inf, delimiter=' ')
        next(reader, None)  # skip the headers
    
        stuff = list(reader)
        N = len(stuff)
        print(N)
        m_x = np.zeros((N,2))
        m_xprime = np.zeros((N,2))
        
        xprime_min = 0.0025989
        yprime_min = 0.0014243
        
        xprime_max = 0.0025989
        yprime_max = 0.0014243
        
        for row in stuff:
            m_x[i,0] = row[0]                
            m_xprime[i,0] = row[1]
            
            if np.abs(m_xprime[i,0]) < xprime_min:
                xprime_min = np.copy(np.abs(m_xprime[i,0]))
            
            if np.abs(m_xprime[i,0]) > xprime_max:
                xprime_max = np.copy(np.abs(m_xprime[i,0]))
                
            m_x[i,1] = row[2]                
            m_xprime[i,1] = row[3]
            
            if np.abs(m_xprime[i,1]) < yprime_min:
                yprime_min = np.copy(np.abs(m_xprime[i,1]))
            
            if np.abs(m_xprime[i,1]) > yprime_max:
                yprime_max = np.copy(np.abs(m_xprime[i,1]))
                
            i = i + 1
            
    print('xprime min: %0.8f'%xprime_min)
    print('xprime max: %0.8f'%xprime_max)
    
    print('yprime min: %0.8f'%yprime_min)
    print('yprime max: %0.8f'%yprime_max)
    
    plt.figure(1)
    ax = plt.subplot(1,1,1)
    
    ax.plot(m_x[:,0],m_xprime[:,0],'.')
    ax.plot(m_x[:,1],m_xprime[:,1],'.')
    
    plt.figure(2)
    ax = plt.subplot(1,1,1)
    
    ax.hist(m_xprime[:,0],30)
    ax.hist(m_xprime[:,1],30)
    
    plt.figure(3)
    ax = plt.subplot(1,1,1)
    
    ax.hist(m_x[:,0],30)
    ax.hist(m_x[:,1],30)
    
    
if __name__ == '__main__':
    
    main()