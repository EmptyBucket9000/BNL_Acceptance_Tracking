# -*- coding: utf-8 -*-
"""
Created on Wed Nov 02 07:02:36 2016

@author: Eric Schmidt
"""

'''
Loads the muon data from the .csv and returns it to 'main.py'. Muon local
coordinates are used.

See README.md for more information.
'''

import numpy as np
import csv
import random
            
def muon(N,file_name,p_magic,x_pos_range,x_prime_range):
    
    # Reads the muon data from the .csv then picks the desired number of muons
    # at random from the returned list.
    
    m,length = readFile(file_name,x_pos_range,x_prime_range)
    lines = random.sample(range(length), N)
    m_x = np.zeros((N,2))               # Muon position
    m_p = np.zeros((N,3))               # Muon momentum
    i = 0
    for line in lines:
        
        # Convert the linear momentum from the file, which is provided as a
        # fraction away from the magic momentum.
        m_p[i,2] = p_magic + (m[line,4]*p_magic)
        m_x[i,0] = m[line,0]                # x-position
        m_p[i,0] = m[line,1]                # x-momentum
        m_x[i,1] = m[line,2]                # y-position
        m_p[i,1] = m[line,3]                # y-momentum
        i = i + 1
        
    return m_x,m_p

def readFile(file_name,x_pos_range,x_prime_range):
    
    # csv: [x,xprime,y,yprime,p_z]
    i = 0

    with open(file_name, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
    
        stuff = list(reader)
        
        # Get the total number of muons available
        length = len(stuff)
        
        # Pre-allocate the muon array
        m = np.zeros((length,5))
        
        # Returns only the muons that exist between within the x-position
        # range defined in 'main.py'.
        
        for row in stuff:
        
            if float(row[0]) >= x_pos_range[0] and \
                float(row[0]) <= x_pos_range[1] and \
                float(row[1]) >= x_prime_range[0] and \
                float(row[1]) <= x_prime_range[1]:

                m[i,0] = row[0]           
                m[i,1] = row[1]
                m[i,2] = row[2]
                m[i,3] = row[3]
                m[i,4] = row[4]
                    
                i = i + 1
        
        # Remove any muon rows that did not meet the range requirements.
        m = m[np.any(m != 0,axis=1)]
        length = len(m)
            
        return m,length