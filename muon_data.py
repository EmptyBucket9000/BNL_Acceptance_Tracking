# -*- coding: utf-8 -*-
"""
Created on Wed Nov 02 07:02:36 2016

@author: Eric Schmidt
"""

'''
See README.md for information.
'''

import numpy as np
import csv
import random
            
def muon(N,file_name,p_magic):
    
    m,length = readFile(file_name)
    lines = random.sample(range(length), N)
    m_x = np.zeros((N,2))
    m_p = np.zeros((N,3))
    i = 0
    for line in lines:
        m_p[i,2] = p_magic + (m[line,4]*p_magic)
        m_x[i,0] = m[line,0]
        m_p[i,0] = m[line,1]
        m_x[i,1] = m[line,2]
        m_p[i,1] = m[line,3]
        i = i + 1
        
    return m_x,m_p

def readFile(file_name):
    
    # csv: [x,xprime,y,yprime,p_z]
    i = 0

    with open(file_name, "rt") as inf:
        reader = csv.reader(inf, delimiter=',')
        next(reader, None)  # skip the headers
    
        stuff = list(reader)
        length = len(stuff)
        m = np.zeros((length,5))
        
        for row in stuff:
            m[i,0] = row[0]           
            m[i,1] = row[1]
            m[i,2] = row[2]
            m[i,3] = row[3]
            m[i,4] = row[4]
                
            i = i + 1
            
        return m,length