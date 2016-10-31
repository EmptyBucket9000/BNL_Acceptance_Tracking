# -*- coding: utf-8 -*-
"""
Created on Tue Oct 18 08:00:03 2016

@author: schmidte
"""

import numpy as np
import matplotlib.pyplot as plt

#N = 100000
#i = 0
#val = np.zeros((N))
#x = np.linspace(0,1,100000)
#
#while i < N:
#    
#    val[i] = np.random.random()
#    i = i + 1
#    
#plt.plot(x,val,'.',markersize=0.1)

def main():

    lambda_min = energy2Wavelength(3*10**9)
    
    # Get the largest wavelength to be used
    lambda_max = energy2Wavelength(0.0005*10**9)
        
    # Range of possible energies
    dist_range = np.linspace(lambda_min,lambda_max,10000) # (m)
    
    # Weights for the random mean position based on the distribution    
    dist_weights1 = (dist_range/lambda_min - 1)*(1/dist_range**2)
    
    # Normalize the weights to give a total probability of 1
    dist_weights = dist_weights1/sum(dist_weights1)
    
    # Randomly choose a mean within allowed range
    l = np.random.choice(dist_range, 1, p=dist_weights)[0]

    print('Wavelength: %e'%l)
    print('Energy: %0.3f'%(1240/(l*10**9)))
    
    plt.plot(dist_range,dist_weights1)


def energy2Wavelength(E):
    
    l = 1240/E
    
    return l
    
    
main()