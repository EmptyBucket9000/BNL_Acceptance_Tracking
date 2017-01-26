# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 09:52:12 2017

@author: Eric Schmidt
"""

'''
Returns the possible values of the muon phase-space position and number of
muons desired from those available. If adding a new arg value, you must
determine what the maximum value of N can be.
'''

import numpy as np

def do(arg):
    
    # Ranges, set all elements to zero for full range or set specific limits
    # to allow limited variability.
    
    if arg == '1':
        
        # All muons
    
        x_pos_range = np.array([0,0])/100
        x_prime_range = np.array([0,0])
    
        y_pos_range = np.array([0,0])/100
        y_prime_range = np.array([0,0])
    
        N = 5171
        
    elif arg == '2':
        
        # For x = [-4.5,-2.5] and all x',y,y'
    
        x_pos_range = np.array([-4.5,-2.5])/100
        x_prime_range = np.array([0,0])
    
        y_pos_range = np.array([0,0])/100
        y_prime_range = np.array([0,0])
    
        N = 223
        
    elif arg == '3':
        
        # For x = [2.5,4.5] and all x',y,y'
    
        x_pos_range = np.array([2.5,4.5])/100
        x_prime_range = np.array([0,0])
    
        y_pos_range = np.array([0,0])/100
        y_prime_range = np.array([0,0])
    
        N = 386
        
    elif arg == '4':
        
        # For x' = [-0.005,-0.002] and all x,y,y'
    
        x_pos_range = np.array([0,0])/100
        x_prime_range = np.array([-0.005,-0.002])
    
        y_pos_range = np.array([0,0])/100
        y_prime_range = np.array([0,0])
    
        N = 356
        
    elif arg == '5':
        
        # For x' = [0.002,0.005] and all x,y,y'
    
        x_pos_range = np.array([0,0])/100
        x_prime_range = np.array([0.002,0.005])
    
        y_pos_range = np.array([0,0])/100
        y_prime_range = np.array([0,0])
    
        N = 404
        
    elif arg == '6':
        
        # For y = [-5,-2.5] and all x,x',y'
    
        x_pos_range = np.array([0,0])/100
        x_prime_range = np.array([0,0])
    
        y_pos_range = np.array([-5,-2.5])/100
        y_prime_range = np.array([0,0])
    
        N = 143
        
    elif arg == '7':
        
        # For y = [2.5,5] and all x,x',y'
    
        x_pos_range = np.array([0,0])/100
        x_prime_range = np.array([0,0])
    
        y_pos_range = np.array([2.5,5])/100
        y_prime_range = np.array([0,0])
    
        N = 217
        
    elif arg == '8':
        
        # For y' = [-0.005,-0.0015] and all x,x',y
    
        x_pos_range = np.array([0,0])/100
        x_prime_range = np.array([0,0])
    
        y_pos_range = np.array([0,0])/100
        y_prime_range = np.array([-0.005,-0.0015])
    
        N = 182
        
    elif arg == '9':
        
        # For y' = [0.0015,0.005] and all x,x',y
    
        x_pos_range = np.array([0,0])/100
        x_prime_range = np.array([0,0])
    
        y_pos_range = np.array([0,0])/100
        y_prime_range = np.array([0.0015,0.005])
    
        N = 200
    
    return x_pos_range,x_prime_range,y_pos_range,y_prime_range,N