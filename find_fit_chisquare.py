# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 07:14:28 2016

@author: Eric Schmidt
"""

'''
Find the best fit via the minimum chi-square value. Return the results to the
calling function.

See README.md for more information.

'''

import numpy as np
from scipy.optimize import minimize
    
def fit(x,data,err,init_guess,fit_type):
    
    # Change fit_type to lower case to ensure consistancy.
    fit_type = fit_type.lower()
    
    # Number of elements in the guess array, this determines what order
    # polynomial will be used.
    l = len(init_guess)
        
    # Creates an array of 'powers', one for each order in the polynomial.
    p = range(l)
    
    # Call the minimize function and return its results.
    minimum = minimize(f,init_guess,args=(x,data,err,fit_type,l,p))
    
    uncertainties = unc(x,minimum.x,minimum.fun,data,err,l,p)
    
    return minimum,uncertainties

def unc(x,var,x2,data,err,l,p):
    
    # Find the uncertainties in the variables via the chi-square +1 method
    
    uncertainties = np.zeros((l))
        
    i = 0
    
    while i < l:
        
        # Make a copy as we will be changing it
        temp_var = np.copy(var)
        
        while True:
            
            j = 0
            fun = 0
            
            # Rebuild the function with the current constants
            while j < l:
                fun = fun + temp_var[j]*x**p[j]
                j = j + 1
            
            # Find the chi-square value
            val = np.sum(((data - fun)/err)**2)
            
            # Check if the current chi-square value is larger than the chi-
            # square minimum + 1, if so, exit the loop, if not, add a small
            # value to the currently focused coefficient in the function.
            if val < x2 + 1:
                temp_var[i] = temp_var[i] + temp_var[i]*0.001
                
            else:

                uncertainties[i] = temp_var[i]
                break
            
        i = i + 1
        
    uncertainties = np.abs(uncertainties - var)
    
    return uncertainties
                
                

def f(init_guess,x,data,err,fit_type,l,p):
    
    if fit_type == "poly" or fit_type == "polynomial":
        
        i = 0
        
        # Set the initial function to zero.
        fun = 0
        
        # For each order in the polynomial, add that order to the function.
        while i < l:
            fun = fun + init_guess[i]*x**p[i]
            i = i + 1
            
    # Final function to be minimized and returned.
    return np.sum(((data - fun)/err)**2)