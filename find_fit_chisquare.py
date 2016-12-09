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
    
    # Call the minimize function and return its results.
    return minimize(f,init_guess,args=(x,data,err,fit_type))

def f(init_guess,x,data,err,fit_type):
    
    if fit_type == "poly" or fit_type == "polynomial":
    
        # Number of elements in the guess array, this determines what order
        # polynomial will be used.
        l = len(init_guess)
        
        # Creates an array of 'powers', one for each order in the polynomial.
        p = range(l)
        
        i = 0
        
        # Set the initial function to zero.
        fun = 0
        
        # For each order in the polynomial, add that order to the function.
        while i < l:
            fun = fun + init_guess[i]*x**p[i]
            i = i + 1
    
    # Final function to be minimized and returned.
    return np.sum(((data - fun)/err)**2)