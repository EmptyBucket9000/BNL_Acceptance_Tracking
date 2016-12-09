# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 07:14:28 2016

@author: Eric Schmidt
"""

import numpy as np
from scipy.optimize import minimize
    
def fit(x,data,err,init_guess,fit_type):
    
    fit_type = fit_type.lower()
    
    return minimize(f,init_guess,args=(x,data,err,fit_type))

def f(init_guess,x,data,err,fit_type):
    
    if fit_type == "poly" or fit_type == "polynomial":
    
        l = len(init_guess)
        p = range(l)
        
        i = 0
        fun = 0
        while i < l:
            fun = fun + init_guess[i]*x**p[i]
            i = i + 1
    
    return np.sum(((data - fun)/err)**2)