# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 09:03:04 2016

@author: schmidte
"""

import numpy as np
import matplotlib.pyplot as plt
import random

SecR = np.pi/4 # (rad) Radians in each section. 1 calorimeter -> 1 section.
N = 1 # Number of Muons in system
L = 0.3 # (m) Length of the calorimeter
D = 0.6 # (m) Distance from the optimum beam to outside of calorimeter
R = 7.112 # (m) Radius of optimum beam position
d = 0.045 # (m) Maximum distance from R to outer muon
c = 2.99792458*10**8 # (m/s) Speed of light
#q = 1 # (eV) Charge of an electron
q = -1.60217662*10**(-19) # (C) Charge of an electron
qv = -1*10**(-9) # (GeV) Charge of an electron
mv = 0.510998928*10**(-3) # (GeV/c**2) Mass of an electron
m = 9.10938356*10**(-31) # (Kg) Mass of an electron
Mv = 105.6583715*10**(-3) # (Gev/c**2) Mass of a muon
M = 1.883531475*10**(-28) # (Kg) Mass of a muon
a_mu = 11659208*10**(-10) # Initial anomalous mag moment for finding p
w_c = 6.7 # (MHz) Muon cyclotron frequency
p = (Mv/np.sqrt(a_mu)) # (GeV/c) Magic momentum of the muon
beta = np.sqrt((p**2)/(p**2 + Mv**2))
gamma = 1/np.sqrt(1-(beta)**2)
#v = beta*c # (m/s) Velocity of the muon
#beta = 1-1/gamma**2
#B = np.abs(gamma*M*v/(q*R)) # (T) Magnetic field
B = np.abs(gamma*Mv*beta*c/(c**2*qv*R)) # (T) Magnetic field
print(B)

def main():
    
    i = 0
    while i < N:
        E = getElectronEnergyAtDecay()
        
        s = getPositionAtDecay()
        
        beta = energy2Beta(E,m)
        
        i = i + 1
    
def getPositionAtDecay():
    
    rPos = random.uniform(R - d, R + d) # (m) Generates radial position
    s = random.uniform(0,SecR) # (m) Generates linear position inside section
    
    return s

def updateElectronPostion(x,E):
    
    beta = energy2Beta(E,m)
    
    return beta
    
def getElectronEnergyAtDecay():
    
    E = 52.8 # (MeV)
#    E = 8.459492263e-12 # (J) Currently set to E_(e,max)
    return E
    
def energy2Beta(E,m):
    
#    p = np.sqrt(E**2 - m**2*c**4)/c
#    v = np.sqrt((p**2)/((p**2)/(c**2) + (m**2)/(c**4)))
    
    p = np.sqrt(E**2 - mv**2)
    beta = np.sqrt((p**2)/(p**2 + mv**2))
    
    return beta
    
main()