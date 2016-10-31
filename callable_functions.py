# -*- coding: utf-8 -*-
"""
Created on Fri Sep 30 08:29:11 2016

@author: Eric Schmidt
"""

import numpy as np
#==============================================================================
# Global Constants
#==============================================================================

c = 2.99792458*10**8 # (m/s) Speed of light

#==============================================================================
# Particle movement and position
#==============================================================================

## Calculate the Lorentz force on the particle

def forceDueFields(v,B,E,q):
    
    F = (q*(E + np.cross(v,B)))
    
    return F
    
## Get the electric field based on position

def getElectricField(x):
    
    E = ([0,0,0])
    
    return E
    
## Get the angle of the particle position
    
def getParticleTheta(x):
    
    # Angle of particle position
    theta = np.arctan2(x[1],x[0])
    
    # Translate to 0-2pi
    if theta < 0:
        theta = 2*np.pi + theta
    
    return theta
    
## Get the particle's distance from origin
    
def getParticleRadialPosition(x):
    
    # Radial distance of particle from origin
    r = np.sqrt(x[0]**2 + x[1]**2)
    
    return r
    
## Check if pair-production occurs
    
def ifPairProduction(E,photon_dt):
    
    # Probability of pair production as a function of photon energy
    P = (0.023248 - 9.56748*10**-9 / E**5 + 5.52086*10**-7 / E**4 - 
         0.0000122004 / E**3 + 0.000139708 / E**2 - 0.00125509 / E)
         
    # Adjust the probability based on actual photon_dt as the above equation
    # assumes photon_dt = 10**-11
    P = P * (photon_dt/10**-11)*20
    
    rn = np.random.random()
    
    if rn < P:
        return True
    
## Create electron and positron and note their existance for future tracking

def doPairProduction(E,particle_count,particle_proc,m,v_norm,x,step_counter):
    
    E_part = E/2                    # (eV) Energy of each particle
    p = np.sqrt(E_part**2 - m**2)   # (eV/c) Momentum of each particle
    beta = momentumScalar2Beta(p,m) # Relativistic beta
    gamma = beta2Gamma(beta)      # Relativistic gamma
    v = (p/(gamma*m))*v_norm*c        # (m/s) Velocity vector of each particle
    
    i = 0
    while i < 2:
        particle_count = particle_count + 1
        particle_proc[particle_count] = np.array([x[0],x[1],x[2],
                                                  v[0],v[1],v[2],(-1)**i,0,
                                                  step_counter])
        i = i + 1
    print(p/(10**9))
    return particle_proc,particle_count

## Check if photon is released between energies k_min and k_max while particle
## is traveling through material
    
def isPhotonReleased(k_min,k_max,X0,v,dt,m):
    
    E_tot = velocity2Energy(v,m)
    v_mag = mag(v)
    d = v_mag*dt
    
    # Determine the average number of photons released in distance d
    N = (d/X0) * ((4/3)*np.log(E_tot/k_min) - (4*(E_tot-k_min))/(3*E_tot) +
        (E_tot**2-k_min**2)/(2*E_tot**2))
        
    # Determine the probability that 1,2,3 photon(s) is(are) released
    P1 = poissonProb(N,1)
    P2 = poissonProb(N,2)
    P3 = poissonProb(N,3)
    P4 = poissonProb(N,4)
    P5 = poissonProb(N,5)
    P6 = poissonProb(N,6)
    rn = np.random.rand()
#    print('P1: %0.5f'%P1)
#    print('P2: %0.5f'%P2)
#    print('P3: %0.5f'%P3)
#    print('P4: %0.5f'%P4)
#    print('P5: %0.5f'%P5)
#    print('P6: %0.5f'%P6)
    if rn < P1:
        return 1
    elif rn > P1 and rn < P1 + P2:
        return 2
    elif rn > P1 + P2 and rn < P1 + P2 + P3:
        return 3
    elif rn > P1 + P2 + P3 and rn < P1 + P2 + P3 + P4:
        return 4
    elif rn > P1 + P2 + P3 + P4 and rn < P1 + P2 + P3 + P4 + P5:
        return 5
    elif rn > P1 + P2 + P3 + P4 + P5 and rn < P1 + P2 + P3 + P4 + P5 + P6:
        return 6
    else:
        return 0
        
## Adjust momentum of particle when photon is released
        
def adjustParticleVelocityFromBremsstrahlung(v,k_min,k_max,m):
    
    # Get the magnitude of the velocity vector
    v_mag = mag(v)
    
    # Normalize the velocity vector
    v_vec_norm = v/v_mag
    
    ''' Randomly get the energy of the released photon between k_min and k_max
    based on the Bremsstrahlung intensity plot '''
    
    # Get the smallest wavelength available equal to the particle momentum
    lambda_min = energy2Wavelength(mag(velocity2Energy(v,m)))
    
    # Get the largest wavelength to be used
    lambda_max = energy2Wavelength(k_min)
        
    # Range of possible energies
    dist_range = np.linspace(lambda_min,lambda_max,10000) # (m)
    
    # Weights for the random mean position based on the distribution    
    dist_weights = (dist_range/lambda_min - 1)*(1/dist_range**2)
    
    # Normalize the weights to give a total probability of 1
    dist_weights = dist_weights/sum(dist_weights)
    
    # Randomly choose a mean within allowed range
    l = np.random.choice(dist_range, 1, p=dist_weights)[0]

    # Get photon energy from wavelength
    photon_energy = wavelength2Energy(l)
    
    # Get the energy of the particle before photon release
    E_tot = velocity2Energy(v,m)
    
    # Get the new particle energy
    E_new = E_tot - photon_energy
    
    # Get the new particle momentum scalar
    p_mag_new = energy2Momentum(E_new,m)

    # Get the new velocity scalar    
    v_mag_new = momentumScalar2Beta(p_mag_new,m)*c
    
    # Get the new velocity vector
    v_new = v_mag_new*v_vec_norm
    
    return v_new,photon_energy
    
## If bremsstrahlung procs, release a photon and update appropriate variables
    
def bremsstrahlung(v,m,k_min,k_max,i,photon_count,
                   min_detectable_energy):
#                   photon_energies
    
    p = beta2Momentum(v/c,m)
    Energy = momentum2Energy(p,m)
    
    # The energy of the photon can't be larger thaphoton_energyn the energy of the particle
    if k_max > Energy:
        k_max = Energy
    
    # Readjust the particle velocity from the photon's release
    v,photon_energy = \
        adjustParticleVelocityFromBremsstrahlung(v,k_min,k_max,m)
    
    # Add the location of the bremsstrahlung event
#    brem_array = np.append(brem_array,i)
#    photon_energies = np.vstack((photon_energies,np.array([photon_energy,i])))
    
    # Update momentum
    p = beta2Momentum(v/c,m)
    
    # Update gamma
    gamma = beta2Gamma(v/c)
    
    # Note that a new high energy photon was released
    if photon_energy > min_detectable_energy:
        photon_count[1] = photon_count[1] + 1
    
    # Note that a new photon was released
    photon_count[0] = photon_count[0] + 1
    
    return p,k_max,v,photon_count,gamma,photon_energy

#==============================================================================
# Contact Functions
#==============================================================================

## Checks if inside an element it can pass through
    
def passthroughElementContact(x,el_rad,el_theta):    

    # Get particle's r and theta positions
    
    theta = getParticleTheta(x)
    r = getParticleRadialPosition(x)

    # Check if particle is inside an electrode

    for row in el_theta:
        if row[0] > theta and \
            row[1] < theta and \
            el_rad[0] < r and \
            el_rad[1] > r:
                
            return True

## Checks if contacts the front panel of the calorimeter

def noPassthroughElementContact(x,cal_rad,cal_theta):
    
    # Get particle's r and theta positions
    
    theta = getParticleTheta(x)
    r = getParticleRadialPosition(x)
    
    # Check if particle is inside a calorimeter box
    
    for row in cal_theta:
        if row[0] > theta and \
            row[1] < theta and \
            cal_rad[0] < r and \
            cal_rad[1] > r:
    
            return True
    
## Checks if reaches inner radial limit

def innerLimit(x,R_i):
    if np.sqrt(x[0]**2 + x[1]**2) < R_i:
        return True
    
## Checks if reaches outer radial limit

def outerLimit(x,R):
    if np.sqrt(x[0]**2 + x[1]**2) > R:
        return True

#==============================================================================
# Misc. functions
#==============================================================================

# Updates some tracked variables when particle is inside matter

def updateInsideMatter(v,dt,steps_inside,d_matter):
                
    steps_inside = steps_inside + 1
    d_matter = d_matter + mag(v)*dt
    
    return steps_inside,d_matter

## Poisson probability

def poissonProb(l,n):
    
    P = (l**n * np.exp(-l))/(np.math.factorial(n))
    
    return P

## Get energy of photon from its wavelength

def wavelength2Energy(l):
    
    E = 1240/l
    
    return E

## Get wavelength of photon from its energy

def energy2Wavelength(E):
    
    l = 1240/E
    
    return l

## Get momentum from total energy and mass

def energy2Momentum(E,m):
    
    if E >= m:
        p = np.sqrt(E**2 - m**2)
    else:
        p = 0
    
    return p

## Get total energy from velocity and mass

def velocity2Energy(v,m):
    
    v_mag = mag(v)
    beta = v_mag/c
    p = beta2Momentum(beta,m)
    E_tot = momentum2Energy(p,m)
    
    return E_tot

## Get total energy from momentum and mass

def momentum2Energy(p,m):

    p_mag = mag(p)
    E_tot = np.sqrt(p_mag**2 + m**2)
 
    return E_tot

## Convert beta to momentum

def beta2Momentum(beta,m):

    v = beta*c    
    v_mag = mag(v)
    beta_mag = v_mag/c
#    beta_mag = mag(beta)
    gamma = beta2Gamma(beta_mag)
    p = gamma*m*beta
    
    return p

## Convert total energy to relativistic gamma

def energy2Gamma(E,m):
    
    gamma = E/m
    
    return gamma
    
## Convert relativistic beta to relativistic gamma
    
def beta2Gamma(beta):

    gamma = 1/np.sqrt(1-mag(beta)**2)
    
    return gamma
    
## Convert relativistic gamma to relativistic beta

def gamma2Beta(g):

    beta = np.sqrt((g**2-1)/g**2)
    
    return beta
    
## Convert particle momentum to relativistic beta vector
    
def momentumScalar2Beta(p,m):

    beta = np.sqrt(p**2/(m**2 + p**2))
    
    return beta
    
def momentum2Beta(p,m):
    
    p = np.array([p[0],p[1],p[2]])
    p_norm = mag(p)
    
    beta = 1/np.sqrt((m/p_norm)**2 + 1)
    v = beta*c
    v = v*(p/p_norm)
    beta = v/c

    return beta
    
# Get the magnitude of some vector
    
def mag(v):
    
    return np.sqrt(np.dot(v,v))

#==============================================================================
# Setting initial conditions
#==============================================================================
    
## Return the muon momentum away from the magic momentum at decay

def getMuonMomentumAtDecay():
    
    # Muon magic momentum
#    p_magic = 3.094*10**9 # (eV/c)
    
    m_p = np.zeros((3))
    
    return m_p
    
## Return the particle momentum at decay
    
def getParticleMomentumAtDecay(x,theta,m_x_momentum,m,p_range,p_n):

    # Currently assumes x' = y' = z' = 0

#    p_array = np.linspace([p_range[0],p_range[1],p_n])

    # Magnitude of the particle momentum
    p_tot = p_range[1] - (p_range[1] - p_range[0])*np.random.rand() # (eV/c)
    
    # Convert to momentum vector based on position
    p = np.array([p_tot*np.sin(theta),-p_tot*np.cos(theta),0])
    
    return p
    
## Return the muon beam distribution width
    
def getParticleSigma(m_sigma_amp,m_sigma_ideal,m_sigma_0,num):
    
    sigma = np.zeros(len(m_sigma_amp))    
    
    i = 0
    while i < 2:
    
        if m_sigma_ideal[i] == 1:
            
            sigma[i] = 0
            
        elif m_sigma_ideal[i] == 0:
        
            # Range of possible particle positions. The subraction is due to the fact
            # that the probability function of a cosine goes to infinity at the
            # cosine amplitudes.
            dist_range = np.linspace(-(m_sigma_amp[i] - m_sigma_amp[i]/30),
                                     m_sigma_amp[i] - m_sigma_amp[i]/30,
                                     num[i]) # (m)
            
            # Weights for the random mean position based on the distribution
            
            dist_weights = (np.pi*np.sqrt(1-(dist_range/m_sigma_amp[i])**2))**(-1)
            
            # Normalize the weights to give a total probability of 1
            dist_weights = dist_weights/sum(dist_weights)
            
            # Randomly choose a mean within allowed range
            sigma_temp = np.random.choice(dist_range, 1, p=dist_weights)
            
            # Add xbar to the value around which it oscillates
            sigma[i] = sigma_temp[0] + m_sigma_0[i]
    
        i = i + 1
        
    return sigma
    
## Return the mean particle position related to current particle being tracked
    
def getParticleXBar(m_xbar_amp,m_xbar_ideal,num,m_xbar_0,ax):
    
    if m_xbar_ideal == 1:
        
        xbar = 0
        
    elif m_xbar_ideal == 0:
    
        # Range of possible particle positions. The subraction is due to the fact
        # that the probability function of a cosine goes to infinity at the
        # cosine amplitudes.
        dist_range = np.linspace(-(m_xbar_amp - m_xbar_amp/30),
                                 m_xbar_amp - m_xbar_amp/30,
                                 num) # (m)
        
        # Weights for the random mean position based on the distribution
        
        if ax == 'x':
            dist_weights = (np.pi*np.sqrt(1-(dist_range/m_xbar_amp)**2))**(-1)
            
        if ax == 'y':
            dist_weights = (np.pi*np.sqrt(1-(dist_range/m_xbar_amp)**2))**(-1)
        
        # Normalize the weights to give a total probability of 1
        dist_weights = dist_weights/sum(dist_weights)
        
        # Randomly choose a mean within allowed range
        xbar = np.random.choice(dist_range, 1, p=dist_weights)
        
        # Add xbar to the value around which it oscillates
        xbar = xbar[0] + m_xbar_0
    
    return xbar
    
## Return the local muon x-position based on distribution function

def getXParticlePositions(xbar,sigma,xmin_init,xmax_init,num,fit):
    
    # Range of possible particle positions
    init_dist_range = np.linspace(xmin_init, xmax_init, num) # (m)
    dist_range = xbar + init_dist_range # (m)
    
    # Weights for the random position based on the distribution 
    
    if sigma != 0:
    
        if fit == "Gaussian":
            dist_weights = (np.sqrt(2*np.pi*sigma**2)**(-1)) * \
                np.exp(-((init_dist_range-xbar)**2)/(2*sigma**2))
    
        # Normalize the weights to give a total probability of 1
        dist_weights = dist_weights/sum(dist_weights)
        
        # Randomly choose a position within allowed range
        particle_position = np.random.choice(dist_range, 1, p=dist_weights)
        
        # Convert from single-element array to float
        particle_position = particle_position[0]
        
    elif sigma == 0:
        
        particle_position = xbar
    
    return particle_position
    
## Return the local muon y-position based on distribution function
    
def getYParticlePositions(ybar,sigma,ymin_init,ymax_init,num,fit):
    
    # See 'getXParticlePositions' for comments
    
#    init_dist_range = np.linspace(ymin_init, ymax_init, num)
#    dist_range = ybar + init_dist_range
#    
#    dist_weights = 1 * init_dist_range
#    dist_weights = dist_weights/sum(dist_weights)
#    particle_position = np.random.choice(dist_range, 1, p=dist_weights)
#    
#    particle_position = particle_position[0]
    
    particle_position = 0
    
    return particle_position