# BNL_Admittance_Tracking

Two coordinate systems are used, global and local:
 
    Global is the entire ring in 3-dimensions cartesian coordinates where the
    origin is the center of the ring and the z-direction indicates 'up' and
    'down', from the perspective of a person standing in the ring, and is
    parallel to the y-direction in the local coordinate system.
     
    Local is a 2-dimensional x,y cartesian coordinate where the origin is the
    beam centroid, the negative x-direction points towards the center of the
    ring and the y-direction is parallel to the z-direction in the global
    system. Local is primarily used in determining particle initial conditions,
    i.e. the conditions describing the muon at decay
 
Variable names formatted 'm_*' are for muons, all others are for particles
 
The muon density function is Gaussian, following:
    (np.sqrt(2*np.pi*sigma**2)**(-1))*np.exp(-((x-xbar)**2)/(2*sigma**2))
 
Magic momentum is set to 3.09435 GeV/c giving a magnetic field of 1.4513 T
for a radius of 7.112 m.
