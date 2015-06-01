"""
Created on Thu Oct 10 23:37:58 2013
"""

# Random Sphere
"""
Generates random uniformly distributed points on the surface of a sphere. 
"""

import numpy as np


#****************************************************
#----------------------------------------------------
# N = number of random points to generate on sphere
# R = radius of sphere
#----------------------------------------------------
def RandomSphere(N, R):
    
    # generates randomly distributed theta,phi points at r = 1
    r = np.ones(N)*R
    # this prevents clumping at poles for theta,phi
    theta0 = np.random.uniform(-1,1,N)    
    theta = np.arccos(theta0)
    phi = np.random.uniform(0,2*np.pi,N)
    
    # converts spherical [r,theta,phi] points to cartesian [x,y,z] points
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)

    return x,y,z
#****************************************************