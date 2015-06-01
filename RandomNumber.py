"""
Created on Sat Sep 21 23:44:35 2013
"""

# Pseudo-Random Number Generator
"""
The following random number methods are below.

-------- Uniform Distribution --------
1. Lehmer modulo generator

-------- Gaussian Distribution --------
2. Polar method generator
3. Box-Muller generator
"""

import numpy as np
from math import sqrt,log,cos,sin,pi


# 1
# Lehmer modulo generator
#*************************************************************
# NOTE: this generates a random uniform distribution over [0,1]
#-----------------------------------------------------------
# seed = initial seed to start the random sequence
# n = how many random numbers to generate in array
#----------------------------------------------------------- 
def randUniform(seed,n):
    # define constants to use in algorithm
    a = 16807.0
    xm = 2147483647.0
    x0 = 2147483711.0    
    
    # array to store random uniform numbers
    rand = []
    
    # loop through modulo routine to generate numbers
    for i in range(n):
        seed = (a*seed) % xm
        rand.append(seed/x0)
        
    return rand
#*************************************************************
    
    
# 2
# Polar method generator
#*************************************************************
# NOTE: this generates a random normal distribution over [0,1]
#       from a random uniform generator
#-----------------------------------------------------------
# n = how many random numbers to generate in array
#----------------------------------------------------------- 
def randNormal1(n): 
    # array to store random normal numbers
    rand = []
    
    count = 0
    # loop to generate random uniform numbers
    # and turn into random normal numbers
    while count < n:
        # get two random uniform numbers to start polar method 
        # NOTE: uses built-in numpy.random.uniform 
        u1 = np.random.uniform(0,1,1)[0]
        u2 = np.random.uniform(0,1,1)[0]
        
        # define two new values to start polar method
        v1 = 2.0*u1 - 1.0
        v2 = 2.0*u2 - 1.0
        
        # define check value
        s = v1**2 + v2**2
        
        # check condition
        if s < 1:
            # define random normal numbers
            rand1 = v1*sqrt((-2.0*log(s))/s)
            rand2 = v2*sqrt((-2.0*log(s))/s)
            
            # append to random array
            rand.append(rand1)
            rand.append(rand2)
            
            # increment counter
            count = count + 1
        
    return rand
#*************************************************************
    
    
# 3
# Box-Muller generator
#*************************************************************
# NOTE: this generates a random normal distribution over [0,1]
#       from a random uniform generator
#-----------------------------------------------------------
# n = how many random numbers to generate in array
#----------------------------------------------------------- 
def randNormal2(n): 
    # array to store random normal numbers
    rand = []
    
    count = 0
    # loop to generate random uniform numbers
    # and turn into random normal numbers
    while count < n:
        # get two random uniform numbers to start method 
        # NOTE: uses built-in numpy.random.uniform 
        u1 = np.random.uniform(0,1,1)[0]
        u2 = np.random.uniform(0,1,1)[0]
        
        # define random normal numbers
        rand1 = sqrt(-2.0*log(u1))*cos(2*pi*u2)
        rand2 = sqrt(-2.0*log(u1))*sin(2*pi*u2)
        
        # append to random array
        rand.append(rand1)
        rand.append(rand2)
        
        # increment counter
        count = count + 1
        
    return rand
#*************************************************************