"""
Created on Fri Feb 28 23:36:46 2014
"""

# Lennard-Jones Potential
"""
Calculates the Lennard-Jones interatomic potential for a system of atoms. 
"""


#****************************************************
#----------------------------------------------------
# x,y,z = arrays of particle positions
# eps = epsilon (depth of potential well)
# sig = sigma (distance of zero potential)
#----------------------------------------------------
def LJ(x,y,z,eps,sig):
    
    N = len(x)
    LJsum = 0
    
    # loop through atoms ij
    for i in range(N-1):
        for j in range(i+1,N):
            
            # calculate distance between atoms i and j
            r = ((x[j]-x[i])**2 + (y[j]-y[i])**2 + (z[j]-z[i])**2)**0.5
            # sum LJ terms to get potential
            LJsum = LJsum + (sig/r)**12 - (sig/r)**6
            
    return 4.0*eps*LJsum
#****************************************************