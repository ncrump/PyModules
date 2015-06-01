"""
Created on Thu Oct 10 23:37:58 2013
"""

# Random Walks
"""
Simulation of a random walk on a lattice of evenly spaced points. 
1. 2D Random Walk
2. 3D Random Walk 
"""

import numpy as np


# 1
# 2D Random Walk
#*******************************************************************
#-------------------------------------------------------------------
# N = number of steps to take on random walk
#------------------------------------------------------------------- 
def RandomWalk2D(N):
    
    # number of steps to take
    Nsteps = range(N)
    
    # define start position at origin of 2-D x,y grid
    pos = [0,0]
    xp = [0]
    yp = [0]
    
    # generate array of random uniform numbers between [0,1)
    rand = np.random.uniform(0,1,N)
        
    # loop through array of random numbers for random walk
    for i in Nsteps:
        # depending on random number go right, left, up, or down by 1 unit
        if 0.00 <= rand[i] < 0.25: pos[0] = pos[0]+1
        if 0.25 <= rand[i] < 0.50: pos[0] = pos[0]-1
        if 0.50 <= rand[i] < 0.75: pos[1] = pos[1]+1
        if 0.75 <= rand[i] < 1.00: pos[1] = pos[1]-1
            
        # append x,y positions
        xp.append(pos[0])
        yp.append(pos[1])
        
    # get final displacement from start position
    dist = (xp[N]**2 + yp[N]**2)**0.5
        
    return xp,yp,dist
#*******************************************************************
    
    
# 2
# 3D Random Walk
#*******************************************************************
#-------------------------------------------------------------------
# N = number of steps to take on random walk
#------------------------------------------------------------------- 
def RandomWalk3D(N):
    
    # number of steps to take
    Nsteps = range(N)
    
    # define start position at origin of 3-D x,y,z cube
    pos = [0,0,0]
    xp = [0]
    yp = [0]
    zp = [0]
    
    # generate array of random uniform numbers between [0,1)
    rand = np.random.uniform(0,1,N)
        
    # loop through array of random numbers for random walk
    for i in Nsteps:
        # depending on random number go +/-x, +/-y, or +/-z by 1 unit
        if 0.00000 <= rand[i] < 1.0/6.0: pos[0] = pos[0]+1
        if 1.0/6.0 <= rand[i] < 2.0/6.0: pos[0] = pos[0]-1
        if 2.0/6.0 <= rand[i] < 3.0/6.0: pos[1] = pos[1]+1
        if 3.0/6.0 <= rand[i] < 4.0/6.0: pos[1] = pos[1]-1
        if 4.0/6.0 <= rand[i] < 5.0/6.0: pos[2] = pos[2]+1
        if 5.0/6.0 <= rand[i] < 6.0/6.0: pos[2] = pos[2]-1
            
        # append x,y,z positions
        xp.append(pos[0])
        yp.append(pos[1])
        zp.append(pos[2])
        
    # get final displacement from start position
    dist = (xp[N]**2 + yp[N]**2 + zp[N]**2)**0.5
        
    return xp,yp,zp,dist
#*******************************************************************