"""
Created on Fri Feb 28 23:36:46 2014
"""

# Cubic Lattice
"""
Generates positions of atoms that form clusters of the following unit cells:
1. Simple Cubic (SC)
2. Body Centered Cubic (BCC)
3. Face Centered Cubic (FCC) 
"""

import numpy as np


# 1
# simple cubic lattice (SC)
#*******************************************************************
#-------------------------------------------------------------------
# d = distance between atoms
# N = number of times to translate cell to form lattice
#-------------------------------------------------------------------
def SC(d, N):   
    
    # get translation step of cube
    steps = np.linspace(d, d*N, N)
    
    # get initial unit cell positions
    x = np.zeros(1)
    y = np.zeros(1)
    z = np.zeros(1)
    
    # translate in x direction
    for xi in steps:
        x = np.append(x,x[0:1]+xi)
        y = np.append(y,y[0:1])
        z = np.append(z,z[0:1])
    
    # translate in y direction
    for yi in steps:
        x = np.append(x,x[0:N+1])
        y = np.append(y,y[0:N+1]+yi)
        z = np.append(z,z[0:N+1])
            
    # translate in z direction
    for zi in steps:
        x = np.append(x,x[0:(N+1)**2])
        y = np.append(y,y[0:(N+1)**2])
        z = np.append(z,z[0:(N+1)**2]+zi)    
    
    return x,y,z
#*******************************************************************
    

# 2  
# body centered cubic lattice (BCC)
#*******************************************************************
#-------------------------------------------------------------------
# d = distance between atoms
# N = number of times to translate cell to form lattice
#-------------------------------------------------------------------
def BCC(d, N):   
        
    # get translation step of cube
    a = d*(2.0**0.5)
    steps = np.linspace(a, a*N, N)
        
    # get initial unit cell positions
    x = np.array([0,1])*0.5*a
    y = np.array([0,1])*0.5*a
    z = np.array([0,1])*0.5*a
    
    # translate in x direction
    for xi in steps:
        x = np.append(x,x[0:2]+xi)
        y = np.append(y,y[0:2])
        z = np.append(z,z[0:2])
    
    # translate in y direction
    for yi in steps:
        x = np.append(x,x[0:2*(N+1)])
        y = np.append(y,y[0:2*(N+1)]+yi)
        z = np.append(z,z[0:2*(N+1)])
            
    # translate in z direction
    for zi in steps:
        x = np.append(x,x[0:2*(N+1)**2])
        y = np.append(y,y[0:2*(N+1)**2])
        z = np.append(z,z[0:2*(N+1)**2]+zi)    
    
    return x,y,z
#*******************************************************************
    

# 3 
# face centered cubic lattice (FCC)
#*******************************************************************
#-------------------------------------------------------------------
# d = distance between atoms
# N = number of times to translate cell to form lattice
#-------------------------------------------------------------------
def FCC(d, N):   
        
    # get translation step of cube
    a = d*(2.0**0.5)
    steps = np.linspace(a, a*N, N)
        
    # get initial unit cell positions
    x = np.array([0,1,0,1])*0.5*a
    y = np.array([0,1,1,0])*0.5*a
    z = np.array([0,0,1,1])*0.5*a
    
    # translate in x direction
    for xi in steps:
        x = np.append(x,x[0:4]+xi)
        y = np.append(y,y[0:4])
        z = np.append(z,z[0:4])
    
    # translate in y direction
    for yi in steps:
        x = np.append(x,x[0:4*(N+1)])
        y = np.append(y,y[0:4*(N+1)]+yi)
        z = np.append(z,z[0:4*(N+1)])
            
    # translate in z direction
    for zi in steps:
        x = np.append(x,x[0:4*(N+1)**2])
        y = np.append(y,y[0:4*(N+1)**2])
        z = np.append(z,z[0:4*(N+1)**2]+zi)    
    
    return x,y,z
#*******************************************************************