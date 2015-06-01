"""
Created on Tue Jun 04 17:40:50 2013
"""

# Interpolation
"""
The following interpolation methods are below.
1. Lagrange method
2. Neville method
3. Lagrange Piecewise method
4. Cubic Spline method
"""

import numpy as np
import sympy as sym


# 1
# Lagrange Interpolation method
#*******************************************************************
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use (degree of polynomial)
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------  
def lagrange(x, y, nPoints, nxIncrement):
    t = sym.Symbol('t')  # define symbolic variable
    n = range(nPoints)   # index for interpolation 
    p = 0
    
    for i in n:
        term = 1
        
        # Lagrange Method
        for j in [k for k in n if k != i]:
            term = term * ((t - x[j]) / (x[i] - x[j]))
            
        p = p + y[i]*term  # build symbolic interpolating polynomial
    
    fInterp = sym.lambdify(t,p)                                  # turn symbolic poly into real poly    
    xx = np.arange(x[0], x[nPoints-1]+nxIncrement, nxIncrement)  # create x interped values
    yy = [fInterp(z) for z in xx]                                # get y interped values
    
    return xx,yy
#*******************************************************************
    
    
# 2
# Neville Interpolation method
#*******************************************************************
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use (degree of polynomial)
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------
def neville(x, y, nPoints, nxIncrement):
    xx = np.arange(x[0], x[nPoints-1]+nxIncrement, nxIncrement)  # create x interped values
    yy = []  # for storing y interped values
    
    for xval in xx:
    
        arry = y[0:nPoints]                # get intial y values
        arr = np.zeros((nPoints,nPoints))  # create Neville matrix
        arr[0] = arry                      # write y values to matrix
        
        xpos = 0
        
        for i in range(nPoints - 1): 
            xpos = xpos + 1
             
            # Neville Method & build matrix
            for j in range(nPoints - 1 - i):
                term = ((xval-x[j+xpos])*arr[i][j] - (xval-x[j])*arr[i][j+1]) / (x[j]-x[j+xpos])
                arr[i+1][j] = term         # build Neville matrix
                   
        yy.append(term) 
    
    return xx,yy
#*******************************************************************
    

# 3
# Lagrange Piecewise Interpolation method
#*******************************************************************
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use (degree of polynomial)
# nLocal = number of local points to use for piecewise
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------
def lagrangePW(x, y, nPoints, nLocal, nxIncrement):
    t = sym.Symbol('t')                       # define symbolic variable
    n = range(0, nPoints-nLocal+1, nLocal-1)  # create x interped values
    
    xx = []  # for storing x interped values
    yy = []  # for storing y interped values
    
    # iterate over data points
    for loc in n:
        p = 0
        pts = range(loc, loc+nLocal)  # define local points to iterate over
    
        # iterate over local points
        for i in pts:
            term = 1
            
            # Lagrange Method
            for j in [k for k in pts if k != i]:
                term = term * ((t - x[j]) / (x[i] - x[j]))
                
            p = p + y[i]*term  # build symbolic interpolating polynomial
        
        fInterp = sym.lambdify(t,p)                                        # turn symbolic poly into real poly
        xx0 = np.arange(x[loc], x[loc+nLocal-1]+nxIncrement, nxIncrement)  # get local interped x values
        yy0 = [fInterp(z) for z in xx0]                                    # get local interped y values
        # write local interped x and y values to array
        # NOTE: duplicate data points exist at locations where piecewise function comes together
        # NOTE: but since they're equal it doesn't matter
        xx = np.concatenate((xx,xx0))
        yy = np.concatenate((yy,yy0))

    return xx,yy
#*******************************************************************


# 4
# Cubic Spline Interpolation method
# NOTE: uses 'natrual spline' condition (set 2nd derivative to zero at endpoints)
#*******************************************************************
#-------------------------------------------------------------------
# x = x data array
# y = y data array
# nPoints = number of data points to use
# nxIncrement = x axis interpolation step size
#-------------------------------------------------------------------
def cubicSpline(x, y, nPoints, nxIncrement):
    # start by solving for double primes to use in cubic spline equation
    # create matrices for Ax=b matrix solver to get values of double primes
    dprimesA = np.zeros((nPoints,nPoints))  # create NxN matrix A for double primes
    dprimesB = np.zeros((nPoints))          # create Nx1 matrix b for double primes
    
    # build matrix A and b to solve for double primes
    for i in range(1, nPoints-1):
        # build matrix A of coefficients
        dprimesA[i][i-1] = x[i]-x[i-1]
        dprimesA[i][ i ] = 2*(x[i]-x[i-1]) + 2*(x[i+1]-x[i])
        dprimesA[i][i+1] = x[i+1]-x[i]
        # build matrix b
        dprimesB[i] = 6*((y[i+1]-y[i])/(x[i+1]-x[i])) - 6*((y[i]-y[i-1])/(x[i]-x[i-1]))
        
    # set endpoint values for spline 'natural condition' 
    dprimesA[0][0] = 1
    dprimesA[nPoints-1][nPoints-1] = 1
    dprimesA[1][0] = 0
    
    # solve for double primes using matrix Ax=b solver to get x
    dprimesX = np.linalg.solve(dprimesA,dprimesB)    
    
    #--------------------------------------------------------------
    
    xx = []  # for storing x interped values
    yy = []  # for storing y interped values
    
    t = sym.Symbol('t')    # define symbolic variable
    n = range(nPoints-1)   # index for interpolation 
    p = 0
    
    for i in n:   
        # Cubic Spline Method
        p0 = y[i]
        p1 = (y[i+1]-y[i])/(x[i+1]-x[i]) - (x[i+1]-x[i])*dprimesX[i+1]/6 - (x[i+1]-x[i])*dprimesX[i]/3
        p2 = dprimesX[i]/2
        p3 = (dprimesX[i+1]-dprimesX[i])/(6*(x[i+1]-x[i]))
        
        # build symbolic cubic function at each point
        p = p0 + p1*(t-x[i]) + p2*(t-x[i])**2 + p3*(t-x[i])**3
        
        fInterp = sym.lambdify(t,p)                             # turn symbolic p into real p
        xx0 = np.arange(x[i], x[i+1]+nxIncrement, nxIncrement)  # get interped values from x(i) to x(i+1)
        yy0 = [fInterp(z) for z in xx0]                         # get interped y values
        # write local interped x and y values to array
        # NOTE: duplicate data points exist at locations where piecewise function comes together
        # NOTE: but since they're equal it doesn't matter
        xx = np.concatenate((xx,xx0))
        yy = np.concatenate((yy,yy0))
                
    return xx,yy
#*******************************************************************  