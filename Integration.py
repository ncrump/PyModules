"""
Created on Thu Jun 13 20:47:04 2013
"""

# Numerical Integration
"""
The following numerical integration methods are below.

-------- Single Integral Routines --------
1. Trapezoid method
2. Simpson method
3. Romberg method
4. Gauss-Legendre method  (interval must be FINITE)
5. Gauss-Laguerre method  (interval is ALWAYS [0,infty])

-------- Double Integral Routines --------
6. Gauss-Legendre 2-D method  (HACK!)
7. Monte Carlo 2D method
"""

import numpy as np
import sympy as sym
from scipy.special import p_roots,l_roots
from math import sin,cos,tan,exp,log,pi


# 1
# Trapezoid method
#*******************************************************************
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of points to evaluate
#-------------------------------------------------------------------  
def trapezoid(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate

    h = (b-a)/n              # interval size (using n rather than n-1 for Python)
    x = np.arange(a,b,h)     # x array defined by interval size
    fx = [f(z) for z in x]   # y array of function values evaluated at x
    
    m = len(x)               # number of points to integrate over
    intSum = 0               # collects sum as integration iterates
    
    # loop over points to sum trapezoids in each interval to get integral
    for i in range(m):
        # check for using correct weight value in summation
        if i == 0 or i == n-1:
            w = 0.5
        else: 
            w = 1
            
        # Trapezoid routine
        intSum = intSum + fx[i]*h*w
        
    return intSum
#*******************************************************************    


# 2
# Simpson method
#*******************************************************************
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of points to evaluate
#------------------------------------------------------------------- 
def simpson(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate
    
    # Simpson requires that n be odd
    # check if n is odd and if not then add 1
    if n % 2 == 0:
        n = n + 1
    
    h = (b-a)/n              # interval size (using n rather than n-1 for Python)
    x = np.arange(a,b,h)     # x array defined by interval size
    fx = [f(z) for z in x]   # y array of function values evaluated at x
    
    m = len(x)               # number of points to integrate over
    intSum = 0               # collects sum as integration iterates
    
    # loop over points to sum parabolas in each interval to get integral
    for i in range(m):
        # check for using correct weight value in summation
        if i == 0 or i == n-1:  # for endpoints use w = 1/3
            w = 1.0/3.0
        elif i % 2 != 0:        # for odd index use w = 4/3
            w = 4.0/3.0
        elif i % 2 == 0:        # for even index use w = 2/3
            w = 2.0/3.0
                        
        # Simpson routine
        intSum = intSum + fx[i]*h*w
        
    return intSum
#*******************************************************************


# 3
# Romberg method
#*******************************************************************
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of interval steps up to 2**(n-1)
#-------------------------------------------------------------------   
def romberg(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate
    
    Rnm = np.zeros((n,n))    # create Romberg matrix
        
    # this part populates the initial column of the Romberg matrix 
    # using modified Trapezoid rule
    # ---------------------------------------------------
    # iterate over interval step size up to n
    for s in range(n):              
        
        h = (b-a)/(2**s)  # interval size
        pSum = 0          # partial sum in integral routine
        
        if s == 0:
            R = 0.5*(b-a)*(f(a)+f(b))  # R(0,0) of Romberg method
            Rnm[s][0] = R              # populate first value of matrix
            
        else: 
            # iterate over interval size up to 2**(n-1)
            for k in range(1, 2**s):
                if k % 2 != 0:                # only odds (k=1,3,5,..2**(n-1))
                    pSum = pSum + f(a + k*h)  # partial sum in integral routine
                  
            R = 0.5*Rnm[s-1][0] + h*pSum      # R(n,0) of Romberg method
            Rnm[s][0] = R                     # populate first column of matrix
    # ---------------------------------------------------
    
    # this part uses previous values in Romberg matrix to get new values
    # each new value is closer to the actual answer
    # --------------------------------------------------- 
    # iterate over rows and columns of Romberg matrix to get new values      
    for col in range(1,n):
        for row in range(col,n):
            romb1 = Rnm[row][col-1]
            romb0 = Rnm[row-1][col-1]
            coeff = 1/((4.0**col)-1)
                        
            # R(n,m) of Romberg method
            Rnm[row][col] = romb1 + coeff*(romb1 - romb0)
    # ---------------------------------------------------
    
    return Rnm[n-1][n-1]
#*******************************************************************   


# 4
# Gauss-Legendre method
#*******************************************************************
# NOTE: the limits of integration must be FINITE
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of points to generate and evaluate
#-------------------------------------------------------------------
def gaussLegendre(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate
    
    # get abscissas (xi), and weights (wi) from p_roots(n) function below
    # scipy.special.p_roots(n) returns abscissas & weights of Legendre polynomials
    abscissa, weight = p_roots(n)  
    intSum = 0
    
    # Gauss-Legendre method is valid only over interval [-1,1]            
    # map input [a,b] to interval [-1,1]
    for i in range(n):
        # evaluate integral with interval transformed by x--> 0.5(b-a)x + 0.5(b+a)
        intSum = intSum + weight[i]*f(0.5*(b-a)*abscissa[i] + 0.5*(b+a))  

    return 0.5*(b-a)*intSum
#*******************************************************************


# 5
# Gauss-Laguerre method
#*******************************************************************
# NOTE: the limits of integration are assumed [0,infty]
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)*exp(-x)')
# n = number of points to generate and evaluate
#-------------------------------------------------------------------
def gaussLaguerre(g, n):
    f = lambda x: eval(g)    # define function to evaluate
    
    # get abscissas (xi), and weights (wi) from l_roots(n) function below
    # scipy.special.l_roots(n) returns abscissas & weights of Laguerre polynomials
    abscissa, weight = l_roots(n)  
    intSum = 0
    
    # Gauss-Laguerre method is valid  over interval [0,infty]
    for i in range(n):
        # evaluate integral
        intSum = intSum + weight[i]*f(abscissa[i])*exp(abscissa[i])
    
    return intSum
#*******************************************************************


# 6
# Gauss-Legendre 2-D method
#*******************************************************************
# NOTE: this one is a HACK !!!!!
# NOTE: have to enter functions as 'sym.function(x,y)'
#-------------------------------------------------------------------
# g = 'function of x,y' entered in quotes (Ex. 'sym.sin(x*y)')
# x1 = integration region x start point (entered as decimal)
# x2 = integration region x end point (entered as decimal)
# y1 = integration region y start point (entered as decimal)
# y2 = integration region y end point (entered as decimal)
# n = number of points to generate and evaluate
#-------------------------------------------------------------------
def gaussLegendre2D(g, x1, x2, y1, y2, n):
    # define symbolic variable for function of y
    y = sym.Symbol('y')
    
    # integrate function over x interval first
    # ------------------------------------------
    # get abscissas (xi), and weights (wi) from p_roots(n) function below
    abscissa, weight = p_roots(n)  
    intSum = 0
    
    # map x interval to [-1,1]
    for i in range(n):
        # transform x interval
        x = 0.5*(x2-x1)*abscissa[i] + 0.5*(x2+x1)
        fxy = eval(g)
        
        # evaluate integral with x interval transformed
        intSum = intSum + weight[i]*fxy  

    # integral functin in y now
    fy = 0.5*(x2-x1)*intSum
    Fy = sym.lambdify(y,fy)  
    
    
    # integrate function over y interval second
    # ------------------------------------------
    intSum = 0
    
    # map y interval to [-1,1]
    for i in range(n):
        # transform y interval
        yi = 0.5*(y2-y1)*abscissa[i] + 0.5*(y2+y1)
        
        # evaluate integral with y interval transformed
        intSum = intSum + weight[i]*Fy(yi)    

    return 0.5*(y2-y1)*intSum    
#*******************************************************************
    
    
# 7
# Monte Carlo 2-D method
#*******************************************************************
# NOTE: uses built-in random number generator numpy.random.uniform
#-------------------------------------------------------------------
# g = 'function of x,y' entered in quotes (Ex. 'sin(x*y)')
# x1 = integration region x start point (entered as decimal)
# x2 = integration region x end point (entered as decimal)
# y1 = integration region y start point (entered as decimal)
# y2 = integration region y end point (entered as decimal)
# N = total number of random points to generate and evaluate
#-------------------------------------------------------------------
def monteCarlo2D(g, x1, x2, y1, y2, N):    
    f = lambda x,y: eval(g)    # define function to evaluate
    n = int(N**0.5)            # number of random points in each interval
    
    # generate random numbers over x,y integration interval
    Rx = np.random.uniform(x1,x2,n)
    Ry = np.random.uniform(y1,y2,n)
    
    # compute weights for average of function
    Wx = (x2-x1)/n
    Wy = (y2-y1)/n
    
    intSum = 0 
    
    # loop over x,y region of random numbers to evaluate function
    for i in Rx:
        for j in Ry:
            intSum = intSum + f(i,j)

    # return average value of function
    return Wx*Wy*intSum
#*******************************************************************