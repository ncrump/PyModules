"""
Created on Thu Jun 13 20:47:04 2013
"""

# Numerical Integration
"""
The following numerical integration methods are below.

-------- Single Integral Routines --------
  Newton-Cotes Methods
    1. Trapezoid method
    2. Simpson method
  Gauss-Quadrature Methods
    3. Romberg method
    4. Gauss-Legendre method  (interval must be FINITE)
    5. Gauss-Laguerre method  (interval is ALWAYS [0,infty])

-------- Double Integral Routines --------
    6. Monte Carlo 2D method
"""

import numpy as np
from scipy.special import p_roots,l_roots


# 1
# Trapezoid method
#*******************************************************************
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of intervals to evaluate
#-------------------------------------------------------------------
def trapezoid(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate

    h  = (b-a)/float(n)      # interval size
    x  = np.arange(a,b+h,h)  # x array defined by interval size
    fx = f(x)                # y array of function values evaluated at x

    # loop over points to sum trapezoids in each interval to get integral
    intSum = fx[0]+fx[n]
    for i in range(1,n):
        # get weight for summation
        intSum += 2*fx[i]
    # get final integral
    trap = (h/2.0)*intSum

    return trap
#*******************************************************************


# 2
# Simpson method
#*******************************************************************
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = number of intervals to evaluate
#-------------------------------------------------------------------
def simpson(g, a, b, n):
    f = lambda x: eval(g)    # define function to evaluate

    # Simpson 1/3 requires that n be even
    # if n not even then add 1
    if n % 2 != 0:
        n += 1

    h  = (b-a)/float(n)      # interval size
    x  = np.arange(a,b+1,h)  # x array defined by interval size
    fx = f(x)                # y array of function values evaluated at x

    # loop over points to sum parabolas in each interval to get integral
    intSum = fx[0]+fx[n]
    for i in range(1,n):
        # get weight for summation
        if i%2 == 0: w = 2
        else:        w = 4
        intSum += w*fx[i]
    # get final integral
    simp = (h/3.0)*intSum

    return simp
#*******************************************************************


# 3
# Romberg method
#*******************************************************************
#-------------------------------------------------------------------
# g = 'function of x' entered in quotes (Ex. 'sin(x)')
# a = integration region x start point (entered as decimal)
# b = integration region x end point (entered as decimal)
# n = order of extrapolation (Ex. n=10)
#-------------------------------------------------------------------
def romberg(g, a, b, n=10):
    f = lambda x: eval(g)    # define function to evaluate

    # initialize Romberg matrix
    R      = np.zeros((n,n))
    R[0,0] = 0.5*(b-a)*(f(a)+f(b))

    # this part populates the initial column of the Romberg matrix
    # using modified Trapezoid rule
    # ---------------------------------------------------
    # iterate over interval step size up to n
    for i in range(1,n):
        pts  = 2**i
        h    = (b-a)/float(pts)  # interval size
        intSum = 0               # partial sum in integral routine
        # iterate over interval size up to 2**(n-1)
        for j in range(1, pts):
            if j % 2 != 0:            # only odds (k=1,3,5,..2**(n-1))
                intSum += f(a + j*h)  # partial sum in integral routine

        R[i,0] = 0.5*R[i-1,0] + h*intSum  # R(n,0) of Romberg method
    # ---------------------------------------------------

    # this part uses previous values in Romberg matrix to get new values
    # each new value is closer to the actual answer
    # ---------------------------------------------------
    # iterate over rows and columns of Romberg matrix to get new values
    for col in range(1,n):
        for row in range(col,n):
            romb1 = R[row,col-1]
            romb0 = R[row-1,col-1]
            coeff = 1.0/(4.0**col-1)

            # R(n,m) of Romberg method
            R[row,col] = romb1 + coeff*(romb1 - romb0)
    # ---------------------------------------------------

    return R[-1,-1]
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
        intSum = intSum + weight[i]*f(abscissa[i])*np.exp(abscissa[i])

    return intSum
#*******************************************************************


# 6
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