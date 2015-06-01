"""
Created on Wed Jul 10 20:41:20 2013
"""

# Ordinary Differential Equations (ODEs)
"""
The following numerical ODE methods are below.

-------- First Order ODEs --------
1. Simple Euler Method
2. Modified Euler Method
3. Improved Euler Method
4. 4th-Order Runge-Kutta Method (RK4)
5. 5th-Order Runge-Kutta-Fehlberg (RKF45) (adaptive step)

-------- Higher Order ODEs --------
6. 5th-Order Runge-Kutta-Fehlberg (RKF45) (adaptive step)
"""

import math as m
import numpy as np


# 1
# Simple Euler Method
#*******************************************************************
#-------------------------------------------------------------------
# f = ODE function (** callable function **)
# x0 = initial condition
# y0 = initial condition
# xf = upper bound of function domain
# h = step size
#------------------------------------------------------------------- 
def simpleEuler(f,x0,y0,xf,h):
    # this uses the derivative at the start of the interval across each step  
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1

    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Simple Euler Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        y[i+1] = y[i] + h*f(x[i],y[i])    # Euler approx to function
        
    # returns solution values of ODE
    return x,y
#*******************************************************************


# 2
# Modified Euler Method
#*******************************************************************
#-------------------------------------------------------------------
# f = ODE function (** callable function **)
# x0 = initial condition
# y0 = initial condition
# xf = upper bound of function domain
# h = step size
#------------------------------------------------------------------- 
def modifiedEuler(f,x0,y0,xf,h):
    # this uses the derivative at the midpoint of the interval across each step
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1
    
    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Modified Euler Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        xmid = x[i] + 0.5*h               # get x midpoint across interval
        ymid = y[i] + 0.5*h*f(x[i],y[i])  # get y at x midpoint
        
        y[i+1] = y[i] + h*f(xmid,ymid)    # Euler approx to function at midpoint
        
    # returns solution values of ODE
    return x,y
#*******************************************************************


# 3
# Improved Euler Method
#*******************************************************************
#-------------------------------------------------------------------
# f = ODE function (** callable function **)
# x0 = initial condition
# y0 = initial condition
# xf = upper bound of function domain
# h = step size
#------------------------------------------------------------------- 
def improvedEuler(f,x0,y0,xf,h):
    # this uses the average derivative value over the interval across each step
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1

    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Improved Euler Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        dy1 = f(x[i],y[i])               # get derivative at interval start
        dy2 = f(x[i] + h, y[i] + h*dy1)  # get derivative at interval end
        dyAvg = 0.5*(dy1 + dy2)          # take average of derivatives
                
        y[i+1] = y[i] + h*(dyAvg)        # Euler approx to function using dyAvg
        
    # returns solution values of ODE
    return x,y
#*******************************************************************


# 4
# 4th-Order Runge-Kutta Method (RK4)
#*******************************************************************
#-------------------------------------------------------------------
# f = ODE function (** callable function **)
# x0 = initial condition
# y0 = initial condition
# xf = upper bound of function domain
# h = step size
#------------------------------------------------------------------- 
def RK4(f,x0,y0,xf,h):
    # this method uses the 4th-order Runge-Kutta Method (RK4)
    
    # get number of points over interval for step size h
    pts = m.floor(abs(xf-x0)/h) + 1

    # create x array of points and y array for storing solution values
    x = np.linspace(x0,xf,pts)
    y = np.zeros(pts)
    y[0] = y0    
    
    # do Runge-Kutta Method to solve ODE over each step h
    for i in np.arange(0,pts-1):
        y1 = f(x[i],y[i])
        y2 = f(x[i] + 0.5*h, y[i] + 0.5*h*y1)
        y3 = f(x[i] + 0.5*h, y[i] + 0.5*h*y2)
        y4 = f(x[i] + h, y[i] + h*y3)
        
        y[i+1] = y[i] + (h/6.0)*(y1 + 2*y2 + 2*y3 + y4)
        
    # returns solution values of ODE
    return x,y
#*******************************************************************


# 5
# 5th-Order Runge-Kutta-Fehlberg Method (RKF45) (adaptive step)
#*******************************************************************
#-------------------------------------------------------------------
# f = ODE function (** callable function **)
# x0 = initial condition
# y0 = initial condition
# xf = upper bound of function domain
# h = initial step size
# eps = relative error tolerance for computing adaptive step h
#------------------------------------------------------------------- 
def RKF45(f,x0,y0,xf,h,eps):
    # this method uses the 5th-order Runge-Kutta-Fehlberg Method (RKF45)

    # create x and y arrays for storing solution values
    x = [x0]
    y = [y0]
    
    # set initial values
    xi = x0
    yi = y0
        
    # RKF45 coefficients
    # for intermediate function f1
    a1, b1, = 0.25, 0.25
    # for intermediate function f2
    a2, b2, c2 = 3.0/8.0, 3.0/32.0, 9.0/32.0
    # for intermediate function f3
    a3, b3, c3, d3 = 12.0/13.0, 1932.0/2197.0, 7200.0/2197.0, 7296.0/2197.0
    # for intermediate function f4
    a4, b4, c4, d4, e4 = 1, 439.0/216.0, 8, 3680.0/513.0, 845.0/4104.0
    # for intermediate function f5
    a5, b5, c5, d5, e5, f5 = 0.5, 8.0/27.0, 2, 3544.0/2565.0, 1859.0/4104.0, 11.0/40.0
    
    # for 5th-order solution
    Y1, Y2, Y3, Y4, Y5 = 16.0/135.0, 6656.0/12825.0, 28561.0/56430.0, 9.0/50.0, 2.0/55.0
    # for error estimate 
    E1, E2, E3, E4, E5 = 1.0/360.0, 128.0/4275.0, 2197.0/75240.0, 1.0/50.0, 2.0/55.0
    
    # RKF45 method to solve ODE over each adaptive step h
    while xi <= xf:
        
        # compute initial intermediate function f0
        f0 = f(xi,yi)
        
        # compute dependent variable arrays and intermediate function arrays f1-f5
        xI = xi + a1*h
        yI = yi + b1*h*f0
        f1 = f(xI,yI)
        
        xI = xi + a2*h
        yI = yi + b2*h*f0 + c2*h*f1
        f2 = f(xI,yI)
        
        xI = xi + a3*h
        yI = yi + b3*h*f0 - c3*h*f1 + d3*h*f2
        f3 = f(xI,yI)
        
        xI = xi + a4*h
        yI = yi + b4*h*f0 - c4*h*f1 + d4*h*f2 - e4*h*f3
        f4 = f(xI,yI)
        
        xI = xi + a5*h
        yI = yi - b5*h*f0 + c5*h*f1 - d5*h*f2 + e5*h*f3 - f5*h*f4
        f5 = f(xI,yI)
        
        # compute solution function and error for adaptive step
        yy = yi + h*(Y1*f0 + Y2*f2 + Y3*f3 - Y4*f4 + Y5*f5)        
        error = h*(E1*f0 - E2*f2 - E3*f3 + E4*f4 + E5*f5)
        RelError = abs(error/yy)
    
        # set adaptive step based on specified input tolerance 
        hNew = h*(eps/RelError)**0.2
        if hNew >= h: 
            xi = xi + h  
            yi = yy
            x.append(xi)        
            y.append(yy)
        
        # try smaller step this time
        h = 0.9*hNew
              
    # returns solution values of ODE
    return x,y
#*******************************************************************
    
    
# 6
# 5th-Order Runge-Kutta-Fehlberg Method (RKF45) (adaptive step)
#*******************************************************************
# NOTE: initial conditions must be consistent with callable function
# NOTE: callable function must return a numpy array
#-------------------------------------------------------------------
# f = ODE function (** callable function - must return numpy array **)
# IC = array of initial solution conditions (** order matters **)
# t0 = lower bound of independent variable
# tf = upper bound of independent variable
# h = initial step size
# eps = relative error tolerance for computing adaptive step h
# stp = option to use adaptive step or h (0 = h, 1 = adaptive)
#-------------------------------------------------------------------
def RKF45HO(f,IC,t0,tf,h,eps,stp=1):
    # this method uses the 5th-order Runge-Kutta-Fehlberg Method (RKF45)  
    
    # create arrays for storing solution values
    N = len(IC)
    t = np.array([t0])
    s = np.reshape(IC,(N,1))    
    
    # set initial values
    ti = t0
    Si = np.array(IC)
            
    # RKF45 method to solve ODE over each adaptive step h
    # intermediate function f1
    a1, b1, = 0.25, 0.25
    # intermediate function f2
    a2, b2, c2 = 3.0/8.0, 3.0/32.0, 9.0/32.0
    # intermediate function f3
    a3, b3, c3, d3 = 12.0/13.0, 1932.0/2197.0, 7200.0/2197.0, 7296.0/2197.0
    # intermediate function f4
    a4, b4, c4, d4, e4 = 1, 439.0/216.0, 8, 3680.0/513.0, 845.0/4104.0
    # intermediate function f5
    a5, b5, c5, d5, e5, f5 = 0.5, 8.0/27.0, 2, 3544.0/2565.0, 1859.0/4104.0, 11.0/40.0

    # for 5th-order solution
    Y1, Y2, Y3, Y4, Y5 = 16.0/135.0, 6656.0/12825.0, 28561.0/56430.0, 9.0/50.0, 2.0/55.0
    # for error estimate
    E1, E2, E3, E4, E5 = 1.0/360.0, 128.0/4275.0, 2197.0/75240.0, 1.0/50.0, 2.0/55.0
    
    # RKF45 method to solve ODE over each adaptive step h
    # NOTE: intermediate functions are vectors represented as arrays
    # NOTE: solution and error functions are also vector functions
    while ti <= tf:
                
        # compute initial intermediate function f0
        f0 = f(ti,Si)
        
        # compute dependent variable arrays and intermediate function arrays f1-f5
        SI = Si + b1*h*f0
        f1 = f(ti,SI)
        
        SI = Si + b2*h*f0 + c2*h*f1
        f2 = f(ti,SI)
        
        SI = Si + b3*h*f0 - c3*h*f1 + d3*h*f2
        f3 = f(ti,SI)
        
        SI = Si + b4*h*f0 - c4*h*f1 + d4*h*f2 - e4*h*f3
        f4 = f(ti,SI)
        
        SI = Si - b5*h*f0 + c5*h*f1 - d5*h*f2 + e5*h*f3 - f5*h*f4
        f5 = f(ti,SI)
        
        # compute solution function and error for adaptive step
        F = Si + h*(Y1*f0 + Y2*f2 + Y3*f3 - Y4*f4 + Y5*f5)
        error = h*(E1*f0 - E2*f2 - E3*f3 + E4*f4 + E5*f5)
        RMS = (sum(error**2) / len(F))**0.5
        
        # option to use adaptive step size
        if stp == 1:        
            # set adaptive step size based on RMS error vs specified max acceptable error
            hNew = h*(eps/RMS)**0.2
            
            # if error is small enough, accept solution and continue
            # otherwise iterate again with adjusted smaller step size
            if hNew >= h: 
                ti = ti + h                 # update independent variable
                Si = F                      # update dependent variables
                rSi = np.reshape(Si,(N,1))  # reshape Si vector for storing
                t = np.append(t,ti)         # store t values to array
                s = np.append(s,(rSi),1)    # store s values to array
    
            # try smaller step size this time
            h = 0.9*hNew 

        # option to use step h instead of adaptive step
        else: 
            ti = ti + h                     # update independent variable
            Si = F                          # update dependent variables
            rSi = np.reshape(Si,(N,1))      # reshape Si vector for storing
            t = np.append(t,ti)             # store t values to array
            s = np.append(s,(rSi),1)        # store s values to array              
              
    # returns array of solution vectors s for independent variable t
    return t,s
#*******************************************************************