"""
Created on Sat Sep 21 23:44:35 2013
"""

# First Derivative Approximations
"""
The following first derivative methods are below.
1. Backward difference method
2. Forward difference method
3. Central difference method
"""


# 1
# backward difference method for first derivative approximation
#*************************************************************
#-----------------------------------------------------------
# f = function of x (**enter as lambda function**)
# x = x point to evaluate derivative
# h = step size
#----------------------------------------------------------- 
def backwardDiff(f,x,h):
    df = (f(x) - f(x-h)) / h
    return df
#*************************************************************


# 2
# forward difference method for first derivative approximation
#*************************************************************
#-----------------------------------------------------------
# f = function of x (**enter as lambda function**)
# x = x point to evaluate derivative
# h = step size
#----------------------------------------------------------- 
def forwardDiff(f,x,h):
    df = (f(x+h) - f(x)) / h
    return df
#*************************************************************


# 3
# central difference method for first derivative approximation
#*************************************************************
#-----------------------------------------------------------
# f = function of x (**enter as lambda function**)
# x = x point to evaluate derivative
# h = step size
#----------------------------------------------------------- 
def centralDiff(f,x,h):
    df = (f(x+h) - f(x-h)) / (2*h)
    return df
#*************************************************************