"""
Created on Thu May 30 21:01:03 2013
"""

# Root Finding Methods
"""
The following root finding methods are below.
1. Bisection method
2. Hybrid Bisection/Newton-Raphson method
3. Hybrid Bisection/Secant method
4. Hybrid Bisection/Muller-Brent method
"""

from math import isnan,isinf,sqrt


#--------------------------------------------------------------------------------------
# called as 'method'('function(x)', xInitial, xFinal, Tolerance, MaxIterations)
# function(x) = input function of x (**as a lambda function**)
# xInitial = initial x bracket value on one side of the root (entered as decimal)
# xFinal = initial x bracket value on the other side of the root (entered as decimal)
# Tolerance = the degree of accuracy in which to compute the root (Ex: 10e-5)
# MaxIterations = max number of iterations to compute if Tolerance has not been met
#--------------------------------------------------------------------------------------


# 1
# Bisection method
#***********************************************************************
def bisection(f, xI, xF, Tol, nMax):
        # ---------------------------------------------------------------------
        # check input bracket condition to ensure root is between bracket points
        # check to make sure function is defined at bracket points too
        # this checks for examples like ln(0) = -inf or ln(0)*ln(4) = nan
        notnum1 = isnan(f(xI)*f(xF))
        notnum2 = isinf(f(xI)*f(xF))
        if f(xI)*f(xF) >= 0:
            print '\nEITHER NO ROOTS OR MULTIPLE ROOTS EXIST BETWEEN BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        if notnum1 == True or notnum2 == True:
            print '\nFUNCTION IS UNDEFINED AT BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        # ---------------------------------------------------------------------

        # initialize variables
        error = 1
        n = 1
        xiMid = 0  # initial midpoint value to store the n-1 value

        # loop until error is less than input tolerance
        while error > Tol:
            xMid = 0.5*(xI+xF)

            # set up main Bisection method:
            # make bracket interval smaller each iteration until root is found
            # check conditions and update bracket points
            if f(xI)*f(xMid) > 0:
                xI = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint

            elif f(xI)*f(xMid) < 0:
                xF = xMid
                error = abs(xMid - xiMid)  # calculate approx error
                n = n + 1
                xiMid = xMid               # store the n-1 midpoint

        print 'Iterations = ', n-1
        print 'Approximate Error =', error

        return xMid
#***********************************************************************


# 2
# hybrid Bisection/Newton-Raphson method
#***********************************************************************
# NOTE: requires input step size h to compute derivative of f
# NOTE: uses central difference method to compute derivative
def newtonRaphson(f, xI, xF, h, Tol, nMax):
        # ---------------------------------------------------------------------
        # check input bracket condition to ensure root is between bracket points
        # check to make sure function is defined at bracket points too
        # this checks for examples like ln(0) = -inf or ln(0)*ln(4) = nan
        notnum1 = isnan(f(xI)*f(xF))
        notnum2 = isinf(f(xI)*f(xF))
        if f(xI)*f(xF) >= 0:
            print '\nEITHER NO ROOTS OR MULTIPLE ROOTS EXIST BETWEEN BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        if notnum1 == True or notnum2 == True:
            print '\nFUNCTION IS UNDEFINED AT BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        # ---------------------------------------------------------------------

        # initialize variables
        error = 1
        n = 1

        # check condition for starting x value
        if f(xI) < f(xF): xn0 = xI
        else: xn0 = xF

        # loop until error is less than input tolerance
        while error > Tol:

            # set up main Newton-Raphson method:
            # use derivative of function as tangent line to get new x point
            # calcuate new x value from f(x) and df(x)

            df = (f(xn0+h) - f(xn0)) / h   # central difference method
            xn1 = xn0 - (f(xn0)/df)

            # if new x point is outside bracket interval then use midpoint instead
            if xn1 < xI or xn1 > xF:
                xn1 = 0.5*(xI+xF)

                # check conditions and update bracket points
                if f(xI)*f(xn1) > 0:
                    xI = xn1

                elif f(xI)*f(xn1) < 0:
                    xF = xn1

            error = abs(xn1 - xn0)  # calculate approx error
            n = n + 1
            xn0 = xn1  # store the n-1 value for x

        print 'Iterations = ', n-1
        print 'Approximate Error =', error

        return xn1
#***********************************************************************


# 3
# hybrid Bisection/Secant method
#***********************************************************************
def secant(f, xI, xF, Tol, nMax):
        # ---------------------------------------------------------------------
        # check input bracket condition to ensure root is between bracket points
        # check to make sure function is defined at bracket points too
        # this checks for examples like ln(0) = -inf or ln(0)*ln(4) = nan
        notnum1 = isnan(f(xI)*f(xF))
        notnum2 = isinf(f(xI)*f(xF))
        if f(xI)*f(xF) >= 0:
            print '\nEITHER NO ROOTS OR MULTIPLE ROOTS EXIST BETWEEN BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        if notnum1 == True or notnum2 == True:
            print '\nFUNCTION IS UNDEFINED AT BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        # ---------------------------------------------------------------------

        # initialize variables
        error = 1
        n = 1

       # initialize starting values
        xn0 = xI
        xn1 = xF

        # loop until error is less than input tolerance
        while error > Tol:

            # set up main Secant method:
            # use secant line as linear approx to function to get new x point
            # calcuate new x value from secant line
            xn2 = xn1 - ((f(xn1)*(xn1 - xn0)) / (f(xn1) - f(xn0)))

            # if new x point is outside bracket  then use midpoint instead
            if xn2 < xI or xn2 > xF:
                xn2 = 0.5*(xI+xF)

                # check conditions and update bracket points
                if f(xI)*f(xn2) > 0:
                    xI = xn2

                elif f(xI)*f(xn2) < 0:
                    xF = xn2

            error = abs(xn2 - xn1)  # calculate approx error
            n = n + 1
            xn0 = xn1  # store the n-1 value for x
            xn1 = xn2  # store the nth value for x

        print 'Iterations = ', n-1
        print 'Approximate Error =', error

        return xn2
#***********************************************************************


# 4
# hybrid Bisection/Muller-Brent method
#***********************************************************************
def mullerBrent(f, xI, xF, Tol, nMax):
        # ---------------------------------------------------------------------
        # check input bracket condition to ensure root is between bracket points
        # check to make sure function is defined at bracket points too
        # this checks for examples like ln(0) = -inf or ln(0)*ln(4) = nan
        notnum1 = isnan(f(xI)*f(xF))
        notnum2 = isinf(f(xI)*f(xF))
        if f(xI)*f(xF) >= 0:
            print '\nEITHER NO ROOTS OR MULTIPLE ROOTS EXIST BETWEEN BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        if notnum1 == True or notnum2 == True:
            print '\nFUNCTION IS UNDEFINED AT BRACKET POINTS.'
            print 'ADJUST BRACKET VALUES.\n'
            return
        # ---------------------------------------------------------------------

        # initialize variables
        error = 1
        n = 1

        # initialize starting values
        xn0 = xI
        xn1 = xF
        # use midpoint for 3rd point for use in quadratic approx to function
        xn2 = 0.5*(xn0+xn1)

        # loop until error is less than input tolerance
        while error > Tol:

            # set up main Muller-Brent method:
            # use 3 points as quadratic approx to function to get new x point
            # first calculate a, b, c from function values at 3 known points
            c = f(xn2)
            bNumer = (((xn0-xn2)**2)*(f(xn1)-f(xn2))) - (((xn1-xn2)**2)*(f(xn0)-f(xn2)))
            aNumer = ((xn1-xn2)*(f(xn0)-f(xn2))) - ((xn0-xn2)*(f(xn1)-f(xn2)))
            Denom = (xn0-xn1)*(xn0-xn2)*(xn1-xn2)
            b = bNumer / Denom
            a = aNumer / Denom

            # then calculate new x value from quadratic eq
            # check b and use alternate quadratic eq to avoid subtractive cancellations
            if b >= 0:
                xn3 = xn2 - (2*c) / (b + sqrt(b**2 - 4*a*c))
            else:
                xn3 = xn2 + (2*c) / (-b + sqrt(b**2 - 4*a*c))

            # if new x point is outside bracket interval then use midpoint instead
            if xn3 < xI or xn3 > xF:
                xn3 = 0.5*(xI+xF)

                # check conditions and update bracket points
                if f(xI)*f(xn3) > 0:
                    xI = xn3

                elif f(xI)*f(xn3) < 0:
                    xF = xn3

            error = abs(xn3 - xn2)  # calculate approx error
            n = n + 1
            xn0 = xn1  # store the n-1 value for x
            xn1 = xn2  # store the nth value for x
            xn2 = xn3  # store the n+1 value for x

        print 'Iterations = ', n-1
        print 'Approximate Error =', error

        return xn3
#***********************************************************************