"""
Created on Thu Jun 20 06:14:19 2013
"""

# Least Squares Fitting
"""
The following least square methods are below.
1. Polynomial fit for coefficients with LINEAR dependence
2. Gauss-Newton fit for NON-LINEAR coefficients with 2 parameters
3. Gauss-Newton fit for NON-LINEAR coefficients with 3 parameters
4. Gauss-Newton fit for NON-LINEAR coefficients with 4 parameters
"""

import numpy as np


# 1
# Polynomial least squares fitting method for linear data
#*******************************************************************
#-------------------------------------------------------------------
# x = x data set
# y = y data set
# n = degree of polynomial fit
#------------------------------------------------------------------- 
def polyFit(x,y,n):
    # fit model is polynomial: y = c0 + c1(x) + c2(x**2) +...+ cn(x**n)
    
    pts = len(x)  # number of poins in data set
    polyval = []  # array to store poly fit values
    
    # if degree of polynomial is greater than number of points print error
    if n+1 >= pts:
        print 'INVALID ENTRY'
        print 'POLYNOMIAL DEGREE CANNOT BE GREATER THAN NUMBER OF DATA POINTS'
    
    # otherwise compute coefficients for polynomial fit
    # solves for coefficients using matrix solve of 'normal equations'
    # turns matrix equation Ax=b into 'normal equation' (A^T)Ax=(A^T)b
    # where b = y values, A = x values^power
    # this needed because system of equations is overdetermined (#equations > #unknowns)
    else:
        matrixA = np.zeros((pts,n+1))             # initialize matrix A
        matrixB = np.zeros((pts,1))               # initialize matrix b
        
        # loop to populate arrays
        for i in range(pts):
            matrixB[i][0] = y[i]                  # b gets y values
            
            for j in range(n+1):
                matrixA[i][j] = (x[i])**j         # A gets x values^j
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)                # A transpose (A^T)
        AtA = np.dot(At,matrixA)                  # matrix multiply (A^T)A
        AtB = np.dot(At,matrixB)                  # matrix multiply (A^T)b
        
        coeff = np.linalg.solve(AtA,AtB)          # solve normal equation (A^T)Ax=(A^T)b for x
        
        
        # we now have the coefficients for the polynomial fit
        # next plug x values back into fit model and compute new y fit values
        varsum = 0     
        # loop to get y fit values and compute variance on fit                           
        for i in range(pts):
            yval = 0
            for j in range(n+1):
                yval = yval + coeff[j]*(x[i])**j  # gets new y value for each x using fit coeffs
                            
            polyval.append(yval)                  # collect new y values
            
            # compute variance from erorr (error=actualY-fitY)
            error = y[i]-yval
            varsum = varsum + error**2
        
        # using variance as measure for 'goodness of fit'
        # computed as variance = (sum of squares of errors) / (nDataPts-PolyDegree-1)
        variance = round(varsum/(pts-n-1), 4)
        
        
        print 'coeffs=',coeff        
        print 'variance=',variance        
        return x,polyval 
#*******************************************************************


# 2
# Gauss-Newton fitting method for nonlinear data and 2 fit parameters
#*******************************************************************
#-------------------------------------------------------------------
# x = x data set
# y = y data set
# f = fit model (**as a lambda function**)
# dfa = derivative of f wrt a (**as a lambda function**)
# dfb = derivative of f wrt b (**as a lambda function**)
# a = guess at first fit parameter
# b = guess at second fit parameter
# n = number of iterations to converge on fit parameters
#------------------------------------------------------------------- 
def gaussNewton2(x,y,f,dfa,dfb,a,b,n):
    # fit model is user defined for 2 parameters
    # solves for coefficients using a matrix solve of modified 'normal equations'
    # turns matrix equation Ax=b into 'normal equation' (A^T)Ax=(A^T)b
    # where b = errors to minimize (actualY-fitY), A = numeric derivative of b
                
    pts = len(x)  # number of points in dataset                                  
    fitVals = []  # stores fit values
    
    # loop over specified number of iterations to update coefficients a,b
    for k in range(n):
        matrixA = np.zeros((pts,2))
        matrixB = np.zeros((pts,1))
   
       # loop to populate arrays
        for i in range(pts):            
            # populate matrix elements 
            matrixB[i][0] = y[i]-f(a,b,x[i])  # matrix of residual errors                  
            matrixA[i][0] = dfa(a,b,x[i])     # 1st column gets derivative wrt a
            matrixA[i][1] = dfb(a,b,x[i])     # 1st column gets derivative wrt b
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)            # A transpose (A^T)
        AtA = np.dot(At,matrixA)              # matrix multiply (A^T)A
        AtB = np.dot(At,matrixB)              # matrix multiply (A^T)b
        
        coeff = np.linalg.solve(AtA,AtB)      # solve normal equation (A^T)Ax=(A^T)b
                      
        # get coefficients
        a = a + coeff[0]                           
        b = b + coeff[1]        
    
    
    # we now have the coefficients for the fit
    # next plug x values back into fit model and compute new y fit values
    varsum = 0     
    
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        val = f(a,b,x[i])
        
        fitVals.append(val)
                                
        # compute variance from erorr (error=actualY-fitY)
        error = y[i] - val
        
        varsum = varsum + error**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    variance = varsum/(pts-2-1)
    
    
    print 'coeffs=',a,b
    print 'variance=',variance
    return x,fitVals
#******************************************************************* 


# 3
# Gauss-Newton fitting method for nonlinear data and 3 fit parameters
#*******************************************************************
#-------------------------------------------------------------------
# x = x data set
# y = y data set
# f = fit model (**as a lambda function**)
# dfa = derivative of f wrt a (**as a lambda function**)
# dfb = derivative of f wrt b (**as a lambda function**)
# dfc = derivative of f wrt c (**as a lambda function**)
# a = guess at first fit parameter
# b = guess at second fit parameter
# c - guess at third fit parameter
# n = number of iterations to converge on fit parameters
#------------------------------------------------------------------- 
def gaussNewton3(x,y,f,dfa,dfb,dfc,a,b,c,n):
    # fit model is user defined for 3 parameters
    # solves for coefficients using a matrix solve of modified 'normal equations'
    # turns matrix equation Ax=b into 'normal equation' (A^T)Ax=(A^T)b
    # where b = errors to minimize (actualY-fitY), A = numeric derivative of b
                
    pts = len(x)  # number of points in dataset                                  
    fitVals = []  # stores fit values
    
    # loop over specified number of iterations to update coefficients a,b,c
    for k in range(n):
        matrixA = np.zeros((pts,3))
        matrixB = np.zeros((pts,1))
   
       # loop to populate arrays
        for i in range(pts):            
            # populate matrix elements 
            matrixB[i][0] = y[i]-f(a,b,c,x[i])  # matrix of residual errors                  
            matrixA[i][0] = dfa(a,b,c,x[i])     # 1st column gets derivative wrt a
            matrixA[i][1] = dfb(a,b,c,x[i])     # 1st column gets derivative wrt b
            matrixA[i][2] = dfc(a,b,c,x[i])     # 1st column gets derivative wrt c
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)              # A transpose (A^T)
        AtA = np.dot(At,matrixA)                # matrix multiply (A^T)A
        AtB = np.dot(At,matrixB)                # matrix multiply (A^T)b
        
        coeff = np.linalg.solve(AtA,AtB)        # solve normal equation (A^T)Ax=(A^T)b
                      
        # get coefficients
        a = a + coeff[0]                           
        b = b + coeff[1]
        c = c + coeff[2]
        
    
    # we now have the coefficients for the fit
    # next plug x values back into fit model and compute new y fit values
    varsum = 0     
    
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        val = f(a,b,c,x[i])
        
        fitVals.append(val)
                                
        # compute variance from erorr (error=actualY-fitY)
        error = y[i] - val
        
        varsum = varsum + error**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    variance = varsum/(pts-3-1)
    
    
    print 'coeffs=',a,b,c
    print 'variance=',variance
    return x,fitVals
#******************************************************************* 


# 4
# Gauss-Newton fitting method for nonlinear data and 4 fit parameters
#*******************************************************************
#-------------------------------------------------------------------
# x = x data set
# y = y data set
# f = fit model (**as a lambda function**)
# dfa = derivative of f wrt a (**as a lambda function**)
# dfb = derivative of f wrt b (**as a lambda function**)
# dfc = derivative of f wrt c (**as a lambda function**)
# dfd = derivative of f wrt d (**as a lambda function**)
# a = guess at first fit parameter
# b = guess at second fit parameter
# c - guess at third fit parameter
# d - guess at fourth fit parameter
# n = number of iterations to converge on fit parameters
#------------------------------------------------------------------- 
def gaussNewton4(x,y,f,dfa,dfb,dfc,dfd,a,b,c,d,n):
    # fit model is user defined for 4 parameters
    # solves for coefficients using a matrix solve of modified 'normal equations'
    # turns matrix equation Ax=b into 'normal equation' (A^T)Ax=(A^T)b
    # where b = errors to minimize (actualY-fitY), A = numeric derivative of b
                
    pts = len(x)  # number of points in dataset                                  
    fitVals = []  # stores fit values
    
    # loop over specified number of iterations to update coefficients a,b,c,d
    for k in range(n):
        matrixA = np.zeros((pts,4))
        matrixB = np.zeros((pts,1))
   
       # loop to populate arrays
        for i in range(pts):            
            # populate matrix elements 
            matrixB[i][0] = y[i]-f(a,b,c,d,x[i])  # matrix of residual errors                  
            matrixA[i][0] = dfa(a,b,c,d,x[i])     # 1st column gets derivative wrt a
            matrixA[i][1] = dfb(a,b,c,d,x[i])     # 1st column gets derivative wrt b
            matrixA[i][2] = dfc(a,b,c,d,x[i])     # 1st column gets derivative wrt c
            matrixA[i][3] = dfd(a,b,c,d,x[i])     # 1st column gets derivative wrt d
        
        # create normal equation (A^T)Ax=(A^T)b 
        At = np.transpose(matrixA)              # A transpose (A^T)
        AtA = np.dot(At,matrixA)                # matrix multiply (A^T)A
        AtB = np.dot(At,matrixB)                # matrix multiply (A^T)b
        
        coeff = np.linalg.solve(AtA,AtB)        # solve normal equation (A^T)Ax=(A^T)b
                      
        # get coefficients
        a = a + coeff[0]                           
        b = b + coeff[1]
        c = c + coeff[2]
        d = d + coeff[3]
        
    
    # we now have the coefficients for the fit
    # next plug x values back into fit model and compute new y fit values
    varsum = 0     
    
    # loop to get y fit values and compute variance on fit                           
    for i in range(pts):
        val = f(a,b,c,d,x[i])
        
        fitVals.append(val)
                                
        # compute variance from erorr (error=actualY-fitY)
        error = y[i] - val
        
        varsum = varsum + error**2
    
    # using variance as measure for 'goodness of fit'
    # computed as variance = (sum of squares of errors) / (nDataPts-nFitCoefficients-1)
    variance = varsum/(pts-4-1)
    
    
    print 'coeffs=',a,b,c,d
    print 'variance=',variance
    return x,fitVals
#*******************************************************************