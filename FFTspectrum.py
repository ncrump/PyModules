"""
Created on Mon Nov 11 10:15:15 2013
"""

# FFT Spectrum
"""
Computes the FFT of a time signal and returns the real one-sided spectrum.
"""

import numpy as np


#****************************************************
#----------------------------------------------------
# t = discrete time points of signal (sec)
# f = discrete values of signal at time points
# rate = sample rate (Hz)
#----------------------------------------------------
def FFTspectrum(t,f,rate):
    
    # get number of points in dataset
    pts = len(t)
    
    # compute FFT and take only real part
    s = abs(np.fft.fft(f))
    # take only one side of spectrum 
    spectrum = s[0:pts/2.0]
    
    # create freq spectrum
    npts = len(spectrum)
    fs = np.linspace(0,rate/2.0,npts)
    
    # return spectrum in freq domain
    return fs, spectrum
#****************************************************