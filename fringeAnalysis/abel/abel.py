import numpy as np
import scipy.fftpack 
fft = scipy.fftpack.fft
ifft = scipy.fftpack.ifft
from scipy.optimize import fmin
from matplotlib import pyplot as plt
from scipy import signal

def gaussian(x, mu, sigma):
    return(
        np.exp(-(x-mu)**2/sigma**2/2.0)/(np.sqrt(2)*np.sqrt(np.pi)*abs(sigma))
        )

def dgaussiandx(x, mu, sigma):
    return(
        -(x-mu)*np.exp(-(x-mu)**2/sigma**2/2.0)/(sigma**2*np.sqrt(2*np.pi*sigma**2))
         )
    
def periodic_gaussian(x, sigma):
    nx = len(x)
    # put the peak in the middle of the array: 
    mu = x[nx/2]
    g = gaussian(x, mu, sigma)
    # need to normalize, the discrete approximation will be a bit off:
    g = g / sum(g) 
    # reorder to split the peak accross the boundaries: 
    return(np.append(g[nx/2:nx], g[0:nx/2]))

def periodic_gaussian_deriv(x, sigma):
    nx = len(x)
    # put the peak in the middle of the array: 
    mu = x[nx/2]
    g = dgaussiandx(x, mu, sigma)
    # need to normalize, the discrete approximation will be a bit off:
    g = g / (-np.sum(x * g))
    # reorder to split the peak accross the boundaries: 
    return(np.append(g[nx/2:nx], g[0:nx/2]))

def abel(x, dfdx):
    nx = len(x)
    # do each half of the signal seperately (they should be the same
    # up to noise)
    integral = np.zeros((2,nx/2), dtype=float)
    for i in xrange(nx/2, nx-1):
        divisor = np.sqrt(x[i:nx]**2 - x[i]**2)
        integrand = dfdx[i:nx] / divisor
        integrand[0] = integrand[1] # deal with the singularity at x=r
        integral[0][i-nx/2] = - np.trapz(integrand, x[i:nx]) / np.pi

    for i in xrange(nx/2, 1, -1):
        divisor = np.sqrt(x[i:0:-1]**2 - x[i]**2)
        integrand = dfdx[i:0:-1] / divisor
        integrand[0] = integrand[1] # deal with the singularity at x=r
        integral[1][-i+nx/2] = - np.trapz(integrand, x[i:0:-1]) / np.pi
    return(integral)


def smooth(x, data,sigma_reg=5):
    # 'regularize'
    sigma_reg = 5 # about half our 'true' sigma 
    dx = x[1]-x[0]
    # smooth it
    fmeas = fft(data) * fft(periodic_gaussian(x, sigma_reg*dx))
    return ifft(fmeas ).real


'''
perform abel inversion on data with possible center discontinuity (value = 0)
'''
def perform_abel_discontinuous(x, data):

    if data.ndim == 1:
        data = np.atleast_2d(data)

    ny, nx = data.shape

    dx = x[1]-x[0]
    r = x[nx/2:nx-1]
    nr = len(r)

    smoothed = np.zeros( (ny,nx) )
    deriv = np.zeros( (ny,nx) )
    f_abel = np.zeros( (ny,nr) )
    sigma_deriv = 2

    for i in range(0,ny):

        w_discont = int(np.where( data[i,nx/2:] > 0)[0][0])

        if w_discont == 0:
            smoothed[i, :] = signal.wiener(data[i,:], mysize=10)
            deriv[i,:]= ifft(fft(smoothed[i,:]) * fft(periodic_gaussian_deriv(x, sigma_deriv*dx))).real
            tmp= abel(x, deriv[i,:])
            f_abel[i,:]= tmp[0]

        else:
            # avoid performing operations over discontinuity
            left = range(0, int(nx/2 - w_discont +1 ))
            right = range( int(nx/2 + w_discont) , nx )

            smoothed[i, left] = signal.wiener(data[i, left], mysize=10)
            smoothed[i, right] = signal.wiener(data[i, right], mysize=10)        
            deriv[i,:]= ifft(fft(smoothed[i,:]) * fft(periodic_gaussian_deriv(x, sigma_deriv*dx))).real
            #deriv[i, left] = ifft(fft(smoothed[i,left]) * fft(periodic_gaussian_deriv(x[left], sigma_deriv*dx))).real
            #deriv[i, right] = ifft(fft(smoothed[i,right]) * fft(periodic_gaussian_deriv(x[right], sigma_deriv*dx))).real  
            
            tmp= abel(x, deriv[i,:])
            f_abel[i,:]= tmp[0]
            f_abel[i, 0: w_discont+1] = 0

    return (r, f_abel, deriv, smoothed)

