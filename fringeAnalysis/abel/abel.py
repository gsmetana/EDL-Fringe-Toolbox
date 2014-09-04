import numpy as np
import scipy.fftpack 
fft = scipy.fftpack.fft
ifft = scipy.fftpack.ifft
from scipy.optimize import fmin
from matplotlib import rc
#rc('text', usetex=True)
from matplotlib import pyplot as plt

'''
\delta n = K \delta rho

\delta rho = \lambda / L * \delta N / K


K_N2 = 2.38 * 10^-4 m^3 /kg
'''
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

def abel(dfdx, x):
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

def get_shift(data):
	indices = []
	nx = len(data)
	try_shift = range(0, nx/10) + range(nx - nx/10,nx) 
	for i in try_shift: 
	    indices.append(np.append(range(i,nx), range(0,i)).astype(int))
	
	even_norm = np.zeros(len(indices), dtype=float)
	for i, cyclic_perm in enumerate(indices):
	    even_norm[i] = np.sum(fft(data[cyclic_perm]).real**2)
	max_even_shift = try_shift[even_norm.argmax()]
	print max_even_shift

	if max_even_shift > nx/2:
		nshift = nx - max_even_shift
	else:
		nshift = max_even_shift # ???
		
	return nshift, indices[even_norm.argmax()]

'''
perform abel inversion on row of processed data
'''
def perform_abel(data, d_glowplug, pix_glowplug):
    #d_glowplug = 0.0051 # [m] from phil
    #pix_glowplug = 182
    ny, nx = data.shape
    dx = d_glowplug / pix_glowplug
    x = dx*np.arange(-nx/2, nx/2)

    smoothed = np.zeros( (ny,nx) )
    for i in range(0,ny):
	    smoothed[i,:] = smooth(x, data[i,:])

    nshift, shift_ind = get_shift(smoothed[ny/2, :]) # pick more intelligently?

    shifted = np.zeros( (ny,nx) )
    for i in range(0,ny):
	    shifted[i,:] = smoothed[i, shift_ind]

    if nshift > 0:
     	xs = x[2*nshift+1:]  # still need the +1 if nx even/odd?
	    #for i in range(0,ny):
        shifted = shifted[:, 2*nshift+1:]
    else:
     	xs = x[2*nshift+1:]  # still need the +1 if nx even/odd?
	    #for i in range(0,ny):
        shifted = shifted[:, 2*nshift+1:]
        
    nxs = len(xs)
    r = xs[nxs/2:nxs -1]
    nr = len(r)

    # calculate some derivatives: 
    sigma_deriv = 5 # this can be small because we've already smoothed a lot
    deriv = np.zeros( (ny,nxs) )
    f_abel = np.zeros( (ny,nr) )
    for i in range(0,ny): 
	    deriv[i,:]= ifft(fft(shifted[i,:]) * fft(periodic_gaussian_deriv(xs, sigma_deriv*dx))).real
	    # do some inverse Abel transforms: 
	    # assume the signal is centered in x 
	    tmp= abel(deriv[i,:], xs)
	    #f_abel[i,:]= (tmp[0] - tmp[1])/2.0
	    f_abel[i,:]= tmp[0]

    return f_abel
'''
plt.figure()
plt.imshow(f_abel)

K = 2.38e-4
lam = 532e-9
delta_rho = f_abel / K * lam / np.pi
plt.figure()
plt.imshow(delta_rho)
print delta_rho
plt.colorbar()
plt.show()

plt.figure(figsize=(11,8.5))
plt.plot(x, data, label='measurement')
plt.plot(x, cyl_deconv_smooth, label='smoothed')
plt.plot(xs, cyl_shifted, label='shifted')
plt.legend(loc=0)
plt.xlabel(r'$y$')
plt.ylabel(r'$F(r)$')

plt.figure(figsize=(11,8.5))
plt.plot(xs, cyl_deriv)
plt.xlabel(r'$y$')
plt.ylabel(r'$\frac{\partial F(r)}{\partial y}$')
plt.legend(loc=0)

plt.figure(figsize=(11,8.5))
plt.plot(r, f_abel[0])
plt.legend(loc=0)
plt.xlabel(r'$r$')
plt.ylabel(r'$f(r)$')
plt.show()
'''


