import numpy as np
from scipy import stats

def unbias_image(image):
    """Remove the background bias from image"""

    n, m = image.shape
    # simply fit a line at boundary
    bias = image[:,0]
    slope, intercept, r_value, p_value, std_err = stats.linregress(range(0,n), bias)
    bias = [i * slope for i in range(n)] + intercept
    # extend to 2D??

    background = np.outer(bias,np.ones((m,)))
    return image - background
