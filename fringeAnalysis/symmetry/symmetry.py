import numpy as np

def center_image_x(image, centerX, centerY, dx):
    
    ny, nx = image.shape

    # have to think carefully about this
    # want to end up with nx odd and have x=0 at centerX

    if nx/2 >= centerX:
        radius = centerX - 1
        centered = image[:, 0:2*radius + 1]

    else:
        radius = nx - centerX
        centered = image[:, nx - 2*radius - 1:nx]

    x = dx*np.arange(-radius, radius + 1)
    y = dx*np.arange(-centerY, ny - centerY)

    return (x,y, centered)

def symmetry_error_x(image):
    ny, nx = image.shape
    left = np.fliplr(image[:,0:nx/2+1] )
    right = image[:,nx/2:]
    return (left - right)/left

