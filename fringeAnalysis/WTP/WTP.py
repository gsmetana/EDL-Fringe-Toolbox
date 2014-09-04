import numpy as np
from skimage import io, exposure
from struct import pack, unpack
from matplotlib import pyplot as plt
import sys,os,os.path

from ._perform_WTP import perform_WTP
from .constants import *

m = -1
n = -1

fringeDatFilename = ''
wrappedDatFilename = ''

def imagefile2dat(imageFilename, rotate = False, overwrite = False):
    """Load an image file and save in format to be read by C code"""

    global m
    global n
    global fringeDatFilename
    global wrappedDatFilename

    # read image file
    orig = io.imread(imageFilename, as_grey=True)

    img = exposure.equalize_adapthist(orig)
    img = exposure.rescale_intensity(img,out_range=(0, 255))
    if rotate:
        img = np.transpose(img)

    n = len(img) 
    m = len(img[0])

    fileroot, ext = os.path.splitext(imageFilename)
    fringeDatFilename = fileroot+'.dat'
    wrappedDatFilename = fileroot+'W.dat'
    
    if os.path.isfile(fringeDatFilename) == False or overwrite == True:
        print 'Writing '+fringeDatFilename
        # write in proper binary format
        data = np.reshape( np.transpose(img), (n*m,1))
        newFile = open (fringeDatFilename, "wb")
        newFile.write(pack(str(n*m)+'B', *data))
        newFile.close()
    else:
        print 'Skipped overwriting '+fringeDatFilename

    return img

def performWTP(to_extend_image = NO, extend_by_pixel = 50, use_FFT = YES,
    wavelet_type=PAUL4, ridge_alg=MAX, starting_scale = 3, 
    ending_scale = 200, numDaughters = 50, Morlet_sigma = 1, overwrite = False):
    """Execute C code to perform WTP"""

    global n
    global m
    global fringeDatFilename
    global wrappedDatFilename

    if n == -1 or m == -1:
        raise NameError(fringeDatFilename+' must be created before processing')

    if os.path.isfile(wrappedDatFilename) == False or overwrite == True:
        scale_step = (ending_scale - starting_scale)/ float(numDaughters)
        perform_WTP(fringeDatFilename, wrappedDatFilename, m,n, to_extend_image, extend_by_pixel, use_FFT,
                wavelet_type, ridge_alg, starting_scale, scale_step, ending_scale, Morlet_sigma)
    
    else:
        print 'Skipped overwriting '+wrappedDatFilename

    return dat2image(wrappedDatFilename)

def dat2image(datFilename):
    """Read image from C code output"""
    global n
    global m
    
    print 'Reading '+datFilename

    if n == -1 or m == -1 :
        raise NameError(datFilename+' must be created before processing')

    # read binary information
    newFile = open(datFilename, "rb")
    bytes_read = newFile.read()
    newFile.close()
    img = unpack(str(n*m)+'f', bytes_read)
    img = np.reshape(img, (m, n))
    return np.transpose(img)
