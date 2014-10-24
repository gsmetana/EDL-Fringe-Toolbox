import numpy as np
from skimage import io, exposure
from struct import pack, unpack
from matplotlib import pyplot as plt
import sys,os,os.path

from ._perform_WFT import perform_WFR

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
    img = img - np.mean(img)
    img = exposure.rescale_intensity(img,out_range=(0, 255))
    img = np.round(img)


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

def performWFR(sigmax, wxl, wxi,wxh, sigmay,wyl, wyi, wyh, thr, overwrite = False):
    """Execute C code to perform WFR"""

    global n
    global m
    global fringeDatFilename
    global wrappedDatFilename

    if n == -1 or m == -1:
        raise NameError(fringeDatFilename+' must be created before processing')

    if os.path.isfile(wrappedDatFilename) == False or overwrite == True:
        perform_WFR(fringeDatFilename, wrappedDatFilename, m,n, sigmax, wxl, wxi,wxh, sigmay,wyl, wyi, wyh, thr)
    
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
    img = unpack(str(n*m)+'d', bytes_read)
    img = np.reshape(img, (m, n))
    return np.transpose(img)

