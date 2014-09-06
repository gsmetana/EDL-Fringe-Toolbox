# -*- coding: utf-8 -*-
# glowplug.py --- Remove glow plug from images

import numpy as np
import matplotlib.pyplot as plt
from skimage.morphology import disk
from skimage import data, transform
from skimage.feature import match_template
from skimage import filter, io, exposure
from scipy import misc, ndimage
from os.path import abspath, dirname, join as pjoin

from fringeAnalysis.data import autolite

def get_autolite_template(glowplug_width, glowplug_height):
    
    template = autolite()
    templateHeight, templateWidth = template.shape
  
    # first scale to width
    height = templateHeight*glowplug_width / templateWidth
    scaledTemplate = transform.resize(template, (height, glowplug_width))
    scaledTemplate = np.round(scaledTemplate)

    # crop bottom
    croppedTemplate = scaledTemplate[:glowplug_height, :]    
    
    return croppedTemplate

def fill_autolite_template(template):
   
    # close bottom edge
    closed = (template == 1) 
    closed[-1, 0:-1] = 1
    
    # fill 
    closed = ndimage.binary_fill_holes(closed)
    
    return np.invert(closed)

def locate_glowplug(img, template, debug=False):

    image = exposure.equalize_adapthist(img)
    image = filter.rank.median(image, disk(3))
    image = filter.rank.gradient(image, disk(3))
    imageEdge = filter.canny(image, sigma=4)

    templateHeight, templateWidth = template.shape

    # even though we have a pretty good template size,
    # check if a slightly different one fits better 

    minSearchWidth = templateWidth -3
    maxSearchWidth = templateWidth +3
    widths = np.arange(minSearchWidth, maxSearchWidth+1)

    fit = np.zeros(maxSearchWidth - minSearchWidth+1)
    xw = np.zeros(maxSearchWidth - minSearchWidth+1)
    yw = np.zeros(maxSearchWidth - minSearchWidth+1)

    for width in widths:

        height = templateHeight*width / templateWidth
        scaledTemplate = transform.resize(template, (height, width))

        result = match_template(imageEdge, scaledTemplate)
        ij = np.unravel_index(np.argmax(result), result.shape)
        x, y = ij[::-1]
        if debug:
            print 'width =', width
            print '    location:',x,',', y,' quality =', np.amax(result)

        fit[width - minSearchWidth] = np.amax(result)
        xw[width - minSearchWidth] = x
        yw[width - minSearchWidth] = y

    bestWidth = widths[np.argmax(fit)]
    bestHeight = templateHeight*bestWidth / templateWidth
    bestX = xw[np.argmax(fit)]
    bestY = yw[np.argmax(fit)]

    scaledTemplate = misc.imresize(template, (bestHeight, bestWidth), interp = 'nearest')

    if debug:
        print 'best width=',bestWidth,' location ', bestX,',', bestY
        result = match_template(imageEdge, scaledTemplate)

        fig, (ax1, ax2, ax3,ax4) = plt.subplots(ncols=4, figsize=(8, 3))
        plt.gray()
        ax1.imshow(image, interpolation='none')
        ax1.set_axis_off()
        ax1.set_title('image')

        ax2.imshow(scaledTemplate, interpolation='none')
        ax2.set_axis_off()
        ax2.set_title('template')

        ax3.imshow(imageEdge, interpolation = 'none')
        ax3.set_axis_off()
        ax3.set_title('image')
        # highlight matched region
        rect = plt.Rectangle((bestX, bestY), bestWidth, bestHeight, edgecolor='r', facecolor='none')
        ax3.add_patch(rect)

        ax4.imshow(result, interpolation='none')
        ax4.set_axis_off()
        ax4.set_title('`match_template`\nresult')
        # highlight matched region
        ax4.autoscale(False)
        ax4.plot(bestX, bestY, 'o', markeredgecolor='r', markerfacecolor='none', markersize=10)

    return (bestX, bestY, bestWidth )


def remove_glowplug(original, templateX, templateY, template):

    templateFill = fill_autolite_template(template)
    templateHeight, templateWidth = template.shape

    result = np.copy(original)
    result[templateY:templateY+templateHeight, templateX:templateX+templateWidth] =  result[templateY:templateY+templateHeight, templateX:templateX+templateWidth]*templateFill
   #print result[bestY+bestHeight:, bestX+1:bestX+bestWidth-1]    
    result[templateY+templateHeight:, templateX:templateX+templateWidth] = 0   

    return result

