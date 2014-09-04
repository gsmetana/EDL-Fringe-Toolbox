# -*- coding: utf-8 -*-
# abel.py --- Perform Abel Inversion


import numpy as np
import matplotlib.pyplot as plt
from skimage.morphology import disk
from skimage import data
from skimage.feature import match_template
from skimage import filter, io, exposure
from scipy import misc
from os.path import abspath, dirname, join as pjoin

def remove_glowplug(original, imageToMask, template, debug = False, outputPlot = False):

    image = exposure.equalize_adapthist(original)
    image = filter.rank.median(image, disk(4))
    image = filter.rank.gradient(image, disk(4))
    imageEdge = filter.canny(image, sigma=3)

   # templateFilename = 'autolite_mask180x300.jpg'
   # templatePath = pjoin(dirname(abspath(__file__)), templateFilename)
   # template = io.imread(templatePath, as_grey=True)
    templateHeight, templateWidth = template.shape

    minSearchHeight = 178
    maxSearchWidth = 186
    widths = np.arange(minSearchWidth, maxSearchWidth+1)



    fit = np.zeros(maxSearchWidth - minSearchWidth+1)
    xw = np.zeros(maxSearchWidth - minSearchWidth+1)
    yw = np.zeros(maxSearchWidth - minSearchWidth+1)

    for width in widths:

        height = templateHeight*width / templateWidth
        scaledTemplate = misc.imresize(template, (height, width) , interp = 'nearest')
        scaledTemplateEdge = filter.canny(scaledTemplate, sigma=1)

        result = match_template(imageEdge, scaledTemplateEdge)
        ij = np.unravel_index(np.argmax(result), result.shape)
        x, y = ij[::-1]
        if debug:
            print 'width =', width
            print x, y, np.amax(result)

        fit[width - minSearchWidth] = np.amax(result)
        xw[width - minSearchWidth] = x
        yw[width - minSearchWidth] = y



    bestWidth = widths[np.argmax(fit)]
    bestHeight = templateHeight*bestWidth / templateWidth
    bestX = xw[np.argmax(fit)]
    bestY = yw[np.argmax(fit)]

    scaledTemplate = misc.imresize(template, (bestHeight, bestWidth), interp = 'nearest')

    if debug:
        print np.amax(fit)
        print bestWidth, bestX, bestY
        scaledTemplateEdge = filter.canny(scaledTemplate, sigma=1)
        result = match_template(imageEdge, scaledTemplateEdge)

        fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, figsize=(8, 3))

        ax1.imshow(scaledTemplate)
        ax1.set_axis_off()
        ax1.set_title('template')

        ax2.imshow(imageEdge)
        ax2.set_axis_off()
        ax2.set_title('image')
        # highlight matched region
        rect = plt.Rectangle((bestX, bestY), bestWidth, bestHeight, edgecolor='r', facecolor='none')
        ax2.add_patch(rect)

        ax3.imshow(result)
        ax3.set_axis_off()
        ax3.set_title('`match_template`\nresult')
        # highlight matched region
        ax3.autoscale(False)
        ax3.plot(bestX, bestY, 'o', markeredgecolor='r', markerfacecolor='none', markersize=10)
        plt.show()

    scaledTemplateMask = (scaledTemplate == 255)
    #show resulting figure with glowplug removed


    result = 1*imageToMask
    result[bestY:bestY+bestHeight, bestX:bestX+bestWidth] = result[bestY:bestY+bestHeight, bestX:bestX+bestWidth]*scaledTemplateMask
   #print result[bestY+bestHeight:, bestX+1:bestX+bestWidth-1]    
    result[bestY+bestHeight:, bestX+1:bestX+bestWidth-1] = 0   

    if outputPlot:
        fig2, (ax1, ax2) = plt.subplots(ncols=2, figsize=(8, 3))

        ax1.imshow(imageToMask, cmap=plt.cm.gray)
        ax1.set_axis_off()
        ax1.set_title('original')

        ax2.imshow(result, cmap=plt.cm.gray)
        ax2.set_axis_off()
        ax2.set_title('result')
        ax2.autoscale(False)

        plt.show()

    return result

