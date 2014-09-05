# -*- coding: utf-8 -*-
# sphere.py --- Remove sphere from images

import numpy as np
import matplotlib.pyplot as plt

from skimage import data, filter, color
from skimage.transform import hough_circle
from skimage.feature import peak_local_max
from skimage import draw
from skimage.util import img_as_ubyte

from skimage import filter, io, exposure

def remove_sphere(image, imageToMask, r_min, r_max, debug = False, outputPlot = False):

    edges = filter.canny(image, sigma=3, low_threshold=10, high_threshold=50)

    print 'hi3'

    # Detect two radii
    hough_radii = np.arange(r_min, r_max, 2)
    hough_res = hough_circle(edges, hough_radii)

    centers = []
    accums = []
    radii = []

    for radius, h in zip(hough_radii, hough_res):
        # For each radius, extract two circles
        peaks = peak_local_max(h, num_peaks=2)
        centers.extend(peaks)
        accums.extend(h[peaks[:, 0], peaks[:, 1]])
        radii.extend([radius, radius])

    idx = np.argsort(accums)[::-1][0]
    center_x, center_y = centers[idx]
    radius = radii[idx]
    cx, cy = draw.circle(center_y, center_x, radius)

    masked = np.copy(imageToMask)
    masked[cy, cx] = 0

    return masked
