""" Glow Plug template

"""

import os as _os

from skimage.io import imread
from fringeAnalysis import data_dir


__all__ = ['load',
           'autolite']


def load(f):
    """Load an image file located in the data directory.

    Parameters
    ----------
    f : string
        File name.

    Returns
    -------
    img : ndarray
        Image loaded from skimage.data_dir.
    """
    return imread(_os.path.join(data_dir, f))


def autolite():
    """ Autolite glowplug

    """
    return load("autolite_mask180x536.jpg")

