#!/usr/bin/env python

import os

from skimage._build import cython

base_path = os.path.abspath(os.path.dirname(__file__))


def configuration(parent_package='', top_path=None):
    from numpy.distutils.misc_util import Configuration, get_numpy_include_dirs

    print base_path

    config = Configuration('mask', parent_package, top_path)
    config.add_data_files('autolite_mask180x536.jpg', 'autolite_mask180x300.jpg')
    return config


if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(maintainer='scikit-image Developers',
          author='scikit-image Developers',
          maintainer_email='scikit-image@googlegroups.com',
          description='Restoration',
          url='https://github.com/scikit-image/scikit-image',
          license='SciPy License (BSD Style)',
          **(configuration(top_path='').todict())
        )
