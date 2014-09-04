EDL-Fringe-Toolbox
==================

Collection of scripts to process interferogram fringes

The documentation in EDL-Fringe-Toolbox/doc aims to be a complete description of the toolbox and the methods it uses. A set of scripts illustrating its use is included in a repo at https://github.com/gsmetana/EDL-Fringe-demos

Source
------
https://github.com/gsmetana/EDL-Fringe-Toolbox/

Installation from source
------------------------
Refer to DEPENDS.txt for a list of dependencies.

The toolbox may be installed globally using

    $ python setup.py install

or locally using

    $ python setup.py install --prefix=${HOME}

If you prefer, you can use it without installing, by simply adding
this path to your PYTHONPATH variable and compiling the extensions:

    $ python setup.py build_ext -i
