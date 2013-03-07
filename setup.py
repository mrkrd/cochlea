#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

extensions = [
    Extension(
        "cochlea.zilany2009._pycat",
        [
            "cochlea/zilany2009/_pycat.pyx",
            "cochlea/zilany2009/catmodel.c",
            "cochlea/zilany2009/complex.c"
        ]
    ),
    Extension(
        "cochlea.holmberg2007._traveling_waves",
        [
            "cochlea/holmberg2007/_traveling_waves.pyx",
        ]
    ),
]


setup(
    name = "cochlea",
    version = "6",
    description = "Collection of inner ear models",
    author = "Marek Rudnicki",
    author_email = "marek.rudnicki@tum.de",
    packages = [
        "cochlea",
        "cochlea.stats",
        "cochlea.pycat",
        "cochlea.traveling_waves",
    ],
    package_data = {
        "cochlea": ["*.csv"]
    },
    include_dirs = [numpy.get_include()],
    ext_modules = cythonize(extensions)
)
