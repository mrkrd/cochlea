#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy

setup(
    name = "cochlea",
    version = "4",
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
        "cochlea": ["data/*.csv", "pars/*.par"]
    },
    include_dirs = [numpy.get_include()],
    ext_modules = cythonize([
        Extension("cochlea.pycat._pycat", ["cochlea/pycat/_pycat.pyx"]),
        Extension("cochlea.traveling_waves/_tw", ["cochlea/traveling_waves/_tw.pyx"]),
    ])
)
