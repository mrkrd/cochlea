#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

import numpy
numpy_include = numpy.get_include()

setup(
    name = "cochlea",
    version = "4",
    description = "Collection of inner ear models",
    author = "Marek Rudnicki",
    author_email = "marek.rudnicki@tum.de",
    packages = [
        "cochlea",
        "cochlea.stats"
    ],
    package_dir = {"cochlea": "src"},
    package_data = {
        "cochlea": ["data/*.txt", "pars/*.par"]
    },
    ext_package = "cochlea",
    ext_modules = cythonize([
        Extension(
            "_pycat",
            ["src/_pycat.pyx", "src/catmodel.c", "src/complex.c"],
            include_dirs=[numpy_include]
        ),
        Extension(
            "_tw",
            ["src/traveling_waves/_tw.pyx"],
            include_dirs=[numpy_include]
        )
    ])
)
