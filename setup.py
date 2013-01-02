#!/usr/bin/env python

from distutils.core import setup
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
    ext_modules = cythonize(
        ["src/_pycat.pyx", "src/traveling_waves/_tw.pyx"],
        include_path=[numpy_include]
    )
)
