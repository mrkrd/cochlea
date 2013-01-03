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
        "cochlea.stats"
    ],
    package_dir = {"cochlea": "src"},
    package_data = {
        "cochlea": ["data/*.txt", "pars/*.par"]
    },
    ext_package = "cochlea",
    include_dirs = [numpy.get_include()],
    ext_modules = cythonize(
        ["src/_pycat.pyx", "src/traveling_waves/_tw.pyx"]
    )
)
