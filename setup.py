#!/usr/bin/env python

from distutils.core import setup
from Cython.Build import cythonize


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
    package_data = {
        "cochlea": ["data/*.txt", "pars/*.par"]
    },
    ext_package = "cochlea",
    ext_modules = cythonize(
        ["cochlea/_pycat.pyx", "cochlea/traveling_waves/_tw.pyx"]
    )
)
