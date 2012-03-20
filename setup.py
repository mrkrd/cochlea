#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
numpy_include = numpy.get_include()

setup(
    name = "cochlea",
    version = "3",
    author = "Marek Rudnicki",
    packages = [
        "cochlea",
        "cochlea.stats"
    ],
    package_dir = {
        "cochlea": "src",
        "cochlea.stats": "src/stats"
    },
    package_data = {"cochlea": ["data/*.txt", "pars/*.par"]},
    ext_package = "cochlea",
    ext_modules = [
        Extension("_pycat",
                  ["src/_pycat.pyx", "src/catmodel.c",
                   "src/complex.c"],
                  include_dirs=[numpy_include]),
        Extension("_tw",
                  ["src/traveling_waves/_tw.pyx"],
                  include_dirs=[numpy_include])
        ],
    cmdclass = {"build_ext": build_ext}
)
