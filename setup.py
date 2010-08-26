#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
numpy_include = numpy.get_include()

setup(
    name = "traveling_waves",
    version = "1.0",
    author = "Marek Rudnicki",
    packages = ["traveling_waves"],
    package_dir = {"traveling_waves": "src/traveling_waves"},
    ext_package = "traveling_waves",
    ext_modules = [
        Extension("_tw",
                  ["src/traveling_waves/_tw.pyx"],
                  include_dirs=[numpy_include])
        ],
    cmdclass = {"build_ext": build_ext}
)
