#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import numpy
numpy_include = numpy.get_include()

setup(
    name = "pycat",
    version = "1.0",
    author = "Marek Rudnicki",
    packages = ["pycat"],
    package_dir = {"pycat": "src/pycat"},
    package_data = {"pycat": ["data/*.txt"]},
    ext_package = "pycat",
    ext_modules = [
        Extension("_pycat",
                  ["src/pycat/_pycat.pyx", "src/pycat/catmodel.c",
                   "src/pycat/complex.c"],
                  include_dirs=[numpy_include])
        ],
    cmdclass = {"build_ext": build_ext}
)
