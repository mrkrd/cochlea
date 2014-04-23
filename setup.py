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
    Extension(
        "cochlea.zilany2013._zilany2013",
        [
            "cochlea/zilany2013/_zilany2013.pyx",
            "cochlea/zilany2013/model_IHC.c",
            "cochlea/zilany2013/model_Synapse.c",
            "cochlea/zilany2013/complex.c"
        ]
    ),
]


setup(
    name="cochlea",
    version="0.6",
    description="Collection of inner ear models",
    author="Marek Rudnicki",
    author_email="marek.rudnicki@tum.de",
    packages=[
        "cochlea",
        "cochlea.stats",
        "cochlea.zilany2009",
        "cochlea.zilany2013",
        "cochlea.holmberg2007",
    ],
    package_data={
        "cochlea": ["*.csv"]
    },
    include_dirs=[numpy.get_include()],
    ext_modules=cythonize(extensions),
    scripts=["scripts/run_zilany2014"],
)
