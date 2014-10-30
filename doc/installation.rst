
.. toctree::
   :maxdepth: 2



Installation
============



Requirements
------------

Those are Python packages necessary to run *cochlea*.  Currently only
Python 2.7 (64-bit) is supported.

- Numpy
- Scipy
- Cython
- Pandas

- Matplotlib (optional, for examples)
- docopt (optional, for the command line scripts)
- thorns_ (optional, for examples and stats)
- matlab_wrapper_ (optional, for the MAP external model)


Note: On Windows you can install a Python distribution such as
Anaconda_ or `Python(x,y)`_ to fulfill most/all of the dependencies.
Make sure that you have 64-bit version of Python.


.. _thorns: https://github.com/mrkrd/thorns
.. _matlab_wrapper: https://github.com/mrkrd/matlab_wrapper
.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _`Python(x,y)`: https://code.google.com/p/pythonxy/




Installation from PyPI
----------------------

This is the easiest and recommended way to install *cochlea*.  On
GNU/Linux/BSD, it will download the source package from Python Package
Index (PyPI) repository, compile and install it on your system.  On
Windows, it will download a pre-compiled binary package and install
it.  Because the package is already compiled, you will not need Cython
nor C compiler on Windows.  Make sure that you have at least `Numpy`,
`Scipy` and `Pandas` on your system, otherwise the installation
program (pip) will try to install them from PyPI.

Type in the console/terminal as root/administrator::

  pip install cochlea

or to install for a single user (no root needed)::

  pip install --user cochlea





Installation from Source Code
-----------------------------



First, clone the repository::

  git clone https://github.com/mrkrd/cochlea.git


Second, inside of the repository. call the installation script::

  python setup.py install


The above is the standard way of installing a Python package from a
source code.  However, if I want follow changes to the source code
closely, I prefer to do the following::

  python setup.py develop --user


which installs just a link to the cloned repository.  In this way, any
change in the source code (git repo) is immediately reflected in the
installation (no re-installation necessary).

It might be necessary to compile C/Cython code::

  make




Installation from Source Code on Windows
----------------------------------------


Generally, the installation steps are the same as on GNU/Linux/BSD,
but you will need to set up a proper environment to compile the
Cython/C extensions.

I followed roughly the instructions from
64BitCythonExtensionsOnWindows_.  First of all, you will need Windows
SDK (.NET 3.5 SP1).  Next, I installed necessary requirements from the
`Unofficial Windows Binaries for Python Extension Packages`_
repository (I'm not sure, if it would work with Anaconda's
Numpy/Scipy/etc).  And finally, you can run a script from the
repository to compile the binary extensions::

  build_windows_dist.bat

Followed by::

  python setup.py install



If you get an error about missing ``_pycat`` while importing cochlea
in Python, most likely the binary sub-modules are not
compiled/installed.


.. _64BitCythonExtensionsOnWindows: https://github.com/cython/cython/wiki/64BitCythonExtensionsOnWindows
.. _`Unofficial Windows Binaries for Python Extension Packages`: http://www.lfd.uci.edu/~gohlke/pythonlibs/
