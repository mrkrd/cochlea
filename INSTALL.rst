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


Note: We recommend Anaconda_ Python distribution for Windows users.
It comes with all required dependencies in one package.


.. _thorns: https://github.com/mrkrd/thorns
.. _matlab_wrapper: https://github.com/mrkrd/matlab_wrapper
.. _Anaconda: https://store.continuum.io/cshop/anaconda/




Installation from PyPI
----------------------

This is the easiest and recommended way to install *cochlea*.

On GNU/Linux/BSD, it will download the source package from Python
Package Index (PyPI) repository, compile and install it on your
system.  On Windows, it will download a pre-compiled binary package
and install it.  Because the package is already compiled, you will not
need Cython nor C compiler on Windows.  Make sure that you have at
least `Numpy`, `Scipy` and `Pandas` on your system, otherwise the
installation program (pip) will try (and probably fail) to install
them from PyPI.

Type in the console/terminal as root/administrator::

  pip install cochlea

or to install for a single user (no root needed)::

  pip install --user cochlea





Installation from Source Code
-----------------------------



First, clone the repository::

  git clone https://github.com/mrkrd/cochlea.git


Second, call the installation script inside of the repository::

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
CythonExtensionsOnWindows_.  First of all, you will need Windows
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


.. _CythonExtensionsOnWindows: https://github.com/cython/cython/wiki/CythonExtensionsOnWindows
.. _`Unofficial Windows Binaries for Python Extension Packages`: http://www.lfd.uci.edu/~gohlke/pythonlibs/
