Inner ear models in Python.




Currently implemented models
============================

Holmbert (2007)

Zilany et al. (2009)

Zilany et al. (2013/2014)






Requirements
============

Python
Numpy
Scipy
Cython
Pandas

matplotlib (optional, for demos)
docopt (optional, for the command line scripts)
mrlib (optional, for demos and stats)




Installation
============

Developer
---------

./build_inplace.sh
python setup.py develop --user


Administrator
-------------

python setup.py install


User
----

python setup.py install --user






Directory structure
===================

cochlea: model implementation, e.g.,

cochlea/holmberg2007
cochlea/zilany2009
cochlea/zilany2013

demos: small demos scripts, a good place to start

tests: unit tests (nose)

scripts: command line interface




Usage
=====

Examples are in the `demos' directory.
