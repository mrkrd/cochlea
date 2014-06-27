cochela
=======

*Inner ear models in Python*




Currently implemented models
----------------------------

- Holmbert (2007)
- Zilany et al. (2009)
- Zilany et al. (2013/2014)


External
--------

- MATLAB Auditory Periphery (Meddis et al.)




Requirements
------------

- Python (2.7)
- Numpy
- Scipy
- Cython
- Pandas

- matplotlib (optional, for examples)
- docopt (optional, for the command line scripts)
- thorns (optional, for examples and stats)




Installation
------------

Developer::

  make
  python setup.py develop --user


Administrator::

  python setup.py install


User::

  python setup.py install --user





Directory structure
-------------------

cochlea: model implementation, e.g.,

cochlea/holmberg2007
cochlea/zilany2009
cochlea/zilany2014

external: external models

examples: small demo scripts, a good place to start

tests: unit tests (nose)

scripts: command line interface




Usage
-----

Please check examples directory.




Other implementations
---------------------

- `Carney Lab`_
- `Matlab Auditory Periphery`_
- DSAM_
- `Brian Hears`_

.. _`Carney Lab`: http://www.urmc.rochester.edu/labs/Carney-Lab/publications/auditory-models.cfm
.. _DSAM: http://dsam.org.uk/
.. _`Matlab Auditory Periphery`: http://www.essexpsychology.macmate.me/HearingLab/modelling.html
.. _`Brian Hears`: http://www.briansimulator.org/docs/hears.html
