cochela
=======

*Inner ear models in Python*


::

                           +-----------+     __|____|_____|____
   .-.     .-.     .-.     |           |-->  _|______|_____|___
  /   \   /   \   /   \ -->|  Cochlea  |-->  ___|____|___|_____
       '-'     '-'         |           |-->  __|____|_____|____
                           +-----------+
            Sound                               Spike Trains



:Name: cochela
:Author: Marek Rudnicki
:Email: marek.rudnicki@tum.de
:URL: https://github.com/mrkrd/cochlea
:License: GNU General Public License v3 or later (GPLv3+)



Description
-----------

*cochlea* is a collection of inner ear models.  They are easily
accessible as Python functions.  The models take sound signal as input
and return spike trains of the auditory nerve fibers.



Models
------

- Holmberg, M. (2007). Speech Encoding in the Human Auditory
  Periphery: Modeling and Quantitative Assessment by Means of
  Automatic Speech Recognition. PhD thesis, Technical University
  Darmstadt.
- Zilany, M. S., Bruce, I. C., Nelson, P. C., &
  Carney, L. H. (2009). A phenomenological model of the synapse
  between the inner hair cell and auditory nerve: long-term adaptation
  with power-law dynamics. The Journal of the Acoustical Society of
  America, 126(5), 2390-2412.
- Zilany, M. S., Bruce, I. C., & Carney, L. H. (2014). Updated
  parameters and expanded simulation options for a model of the
  auditory periphery. The Journal of the Acoustical Society of
  America, 135(1), 283-286.
- MATLAB Auditory Periphery by Meddis et al. (external model, not
  implemented in the package, but easily accessible through
  matlab_wrapper_).


.. _matlab_wrapper: https://github.com/mrkrd/matlab_wrapper


Usage
-----

Initialize the modules::

  import cochlea
  import thorns as th
  import thorns.waves as wv


Generate sound::

  fs = 100e3
  sound = wv.ramped_tone(
      fs=fs,
      freq=1000,
      duration=0.1,
      dbspl=50
  )


Run the model (responses of 200 cat HSR fibers)::

  anf_trains = cochlea.run_zilany2014(
      sound,
      fs,
      anf_num=(200,0,0),
      cf=1000,
      seed=0,
      species='cat'
  )


Plot the results::

  th.plot_raster(anf_trains)
  th.show()



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
- matlab_wrapper (optional, for the MAP external model)



Installation
------------

Quick install::

  pip install cochlea


Developer::

  make
  python setup.py develop --user


Administrator::

  python setup.py install


User::

  python setup.py install --user






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





Acknowledgments
---------------
