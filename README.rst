cochela
=======

*Inner ear models in Python*


:Name: cochela
:Author: Marek Rudnicki
:Email: marek.rudnicki@tum.de
:URL: https://github.com/mrkrd/cochlea
:License: GNU General Public License v3 or later (GPLv3+)



Description
-----------

*cochlea* is a collection of inner ear models.  The models are easily
accessible as Python functions.  They take sound signal as input and
return spike trains of the auditory nerve fibers::



                           +-----------+     __|______|______|____
   .-.     .-.     .-.     |           |-->  _|________|______|___
  /   \   /   \   /   \ -->|  Cochlea  |-->  ___|______|____|_____
       '-'     '-'         |           |-->  __|______|______|____
                           +-----------+
            Sound                               Spike Trains



The implemented models are state-of-the-art biophysical models and
give realistic approximation of the auditory nerve activity.

They are implemented using the original code from their authors as
much as possible.  We wanted to have implementations that gives the
same results as the original models.  Whenever possible, it was
verified with unit testing (tests directory).

The models are also fast.  It is easy to generate responses of
hundreds or even thousands of auditory nerve fibers (ANFs).  It is
possible, for example, to generate responses of the whole human
auditory nerve (around 30,000 ANFs).  We usually tested the models
with sounds up to 1 second in duration.

I developed *cochlea* during my PhD in the group of `Werner Hemmert`_
at TUM.  It went through several versions and rewrites.  Now, it is
quite stable and we decided to release it for the community.

*cochlea* could be used by:

- researchers doing computational neuroscience (combines well with
  NEURON_ and Brian_),
- experimenters that require realistic inputs spike trains,
- researchers improving inner ear models,
- everyone interested in understanding how hearing works.


.. _`Werner Hemmert`: http://www.imetum.tum.de/research/bai/home/?L=1
.. _NEURON: http://www.neuron.yale.edu/neuron/
.. _Brian: http://briansimulator.org/


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



Spike Train Format
------------------

Spike train data format is based on a standard DataFrame_ format from
an excellent pandas_ library.  Spike trains and their meta data are
stored in DataFrame, where each row corresponds to a single neuron:

=====  ========  ====  ====  =================================================
index  duration  type    cf                                             spikes
=====  ========  ====  ====  =================================================
0          0.15   hsr  8000  [0.00243, 0.00414, 0.00715, 0.01089, 0.01358, ...
1          0.15   hsr  8000  [0.00325, 0.01234, 0.0203, 0.02295, 0.0268, 0....
2          0.15   hsr  8000  [0.00277, 0.00594, 0.01104, 0.01387, 0.0234, 0...
3          0.15   hsr  8000  [0.00311, 0.00563, 0.00971, 0.0133, 0.0177, 0....
4          0.15   hsr  8000  [0.00283, 0.00469, 0.00929, 0.01099, 0.01779, ...
5          0.15   hsr  8000  [0.00352, 0.00781, 0.01138, 0.02166, 0.02575, ...
6          0.15   hsr  8000  [0.00395, 0.00651, 0.00984, 0.0157, 0.02209, 0...
7          0.15   hsr  8000  [0.00385, 0.009, 0.01537, 0.02114, 0.02377, 0....
=====  ========  ====  ====  =================================================

The column 'spikes' is the most important and stores an array with
spike times (time stamps) in seconds for every action potential.  The
column 'duration' is the duration of the sound.  The column 'cf' is
the characteristic frequency (CF) of the fiber.  The column 'type'
tells us what auditory nerve fiber generated the spike train.  'hsr'
is for high-spontaneous rate fiber, 'msr' and 'lsr' for medium- and
low-spontaneous rate fibers.

Advantages of the format:

- easy addition of new meta data,
- efficient grouping and filtering of trains using _DataFrame
  functionality,
- export to MATLAB struct array through mat files::

    scipy.io.savemat(
        "spikes.mat",
        {'spike_trains': spike_trains.to_records()}
    )

Please, check thorns_ for more information and functions to manipulate
spike trains.


.. _DataFrame: http://pandas.pydata.org/pandas-docs/stable/generated/pandas.DataFrame.html
.. _pandas: http://pandas.pydata.org/
.. _thorns: https://github.com/mrkrd/thorns


Requirements
------------

- Python (2.7)
- Numpy
- Scipy
- Cython
- Pandas

- Matplotlib (optional, for examples)
- docopt (optional, for the command line scripts)
- thorns (optional, for examples and stats)
- matlab_wrapper (optional, for the MAP external model)


Note: On Windows you can install a Python distribution such as
Anaconda_ or `Python(x,y)`_ to fulfill most of the dependencies.


.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _`Python(x,y)`: https://code.google.com/p/pythonxy/


Installation
------------

Quick install::

  pip install cochlea


Note: Windows binary packages (wheels) are coming soon.  At the moment
you need `cython and a working compiler`_ to install the package.

.. _`cython and a working compiler`: https://github.com/cython/cython/wiki/64BitCythonExtensionsOnWindows



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


I would like to thank `Werner Hemmert`_ for supervising my PhD.

Thanks goes to Marcus Holmberg, who developed the traveling wave based
model.

I would also like to thank Muhammad S.A. Zilany, Ian C. Bruce and
Laurel H. Carney for developing inner ear models and allowing us to
use their code in *cochlea*.

I also grateful to Ray Meddis for support with the Matlab Auditory
Periphery model.

Please, when you use any of the models, always cite the original
publications describing the model.


.. _`Werner Hemmert`: http://www.imetum.tum.de/research/bai/home/?L=1
