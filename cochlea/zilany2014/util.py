"""This ffGn implementation based on the oryginal MATLAB code.

TODO: Unit test against the oryginal.

"""

from __future__ import division, print_function, absolute_import

import numpy as np
from numpy.random import randn
from scipy.signal import resample
from numpy.fft import fft, ifft


def ffGn(N, tdres, Hinput, noiseType, mu, sigma=1, random_debug=None):

    assert (N > 0)
    assert (tdres < 1)
    assert (Hinput >= 0) and (Hinput <= 2)


    # Here we change the meaning of `noiseType', if it's 0, then we
    # return no noise at all.  If necessary, the seed can be set
    # outside by calling np.radnom.seed().
    if noiseType == 0:
        return np.zeros(N)


    # Downsampling No. of points to match with those of Scott jackson (tau 1e-1)
    resamp = int(np.ceil(1e-1 / tdres))
    nop = N
    N = int(np.ceil(N / resamp) + 1)
    if N < 10:
        N = 10

    # Determine whether fGn or fBn should be produced.
    if Hinput <= 1:
        H = Hinput
        fBn = 0
    else:
        H = Hinput - 1
        fBn = 1



    # Calculate the fGn.
    if H == 0.5:
        # If H=0.5, then fGn is equivalent to white Gaussian noise.
        if random_debug is None:
            y = randn(N)
        else:
            y = random_debug
    else:
        # TODO: make variables persistant
        Nfft = int(2 ** np.ceil(np.log2(2*(N-1))))
        NfftHalf = np.round(Nfft / 2)

        k = np.concatenate( (np.arange(0,NfftHalf), np.arange(NfftHalf,0,-1)) )
        Zmag = 0.5 * ( (k+1)**(2*H) -2*k**(2*H) + np.abs(k-1)**(2*H) )

        Zmag = np.real(fft(Zmag))
        assert np.all(Zmag >= 0)

        Zmag = np.sqrt(Zmag)

        if random_debug is None:
            Z = Zmag * (randn(Nfft) + 1j*randn(Nfft))
        else:
            Z = Zmag * (random_debug + 1j*random_debug)

        y = np.real(ifft(Z)) * np.sqrt(Nfft)

        y = y[0:N]

        # Convert the fGn to fBn, if necessary.
        if fBn == 1:
            y = np.cumsum(y)


        # Resampling to match with the AN model
        y = resample(y, resamp*len(y))


        if mu < 0.5:
            sigma = 3
        elif mu < 18:
            sigma = 30          # 7 when added after powerlaw
        else:
            sigma = 200         # 40 when added after powerlaw


        y = y*sigma

        return y[0:nop]



def calc_cfs(cf, species):
    if np.isscalar(cf):
        cfs = [float(cf)]

    elif isinstance(cf, tuple) and ('cat' in species):
        # Based on GenerateGreenwood_CFList() from DSAM
        # Liberman (1982)
        aA = 456
        k = 0.8
        a = 2.1

        freq_min, freq_max, freq_num = cf

        xmin = np.log10(freq_min / aA + k) / a
        xmax = np.log10(freq_max / aA + k) / a

        x_map = np.linspace(xmin, xmax, freq_num)
        cfs = aA * ( 10**( a*x_map ) - k)

    elif isinstance(cf, tuple) and ('human' in species):
        # Based on GenerateGreenwood_CFList() from DSAM
        # Liberman (1982)
        aA = 165.4
        k = 0.88
        a = 2.1

        freq_min, freq_max, freq_num = cf

        xmin = np.log10(freq_min / aA + k) / a
        xmax = np.log10(freq_max / aA + k) / a

        x_map = np.linspace(xmin, xmax, freq_num)
        cfs = aA * ( 10**( a*x_map ) - k)

    elif isinstance(cf, list) or isinstance(cf, np.ndarray):
        cfs = cf

    else:
        raise RuntimeError("CF must be a scalar, a tuple or a list.")

    return cfs
