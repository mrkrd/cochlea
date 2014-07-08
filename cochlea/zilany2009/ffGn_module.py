"""Copyright 2003-2005 B. Scott Jackson
Copyright 2008 M. S. A. Zilany
Copyright 2009-2014 Marek Rudnicki

This file is part of cochlea.

cochlea is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

cochlea is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with cochlea.  If not, see <http://www.gnu.org/licenses/>.


Notes
=====

This implementation is based on the oryginal MATLAB implementation.

TODO: Unit test against the oryginal.

"""

from __future__ import division, print_function, absolute_import

__author__ = "Marek Rudnicki"


import numpy as np
from numpy.random import randn
from scipy.signal import resample
from numpy.fft import fft, ifft

# sigma has been removed as it is never read
def ffGn(N, tdres, Hinput, mu):


    assert (N > 0)
    assert (tdres < 1)
    assert (Hinput >= 0) and (Hinput <= 2)


    # Downsampling No. of points to match with those of Scott jackson (tau 1e-1)
    resamp = np.ceil(1e-1 / tdres)
    nop = N
    N = np.ceil(N / resamp) + 1
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
        y = randn(N)
    else:
        # TODO: make variables persistant
        Nfft = 2 ** np.ceil(np.log2(2*(N-1)))
        NfftHalf = np.round(Nfft / 2)

        k = np.concatenate( (np.arange(0,NfftHalf), np.arange(NfftHalf,0,-1)) )
        Zmag = 0.5 * ( (k+1)**(2*H) -2*k**(2*H) + np.abs(k-1)**(2*H) )

        Zmag = np.real(fft(Zmag))
        assert np.all(Zmag >= 0)

        Zmag = np.sqrt(Zmag)

        Z = Zmag * (randn(Nfft) + 1j*randn(Nfft))

        y = np.real(ifft(Z)) * np.sqrt(Nfft)

        y = y[0:N]

        # Convert the fGn to fBn, if necessary.
        if fBn == 1:
            y = np.cumsum(y)


        # Resampling to match with the AN model
        y = resample(y, resamp*len(y))


        if mu < 0.5:
            sigma = 5
        elif mu < 18:
            sigma = 50          # 7 when added after powerlaw
        else:
            sigma = 200         # 40 when added after powerlaw


        y = y*sigma

        return y[0:nop]


def main():
    y = ffGn(10, 1e-1, 0.2, 1, debug=True)
    print(y)


if __name__ == "__main__":
    main()
