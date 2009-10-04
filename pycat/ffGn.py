from __future__ import division

import numpy as np
from numpy.random import randn
from scipy.signal import resample
from numpy.fft import fft, ifft

def ffGn(N, Hinput, mu, sigma):

    # TODO: assert proper input

    # Downsampling No. of points to match with those of Scott jackson (tau 1e-1)
    tdres = 1e-4
    nop = N
    N = np.ceil(N * tdres / 1e-1)
    if N < 10:
        N = 10

    # Determine whether fGn or fBn should be produced.
    if Hinput <= 1:
        H = Hinput
        fBn = 0
    else:
        H = Hinput - 1
        fBn = 1


    # TODO: check randn('state')

    # Calculate the fGn.
    if H == 0.5:
        # If H=0.5, then fGn is equivalent to white Gaussian noise.
        y = randn(N)
    else:
        # TODO: make variables persistant
        Nfft = 2 ** np.ceil( np.log2(2*(N-1)) )
        NfftHalf = np.round(Nfft / 2)

        k = np.concatenate( (np.arange(0,NfftHalf), np.arange(NfftHalf,0,-1)) )
        Zmag = 0.5 * ( (k+1)**(2*H) -2*k**(2*H) + np.abs(k-1)**(2*H) )

        Zmag = np.real(fft(Zmag))
        assert np.all(Zmag >= 0)

        Zmag = np.sqrt(Zmag)


        Z = Zmag * ( randn(Nfft) + 1j*randn(Nfft) )

        y = np.real(ifft(Z)) * np.sqrt(Nfft)

        y = y[0:N]

        # Convert the fGn to fBn, if necessary.
        if fBn == 1:
            y = np.cumsum(y)


        # Resampling to match with the AN model
        y = resample(y, 1e-1/tdres)


        if mu < 0.5:
            sigma = 5
        elif mu < 18:
            sigma = 50          # 7 when added after powerlaw
        else:
            sigma = 200         # 40 when added after powerlaw


        y = y*sigma

        return y[0:nop]




def main():
    y = ffGn(10, 0.2, 1, 1)
    print y


if __name__ == "__main__":
    main()
