# Author: Marek Rudnicki
# Time-stamp: <2009-10-07 18:27:24 marek>
#
# Description: Model of auditory periphery of: Zilany, M.S.A., Bruce,
# I.C., Nelson, P.C., and Carney, L.H. (manuscript in preparation) 2009

from __future__ import division

import numpy as np
import matplotlib.pyplot as plt


class Carney2009(object):
    def __init__(self, hsr=100, msr=100, lsr=100,
                 freq=1000, animal='cat', syn_implnt='approx'):
        """
        hsr, msr, lsr: number of HSR/MSR/LSR fibers

        freq: CF

        animal: must be cat

        syn_implnt: 'approx'/'acctual' implementation of the power-law
        """
        self.set_freq(freq)


    def set_freq(self, freq):
        if isinstance(freq, int):
            freq = float(freq)
        assert (isinstance(freq, tuple) or
                isinstance(freq, float))

        if isinstance(freq, float):
            self._freq_map = [freq]
        elif isinstance(freq, tuple):
            # Based on GenerateGreenwood_CFList() from DSAM
            aA = 456.0
            k = 0.8
            a = 2.1

            freq_min, freq_max, freq_num = freq

            xmin = np.log10( freq_min / aA + k) / a
            xmax = np.log10( freq_max / aA + k) / a

            x_map = np.linspace(xmin, xmax, freq_num)
            self._freq_map = aA * ( 10**( a*x_map ) - k)
        else:
            assert False


def main():
    ear = Carney2009()

    ear.set_freq( (1000, 5000, 5) )

    print ear._freq_map


if __name__ == "__main__":
    main()

