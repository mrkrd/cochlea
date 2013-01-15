#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

def find_zero(func, kwargs, x1, x2, xtol):
    """We assume that the function is increasing and

    func(xrange[0]) < 0
    func(xrange[1]) > 0

    We perform binary search to find func's zero.

    """

    assert x1 < x2

    if func(x1, **kwargs) > 0:
        return np.nan

    if func(x2, **kwargs) < 0:
        return np.nan

    while x2 - x1 > xtol:
        x = (x1 + x2) / 2
        y = func(x, **kwargs)

        if y > 0:
            x2 = x
        elif y < 0:
            x1 = x
        else:
            x1 = x2 = x

    x = (x1 + x2) / 2

    return x


def main():
    pass


if __name__ == "__main__":
    main()
