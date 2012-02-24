#!/usr/bin/env python

from __future__ import division

__author__ = "Marek Rudnicki"

import numpy as np

def find_threshold(func, args=(), init_range=(0, 10), desired_range=1):
    a, b = init_range

    assert a < b
    assert func(a, *args) < 0
    assert func(b, *args) > 0


    while b - a > desired_range:
        x = (a + b) / 2
        y = func(x, *args)

        if y > 0:
            b = x
        elif y < 0:
            a = x
        else:
            a = b = x

    return (a + b) / 2


def main():
    pass

if __name__ == "__main__":
    main()
