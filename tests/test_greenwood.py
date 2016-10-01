#!/usr/bin/env python

from __future__ import division, print_function, absolute_import

import pytest

import numpy as np
from numpy.testing import assert_almost_equal, assert_equal
import cochlea


def test_invertable():
    """Test if the greenwood function is invertable.

    """
    freq = np.linspace(20, 20e3, 100)
    x = cochlea.greenwood_inverse(freq, 'human')
    freq_greenwood = cochlea.greenwood(x, 'human')

    assert_almost_equal(freq, freq_greenwood)


def test_function():
    """Compare function results with correct value.

    """
    freq_0 = cochlea.greenwood(0, 'human')

    freq_0_target = 165.4 * (1 - 0.88)

    assert_equal(freq_0, freq_0_target)


def test_parameters():
    """Test if providing parameters results in the same as providing species.

    """
    x = np.linspace(0, 35e-3, 100)

    freq = cochlea.greenwood(x, A=165.4, a=60, k=0.88)
    freq_human = cochlea.greenwood(x, 'human')

    x = cochlea.greenwood_inverse(freq, A=165.4, a=60, k=0.88)
    x_human = cochlea.greenwood_inverse(freq, 'human')

    assert_equal(freq, freq_human)
    assert_equal(x, x_human)


def test_too_long():
    with pytest.raises(ValueError):
        cochlea.greenwood(40e-3, 'human')
