#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""External models.

"""

from __future__ import division, absolute_import, print_function

__author__ = "Marek Rudnicki"


import warnings

try:
    from . matlab_auditory_periphery import run_matlab_auditory_periphery
except ImportError:
    warnings.warn("Error importing run_matlab_auditory_periphery().  Probably matlab_wrapper not installed.")
