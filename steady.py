#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 15:38:29 2018

@author: robert
"""

import numpy as np
from scipy.optimize import newton


def get_roll_sensitivity(W_s, h_2, K_f, K_r):
    """ Calculate the roll sensitivity of the car. """
    phi = -W_s*h_2 / (K_f+K_r-W_s*h_2)
    return phi


def get_front_lat_weight_transfer(W_s)
