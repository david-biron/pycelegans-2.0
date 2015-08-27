#!/usr/bin/env python

"""
This module checks the output of PyCelegans by calculating the average length of the worm midline over a sliding window, and rejecting frames whose length deviates significantly from this value.
 
"""

import argparse
import os
import sys

import numpy as np

from calc import moving

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

SPAN = 10*60*10
MAX_VAL = 0.05
MAX_DEV = 0.

def checklength(nosedist, span=SPAN, max_dev=MAX_DEV, max_val=MAX_VAL):
    """
    Remove frames whose lengths do not fit the specified criteria: 
    1. Must fall within the specified number of standard deviations from the mean, calculated over a moving window 
    => set max_dev but not max_val
    2. Must fall within a constant value of the mean, calculated over a moving window
    => set max_val but not max_dev
    
    args:
        nosedist (dict): dict of numpy arrays of distance from nose to tail
    
    kwargs:
        span (int): width of window
        
        max_dev (float): maximum allowed difference between Z-score of length
        
        max_val (float): maximum allowed difference between length and average
    
    returns:
        dict: dict of booleans specifying frames that pass or fail the test
    
    """
    frames = sorted(nosedist.keys())
    span = min(len(nosedist), span)
    
    # Return the length of the worm for each frame in the movie
    vals = np.asarray([nosedist[i].max() for i in frames])
    
    # Calculate the average length over the specified window
    avgs = np.asarray(moving(vals, func=np.mean, span=span))
    
    # Find frames whose lengths do not meet the specified criteria
    if max_dev <= 0:
        cmps = dict(zip(frames, (vals - avgs) / avgs))
        func = lambda x: np.abs(x) > max_val
    else:
        stds = np.asarray(moving(vals, func=np.std, span=span))
        cmps = dict(zip(frames, (vals - avgs) / stds))
        func = lambda x: np.abs(x) > max_dev
    
    # Return dictionary with boolean stating where length fails test
    return {k : func(v) for k, v in cmps.iteritems()}

def main(nosedist, span=SPAN, max_dev=MAX_DEV, max_val=MAX_VAL):
    """
    Find frames where the length fails to meet specified criteria
    
    args:
        nosedist (dict): dict of numpy arrays of distance from nose to tail
    
    kwargs:
        span (int): width of window
        
        max_dev (float): maximum allowed difference between Z-score of length
        
        max_val (float): maximum allowed difference between length and average
    
    returns:
        dict: dict of booleans specifying frames that pass or fail the test
    
    """
    return checklength(nosedist, span, max_dev, max_val)