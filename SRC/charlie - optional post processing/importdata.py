#!/usr/bin/env python

"""
Read output of PyCelegans used for further analysis: convert each midline to 
a numpy array and calculate the distance along the midline from nose to tail.
   
"""

import os
import sys

import numpy as np
from numpy.linalg import norm

from readwrite import readtxt, writetxt
import parseargs

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

FRAME_RANGE = [0, float("inf")]
VERSION = 0

def readnosetail(path_name, frame_range=FRAME_RANGE, version=VERSION):
    """
    Read worm midline (assuming nose-to-tail), converting to numpy array.
    
    args:
        path_name (str): name of directory containing properties .txt files
    
    kwargs:
        frame_range (tuple): range of frames to read in given as [first, last)
        
        version (int): version of PyCelegans process.py code output to use (0: 
        Nick Labello's midline, 1: Marc Goessling's midline)
    
    returns:
        nosetail (dict): dict of numpy arrays of coordinates along midline
    
    """
    
    # Read in the specified midline
    if version == 0:
        parameter = 'midline'
    elif version == 1:
        parameter = 'marcmidline'
    results = readtxt(path_name, parameter, frame_range)
    
    # Convert each datapoint to an array and remove error frames
    results = {k : np.asarray(v) for k, v in results.iteritems()}
    results = {k : v for k, v in results.iteritems() if not np.any(v < 0)}
    if version == 1:
        # By default, Marc's midline runs from tail to nose
        results = {k : v[::-1] for k, v in results.iteritems()}
    return results

def getnosedist(nosetail):
    """
    Calculate distance of each point along worm from nosetip.
    
    args:
        nosetail (dict): dict of numpy arrays of coordinates along midline
    
    returns:
        nosedist (dict): dict of numpy arrays of distance from nose to tail
        
        nosetail (dict): same as input, but with same set of keys as nosedist
    
    """
    nosedist = {}
    for k, v in nosetail.iteritems():
        try:
            # Cumulative sum of difference between successive points
            d = [norm(x) for x in (v[1:] - v[:-1])]
            d.insert(0, 0.0)
            nosedist[k] = np.cumsum(d)
        except:
            pass
    
    # Both nosedist and nosetail should have identical sets of frames
    return nosedist, {k : nosetail[k] for k in nosedist.keys()}

def main(path_name, frame_range=FRAME_RANGE, version=VERSION):
    """
    Import midline and calculate distance from nose to tail.
    
    args:
        path_name (str): name of directory containing properties .txt files
    
    kwargs:
        frame_range (tuple): range of frames to read in given as [first, last)
        
        version (int): version of PyCelegans process.py code output to use (0: 
        Nick Labello's midline, 1: Marc Goessling's midline)
    
    returns:
        nosedist (dict): dict of numpy arrays of distance from nose to tail
        
        nosetail (dict): same as input, but with same set of keys as nosedist
    
    """
    nosetail = readnosetail(path_name, frame_range, version)
    nosedist, nosetail = getnosedist(nosetail)
    return nosedist, nosetail