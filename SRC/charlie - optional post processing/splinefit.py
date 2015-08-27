#!/usr/bin/env python

"""
Construct a uniformly-sampled B-spline of the worm, given the set of points running along the worm midline from nosetip to tailtip and the distance between each of those points.

"""

import argparse
import os
import sys

import numpy as np
from scipy.interpolate import splprep, splev

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

SMOOTHING = 25
NUM_POINTS = 100
SPACING = 0.
LENGTH = 0.

def makebspline(us, vs, s=SMOOTHING, per=0):
    """ Fit the a series of (x, y) coordinates to a smoothing B-spline.
    
    args:
        us (array): euclidean distance between points
        
        vs (array): array of (x, y) points
    
    kwargs:
        s (int): smoothing
        
        per (int): set non-zero for a closed spline, i.e., periodic boundary
        conditions (default = 0)
    
    returns:
        dict: dict of spline fit (knots, coefs, degree)
    
    """
    
    tcks = {}
    for k, v in vs.iteritems():
        try:
            u = us[k]
            tcku = splprep(zip(*v), u=u, s=s, per=per, quiet=1)
            tcks[k] = tcku[0]
        except:
            pass
    return tcks

def evalbspline(tcks, num_points=NUM_POINTS, 
                 spacing=SPACING, lengths=None):
    """
    Evaluate the spline at discrete points.    
    There are three options for constructing the spline: 
    1) Same number of points, variable lengths, variable spacing 
    => must specify num_points (default) and each length 
    2) Variable number of points, variable lengths, same spacing 
    => must specify spacing and each length 
    3) Same number of points, same lengths, same spacing 
    => must specify all same length and either spacing or num_points
    
    args:
        tcks (dict): dict of spline fit (knots, coefs, degree)
    
    kwargs:
        num_points (int): number of points in final midline
        
        spacing (int): spacing between points
        
        lengths (dict): dict of desired lengths of each midline
    
    returns:
        nosedist (dict): dict of uniformly-sampled distances from nose to tail
        
        nosetail (dict): dict of midlines evaluated using polynomial fit at 
        each point in the output nosedist
    
    """
    us, vs, ws = {}, {}, {}
    for k, tck in tcks.iteritems():
        try:
            # Construct spline from uniformly sampled parameterization
            start, stop = 0., lengths[k]
            if spacing > 0: # equal spacing between points
                u = np.arange(start, stop, spacing)
            else:           # uniform number of points
                u = np.linspace(start, stop, num_points)
            
            # Evaluate the spline and its first derivatives
            v = splev(u, tck, der=0)
            w = splev(u, tck, der=1)
            
            us[k] = u
            vs[k] = v
            ws[k] = w
        except:
            pass
    return us, vs, ws

def main(us, vs, smoothing=SMOOTHING, num_points=NUM_POINTS, 
         spacing=SPACING, length=LENGTH, per=0):
    """ Fit the worm to a smoothing B-spline.
    
    args:
        us (array): euclidean distance between points
        
        vs (array): array of (x, y) points
    
    kwargs:
        smoothing (int): smoothing
        
        num_points (int): number of points in final midline
        
        spacing (int): spacing between points
        
        lengths (dict): dict of desired lengths of each midline
        
        per (int): set non-zero for a closed spline, i.e., periodic b.c.
    
    returns:
        nosedist (dict): dict of uniformly-sampled distances from nose to tail
        
        nosetail (dict): dict of midlines evaluated using polynomial fit at 
        each point in the output nosedist
    
    """
    tcks = makebspline(us, vs, smoothing, per)
    if length > 0:
        lengths = {k : length for k, u in us.iteritems()}
    else:
        lengths = {k : max(u) for k, u in us.iteritems()}
    return evalbspline(tcks, num_points, spacing, lengths)