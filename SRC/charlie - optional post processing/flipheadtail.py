#!/usr/bin/env python

"""
This module checks the output of PyCelegans by determining the frames where the nose has been incorrectly identified. Each worm midline is fit to a 
polynomial, and the coefficients from each frame are compared to the
average coefficients (calculated over a sliding window) to find frames
where the posture deviates significantly from the local value. This process 
is performed iteratively to account for outliers, but still relies upon the 
assumption that the majority of the frames within the window are correct. However, if the span is just one frame, only calculate the distance between points along the worm between consecutive frames, identifying flips where this value is smaller where one of the frames has been flipped.
  
"""

import argparse
import os
import sys

import numpy as np
from numpy.polynomial.polynomial import polyfit, polyval

from calc import *
import splinefit

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

DEGREE = 10
NUM_ITERS = 3
SPAN = 10*2
MIN_DEV = 2.0
MAX_DEV = 3.0
SMOOTHING = splinefit.SMOOTHING
NUM_POINTS = splinefit.NUM_POINTS

def makepolyfit(nosedist, nosetail, degree=DEGREE):
    """ Fit the midline to a polynomial.
    
    args:
        nosedist (dict): dict of numpy arrays of distance from nose to tail
        
        nosetail (dict): dict of numpy arrays of coordinates along midline
    
    kwargs:
        degree (int): degree of polynomial fit
    
    returns:
        dict: dict of coefficients, each of which is m X 2 array for m degrees
    
    """
    coefs = {}
    for k, v in nosetail.iteritems():
        u = nosedist[k]
        try:
            # First, smooth the data
            x, y = zip(*v)
            x = fullmovingavg(x, span=5)
            y = fullmovingavg(y, span=5)
            
            # Then fit to requested type of polynomial
            cx = polyfit(u, x, degree)
            cy = polyfit(u, y, degree)
            
            # Return coefficients, converted to regular polynomial basis
            coefs[k] = (cx, cy)
        except:
            pass
    return coefs

def evalpolyfit(coefs, num_points=100, spacing=None, lengths=None):
    """
    Evaluate the polynomial fit at discrete points.
    There are four options for constructing the spline:
    1) Same number of points, variable lengths, variable spacing
    => must specify only num_points (default)
    2) Variable number of points, variable lengths, same spacing
    => must specify only spacing
    3) Same number of points, same lengths, same spacing
    => must specify length and either spacing or num_points
    
    args:
        coefs (dict): dict of coefficients (output of makepolyfit)
    
    kwargs:
        num_points (int): number of points in final midline
        
        spacing (int): spacing between points
        
        lengths (dict): dict of desired lengths of each midline
    
    returns:
        nosedist (dict): dict of uniformly-sampled distances from nose to tail
        
        nosetail (dict): dict of midlines evaluated using polynomial fit at 
        each point in the output nosedist
    
    """
    nosedist, nosetail = {}, {}
    for k, (cx, cy) in coefs.iteritems():
        try:
            # Determine spacing, number of points and length based on inputs
            start, stop = 0., lengths[i]
            if spacing:
                uu = np.arange(start, stop, spacing)
            else:
                uu = np.linspace(start, stop, num_points)
            
            # Construct polynomial from its coefficients
            vx = polyval(uu, cx)
            vy = polyval(uu, cy)
            vv = np.asarray(zip(vx, vy))
            
            nosedist[k] = uu
            nosetail[k] = vv
        except:
            pass
    return nosedist, nosetail

def checkflips(coefs_for, coefs_rev, span=SPAN, 
                min_dev=MIN_DEV, max_dev=MAX_DEV):
    """
    Check for nose-tail flips by comparing coefficients corresponding to worm
    oriented in forward and reverse direction. When a value exceeeds the 
    maximum specified deviation from the local window, flip the coefficients.
    
    args:
        coefs_for (dict): dict of coefficients corresponding to nose-to-tail
        
        coefs_rev (dict): dict of coefficients corresponding to tail-to-nose
    
    kwargs:
        span (int): number of points in final midline
        
        min_dev (float): spacing between points
        
        max_dev (float): dict of desired lengths of each midline
    
    returns:
        coefs_for (dict): adjusted values of coefs_for
        
        coefs_rev (dict): adjusted values of coefs_rev
        
        is_flipped (dict): dict of strings with values of 
        'False': no flip (value < min_dev) 
        'True': yes flip (value > max_dev) 
        'Indeterminate': could not determine (min_dev < value < max_dev)
    
    """
    frames = sorted(set(coefs_for.keys()) & set(coefs_rev.keys()))
    num_frames = len(frames)
    span = min(num_frames, span)
    vals_for = np.asarray([coefs_for[i] for i in frames])
    
    # Calculate the deviations of the coefficients over the specified window
    num_dims, num_coefs = np.shape(coefs_for[frames[0]])
    num_coefs = min(num_coefs, 10)
    all_for = np.zeros((num_frames, num_dims, num_coefs))
    for d in range(num_dims):
        for c in range(num_coefs):
            v = vals_for[:,d,c]
            a = np.asarray(moving(v, func=np.mean, span=span))
            all_for[:,d,c] = np.abs(v - a)
    
    # Add together coefficients for each frame
    sum_for = np.zeros(num_frames)
    for i in range(num_frames):
        x = [np.mean(all_for[i,d,:]) ** 2 for d in range(num_dims)]
        sum_for[i] = np.sqrt(np.sum(x))
    
    # Flip all coefficients with large z-scores
    is_flipped = dict.fromkeys(frames, 'False')
    calc_z = lambda x : np.abs(np.mean(x) - x) / np.std(x)
    z_for = calc_z(sum_for)
    for i, v in enumerate(z_for):
        j = frames[i]
        if v > max_dev:
            coefs_for[j], coefs_rev[j] = coefs_rev[j], coefs_for[j]
            is_flipped[j] = 'True'
        elif v >= min_dev and v <= max_dev:
            is_flipped[j] = 'Indeterminate'
    return coefs_for, coefs_rev, is_flipped

def main(nosedist, nosetail, degree=DEGREE, num_iters=NUM_ITERS, 
         span=SPAN, min_dev=MIN_DEV, max_dev=MAX_DEV, 
         smoothing=SMOOTHING, num_points=NUM_POINTS):
    """
    Find frames where head/tail identification has failed by making a 
    polynomial representation of the worm and comparing the value of the 
    coefficients frame-by-frame to the average in a sliding window. If the size of the sliding window is set to just one frame however, do not fit to a polynomial but instead calculate the Euclidean distance between points along splined worms in consecutive frames, counting a flip if the norm between 'forward' and 'reverse' is smaller than between both 'forward'.
    
    args:
        nosedist (dict): dict of numpy arrays of distance from nose to tail
        
        nosetail (dict): dict of numpy arrays of coordinates along midline
    
    kwargs:
        degree (int): degree of polynomial fit
        
        num_iters (int): number of times to repeat the head/tail flip check
        
        span (int): number of points in final midline
        
        min_dev (float): spacing between points
        
        max_dev (float): dict of desired lengths of each midline
    
    returns:
        is_flipped (dict): dict of strings with values of 
        'False': no flip (value < min_dev)
        'True': yes flip (value > max_dev)
        'Indeterminate': could not determine (min_dev < value < max_dev)
    
    """
    
    if span > 1:
        # Make the polynomial fit in the forward and reverse directions
        coefs_for = makepolyfit(nosedist, nosetail, degree=degree)
        taildist = dict.fromkeys(nosedist.keys())
        for k, v in nosedist.iteritems():
            taildist[k] = (max(v) - v)[::-1]
        coefs_rev = makepolyfit(taildist, nosetail, degree=degree)
                
        # Check iteratively for frames that have been flipped (try to correct 
        # for large neighborhoods where nose has been incorrectly identified)
        for i in range(num_iters):
            coefs_for, coefs_rev = checkflips(coefs_for, coefs_rev, span, 
                                               min_dev, max_dev)[:2]
        is_flipped = checkflips(coefs_for, coefs_rev, 
                                span, min_dev, max_dev)[-1]
    else:
        # Create a splined worm with an even number of points for each frame
        nosetail = splinefit.main(nosedist, nosetail, 
                                  smoothing, num_points)[1]
        
        # Just check by comparing the current worm to the flipped last worm
        frames = sorted(nosetail.keys())
        is_flipped = dict.fromkeys(frames, 'False')
        
        x1, y1 = nosetail[frames[0]]
        for i in frames[1:]:
            x2, y2 = nosetail[i]
            
            # Calculate the Euclidean distance between corresponding points
            f = np.sum(np.sqrt((x1 - x2[::+1]) ** 2 + (y1 - y2[::+1]) ** 2))
            r = np.sum(np.sqrt((x1 - x2[::-1]) ** 2 + (y1 - y2[::-1]) ** 2))
            
            # Flip if the norm between current forward worm and previous
            # reverse worm is smaller than between both forward worms
            if f > r:
                is_flipped[i] = 'True'
                x1, y1 = x2[::-1], y2[::-1]
            else:
                x1, y1 = x2[::+1], y2[::+1]
    
    return is_flipped