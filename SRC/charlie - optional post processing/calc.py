#!/usr/bin/env python

"""
This module contains miscellaneous functions for time series data analysis.

"""

import numpy as np
from itertools import islice

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

def window(seq, n=2):    
    """
    Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    From itertools examples.
    
    args:
        seq (iterable): data sequence
    
    kwargs:
        n (int): width of window (default = 2)
    
    returns:
        iter: sliding window
    
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result    
    for e in it:
        result = result[1:] + (e, )
        yield result

def moving(data, func=np.mean, span=5):
    """
    Calculate the value of the function in a moving (sliding) window.
    
    args:
        data (iterable): data sequence (1D)
    
    kwargs:
        func: name of function (default = np.mean)
        
        span (int): width of window
    
    returns:
        list: output of func applied to data over the sliding window is a list 
        of the same size as the input data, where the ends are padded with the 
        first and last values of the input data
    
    """
    # The span should be odd to avoid shifting the data
    span += 1 if span % 2 == 0 else 0
    
    # First check that the span is smaller than the dataset
    if span < len(data):
        # Calculate average over sliding window, padding the front end
        # with the first value and the back end with the last value
        result = [func([data[i] for i in w])
                  for w in window(range(len(data)), span)]
        [result.insert(0, result[0]) for i in range(span / 2)]
        result.extend([result[-1]] * (len(data) - len(result)))
    else:
        result = [func(data), ] * len(data)
    return result

def fullmovingavg(data, span=5):
    """ Moving average, padded on ends with original data.
        
    args:
        data (iterable): data sequence (1D)
    
    kwargs:
        span (int): width of window
    
    returns:
        list: output of func applied to data over the sliding window is a list  
        of the same size as the input data, where instead of padding the ends 
        of the output, the average is calculated but over a smaller window 
        until reaching the ends of the sequence
    
    """
    span += 1 - span % 2
    weightings = np.repeat(1.0, span) / span
    middle = np.convolve(data, weightings, mode='valid')
    half_span = span / 2
    return np.hstack([data[:half_span], middle, data[-half_span:]])