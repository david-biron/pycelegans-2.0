#!/usr/bin/env python

"""
This module contains functions for reading and writing PyCelegans output data.

"""

import itertools
import os

from scipy import io

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

def readtxt(path_name, parameter, frame_range=(0, float("inf")), 
            data_type=float):
    """
    Read data from text file, e.g., .txt output of collectoutput.py.
    The expected input format is that each line contains data separated by 
    commas, with the first value corresponding to the frame number.
    The output is a dictionary with keys corresponding to frame numbers
    and values corresponding to the data contained in each line.
    
    args:
        path_name (str): name of directory containing properties .txt files
        
        parameter (str): name of parameter (file name minus extension)
    
    kwargs:
        frame_range (tuple): range of frames to read in given as [first, last)
        will read in all data from file)
        
        data_type: expected data type
    
    returns:
        dict: dict of values; keys are frame numbers
    
    """
    data_dict = {}
    with open(os.path.join(path_name, parameter + '.txt'), 'r') as f:
        for line in f:
            # Split comma-separated data from each line
            strs = [s for s in line.strip().split(',') if s]
            if len(strs) < 1:
                continue
            
            # Split the strings apart into separate values
            xs = (x.split() for x in strs if x)
            x0 = xs.next()
            
            # First value is the frame number
            k = int(x0.pop(0))
            if k < frame_range[0]:
                continue
            elif k >= frame_range[1]:
                break
            
            # Remaining values are the data
            v = []
            for ys in iter(x for x in itertools.chain([x0], xs) if x):
                u = []
                for y in ys:
                    try:
                        z = data_type(y)    # By default, assume numeric value
                    except ValueError:
                        z = y               # If error, return original string
                    u.append(z)
                u = u if len(u) > 1 else u[0]
                v.append(u)
            v = v if len(v) > 1 else v[0]
            
            data_dict[k] = v
    return data_dict

def writetxt(path_name, parameter, data_dict, frames=[], data_format="%0.1f"):
    """
    Write data to text file. The exact inverse of the read function.
    The inputs are the parameter name (while determines the file name)
    and the data dictionary, and a string describing how to convert
    the data values to text (e.g., %0.4f for float, %r for boolean).
    
    args:
        path_name (str): name of directory containing properties .txt files
        
        parameter (str): name of parameter (file name minus extension)
        
        data_dict (dict): dict of values; keys are frame numbers
        
    kwargs:
        frames (list): list of frames to write to file. If the data_dict does 
        not have a key corresponding to each frame in this list then an 
        error (-1) will be written. If this input is empty, all values 
        from data_dict will be written
        
        data_format: how to format each value to string
    
    """
    os.path.isdir(path_name) or os.mkdir(path_name)
    frames = frames or data_dict.keys()
    with open(os.path.join(path_name, parameter + '.txt'), 'w') as f:
        for i in sorted(frames):
            # Format and write each line to CSV text file
            line = "%d" % i
            if data_dict.has_key(i):
                try:
                    v = data_dict[i]
                    line += ',\t' + data_format % v
                except:
                    try:
                        for x in v:
                            line += ',\t' + data_format % x
                    except:
                        try:
                            for x in v:
                                for y in x:
                                    line += '\t' + data_format % y
                                line += ','
                        except:
                            line += ',\t-1'
            else:
                line += ',\t-1'
            line += '\n'
            f.write(line)

def reformat(data_dict):
    """
    Reformat data dictionary for saving to Matlab. Convert all dictionaries
    to lists, and construct a new value called frames that corresponds to
    the full set of keys for each variable.
    
    args:
        data_dict: dict of dicts
    
    returns:
        dict: dict of lists
    
    """
    
    # First find unique, sorted list of frames
    frames = set()
    for k, v in data_dict.iteritems():
        if isinstance(v, dict):
            if not frames:
                frames = set(data_dict[k].keys())
            else:
                frames &= set(data_dict[k].keys())
    frames = sorted(frames)
    
    # Get values corresponding to each frame, in order
    data_dict2 = {'frames' : [float(v) for v in frames]}
    for k, v in data_dict.iteritems():
        if isinstance(v, dict):
            data_dict2[k] = [v[i] for i in frames]
        else:
            data_dict2[k] = v
    return data_dict2

def savemat(file_name, data_dict):
    """
    Reformat so that dictionaries become lists and save to .mat file. 
    Input is dictionary of values, i.e., {'midline':..., 'nosetip':..., }
    
    args:
        file_name (str): name of .mat file to save data
        
        data_dict (dict): dict of dicts
    
    """
    io.savemat(file_name, reformat(data_dict), oned_as='column')