#!/usr/bin/env python

"""
This is the main script for post-processing. It can run in two modes:
    1. check: run various checks on the output of PyCelegans
    2. spline: spline the worm midline and export to Matlab

"""

import os
import sys
import time

import numpy as np
from scipy import integrate

import parseargs
import importdata
import checklength
import flipheadtail
import splinefit

from readwrite import *

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

def main(inputs):
    """
    Main function for accessing post-process functions.
    
    args:
        inputs: command line inputs
    
    """
    a = parseargs.main(sys.argv[1:])
    a.properties = os.path.abspath(a.properties)

    # Read in midline and parameterize it
    importdata_inputs = (a.properties, a.frame_range, a.version)
    nosedist, nosetail = importdata.main(*importdata_inputs)
    
    # Get full set of frames (without error frames removed)
    is_loop = readtxt(a.properties, "is_loop", a.frame_range, str)
    all_frames = is_loop.keys()
    
    if a.action == 'check':
        # Remove looped frames (if using Nick's midline)
        if a.version == 0:
            for k, v in is_loop.iteritems():
                if nosedist.has_key(k) and v == 'True':
                    del nosedist[k], nosetail[k]
        
        # Find frames with bad lengths and write to .txt file
        checklength_inputs = (nosedist, a.len_span, 
                              a.len_max_dev, a.len_max_val)
        bad_length = checklength.main(*checklength_inputs)
        writetxt(a.properties, "bad_length", bad_length, all_frames, "%r")
        
        # Remove bad length frames before continuing
        for k, v in bad_length.iteritems():
            if v:
                del nosedist[k], nosetail[k]
        
        # Find frames that should be flipped and write to .txt file
        flipheadtail_inputs = (nosedist, nosetail, 
                               a.flip_degree, a.flip_iters, a.flip_span, 
                               a.flip_min_dev, a.flip_max_dev)
        is_flipped = flipheadtail.main(*flipheadtail_inputs)
        writetxt(a.properties, "is_flipped", is_flipped, all_frames, "%s")
        
    elif a.action == 'spline':
        if a.output == "":
            a.output = os.path.join(os.path.dirname(a.properties), 'processed')
        os.path.isdir(a.output) or os.mkdir(a.output)
        
        if a.use_checks:
            # Read in all checks from file (just skip if missing files)
            old_frames = nosedist.keys()
            checks = dict.fromkeys(['is_loop', 'bad_length', 'is_flipped'])
            
            # Initialize all values to 'False'
            for k in checks.keys():
                try:
                    v = readtxt(a.properties, k, a.frame_range, str)
                except IOError:
                    v = dict.fromkeys(old_frames, 'False')
                checks[k] = v
            
            # Do not use is_loop for Marc's midline
            if a.version == 1:
                checks['is_loop'] = dict.fromkeys(old_frames, 'False')
            
            # Do not flip unless thus requested
            if not a.reorient:
                checks['is_flipped'] = dict.fromkeys(old_frames, 'False')
            
            for k in old_frames:
                try:
                    # Remove loops, bad lengths and unknown orientations
                    if (checks['is_loop'][k] == 'True' or 
                        checks['bad_length'][k] == 'True' or
                        checks['is_flipped'][k] == 'Indeterminate'):
                        del nosedist[k], nosetail[k]
                    # Flip frames with bad orientation
                    elif checks['is_flipped'][k] == 'True':
                        nosedist[k] = (max(nosedist[k]) - nosedist[k])[::-1]
                        nosetail[k] = nosetail[k][::-1]
                except KeyError:
                    # Also remove frames for which checks are unavailable
                    del nosedist[k], nosetail[k]
        
        # Spline the worm and write the .txt files to a new directory
        up_rate = 5
        splinefit_inputs = (nosedist, nosetail, up_rate*a.smoothing, 
                            up_rate*a.num_points, a.spacing, a.length)
        nosedist, nosetail, noseder = splinefit.main(*splinefit_inputs)
        
        # Calculate the length by integrating over the midline
        wormlength = {}
        for k, u in nosedist.iteritems():
            dx, dy = noseder[k]
            z = integrate.trapz(np.sqrt(dx ** 2 + dy ** 2), u)
            wormlength[k] = z
        
        # Reformat the splined midline before saving to file
        nosetail = {k : np.asarray(zip(*v))[::up_rate] 
                    for k, v in nosetail.iteritems()}
        
        writetxt(a.output, "wormlength", wormlength, all_frames, "%0.2f")
        writetxt(a.output, "nosetail", nosetail, all_frames, "%0.2f")
        
        # Create dictionary for saving to .mat file
        data_dict = {'nosetail' : nosetail, 'wormlength' : wormlength}
        
        if a.keep_sides:
            # Read information from side paths
            sideonepath = readtxt(a.properties, 'sideonepath', a.frame_range)
            sidetwopath = readtxt(a.properties, 'sidetwopath', a.frame_range)
            
            sidepath = dict.fromkeys(nosedist.keys())
            for k in sidepath.keys():
                v = sideonepath[k][:-1]
                v.extend(sidetwopath[k][::-1])
                if v[0] != v[-1]:
                    v.append(v[0])
                sidepath[k] = np.asarray(v)
            sidedist, sidepath = importdata.getnosedist(sidepath)
            
            splinefit_inputs = (sidedist, sidepath, a.smoothing, 
                                up_rate*2*a.num_points, up_rate*2*a.spacing, 
                                0., 1)
            sidedist, sidepath, sideder = splinefit.main(*splinefit_inputs)
            
            # Calculate the area by integrating around the perimeter
            wormarea = {}
            for k, u in sidedist.iteritems():
                x = sidepath[k][0]
                dy = sideder[k][1]
                z = np.abs(integrate.trapz(x * dy, u))
                wormarea[k] = z
            
            # Reformat the splined perimeter before saving to file
            sidepath = {k : np.asarray(zip(*v))[::up_rate] 
                        for k, v in sidepath.iteritems()}
            
            writetxt(a.output, "wormarea", wormarea, all_frames, "%0.2f")
            writetxt(a.output, "sidepath", sidepath, all_frames, "%0.2f")
            
            data_dict['wormarea'] = wormarea
            data_dict['sidepath'] = sidepath
        
        # Save to .mat file
        savemat(a.output + '.mat', data_dict)
    
    # Print error check on all files in the directory
    fs, ds = [], []
    for f in os.listdir(a.properties):
        if os.path.splitext(f)[1] == '.txt':
            fs.append(f)
            ds.append(a.properties)
    if a.action == 'spline':
        for f in os.listdir(a.output):
            if os.path.splitext(f)[1] == '.txt':
                fs.append(f)
                ds.append(a.output)
    
    w1 = max([len(f) for f in fs])
    if a.action == 'check':
        print 'Frames that passed checks:'
        for p, b in zip(['is_loop', 'bad_length', 'is_flipped'], 
                        [('False', ), ('False', ), ('True', 'False')]):
            try:
                imported = readtxt(a.properties, p, a.frame_range, str)
                n = len(imported)
                m = len([v for v in imported.values() if v in b])
                r1 = '%d/%d' % (m, n)
                r2 = '%0.2f%%' % (100 * float(m) / n) if n > 0 else ""
                print '{:{w1}} {:{w2}} {}'.format(p, r1, r2, w1=w1, w2=15)
            except:
                pass
    elif a.action == 'spline':
        print 'Frames without errors:'
        for f in os.listdir(a.output):
            p = os.path.splitext(f)[0]
            try:
                imported = readtxt(a.output, p, a.frame_range, str)
                n = len(imported)
                m = len([v for v in imported.values() if v != "-1"])
                r1 = '%d/%d' % (m, n)
                r2 = '%0.2f%%' % (100 * float(m) / n) if n > 0 else ""
                print '{:{w1}} {:{w2}} {}'.format(p, r1, r2, w1=w1, w2=15)
            except:
                pass

if __name__ == "__main__":
    main(sys.argv[1:])