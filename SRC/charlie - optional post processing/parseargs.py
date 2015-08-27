#!/usr/bin/env python

"""
This module creates the argument parser for accepting command line inputs.

"""

import argparse
import os
import sys

import importdata
import checklength
import flipheadtail
import splinefit

__author__ = "Charlie Wright"
__email__ = "charles.s.wright@gmail.com"

def addimportdata(parser):
    """ Add arguments to parser for importdata.py """
    parser.add_argument("-p", "--properties", 
                        default=os.path.abspath("./properties"),
                        help="""properties directory; default = ./properties""")
    parser.add_argument("-f", "--frame_range", type=int, nargs=2, 
                        default=importdata.FRAME_RANGE, 
                        help="""frame range, given as [first frame, last 
                        frame); default = '0 inf'""")
    parser.add_argument("-v", "--version", type=int, choices=(0, 1), 
                        default=importdata.VERSION, 
                        help="""version of process code used to find midline 
                        (0 = Nick's, 1 = Marc's); default = %d"""
                        % importdata.VERSION)

def addchecklength(parser):
    """ Add arguments to parser for checklength.py """
    parser.add_argument("-ls", "--len_span", type=int, 
                        default=checklength.SPAN, 
                        help="""size of window used to check lengths, over which the average worm length should be constant; default = %d""" % checklength.SPAN)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-lv", "--len_max_val", type=float, 
                       default=checklength.MAX_VAL, 
                       help="""maximum allowed value of length outside of average value in the moving window set by the span, set as a percentage; set this value but not len_max_dev to use an absolute value of the threshold instead of comparing to a deviation about the average; default = %0.2f""" % checklength.MAX_VAL)
    group.add_argument("-ld", "--len_max_dev", type=float, 
                       default=checklength.MAX_DEV, 
                       help="""maximum allowed deviation of length from average value of the standard score in the moving window set by the span; default = %0.2f""" % checklength.MAX_DEV)

def addflipheadtail(parser):
    """ Add arguments to parser for flipheadtail.py """
    parser.add_argument("-fd", "--flip_degree", type=int, 
                        default=flipheadtail.DEGREE, 
                        help="""degree used in polynomial fit to determine correct orientation; default = %d"""
                        % flipheadtail.DEGREE)
    parser.add_argument("-fi", "--flip_iters", type=int, 
                        default=flipheadtail.NUM_ITERS, 
                        help="""number of times to run flip-check procedure; default = %d""" % flipheadtail.NUM_ITERS)
    parser.add_argument("-fs", "--flip_span", type=int, 
                        default=flipheadtail.SPAN, 
                        help="""size of window used to check flips, over which the posture of the worm should be approximately constant; default = %d""" % flipheadtail.SPAN)
    parser.add_argument("-fl", "--flip_min_dev", type=float, 
                        default=flipheadtail.MIN_DEV, 
                        help="""minimum deviation of coefficients used to fit worm midline in moving window, below which the orientation is considered correct; default = %0.2f"""
                        % flipheadtail.MIN_DEV)
    parser.add_argument("-fu", "--flip_max_dev", type=float, 
                        default=flipheadtail.MAX_DEV, 
                        help="""maximum deviation of coefficients used to fit worm midline in moving window, above which the orientation is considered flipped; default = %0.2f"""
                        % flipheadtail.MAX_DEV)

def addsplinefit(parser):
    """ Add arguments to parser for splinefit.py """    
    parser.add_argument("-c", "--use_checks", action='store_true', 
                        help="""include this flag to use info from various post-processing checks when calculating splines""")
    parser.add_argument("-r", "--reorient", action='store_true', 
                        help="""include this flag to reorient tail-to-head any worms marked as incorrectly flipped; this flag will be ignored if the '--use_checks' flag is not used""")
    parser.add_argument("-k", "--keep_sides", action='store_true', 
                        help="""include this flag to read the sidepaths and save a single splined worm outline to file""")
    parser.add_argument("-o", "--output", default="", 
                        help="""name of directory and .mat file to save splined data; default = ./processed""")
    parser.add_argument("-s", "--smoothing", type=int, 
                        default=splinefit.SMOOTHING, 
                        help="""level of smoothing in creating splined version of worm midline; should be of the order of the number of points in the midline; default = %d""" % splinefit.SMOOTHING)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("-n", "--num_points", type=int, 
                       default=splinefit.NUM_POINTS, 
                       help="""number of points in the each splined worm; default = %d""" % splinefit.NUM_POINTS)
    group.add_argument("-sp", "--spacing", type=float, 
                       default=splinefit.SPACING, 
                       help="""distance between consecutive (x, y) points in each splined worm; this value is exclusive of the length and specifying it may result in splines with unequal numbers of points; default = %0.2f""" % splinefit.SPACING)
    parser.add_argument("-l", "--length", type=float, 
                        default=splinefit.LENGTH, 
                        help="""desired length of each splined midline; if specified all splined worms will be approximately the same length (depending on whether a uniform number of points or spacing is given); default = %0.2f""" % splinefit.LENGTH)

def main(inputs):
    p = argparse.ArgumentParser(description="Post-process PyCelegans output.")
    s = p.add_subparsers(help="sub-command help", dest="action")
    a = s.add_parser("check", help="check help")
    b = s.add_parser("spline", help="spline help")
    addimportdata(a)
    addimportdata(b)
    addchecklength(a)
    addflipheadtail(a)
    addsplinefit(b)
    return p.parse_args(inputs)

if __name__ == "__main__":
    print main(sys.argv[1:])