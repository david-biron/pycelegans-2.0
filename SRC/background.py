#!/usr/bin/env python
########################################################################
#  Author: Nicholas Labello, labello@uchicago.edu
#  The University of Chicago Research Computing Center
########################################################################
'''
==================
background
==================
This executable builds a composite background image by taking the 
brightest pixel from a series of images.  By compositing multiple
images the worm, which blocks light and always darkens a pixel,
is effectively removed leaving a background which can be
subtracted from individual images.  Subtracting a background image
results in much cleaner segmentation and subsequent 
identification of worm properties (head, tail, sides, etc.). 

usage: background.py [-h] [--imgdir IMGDIR] [--Nimg NIMG]

optional arguments:
  -h, --help       show this help message and exit
  --imgdir IMGDIR  directory where input images are stored
  --Nimg NIMG      Number of images to use to construct the background. The images
                    are spaced as far apart as possible.  For example, if 100 images are available
                    and Nimg=10, then image 1,11,21,31...91 would be used.
  --windowsize WINDOWSIZE
                        Width of window used for background image generation.
                        This is the number of sequential frames that will use
                        the same background image. For example, if there are
                        10,000 frames, and --windowsize=1000, ten background
                        images will be generated. (default: 1)
'''

import argparse
import os
import numpy
import scipy.misc
import time

import libutil

def main(WINDOWSIZE,NIMG,IMGDIR,startindex,endindex):
    FileNames = libutil.get_file_names_for_framerange(
     startindex,endindex,IMGDIR,filter='.jpg')
    totalimages = len(FileNames)
    lastchunk = totalimages % WINDOWSIZE

    print FileNames
    bg_image_windows = []
    window_start_index = 0
    window_end_index = WINDOWSIZE

# non-overlapping pass through frames
    while window_end_index < totalimages: 
        print window_start_index,window_end_index
        bg_image_windows.append( 
         {'windowstart':window_start_index,'windowend':window_end_index} ) 
        window_start_index += WINDOWSIZE
        window_end_index += WINDOWSIZE
        
# create overlap
    window_start_index = WINDOWSIZE/2 
    window_end_index = window_start_index + WINDOWSIZE
    while window_end_index < totalimages:
        bg_image_windows.append( 
         {'windowstart':window_start_index,'windowend':window_end_index} )
        window_start_index += WINDOWSIZE
        window_end_index += WINDOWSIZE

# handle the last few frames that are left when the number of images is not evenly
# divisible by WINDOWSIZE 
    if lastchunk != 0:
        window_start_index = totalimages - lastchunk
        window_end_index = totalimages - 1
        bg_image_windows.append( 
         {'windowstart':window_start_index,'windowend':window_end_index} )

    print bg_image_windows
# create the background images
    for item in bg_image_windows:
        startindex = item['windowstart']
        endindex = item['windowend']
        build_startend(FileNames,IMGDIR,startindex,endindex)
    
def build_startend(FileNames,IMGDIR,startindex,endindex):
    nametag = 'window_' + str(startindex).zfill(7) + '_' + str(endindex).zfill(7)
    stepsize = int( (endindex-startindex) / NIMG ) 
#  Check to be sure stepsize != 0
    if stepsize == 0:
        print 'skipping window %i - %i, not enough images' % (startindex,endindex)
        return 1
    FilesToUse = FileNames[startindex:endindex:stepsize]
#  Always use the last image in the range
    print endindex
    print len(FileNames)
    print FileNames[-1]
    FilesToUse.append(FileNames[endindex])

    build_image(FilesToUse,nametag)

def build_image(FileNames,tag):
#  eventually add median intensity image method here?  average of top 80%?
#  does anything but max matter or help?  Currently I don't think so.
    build_maxintensity_image(FileNames,tag)

def build_maxintensity_image(FileNames,tag ):
    print "\n Building Max Intensity Image %s from " % tag
    max_image = get_image(FileNames[0])
    for file in FileNames:
        new_image = get_image(file)
        max_image = numpy.maximum(max_image,new_image)
        print (file)
    scipy.misc.imsave('background_max_' + tag +'.jpg',numpy.asarray(max_image,numpy.uint8))
    
def build_medianintensity_image(FileNames):
    pass

def get_image(file_pathname):
    image = scipy.misc.imread(file_pathname,flatten=True)
    return image

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
     '''Builds a composite background image by taking the 
     brightest pixel from a series of images.''',
     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--imgdir',action='store',dest='IMGDIR',default='./input',
     help='directory where input images are stored')

    parser.add_argument('--windowsize',action='store',dest='WINDOWSIZE',default=1200,
     help='''Width of window used for background image generation.  This is the number of 
     sequential frames that will use the same background image. For example, if 
     there are 10,000 frames, and --windowsize=1000, ten background images will be generated.''')

    parser.add_argument('--Nimg',action='store',dest='NIMG',default=10,
     help='''Number of images to use to construct the background. The images
     are spaced as far apart as possible.  For example, if 100 images are available
     and Nimg=10, then image 1,11,21,31...91 would be used.''')

    parser.add_argument('--framerange',action='store', nargs=2, dest='FRAMERANGE',required=False,
      default=[-1,9999999],
      help='''Follow with a starting frame and ending frame number to restrict processing.
      --framerange 5000 6000 ''')

    args = parser.parse_args()
    IMGDIR = args.IMGDIR
    NIMG   = int(args.NIMG)
    WINDOWSIZE = int(args.WINDOWSIZE)
    startindex,endindex = int(args.FRAMERANGE[0]),int(args.FRAMERANGE[1]) 
    main(WINDOWSIZE,NIMG,IMGDIR,startindex,endindex)
