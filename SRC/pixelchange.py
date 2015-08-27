#!/usr/bin/env python
########################################################################
#  Author: Nicholas Labello, labello@uchicago.edu
#  The University of Chicago Research Computing Center
########################################################################
'''
==================
pixelchange
==================
usage: pixelchange.py [-h] [--imgdir IMGDIR] [--scratchdir SCRATCHDIR]
                      [--framerange FRAMERANGE FRAMERANGE]
                      [--changethresh CHANGETHRESH] [--framedelta FRAMEDELTA]
                      [--stepsize STEPSIZE] [--np NPROCS]

This executable creates a vector that describes the number of pixels that
change frame to frame.

optional arguments:
  -h, --help            show this help message and exit
  --imgdir IMGDIR       directory where input images are stored (default:
                        ./input)
  --scratchdir SCRATCHDIR
                        directory where serialized output is stored while
                        awaiting post-processing (default: ./scratch)
  --framerange FRAMERANGE FRAMERANGE
                        Follow with a starting frame and ending frame number
                        to restrict processing. --framerange 5000 6000
                        (default: [-1, 9999999999])
  --changethresh CHANGETHRESH
                        8 bit gray scale pixels are represented by an integer
                        between 0 and 255, where 0 is completely black and 255
                        is bright white. --changethresh is the change in pixel
                        intensity that signifies a pixel has changed. This
                        argument is COMMA SEPARATED ARRAY, e.g, --changethresh
                        10,30,50,70 (default: 50)
  --framedelta FRAMEDELTA
                        How many frames back should the current image be
                        compared to? For example, if the current frame is 100,
                        and framedelta=5, frame #100 will be compared to frame
                        #95 to determine its pixel change value. (default: 1)
  --stepsize STEPSIZE   Only compute pixel change of every Nth frame. Set to
                        '1' to avoid skipping frames. (default: 4)
  --np NPROCS           number of processors (default: 12)

'''
import argparse
import os
import numpy
import scipy.misc
import time
import sys

from multiprocessing import Pool

import libutil
import celegansmetadata

import pylab

def main(STARTFRAME,ENDFRAME):

    t0 = time.time()

    if STARTFRAME <= FRAMEDELTA:    # handle case where no prior frames exist to compare
        STARTFRAME = FRAMEDELTA + 1

    cstartframe = STARTFRAME - FRAMEDELTA
    cendframe =     ENDFRAME - FRAMEDELTA

    FileNames = libutil.get_file_names_for_framerange(STARTFRAME,ENDFRAME,IMGDIR,filter='.jpg')
    CompareFileNames = libutil.get_file_names_for_framerange(cstartframe,cendframe,IMGDIR,filter='.jpg')

    RefAndCompareFileNames = zip(FileNames,CompareFileNames)

    if STEPSIZE >1:
        RefAndCompareFileNames = RefAndCompareFileNames[::STEPSIZE]

    ChunksOfRefAndCompareFileNames = libutil.get_data_chunks(RefAndCompareFileNames,NPROCS)

    p = Pool(NPROCS)
    x = p.map(do_pixel_change_on_batch,ChunksOfRefAndCompareFileNames)

    aggregate_output()

    t1 = time.time()
    print "Time to run on %i cores: %5.2f" % (NPROCS,t1-t0)

#    try:
#        make_plot(alloutput)
#    except:
#        print "Could not make plot, missing libraries or graphics configuration."

def make_plot(alloutput):
    X = [i[0] for i in alloutput]    
    Y = [i[2] for i in alloutput]    
    pylab.plot(X,Y)
    pylab.savefig('pixelchange.%i.%i.%i.png' % (STARTFRAME,ENDFRAME,STEPSIZE) )

def write_output(alloutput,thresh):
    x = open('pixelchange.thresh' + str(thresh) + '.txt','w')
    for item in alloutput:
        x.write( "%8i, %8i, %8i \n" % (item[0],item[1],item[2]) )

def aggregate_output():
    for thresh in ChangeThreshes:
        alloutput = []
        tempoutfiles = libutil.get_file_names(SCRATCHDIR,'.thresh'+str(thresh)+'.pixtemp')
        for file in tempoutfiles:
            lines = open(file,'r').readlines()
            for line in lines:
                words = line.split(',')
                alloutput.append( ( int(words[0]),int(words[1]),int(words[2]) ) )
        alloutput.sort()
        write_output(alloutput,thresh) 

def do_pixel_change_on_batch(RefAndCompare):
    reffilename,refnameprefix,refseqnumber = celegansmetadata.getjpg_metadata(RefAndCompare[0][0])
    nametag = str(refseqnumber)
    for thresh in ChangeThreshes:
        ScratchFile = open(SCRATCHDIR+"/" + nametag + ".thresh" + str(thresh) + ".pixtemp",'w')
        for item in RefAndCompare: 
            refseqnumber,compseqnumber,Npixels_changed = get_Npixels_changed(
             item[0],item[1],thresh)
            ScratchFile.write("%8i,%8i,%8i \n" % (refseqnumber,compseqnumber,Npixels_changed ))

def get_Npixels_changed(reffile_pathname,compfile_pathname,changethresh):
    reffilename,refnameprefix,refseqnumber = celegansmetadata.getjpg_metadata(reffile_pathname)
    compfilename,compnameprefix,compseqnumber = celegansmetadata.getjpg_metadata(compfile_pathname)
    refimage = get_image(reffile_pathname)
    compimage = get_image(compfile_pathname)
    pixchange = (refimage - compimage) >= changethresh 
    Npixels_changed = numpy.count_nonzero(pixchange)
    return (refseqnumber,compseqnumber,Npixels_changed)

def get_image(file_pathname):
    image = scipy.misc.imread(file_pathname,flatten=True)
    return image

#############################################################################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
     '''This executable creates a vector that describes the number of pixels that
        change frame to frame.''',
     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--imgdir',action='store',dest='IMGDIR',default='./input',
     help='directory where input images are stored')

    parser.add_argument('--scratchdir',action='store',dest='SCRATCHDIR',default='./scratch',
     help='directory where serialized output is stored while awaiting post-processing')

    parser.add_argument('--framerange',action='store', nargs=2, dest='FRAMERANGE',required=False,
      default=[-1,9999999999],
      help='''Follow with a starting frame and ending frame number to restrict processing.
      --framerange 5000 6000 ''')

    parser.add_argument('--changethresh',action='store',dest='CHANGETHRESH',default="50",required=False,
     help='''8 bit gray scale pixels are represented by an integer between 0 and 255, where
     0 is completely black and 255 is bright white.  --changethresh is the change in pixel intensity
     that signifies a pixel has changed.  This argument is COMMA SEPARATED ARRAY, e.g, 
     --changethresh 10,30,50,70 ''')

    parser.add_argument('--framedelta',action='store',dest='FRAMEDELTA',default=1,
     help='''How many frames back should the current image be compared to?  For example,
     if the current frame is 100, and framedelta=5, frame #100 will be compared to frame #95 to
     determine its pixel change value.''')

    parser.add_argument('--stepsize',action='store',dest='STEPSIZE',default=4,
     help='''Only compute pixel change of every Nth frame. Set to '1' to avoid skipping frames.''')

    parser.add_argument('--np',type=int,action='store',dest='NPROCS',default=12,
     help='number of processors')

    args = parser.parse_args()
    IMGDIR = args.IMGDIR
    FRAMERANGE = args.FRAMERANGE
    STARTFRAME = int(FRAMERANGE[0])
    ENDFRAME = int(FRAMERANGE[1])
    STEPSIZE = int(args.STEPSIZE)
    CHANGETHRESH = args.CHANGETHRESH
    FRAMEDELTA = int(args.FRAMEDELTA)
    SCRATCHDIR = args.SCRATCHDIR
    NPROCS = int(args.NPROCS)

    ChangeThreshes = []
    for item in CHANGETHRESH.split(','):
        ChangeThreshes.append(int(item))

    libutil.makedir(SCRATCHDIR)
    main(STARTFRAME,ENDFRAME)
