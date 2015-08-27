#!/usr/bin/env python
########################################################################
#  Author: Nicholas Labello
#  The University of Chicago Research Computing Center
#
#  More to come...
########################################################################

import argparse
from multiprocessing import Pool
import os
import sys
import time

import scipy.misc

import workflow
import libutil

def localmain(IMGDIR,NPROCS):
#  Run the local version without MPI.  Limited to one node.
    FileNames = filter_file_names()
    t1 = time.time()
    
    if NPROCS == 1:
        process_jpgs(FileNames)
    else:
        DataChunks = libutil.get_data_chunks(FileNames,NPROCS)
#        print "LEN: ",len(DataChunks)
        p = Pool(NPROCS)
        x = p.map(process_jpgs,DataChunks)
    t2 = time.time()
    print "Processed images in %4.2f seconds" % (t2-t1)

def mpimain(IMGDIR):
    t1 = time.time()
    comm = MPI.COMM_WORLD
    nprocs = comm.Get_size()
    rank = comm.Get_rank()

    print "Process %i has launched sucessfully" % rank

    if rank == 0:
        FileNames = filter_file_names()
        DataChunks = libutil.get_data_chunks(FileNames,nprocs)
        NDataChunks = len(DataChunks)
        print "Data divided into %i chunks" % NDataChunks
    else:
        FileNames = None
        DataChunks = None

    comm.Barrier()
    
    FileNames = comm.bcast(FileNames,root=0)   
    DataChunks = comm.bcast(DataChunks,root=0)
     
    for r in xrange(0,len(DataChunks)):
        if r == rank:
            process_jpgs(DataChunks[r])
            t2 = time.time()
            nlocalimages = len(DataChunks[r])
            print "Rank %i Processed %i images in %4.2f seconds" % (rank,nlocalimages,t2-t1)

    
def process_jpgs(FileNames):
    for file_pathname in FileNames:
        if FILTER in file_pathname:
            try:
                workflow.process_frame(file_pathname,max_image,THRESH,THRESH_DIFF,NPOINTS,ANGLE_STEP_SIZE,SCRATCHDIR,DOVIZ,VIZDIR,DETMETHOD,WORMWIDTH,WORMLENGTH)
                print file_pathname
            except: # write a dummy pickle file
                filename = file_pathname.split('/')[-1].strip()
                fileprefix = filename.strip('.jpg')
                pickl = open( os.path.join(SCRATCHDIR,fileprefix)+'.pickle','w')
                pickl.write("Error processing this image\n")
                print "Exception on (SKIP): ",file_pathname

                ###raise

def filter_file_names():
    FileNames = libutil.get_file_names_for_framerange(STARTFRAME,ENDFRAME,IMGDIR,filter='.jpg')
    FilteredFileNames = []
    for fn in FileNames:
        if FILTER in fn:
            FilteredFileNames.append(fn)

    FilteredFileNames.sort()
    return FilteredFileNames 



########################################################################
if __name__ == "__main__":

    parser = argparse.ArgumentParser(description=
     '''Process images of C. elegans. 
     Non-MPI Usage: processcelegans.py [Directory of AJP Files]
         MPI Usage: mpirun -np [NP] processcelegans.py [Directory of AJP Files]''',
     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--np',type=int,action='store',dest='NPROCS',default=1,
     help='number of processors (ignored if run with mpirun)')

    parser.add_argument('--imgdir',action='store',dest='IMGDIR',required=True,
     help='directory where input images are stored')

    parser.add_argument('--vizdir',action='store',dest='VIZDIR',default='./vizdir',
     help='directory where color markup images are stored')

    parser.add_argument('--scratchdir',action='store',dest='SCRATCHDIR',default='./scratch',
     help='directory where serialized output is stored while awaiting post-processing')

    parser.add_argument('--bgimage',action='store',dest='BGIMAGE',default='./background.jpg',
     help='background image for subtraction')

    parser.add_argument('--threshdiff',action='store',dest='THRESH_DIFF',default=35,
     help='''Difference between background image intensity and current intensity that
      a pixel must achieve to be kept. The code that does the segmentation is

        DiffImage = ( (max_image - OrigImage) > THRESH_DIFF ) * OrigImage 
      ''')

    parser.add_argument('--thresh',action='store',dest='THRESH',default=130,
     help='''Threshold to identify worm during secondary segmentation
      of the DiffImage.  Value can range betwen 1 and 255 assuming 8-bit grayscale.  
      Corresponds directly to pixel intensity.
      ''')
    parser.add_argument('--npoints',action='store',dest='NPOINTS',default=120,
    help='Number of points along worm midline and sides') 

    parser.add_argument('--angle-step-size',action='store',dest='ANGLE_STEP_SIZE',default=5,
    help='''Number of npoints away from centroid used from angle measurement. 
     npoints/20-40 is a good starting point.''')

    parser.add_argument('--doviz',action='store',dest='DOVIZ',required=False,default='no',
     help='run the visualization code, yes/no.')

    parser.add_argument('--filter',action='store',dest='FILTER',required=False,default='.jpg',
     help='''an optional filter string that must be matched
      by the file name of the images in IMGDIR in order to be
      processed.  e.g., "9.jpg" ''')

    parser.add_argument('--framerange',action='store', nargs=2, dest='FRAMERANGE',required=False,
      default=[-1,9999999],
      help='''Follow with a starting frame and ending frame number to restrict processing.
      --framerange 5000 6000 ''')
    
    parser.add_argument('--detmethod',action='store',dest='DETMETHOD',required=False,default=0,
     help='detection mode, 0 = both methods, 1 = method 1, 2 = method 2')
    
    parser.add_argument('--wormwidth',action='store',dest='WORMWIDTH',required=False,default=20,
     help='worm width in pixel')
    
    parser.add_argument('--wormlength',action='store',dest='WORMLENGTH',required=False,default=1,
     help='worm length as multiple of default length')
    
    args = parser.parse_args()
    NPROCS = int(args.NPROCS)
    IMGDIR = args.IMGDIR
    DOVIZ = args.DOVIZ
    VIZDIR = args.VIZDIR
    SCRATCHDIR = args.SCRATCHDIR
    BGIMAGE = args.BGIMAGE
    THRESH = int(args.THRESH)
    THRESH_DIFF = int(args.THRESH_DIFF)
    NPOINTS = int(args.NPOINTS)
    ANGLE_STEP_SIZE = int(args.ANGLE_STEP_SIZE)
    FILTER = args.FILTER
    STARTFRAME,ENDFRAME = int(args.FRAMERANGE[0]),int(args.FRAMERANGE[1])
    DETMETHOD = int(args.DETMETHOD)
    WORMWIDTH = int(args.WORMWIDTH)
    WORMLENGTH = float(args.WORMLENGTH)

    if STARTFRAME >= ENDFRAME:
        print "Start frame %i >= than End frame %i.  Please correct and re-run." % (STARTFRAME,ENDFRAME)
        sys.exit()

    libutil.makedir(SCRATCHDIR)
    libutil.makedir(VIZDIR)

    max_image = scipy.misc.imread(BGIMAGE)
    try: 
        from mpi4py import MPI
        comm = MPI.COMM_WORLD

        if comm.Get_size() == 1:
            localmain(IMGDIR,NPROCS)
        else:
            mpimain(IMGDIR)
      
    except ImportError:
#        cProfile.run( 'localmain(FileNames,NPROCS) )
#        import pycallgraph
#        pycallgraph.start_trace()
        localmain(IMGDIR,NPROCS)
#        pycallgraph.make_dot_graph('callgraph.png')

