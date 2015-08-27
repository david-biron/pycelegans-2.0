#!/usr/bin/env python
########################################################################
#  Author: Nicholas Labello, labello@uchicago.edu
#  The University of Chicago Research Computing Center
#
#  Input:  root directory with AJP files
#  Output: corresponding directories of images (jpegs)
#
#  python preprocessajp.py
########################################################################
'''
==============
preprocessajp
==============

This executable splits AJP movies into individual jpeg images.  Given a root
directory preprocessajp will traverse down the file system recursively to find
all .ajp files.  Images are named with the movie prefix and a unique frame 
number.

usage: preprocessajp.py [-h] [--np NPROCS] [--imgdir IMGDIR] --ajpdir AJPDIR
                        [--prefix PREFIX] [--nzerom NZEROM] [--nzeroa NZEROA]

Convert AJP movies to images. 
 Usage: preprocessajp.py [Directory of AJP Files]

optional arguments:
  -h, --help       show this help message and exit
  --np NPROCS      number of processors
  --imgdir IMGDIR  directory where output images are stored
  --ajpdir AJPDIR  directory where input AJP movies are found
  --prefix PREFIX  prefix added to jpg file name
  --nzerom NZEROM  width of zero padding on the frame numbers
  --nzeroa NZEROA  width of zero padding on the total frame number, continuous from first movie to last

============== INFO PROVIDED ON FILE FORMAT ==================
The format is:

4 byte int: number of images

Loop number of images

For each image:

4 byte int: number of bytes
Byte array[number of bytes]
Save byte array as jpeg image file

End loop

There is a possibility you will need to swap bytes for the 4 byte ints.
Some environments expect them in reverse order.  If you don't get 100 
[= the number of frames in the sample he sent]
for the first number, try swapping the bytes until you do.

============== INFO PROVIDED ON FILE FORMAT ==================

Reference to byte order swap is big endian v. little endian format.  All systems involved
seem to be big endian.
'''

import argparse
from multiprocessing import Pool
import sys
import struct

import os

def main(NPROCS,IMGDIR,PREFIX,NZEROM,NZEROA,AJPDIR):

    try:
        os.mkdir(IMGDIR)
    except:
        print ('\n\n******************************************************\n' )
        print ('Please ensure that a valid directory path is used,') 
        print ('and that the directory does not already exist.\n') 
        print ('\n******************************************************* \n' )
        raise

    AJPFileNames = get_AJPFileNames(AJPDIR)
    print "Output Log, please verify accuracy: "
    p = Pool(NPROCS)
    x = p.map(ajp_list_to_frames,AJPFileNames)

    sanity_check()

def get_AJPFileNames(AJPDIR):
    AJPFileNames = []
    for root,dirs,files in os.walk(AJPDIR):
        for file in files:
            if '.ajp' in file:
                AJPFileNames.append([os.path.join(root,file)]) 
    return AJPFileNames

def ajp_list_to_frames(AJPList):
    for item in AJPList:
        ajp_to_frames(item)

def ajp_to_frames(ajp_file):
    BYTES_IN_NIMAGES_INTEGER = 4
    BYTES_IN_NBYTES_PER_IMAGE_INTEGER = 4

    data = open(ajp_file,'rb')

    # Read the 4 byte int that gives number of images in the ajp file
    Nimg = struct.unpack('>L',data.read(BYTES_IN_NIMAGES_INTEGER))[0] 

    AllFrames = []

    # Loop over the number of images
    for i in xrange(0,Nimg):
        # Read the 4 byte int that gives number of bytes in an image
        Nbytes_data = data.read(BYTES_IN_NBYTES_PER_IMAGE_INTEGER)
        # convert to an integer
        Nbytes = struct.unpack('>L',Nbytes_data)[0]
        # read the appropriate number of bytes
        AllFrames.append( data.read(Nbytes) )

#  get name of ajp file
    ajp_name = ajp_file.split('/')[-1].strip('.ajp')
    movie_num = int(ajp_name.split('_')[-1].strip('M'))
    frame_number = Nimg * (movie_num - 1)
    first_frame_number = frame_number
#  write each image
    for frame in AllFrames:
  
        image_name = (ajp_name + PREFIX + str(frame_number - first_frame_number).zfill(NZEROM) + 
                      '_' + str(frame_number).zfill(NZEROA) + 
                      '.jpg')
        image_file_name = os.path.join(IMGDIR,image_name)
        imagefile = open(image_file_name,'wb')
        imagefile.write(frame)
        imagefile.close()
        frame_number += 1

    print "AJP M%i: %i (images), frames %i - %i "% (movie_num,Nimg,first_frame_number,frame_number-1)

def sanity_check():
#  count files in output directory
    i = len(os.listdir(IMGDIR))
    print "Sanity check: There are %i files in %s" %(i,IMGDIR)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
     '''Convert AJP movies to images. 
     Usage: preprocessajp.py [Directory of AJP Files]''',
     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('--np',type=int,action='store',dest='NPROCS',default=4,
     help='number of processors')

    parser.add_argument('--imgdir',action='store',dest='IMGDIR',default='./input',
     help='directory where output images are stored')

    parser.add_argument('--ajpdir',action='store',dest='AJPDIR',default='./ajpdir',required=True,
     help='directory where input AJP movies are found')

    parser.add_argument('--prefix',action='store',dest='PREFIX',default='_frame',
     help='prefix added to jpg file name')

    parser.add_argument('--nzerom',action='store',dest='NZEROM',default=5,
     help='width of zero padding on the frame numbers')

    parser.add_argument('--nzeroa',action='store',dest='NZEROA',default=7,
     help='width of zero padding on the total frame number, continuous from first movie to last')

    args = parser.parse_args()
    NPROCS = args.NPROCS
    IMGDIR = args.IMGDIR
    PREFIX = args.PREFIX
    NZEROM  = args.NZEROM
    NZEROA  = args.NZEROA
    AJPDIR = args.AJPDIR

    main(NPROCS,IMGDIR,PREFIX,NZEROM,NZEROA,AJPDIR)
