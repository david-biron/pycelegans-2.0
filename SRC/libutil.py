########################################################################
#  Author: Nicholas Labello, labello@uchicago.edu
#  The University of Chicago Research Computing Center
########################################################################
'''
=========
util
=========
util contains functions that operate on the file system, read data
from the file system, and operate on the aggregate data.  These
functions are outside of the scope of the operations applied to the 
C. elegans image and have more to do with the practical aspects of 
building a cohesive program.  
'''

import os

def makedir(dirname):
    '''
    Make a directory on the file system.
    '''
    try:
        os.mkdir(dirname)
    except OSError:
        pass # directory already exists

def get_file_names_for_framerange(startframe,endframe,dir,filter='.jpg'):
    '''
    Assuming that the Files sort sequentially, return only a range of files
    to be processed
    '''

    filenames = get_file_names(dir,filter)
    Nfiles = len(filenames)

    if startframe < 0: 
        startframe = 0
    if endframe > Nfiles:
        endframe = Nfiles

    rangefilenames = filenames[startframe:endframe]
    return rangefilenames

def get_file_names(dir,filter='.jpg'):
    '''
    Get all file names that contain the string indicated by 
    *filter*.  This function will recursively search down the 
    directory tree starting at *dir*.   
    '''

    FileNames = []
    for root,dirs,files in os.walk(dir):
        for file in files:
            if filter in file:
                FileNames.append(os.path.join(root,file))
    FileNames.sort()
    return FileNames

def get_data_chunks(FileNames,nchunks):
    '''
    Split a list into *nchunks*, return list of lists.
    '''

    NFrames = len(FileNames)
    if (NFrames % nchunks == 0):
        CHUNK_SIZE = (NFrames / nchunks)
    else:
        CHUNK_SIZE = int((NFrames / nchunks)) + 1
    try:
        DataChunks = split_list(FileNames,CHUNK_SIZE)
    except ValueError:
        print FileNames
        print "Number of Files ",len(FileNames)
        print "CHUNK_SIZE ",CHUNK_SIZE
        raise

    return DataChunks

def split_list(l,n=1):
    '''
    split a list l into many lists of size n
    the last item will be less than n if len(l) not evenly divisible
    '''
    return [l[i:i+n] for i in range(0, len(l), n)]
