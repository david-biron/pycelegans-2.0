
#=======================================================================================================================
# Coarse grey-level image.
#=======================================================================================================================

#cython: boundscheck=False
#cython: wraparound=False
#cython: cdivision=True


import numpy as np
cimport numpy as cnp

from libc.math cimport floor, sqrt, ceil
from libcpp cimport bool



def coarse_image(cnp.ndarray[dtype=cnp.uint8_t, ndim=2, mode='c'] image, int block_size):
    """
    Creates a coarse version of a gray-scale image.
    
    Parameters
    ----------
    image : array of shape (nrow, ncol)
        Input image.
    block_size : int
        Size of coarse blocks.
        
    Returns
    -------
    out : array of shape (nrow//block_size, ncol//block_size)
        Coarse image.
    """
    
    # get coarse shape
    cdef int nrow = image.shape[0] // block_size
    cdef int ncol = image.shape[1] // block_size
    
    # sum over blocks
    cdef cnp.ndarray[dtype=int, ndim=2, mode='c'] block_sums = np.zeros((nrow, ncol), dtype=np.int32)
    cdef int row, r, col, c
    for row in range(nrow):
        for r in range(block_size):
            for col in range(ncol):
                for c in range(block_size):
                    block_sums[row, col] += image[row*block_size+r, col*block_size+c]
            
    # divide by block area
    cdef cnp.ndarray[dtype=cnp.uint8_t, ndim=2, mode='c'] out = np.empty((nrow, ncol), dtype=np.uint8)
    cdef int block_area = block_size*block_size
    for row in range(nrow):
        for col in range(ncol):
            out[row, col] = block_sums[row, col] // block_area
            
    return out